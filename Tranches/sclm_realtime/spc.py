"""Lecture des flux .spc de FLIMApp en mode Realtime : chaque photon dans sa passe.

Le .spc est le flux FIFO d'une entrée de la SPC-QC-104, traduit au format FIFO_150 de la
SPC-150N par FLIMApp (src/spc/source_qc.jl) : un en-tête de 4 octets, puis des
enregistrements de 32 bits, petit-boutistes :

    bits  0-11  macrotemps (compteur 12 bits, déborde tous les 4096 tics de 2,048 ns)
    bits 12-15  code de routage (photon) ou marqueurs M0-M3 (marqueur)
    bits 16-27  ADC : 4095 − canal de temps (rééchantillonné sur 12,5 ns par FLIMApp)
    bit 28 MARK, 29 GAP, 30 MTOV, 31 INVALID

- INVALID + MTOV sans MARK : plusieurs débordements, leur nombre dans les bits 0-27 ;
- sinon MTOV : un débordement juste avant l'événement ;
- MARK + INVALID : marqueurs ; ni MARK ni INVALID : photon.

Les passes, comme FLIMApp (src/spc/passes.jl, `fin_par_m3 = true`) : une passe commence à
un M0 (front montant du signal de passe) et finit au M3 suivant (son front descendant) ;
un photon pile sur M0 est dans la passe, pile sur M3 il n'y est plus. Le code de routage
dit ce que faisait la NI : pendant le scan, 1 (`CODE_SANS_ROI`, sans ROI) ou l'indice de la
ROI visitée ; 0 (`CODE_HORS_ROI`) pendant les pauses et, avec le PWM du laser opto, pendant
ses phases opto — FLIMApp jette ces photons-là, on les garde ici à part.

Deux lectures du fichier : les marqueurs d'abord (les bornes des passes), puis les photons,
rangés par blocs de mots sans jamais tout garder en mémoire.
"""
import numpy as np

PERIODE_MT = 4096
CHUNK = 1 << 23                      # 8 M enregistrements par bloc (32 Mo)

# genres de photons rangés par passe
ACQ, OPTO, PAUSE = 0, 1, 2           # code ≠ 0 dans la passe ; code 0 dans la passe ; code 0 dans la pause qui suit
N_GENRES = 3


def _blocs(path, chunk=CHUNK):
    """(t, est_photon, est_marqueur, champ 12-15, adc, gap) par bloc ; t en tics depuis le début."""
    base = 0
    with open(path, "rb") as f:
        f.read(4)
        while True:
            w = np.fromfile(f, dtype="<u4", count=chunk)
            if not len(w):
                break
            inv = (w >> 31) & 1
            mtov = (w >> 30) & 1
            mark = (w >> 28) & 1
            multi = (inv == 1) & (mtov == 1) & (mark == 0)
            inc = np.where(multi, (w & 0x0FFFFFFF).astype(np.int64) * PERIODE_MT,
                           mtov.astype(np.int64) * PERIODE_MT)
            cum = base + np.cumsum(inc)
            t = cum + (w & 0xFFF).astype(np.int64)
            base = int(cum[-1])
            photon = ~multi & (mark == 0) & (inv == 0)
            marqueur = ~multi & (mark == 1) & (inv == 1)
            yield (t, photon, marqueur, ((w >> 12) & 0xF).astype(np.uint8),
                   ((w >> 16) & 0xFFF).astype(np.int32), int(((w >> 29) & 1).sum()))


def marqueurs(path):
    """Temps (tics) des M0 et des M3, nombre de photons, de GAP et d'enregistrements."""
    m0, m3, n_ph, n_gap = [], [], 0, 0
    for t, ph, mk, champ, _, gap in _blocs(path):
        m0.append(t[mk & ((champ & 1) == 1)])
        m3.append(t[mk & (((champ >> 3) & 1) == 1)])
        n_ph += int(ph.sum())
        n_gap += gap
    return np.concatenate(m0), np.concatenate(m3), n_ph, n_gap


def bornes(m0, m3):
    """Passes (début, fin) : chaque M0 et le premier M3 qui le suit, avant le M0 suivant."""
    j = np.searchsorted(m3, m0, side="right")
    ok = j < len(m3)
    fin = np.full(len(m0), -1, np.int64)
    fin[ok] = m3[j[ok]]
    suivant = np.r_[m0[1:], np.iinfo(np.int64).max]
    ok &= fin <= suivant
    return m0[ok], fin[ok]


def ranger(path, debut, fin, stim, tics_echantillon, n_phase=200, n_canaux=256):
    """Photons de chaque passe, par genre et par canal de temps ; profils dans la période du PWM.

    `stim` : booléen par passe (laser opto commandé pendant le scan). `tics_echantillon` :
    durée d'un échantillon de la NI en tics de la carte (les phases du PWM en sont faites).

    Rend H (passes, 3, canaux) : ACQ (code ≠ 0 dans la passe), OPTO (code 0 dans la passe),
    PAUSE (code 0 entre la fin de la passe et le début de la suivante) ; profil (2, 2, n_phase) :
    [stim ?, code 0 ?, phase] photons selon leur temps depuis le M0 de leur passe, modulo deux
    échantillons (une période du PWM, n_phase tranches) ; codes (passes, 16) : photons de la
    passe par code de routage (en mode ROI, le code est l'indice de la ROI visitée) ; et les
    compteurs du reste.
    """
    n = len(debut)
    groupe = 4096 // n_canaux
    H = np.zeros(n * N_GENRES * n_canaux, np.int64)
    prof = np.zeros(2 * 2 * n_phase, np.int64)
    codes = np.zeros(n * 16, np.int64)
    autres = dict(code_hors_passe=0, avant_premiere=0, apres_derniere=0, code_inattendu=0)
    stim = np.asarray(stim, bool)
    per = 2 * tics_echantillon
    for t, ph, _, champ, adc, _ in _blocs(path):
        t, code, adc = t[ph], champ[ph], adc[ph]
        canal = (4095 - adc) // groupe
        i = np.searchsorted(debut, t, side="right") - 1          # dernier M0 ≤ t
        avant = i < 0
        ii = np.maximum(i, 0)
        dans = ~avant & (t < fin[ii])
        apres = ~avant & ~dans & (ii == n - 1)
        pause = ~avant & ~dans & ~apres                           # entre fin[i] et debut[i+1]
        autres["avant_premiere"] += int(avant.sum())
        autres["apres_derniere"] += int(apres.sum())
        autres["code_hors_passe"] += int((pause & (code != 0)).sum())
        autres["code_inattendu"] += int((dans & (code > 1)).sum())
        genre = np.full(len(t), -1, np.int64)
        genre[dans & (code != 0)] = ACQ
        genre[dans & (code == 0)] = OPTO
        genre[pause & (code == 0)] = PAUSE
        g = genre >= 0
        idx = (ii[g] * N_GENRES + genre[g]) * n_canaux + canal[g]
        H += np.bincount(idx, minlength=len(H))
        codes += np.bincount(ii[dans] * 16 + code[dans].astype(np.int64), minlength=len(codes))
        # profil dans la période du PWM (photons des passes)
        phase = ((t[dans] - debut[ii[dans]]) % per) * n_phase // per
        k = (stim[ii[dans]].astype(np.int64) * 2 + (code[dans] == 0)) * n_phase + phase.astype(np.int64)
        prof += np.bincount(k, minlength=len(prof))
    return H.reshape(n, N_GENRES, n_canaux), prof.reshape(2, 2, n_phase), codes.reshape(n, 16), autres
