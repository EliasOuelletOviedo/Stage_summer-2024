"""Lecture des flux .spc de FLIMApp (Realtime en images) : chaque photon dans son pixel,
son image et son canal de temps.

Le .spc est le flux FIFO d'une carte virtuelle de la SPC-QC-104, traduit au format
FIFO_150 de la SPC-150N par FLIMApp (src/spc/source_qc.jl) : un en-tête de 4 octets,
puis des enregistrements de 32 bits, petit-boutistes :

    bits  0-11  macrotemps (compteur 12 bits, déborde tous les 4096 tics)
    bits 12-15  routage (photon) ou marqueurs M0-M3 (marqueur)
    bits 16-27  ADC : 4095 - canal de temps (le temps croît quand l'ADC décroît)
    bit 28 MARK, 29 GAP, 30 MTOV, 31 INVALID

- INVALID + MTOV sans MARK : plusieurs débordements, leur nombre dans les bits 0-27 ;
- sinon MTOV : un débordement juste avant l'événement ;
- MARK + INVALID : marqueurs (M1 = début de ligne, M2 = début de trame) ;
- ni MARK ni INVALID : photon.

Le rangement reproduit exactement celui de FLIMApp (src/spc/rangement.jl), pour que les
déclins par ROI soient ceux que FLIMApp a analysés en direct :

- une ligne commence à son marqueur M1 ; la première ligne qui suit un M2 est la
  ligne 0 d'une trame ; à temps égal, trame avant ligne, ligne avant photon ;
- colonne = floor(Δt_tics × tic_01 / pix_01) − decalage_pixels, en dixièmes de ns
  **entiers** comme FLIMApp : tic_01 = round(tic × 1e10) = 20, soit un pixel de
  25 tics = 51,2 ns (et non 50 ns) ;
- ligne de l'image = rang dans la trame − lignes ignorées en haut (32 pour 1080
  lignes par trame).

Tout est vectorisé par blocs de mots : un seul passage sur le fichier, sans garder les
photons en mémoire. Les événements du dernier tic d'un bloc sont reportés au bloc
suivant (un marqueur de même tic peut tomber de l'autre côté de la coupure).
"""
from dataclasses import dataclass

import numpy as np

PERIODE_MT = 4096
CHUNK_WORDS = 1 << 23          # 8 M enregistrements par bloc (32 Mo lus à la fois)


@dataclass
class Geometrie:
    """Géométrie de l'image, lue dans <série>_acquisition.ini ([geometrie])."""
    tic_s: float
    temps_pixel_ns: float = 50.0
    pixels_par_ligne: int = 1024
    decalage_pixels: int = 21
    lignes_image: int = 1024
    decalage_lignes: int = 32
    canaux_adc: int = 4096

    @property
    def tic_01(self):
        return int(round(self.tic_s * 1e10))

    @property
    def pix_01(self):
        return int(round(self.temps_pixel_ns * 10))


def lire_ini(path):
    """Fichier .ini de FLIMApp -> {section: {clé: texte}} (les commentaires ';' ignorés)."""
    out, sec = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            sec = line[1:-1]
            out[sec] = {}
        elif "=" in line and sec is not None:
            k, v = line.split("=", 1)
            out[sec][k.strip()] = v.strip()
    return out


def reglage_scanner(texte, lignes_trame):
    """(lignes de l'image, lignes ignorées en haut) du réglage du scanner qui a ce nombre
    de lignes par trame, d'après `reglages_scanner = [(1080, 1024, 32), ...]`."""
    import ast
    for lt, ly, dl in ast.literal_eval(texte):
        if lt == lignes_trame:
            return ly, dl
    return None


def geometrie(ini_path, lignes_trame=1080):
    ini = lire_ini(ini_path)
    a, g = ini["acquisition"], ini["geometrie"]
    ly, dl = int(g["lignes_par_image"]), int(g["decalage_lignes"])
    if ly == 0:
        r = reglage_scanner(g["reglages_scanner"], lignes_trame)
        if r is None:
            raise ValueError(f"aucun réglage du scanner à {lignes_trame} lignes par trame")
        ly, dl = r
    return Geometrie(tic_s=float(a["tic_s"]), temps_pixel_ns=float(g["temps_pixel_ns"]),
                     pixels_par_ligne=int(g["pixels_par_ligne"]),
                     decalage_pixels=int(g["decalage_pixels"]),
                     lignes_image=ly, decalage_lignes=dl)


# ------------------------------------------------------------------ décodage brut
def _decoder(w, base0):
    """Enregistrements -> (t, genre, adc, haut), base de débordement finale.

    genre : 0 = rien (débordements seuls, photon rejeté, inattendu), 1 = photon,
    2 = marqueur. `t` en tics depuis le début du flux.
    """
    inv = (w >> 31) & 1
    mtov = (w >> 30) & 1
    mark = (w >> 28) & 1
    multi = (inv == 1) & (mtov == 1) & (mark == 0)
    inc = np.where(multi, (w & 0x0FFFFFFF).astype(np.int64) * PERIODE_MT,
                   mtov.astype(np.int64) * PERIODE_MT)
    base = base0 + np.cumsum(inc)
    t = base + (w & 0xFFF).astype(np.int64)
    genre = np.zeros(len(w), np.int8)
    genre[~multi & (mark == 0) & (inv == 0)] = 1
    genre[~multi & (mark == 1) & (inv == 1)] = 2
    adc = ((w >> 16) & 0xFFF).astype(np.int16)
    haut = ((w >> 12) & 0xF).astype(np.int8)
    return t, genre, adc, haut, (int(base[-1]) if len(base) else base0)


def marqueurs(spc_path, chunk=CHUNK_WORDS):
    """Temps (tics) des débuts de ligne (M1) et de trame (M2), et nombre de photons."""
    lignes, trames, n_ph = [], [], 0
    base = 0
    with open(spc_path, "rb") as f:
        f.read(4)
        while True:
            w = np.fromfile(f, dtype="<u4", count=chunk)
            if not len(w):
                break
            t, genre, _, haut, base = _decoder(w, base)
            m = genre == 2
            lignes.append(t[m & ((haut >> 1) & 1 == 1)])
            trames.append(t[m & ((haut >> 2) & 1 == 1)])
            n_ph += int((genre == 1).sum())
    return np.concatenate(lignes), np.concatenate(trames), n_ph


# ------------------------------------------------------------------ rangement
class Rangement:
    """Parcourt un .spc et rend, bloc de mots par bloc de mots, les photons rangés :
    trame (0 = la première après le premier M2), ligne y, colonne x, canal de temps
    (0..canaux-1, temps croissant), temps t (tics). Seuls les photons qui tombent dans
    l'image sont rendus ; `trames_t` donne les temps (tics) des débuts de trame du bloc.

    Usage :  for ph in Rangement(path, geo, canaux=256): ph["trame"], ph["y"], ...
    """

    def __init__(self, spc_path, geo: Geometrie, canaux=256, chunk=CHUNK_WORDS):
        self.path, self.geo, self.canaux, self.chunk = spc_path, geo, canaux, chunk
        self.groupe = geo.canaux_adc // canaux

    def __iter__(self):
        g = self.geo
        base = 0
        report = None                     # événements du dernier tic, reportés
        n_trames = 0                      # trames commencées avant ce bloc
        rang_report = -1                  # lignes déjà vues dans la trame en cours
        ligne_prec = None                 # (t, trame, rang) de la dernière ligne vue
        with open(self.path, "rb") as f:
            f.read(4)
            fin = False
            while not fin:
                w = np.fromfile(f, dtype="<u4", count=self.chunk)
                fin = len(w) == 0
                if not fin:
                    t, genre, adc, haut, base = _decoder(w, base)
                    keep = genre > 0
                    t, genre, adc, haut = t[keep], genre[keep], adc[keep], haut[keep]
                else:
                    t = np.zeros(0, np.int64); genre = np.zeros(0, np.int8)
                    adc = np.zeros(0, np.int16); haut = np.zeros(0, np.int8)
                if report is not None:
                    t, genre, adc, haut = (np.concatenate([a, b]) for a, b in
                                           zip(report, (t, genre, adc, haut)))
                    report = None
                if not len(t):
                    continue
                if not fin:                # le dernier tic attend le bloc suivant
                    cut = np.searchsorted(t, t[-1], side="left")
                    report = (t[cut:], genre[cut:], adc[cut:], haut[cut:])
                    t, genre, adc, haut = t[:cut], genre[:cut], adc[:cut], haut[:cut]
                    if not len(t):
                        continue

                m = genre == 2
                L = t[m & ((haut >> 1) & 1 == 1)]
                F = t[m & ((haut >> 2) & 1 == 1)]
                # trame de chaque ligne : trames de temps <= ligne (à égalité, la trame d'abord)
                tr_L = n_trames + np.searchsorted(F, L, side="right") - 1
                # rang dans la trame : index de la ligne moins celui de la 1re ligne de sa trame
                idx = np.arange(len(L))
                premier = np.searchsorted(tr_L, tr_L, side="left")
                rang_L = idx - premier
                cur = n_trames - 1
                if rang_report >= 0:
                    rang_L = np.where(tr_L == cur, rang_L + rang_report + 1, rang_L)
                rang_L = np.where(tr_L >= 0, rang_L, -1)

                ph = genre == 1
                tp, adcp = t[ph], adc[ph]
                k = np.searchsorted(L, tp, side="right") - 1     # ligne <= photon
                avant = k < 0
                if ligne_prec is not None:
                    tl = np.where(avant, ligne_prec[0], L[np.maximum(k, 0)] if len(L) else 0)
                    trp = np.where(avant, ligne_prec[1], tr_L[np.maximum(k, 0)] if len(L) else 0)
                    rgp = np.where(avant, ligne_prec[2], rang_L[np.maximum(k, 0)] if len(L) else 0)
                    ok = np.ones(len(tp), bool)
                else:
                    ok = ~avant
                    kk = np.maximum(k, 0)
                    tl = L[kk] if len(L) else np.zeros(len(tp), np.int64)
                    trp = tr_L[kk] if len(L) else np.full(len(tp), -1)
                    rgp = rang_L[kk] if len(L) else np.full(len(tp), -1)
                # Une ligne vue après un M2 mais avant le M2 suivant appartient à la trame
                # courante ; une ligne avant le tout premier M2 a trp = -1 : ignorée.
                y = rgp - g.decalage_lignes
                x = (tp - tl) * g.tic_01 // g.pix_01 - g.decalage_pixels
                ok &= (trp >= 0) & (rgp >= 0) & (y >= 0) & (y < g.lignes_image) \
                    & (x >= 0) & (x < g.pixels_par_ligne)
                canal = (g.canaux_adc - 1 - adcp[ok].astype(np.int32)) // self.groupe

                # état reporté au bloc suivant
                n_trames += len(F)
                if len(L):
                    ligne_prec = (int(L[-1]), int(tr_L[-1]), int(rang_L[-1]))
                    rang_report = int(rang_L[-1]) if tr_L[-1] == n_trames - 1 else -1
                elif len(F):
                    rang_report = -1        # nouvelle trame, pas encore de ligne
                yield {"trame": trp[ok].astype(np.int64), "y": y[ok].astype(np.int32),
                       "x": x[ok].astype(np.int32), "canal": canal.astype(np.int32),
                       "t": tp[ok], "trames_t": F}
