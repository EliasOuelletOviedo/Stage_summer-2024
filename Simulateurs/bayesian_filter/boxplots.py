#!/usr/bin/env python3
"""
boxplots.py — quatre techniques de lecture de [Cl⁻], calibrées puis testées sur les .sdt réels de SCLM_mCer3 (2026-10-16).

  1. les 100 acquisitions de 1 s par concentration sont partagées : numéros impairs (c001, c003, …) pour calibrer,
     pairs pour tester — entrelacées dans le temps, aucun fichier dans les deux ;
  2. calibration sur la somme des fichiers impairs : le modèle du document (step0, sclm_hybride.py) pour le FLIM du
     système, le ratio et l'hybride ; la courbe τ(c) d'une exponentielle unique pour le temps de vie unique ;
     l'hybride ne suppose pas η = 1 : un contexte de 100 × N photons bleus au même état (l'expérience mise en commun)
     fixe η̂ = R / R₁(ĉ_FLIM), puis l'échantillon est lu à η = η̂ ;
  3. échantillons : 1 ou plusieurs fichiers de test tirés au hasard, amincis (binomial) à N photons bleus attendus dans
     la fenêtre (le jaune suit avec la même part : ratio réel). Amincir un comptage de Poisson donne un comptage de
     Poisson : chaque échantillon a la statistique exacte d'une acquisition plus courte, et la variabilité réelle
     d'un fichier à l'autre (décalage d'horloge, dérive) ;
  4. écart = c mesuré − c réel (mM), par technique, concentration et nombre de photons ; un graphique par nombre de photons.

Sorties dans sorties_boxplots/ : boxplots_<N>.png (écart), boxplots_abs_<N>.png (|écart|, échelle log), ecarts_resume.csv, eta_contexte.csv, fig_resume.png, calibration_moitie.{csv,json},
calibration_temps_de_vie.csv, repartition_fichiers.csv, journal.txt.

Avec --eta : chaque échantillon a son propre η, tiré uniformément dans [boxplots] eta_variable (0,9 à 1,1) : ses
photons jaunes, et ceux de son contexte (la même expérience), sont multipliés par η. Le ratio, lu à η = 1, en porte
toute l'erreur ; l'hybride apprend η de son contexte ; les deux lectures par la forme n'en dépendent pas.
Sorties alors dans sorties_boxplots_eta/.

Lancer depuis ce dossier :  python3 boxplots.py  [--eta] [--rapide]      (~5 min ; --rapide : 1 000 échantillons)
"""
import glob, os, re, sys, time
import numpy as np
import pandas as pd
import sdtfile

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ICI)
from src.modele import charger_params, Senseur
from src.estimateurs import Estimateurs
from src.simulation import ETIQUETTES, CONC
from step0_calibration import calibrer
from src import figures

RAPIDE = '--rapide' in sys.argv
ETA_VARIABLE = '--eta' in sys.argv
JOURNAL = []
TECHNIQUES = ['temps de vie unique', 'FLIM du système', 'ratio', 'hybride FRET-FLIM']


def dire(*a):
    s = ' '.join(str(x) for x in a); print(s); JOURNAL.append(s)


def lire_sdt(chemin):
    """Un fichier à deux blocs (un par carte) : voie 1 = bleu, voie 2 = jaune, 256 canaux chacun (cf. step17)."""
    d = sdtfile.SdtFile(chemin).data
    return np.asarray(d[0], float).ravel(), np.asarray(d[1], float).ravel()


def charger_fichiers(dossier):
    """{c: [(numéro, bleu, jaune), …]} pour les 8 concentrations ; « unknown » et irf.sdt écartés."""
    out = {}
    for c, e in zip(CONC, ETIQUETTES):
        nom = f"{c:g}".replace('.', ',')
        fs = sorted(glob.glob(os.path.join(dossier, f"{nom}_c*.sdt")))
        out[c] = [(int(re.search(r'_c(\d+)\.sdt$', f).group(1)), *lire_sdt(f)) for f in fs]
    return out


def main():
    t0 = time.time()
    p = charger_params(); q = p['boxplots']
    sortie = os.path.join(ICI, q['sorties_eta'] if ETA_VARIABLE else q['sorties']); os.makedirs(sortie, exist_ok=True)
    rng = np.random.default_rng(q['graine'])
    n_ech = 1000 if RAPIDE else q['echantillons']

    # ---------------------------------------------------------------- 1. partage des fichiers
    fichiers = charger_fichiers(os.path.join(ICI, q['sdt']))
    ref = pd.read_csv(os.path.join(ICI, p['donnees']['dossier'], p['donnees']['titration']))
    ref.columns = [c.strip() for c in ref.columns]
    a, b = np.where(ref.dans_fenetre_voie1.astype(str).str.strip().str.lower() == 'true')[0][[0, -1]]
    fen = slice(a, b + 1)
    rep, cal_csv, test = [], ref.copy(), {}
    for c, e in zip(CONC, ETIQUETTES):
        tot_B = sum(f[1] for f in fichiers[c]); tot_Y = sum(f[2] for f in fichiers[c])
        assert np.array_equal(tot_B, ref[f'voie1_brut_{e}'].values) and np.array_equal(tot_Y, ref[f'voie2_brut_{e}'].values), \
            f"{c} mM : la somme des .sdt ne redonne pas irf_bruts_mCer3.csv"
        cal = [f for f in fichiers[c] if f[0] % 2 == 1]; tst = [f for f in fichiers[c] if f[0] % 2 == 0]
        cal_csv[f'voie1_brut_{e}'] = sum(f[1] for f in cal); cal_csv[f'voie2_brut_{e}'] = sum(f[2] for f in cal)
        test[c] = (np.array([f[1][fen] for f in tst]), np.array([f[2][fen] for f in tst]))
        rep += [dict(c_mM=c, fichier=f[0], role='calibration' if f[0] % 2 else 'test',
                     photons_bleus_fenetre=f[1][fen].sum(), photons_jaunes_fenetre=f[2][fen].sum()) for f in fichiers[c]]
        dire(f"  {c:5g} mM : {len(cal)} fichiers de calibration, {len(tst)} de test ; photons bleus par fichier de test "
             f"{test[c][0].sum(1).min():.0f} à {test[c][0].sum(1).max():.0f}")
    pd.DataFrame(rep).to_csv(os.path.join(sortie, 'repartition_fichiers.csv'), index=False)
    csv_cal = os.path.join(sortie, 'calibration_moitie.csv'); cal_csv.to_csv(csv_cal, index=False)

    # ---------------------------------------------------------------- 2. calibration sur la moitié impaire
    dire("\ncalibration (modèle du document) sur les fichiers impairs :")
    cal = calibrer(p, csv_cal, os.path.join(sortie, 'calibration_moitie.json'))
    k = cal['constantes']
    dire("  τ1 {t1:.4f}  τ2 {t2:.4f}  τD {tD:.4f}  τA {tA:.4f} ns | λ {lam:.3f}  G {G:.3f}  g2 {G2r:.3f}  κ {k1:.3f}".format(**k)
         + f" | x : K {cal['isothermes']['x']['K']:.1f} mM | déviance/ddl {cal['deviance'] / cal['ddl']:.3f}")
    S = Senseur(cal)
    fond = (np.mean([t['part_fond_bleu'] for t in cal['titration']]), np.mean([t['part_fond_jaune'] for t in cal['titration']]))
    cg = np.concatenate([[0.0], np.geomspace(0.01, 300.0, q['n_c'] - 1)])
    tg = np.geomspace(q['tau_min_ns'], q['tau_max_ns'], q['n_tau'])
    Bc = np.array([cal_csv[f'voie1_brut_{e}'].values[fen] for e in ETIQUETTES], float)
    Yc = np.array([cal_csv[f'voie2_brut_{e}'].values[fen] for e in ETIQUETTES], float)
    E = Estimateurs(S, fond, cg, tg, Bc, Yc, CONC)
    pd.DataFrame(dict(c_mM=CONC, tau_mono_ns=E.tau_cal)).to_csv(os.path.join(sortie, 'calibration_temps_de_vie.csv'), index=False)
    dire(f"  fond : {100 * fond[0]:.2f} % du bleu, {100 * fond[1]:.2f} % du jaune")
    dire("  temps de vie unique τ(c) : " + ", ".join(f"{c:g} mM {t:.3f} ns" for c, t in zip(CONC, E.tau_cal)))
    # les 4 techniques relisent les sommes complètes de calibration (biais « dans l'échantillon ») et de test
    Bs = np.array([test[c][0].sum(0) for c in CONC]); Ys = np.array([test[c][1].sum(0) for c in CONC])
    for nom, (B, Y) in (('calibration', (Bc, Yc)), ('test', (Bs, Ys))):
        eta, _ = E.eta_contexte(B, Y)                                  # la somme sert de son propre contexte
        r = E.toutes(B, Y, eta)
        dire(f"  sommes de {nom:<11s} : η du contexte " + " ".join(f"{e:.3f}" for e in eta))
        dire(f"  sommes de {nom:<11s} relues (écart, mM) : " + " | ".join(
            f"{t} " + " ".join(f"{x - c:+.1f}" for x, c in zip(r[t], CONC)) for t in TECHNIQUES))

    # ---------------------------------------------------------------- 3. échantillons et écarts
    lignes = []; lignes_eta = []; stats = {}; stats_abs = {}
    for N in q['photons_bleus']:
        dire(f"\n{N} photons bleus, {n_ech} échantillons par concentration")
        for c in CONC:
            Bt, Yt = test[c]; F = len(Bt)
            k_f = int(np.ceil(N / (q['part_gardee_max'] * Bt.sum(1).min())))
            if k_f > F:
                raise ValueError(f"{N} photons à {c} mM : il faudrait {k_f} fichiers de test, il y en a {F}")
            choix = np.argsort(rng.random((n_ech, F)), axis=1)[:, :k_f]
            sB, sY = Bt[choix].sum(1), Yt[choix].sum(1)
            part = N / sB.sum(1, keepdims=True)
            # η vrai de chaque échantillon : il multiplie le jaune (part × η ≤ part_gardee_max × 1,1 < 1)
            eta_vrai = rng.uniform(*q['eta_variable'], n_ech) if ETA_VARIABLE else np.ones(n_ech)
            B = rng.binomial(sB.astype(np.int64), part).astype(float)
            Y = rng.binomial(sY.astype(np.int64), part * eta_vrai[:, None]).astype(float)
            # contexte de l'expérience : facteur_contexte × N photons bleus au même état et au même η, Poisson autour de la
            # forme de l'ensemble des fichiers de test ; il fixe η pour l'hybride
            pB, pY = Bt.sum(0), Yt.sum(0); s_ctx = q['facteur_contexte'] * N / pB.sum()
            eta, c_ctx = E.eta_contexte(rng.poisson(s_ctx * pB, (n_ech, len(pB))).astype(float),
                                        rng.poisson(s_ctx * pY * eta_vrai[:, None]).astype(float))
            qe = np.percentile(eta, [2.5, 50, 97.5]); ee = 100 * (eta / eta_vrai - 1)
            qr = np.percentile(ee, [2.5, 50, 97.5])
            lignes_eta.append(dict(photons_bleus=N, photons_bleus_contexte=q['facteur_contexte'] * N, c_mM=c,
                                   eta_q2_5=qe[0], eta_mediane=qe[1], eta_q97_5=qe[2], eta_sd_pct=100 * np.std(eta) / np.mean(eta),
                                   erreur_eta_pct_q2_5=qr[0], erreur_eta_pct_mediane=qr[1], erreur_eta_pct_q97_5=qr[2],
                                   c_contexte_mediane=np.median(c_ctx)))
            r = E.toutes(B, Y, eta)
            for t in TECHNIQUES:
                e = r[t] - c
                qs = np.percentile(e, [2.5, 25, 50, 75, 97.5])
                qa = np.percentile(np.maximum(np.abs(e), figures.PLANCHER_ABS), [2.5, 25, 50, 75, 97.5])
                stats[(N, c, t)] = qs; stats_abs[(N, c, t)] = qa
                lignes.append(dict(photons_bleus=N, c_mM=c, technique=t, echantillons=n_ech, fichiers_par_echantillon=k_f,
                                   photons_bleus_moyens=B.sum(1).mean(), photons_jaunes_moyens=Y.sum(1).mean(),
                                   moyenne=e.mean(), ecart_type=e.std(), rmse=np.sqrt((e ** 2).mean()),
                                   q2_5=qs[0], q25=qs[1], mediane=qs[2], q75=qs[3], q97_5=qs[4],
                                   abs_q2_5=qa[0], abs_q25=qa[1], abs_mediane=qa[2], abs_q75=qa[3], abs_q97_5=qa[4],
                                   part_au_bord_pct=100 * np.mean((r[t] <= 0) | (r[t] >= cg[-1]))))
            dire(f"  {c:5g} mM ({k_f} fichier{'s' if k_f > 1 else ''}) médiane [IQR] : " + " | ".join(
                f"{t} {stats[(N, c, t)][2]:+.2f} [{stats[(N, c, t)][3] - stats[(N, c, t)][1]:.2f}]" for t in TECHNIQUES)
                + f" | η̂/η − 1 {qr[1]:+.2f} % [{qr[0]:+.2f} – {qr[2]:+.2f}]")
        for nom, st, absolu in (('boxplots', stats, False), ('boxplots_abs', stats_abs, True)):
            figures.boxplots(st, N, CONC, TECHNIQUES, n_ech, os.path.join(sortie, f'{nom}_{N}.png'),
                             eta=q['eta_variable'] if ETA_VARIABLE else None, absolu=absolu)
    D = pd.DataFrame(lignes); D.to_csv(os.path.join(sortie, 'ecarts_resume.csv'), index=False, float_format='%.6g')
    pd.DataFrame(lignes_eta).to_csv(os.path.join(sortie, 'eta_contexte.csv'), index=False, float_format='%.6g')
    figures.resume_boxplots(D, TECHNIQUES, os.path.join(sortie, 'fig_resume.png'))
    dire(f"\nterminé en {time.time() - t0:.0f} s")
    with open(os.path.join(sortie, 'journal.txt'), 'w') as f:
        f.write('\n'.join(JOURNAL) + '\n')


if __name__ == '__main__':
    main()
