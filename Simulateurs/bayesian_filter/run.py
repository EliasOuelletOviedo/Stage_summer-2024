#!/usr/bin/env python3
"""
run.py — la course : le filtre bayésien hybride (ratio + FLIM) sur des données simulées d'après une titration
(SCLM_mCer3 par défaut ; SCLM régulier avec --params params_SCLM.toml).

  1. vérifications : les exemples chiffrés des sections 6, 7, 9, 11 et 12 du document ;
  2. information   : bornes de Cramér-Rao (section 10) et leur vérification par 400 tirages ;
  3. exemple       : le soma de la section 15 (η = 1,25 inconnu, 5 000 photons/trame), filtre et lisseur ;
  4. temps réel    : le tableau de la section 14 (1 ROI, 10 ROI partageant η, η connu ; 1 000 à 10 000 photons) ;
  5. amincissement : le filtre sur des trames tirées des histogrammes MESURÉS (fond et défauts du modèle réels),
                     puis la titration complète relue (le plancher de la calibration) ;
  6. robustesse    : les biais de la section 16 à 10 mM, plus le décalage d'horloge et le fond des données ;
  7. inconnu       : l'échantillon « inconnu » mesuré le même jour, lu par le ratio, le FLIM et le filtre ;
  8. η et le modèle de transition : ce que la dynamique supposée de c fait à l'estimation de η.
Sorties dans le dossier [donnees] sorties (sorties/, sorties_SCLM/) : un CSV par partie, les figures, journal.txt.

Lancer depuis ce dossier (après step0_calibration.py) :  python3 run.py  [--params f.toml] [--rapide]   (~2,5 min ; --rapide ~25 s)
"""
import os, sys, time
import numpy as np
import pandas as pd

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ICI)
from src.modele import charger_params, chemin_params, charger_calibration, Senseur, x_vers_c
from src.filtre import (Grille, grille_c, grille_eta, portes_log, transition, filtre_lisseur, FiltreTempsReel,
                        resume_c, resume_eta)
from src.simulation import (trame_modele, trajectoire_exemple, trajectoires_roi, Mesures, inconnu, fond_moyen)
from src.information import bornes, sigma_eta
from src import figures

P = charger_params(chemin_params(sys.argv))
SORTIES = os.path.join(ICI, P['donnees']['sorties'])
RAPIDE = '--rapide' in sys.argv
JOURNAL = []


def dire(*a):
    s = ' '.join(str(x) for x in a); print(s); JOURNAL.append(s)


def ecrire(nom, lignes):
    pd.DataFrame(lignes).to_csv(os.path.join(SORTIES, nom + '.csv'), index=False, float_format='%.6g')


def noyaux(g, f):
    return transition(np.log(g.c), f['s'], f['s_saut'], f['w']), transition(np.log(g.eta), f['s_eta_lisseur'])


# ================================================================================================ 1. vérifications
def verifications(S, g):
    v = []
    def ligne(section, quantite, doc, calcule):
        v.append(dict(section=section, quantite=quantite, document=doc, calcule=calcule))
    iso = S.iso['x']
    ligne(6, 'c pour x = 0,475 (mM)', 10.0, x_vers_c(0.475, iso))
    for c, r in zip((0, 5, 10, 20, 40, 131), (3.61, 3.21, 2.92, 2.50, 2.02, 1.33)):
        ligne(7, f'R1({c} mM)', r, float(S.ratio_calibration(c)))
    NB, NY = 2560.0, 7440.0
    for c, doc in ((5, -5697.8), (10, -5688.3), (20, -5710.7)):
        R = float(S.ratio_calibration(c)); pY = R / (1 + R)
        ligne(9, f'ln L({c} mM), N_B = 2560, N_Y = 7440', doc, NY * np.log(pY) + NB * np.log(1 - pY))
    cf = np.geomspace(0.5, 300, 20001); R = S.ratio_calibration(cf); pY = R / (1 + R)
    ligne(9, 'maximum de vraisemblance (mM)', 10.2, cf[np.argmax(NY * np.log(pY) + NB * np.log(1 - pY))])
    for e, doc in ((0.9, 4.8), (1.0, 10.2), (1.1, 16.0), (1.2, 22.4), (1.3, 29.2)):
        L = NY * np.log(e * R / (1 + e * R)) + NB * np.log(1 / (1 + e * R))
        ligne(11, f'c indiqué par le ratio si η = {e} (mM)', doc, cf[np.argmax(L)])
    # section 12 : a posteriori sur la grille fine (a priori uniforme sur la grille), η = 1, totaux seuls
    L = g.ratio(NB, NY, 1.0)[:, 0]; p = np.exp(L - L.max()); p /= p.sum()
    m, lo, hi = resume_c(p, g.c)
    ligne(12, 'moyenne a posteriori (mM)', 10.1, float(m[0]))
    ligne(12, 'borne basse à 95 % (mM)', 7.6, float(lo[0])); ligne(12, 'borne haute à 95 % (mM)', 12.7, float(hi[0]))
    for x in v:
        dire(f"  § {x['section']:2d} {x['quantite']:<46s} document {x['document']:>9.2f}   calculé {x['calcule']:>9.2f}")
    ecrire('verifications', v)


# ================================================================================================ 2. information
def information(S, g, p, rng):
    q = p['information']; N = q['photons']; lignes = []
    for c in q['concentrations_mM']:
        lignes.append(dict(c_mM=c, **bornes(S, c, N)))
    ecrire('cramer_rao', lignes)
    b = bornes(S, 10.0, N)
    doc = {'forme bleu seule': 10.2, 'forme jaune seule': 39.9, 'formes 2 canaux (eta inconnu)': 9.9,
           'ratio (eta connu)': 1.28, 'ratio + formes (eta connu)': 1.27}
    for k, x in b.items():
        dire(f"  σ_c à 10 mM, {N:.0f} photons, {k:<30s} {x:6.2f} mM (document {doc[k]:5.2f}) ; "
             f"photons pour 1 mM ≈ {N * x ** 2:,.0f}".replace(',', ' '))
    se = sigma_eta(S, 10.0, 1e6)
    dire(f"  sd(ln η) pour 10⁶ photons mis en commun à 10 mM : {100 * se:.2f} % (document 1,8 %)")
    # vérification par simulation (400 tirages à 10 mM) : ratio (η connu) et durée de vie seule
    c0, M = q['c_verif_mM'], (100 if RAPIDE else q['tirages'])
    r, f = [], []
    for _ in range(M):
        B, Y = trame_modele(S, c0, N, rng)
        r.append(g.lire_ratio(B, Y)); f.append(g.lire_flim(B, Y))
    mc = dict(tirages=M, c_mM=c0, photons=N, ratio_moyenne=np.mean(r), ratio_sd=np.std(r),
              flim_moyenne=np.mean(f), flim_sd=np.std(f), flim_mediane=np.median(f),
              borne_ratio=b['ratio (eta connu)'], borne_flim=b['formes 2 canaux (eta inconnu)'], sd_ln_eta_1e6=se)
    ecrire('cramer_rao_verification', [mc])
    dire(f"  {M} tirages à {c0:g} mM : ratio {mc['ratio_moyenne']:.2f} ± {mc['ratio_sd']:.2f} mM (document ± 1,25) ; "
         f"durée de vie seule {mc['flim_moyenne']:.1f} ± {mc['flim_sd']:.1f} mM (document 12,8, ≈ 11)")
    return lignes


# ================================================================================================ 3. exemple (section 15)
METHODES_EXEMPLE = (('ratio avec la calibration (η = 1)', 'c_ratio_eta1'), ('durée de vie seule, trame par trame', 'c_flim'),
                    ('filtre hybride, temps réel', 'c_filtre'), ('filtre hybride, lissé après coup', 'c_lisse'))


def course_exemple(S, g, p, rng):
    """Une course de la section 15 : trames, filtre et lisseur sur la grille complète. Renvoie (table par trame, H, durée)."""
    e, f = p['exemple'], p['filtre']; T = e['trames']
    C = trajectoire_exemple(e)
    A = 1 + e['derive_A'] * np.arange(T) / (T - 1)            # photoblanchiment : le filtre ne le connaît pas
    H = [trame_modele(S, C[k], e['photons'] * A[k], rng, e['eta']) for k in range(T)]
    t0 = time.perf_counter()
    Ls = np.array([g.loglik(B, Y) for B, Y in H])
    Kc, Ke = noyaux(g, f)
    alpha, post = filtre_lisseur(Ls, Kc, Ke)
    duree = time.perf_counter() - t0
    out = []
    for k in range(T):
        B, Y = H[k]
        fm, flo, fhi = resume_c(alpha[k].sum(1), g.c); lm, llo, lhi = resume_c(post[k].sum(1), g.c)
        ef, ef_sd = resume_eta(alpha[k].sum(0), g.eta); el, el_sd = resume_eta(post[k].sum(0), g.eta)
        out.append(dict(trame=k + 1, c_vrai=C[k], N_B=B.sum(), N_Y=Y.sum(), R=Y.sum() / B.sum(),
                        c_ratio_eta1=g.lire_ratio(B, Y), c_flim=g.lire_flim(B, Y),
                        c_filtre=fm[0], filtre_bas=flo[0], filtre_haut=fhi[0], eta_filtre=ef, eta_filtre_sd=ef_sd,
                        c_lisse=lm[0], lisse_bas=llo[0], lisse_haut=lhi[0], eta_lisse=el, eta_lisse_sd=el_sd))
    return pd.DataFrame(out), H, duree


def bilan_exemple(D):
    b = {}
    for nom, col in METHODES_EXEMPLE:
        err = D[col] - D.c_vrai
        b[f'{col}_biais'] = err.mean(); b[f'{col}_erreur'] = np.sqrt((err ** 2).mean())
    b['couverture_filtre'] = np.mean((D.filtre_bas <= D.c_vrai) & (D.c_vrai <= D.filtre_haut))
    b['couverture_lisse'] = np.mean((D.lisse_bas <= D.c_vrai) & (D.c_vrai <= D.lisse_haut))
    b['eta_filtre_final'] = D.eta_filtre.iloc[-1]; b['eta_lisse'] = D.eta_lisse.iloc[0]
    return b


def exemple(S, g, p, rng):
    e = p['exemple']; T = e['trames']
    D, H, duree = course_exemple(S, g, p, rng)
    ecrire('exemple', D.to_dict('records'))
    B, Y = H[0]
    dire(f"  trame 1 : N_B = {B.sum():.0f}, N_Y = {Y.sum():.0f}, R = {Y.sum() / B.sum():.2f} (R₁(0) = {g.R1_zero:.2f}) ; "
         f"ratio à η = 1 → {g.lire_ratio(B, Y):.1f} mM ; FLIM seul → {g.lire_flim(B, Y):.1f} mM")
    L1 = g.ratio(B.sum(), Y.sum())
    dire("  ligne du ratio : " + " ; ".join(f"η {g.eta[j]:.2f} → {g.c[np.argmax(L1[:, j])]:.1f} mM"
                                            for j in [np.argmin(abs(g.eta - x)) for x in (1.1, 1.2, 1.25, 1.3, 1.4)]))
    for k in (1, 10, 50, 70, 100, 240):
        r = D.iloc[k - 1]
        dire(f"  trame {k:3d} : vrai {r.c_vrai:5.1f} | filtre {r.c_filtre:5.1f} [{r.filtre_bas:5.1f} – {r.filtre_haut:5.1f}]"
             f" η {r.eta_filtre:.3f} ± {100 * r.eta_filtre_sd:4.1f} % | lissé {r.c_lisse:5.1f} [{r.lisse_bas:5.1f} – {r.lisse_haut:5.1f}]")
    b1 = bilan_exemple(D)
    dire(f"  répétition 1 (la trace ci-dessus) : η final temps réel {D.eta_filtre.iloc[-1]:.3f} ± {100 * D.eta_filtre_sd.iloc[-1]:.1f} %,"
         f" lissé {b1['eta_lisse']:.3f} (vrai {e['eta']}) ; {1e3 * duree / T:.1f} ms par trame (grille complète, avec lissage)")
    reps = [b1] + [bilan_exemple(course_exemple(S, g, p, rng)[0]) for _ in range(1, 4 if RAPIDE else e['repetitions'])]
    R = pd.DataFrame(reps); R.insert(0, 'repetition', np.arange(1, len(R) + 1)); ecrire('exemple_repetitions', R.to_dict('records'))
    doc = dict(c_ratio_eta1=(-11.5, 11.8), c_flim=(np.nan, 41), c_filtre=(0.3, 1.4), c_lisse=(-0.5, 1.0))
    dire(f"  bilan, répétition 1 | moyenne ± sd sur {len(R)} répétitions | (document)")
    for nom, col in METHODES_EXEMPLE:
        dire(f"  {nom:<38s} biais {b1[col + '_biais']:+6.2f} | {R[col + '_biais'].mean():+6.2f} ± {R[col + '_biais'].std():4.2f} "
             f"({'—' if np.isnan(doc[col][0]) else f'{doc[col][0]:+.1f}'}) ; erreur typique {b1[col + '_erreur']:5.2f} | {R[col + '_erreur'].mean():5.2f} ± "
             f"{R[col + '_erreur'].std():4.2f} ({doc[col][1]:.1f}) mM")
    dire(f"  η final temps réel {R.eta_filtre_final.mean():.3f} ± {R.eta_filtre_final.std():.3f}, lissé {R.eta_lisse.mean():.3f} ;"
         f" couverture à 95 % : filtre {100 * R.couverture_filtre.mean():.0f} %, lissé {100 * R.couverture_lisse.mean():.0f} %"
         f" (minimum {100 * R.couverture_lisse.min():.0f} %)")
    return D


# ================================================================================================ 4. temps réel (section 14)
def temps_reel(S, g_portes, g_cases, p, rng):
    q, f = p['temps_reel'], p['filtre']; T = q['trames']; M = q['n_roi']; deb = q['trame_debut_erreur']
    Kc = transition(np.log(g_portes.c), f['s'], f['s_saut'], f['w'])
    reps = 1 if RAPIDE else q['repetitions']; lignes = []
    for N in q['photons']:
        for r in range(reps):
            C = trajectoires_roi(rng, M, T)
            H = [[trame_modele(S, C[m, k], N, rng, q['eta']) for m in range(M)] for k in range(T)]
            variantes = {'1 ROI': [FiltreTempsReel(g_portes, 1, Kc) for _ in range(M)],
                         f'{M} ROI partageant η': [FiltreTempsReel(g_portes, M, Kc)],
                         'η connu d\'avance': [FiltreTempsReel(g_cases, M, Kc, eta_fixe=q['eta'])]}
            for nom, filtres in variantes.items():
                est = np.zeros((M, T)); duree = 0.0; etas = []
                for k in range(T):
                    B = np.array([H[k][m][0] for m in range(M)]); Y = np.array([H[k][m][1] for m in range(M)])
                    t0 = time.perf_counter()
                    if len(filtres) == M:
                        res = [fl.mettre_a_jour(B[m:m + 1], Y[m:m + 1]) for m, fl in enumerate(filtres)]
                        est[:, k] = [x['c'][0] for x in res]; etas = [(x['eta'], x['eta_sd']) for x in res]
                    else:
                        x = filtres[0].mettre_a_jour(B, Y); est[:, k] = x['c']; etas = [(x['eta'], x['eta_sd'])]
                    duree += time.perf_counter() - t0
                rmse = np.sqrt(np.mean((est[:, deb:] - C[:, deb:]) ** 2))
                lignes.append(dict(photons=N, repetition=r + 1, variante=nom, rmse_mM=rmse,
                                   eta_moyen=np.mean([x[0] for x in etas]), eta_sd_pct=100 * np.mean([x[1] for x in etas]),
                                   ms_par_roi_trame=1e3 * duree / T / M))
    ecrire('temps_reel', lignes)
    D = pd.DataFrame(lignes).groupby(['photons', 'variante'], sort=False).mean(numeric_only=True).reset_index()
    doc = {(1000, 0): 2.6, (1000, 1): 2.5, (1000, 2): 1.7, (3000, 0): 1.5, (3000, 1): 1.4, (3000, 2): 1.1,
           (10000, 0): 1.5, (10000, 1): 0.9, (10000, 2): 0.8}
    for _, x in D.iterrows():
        j = list(D.variante.unique()).index(x.variante)
        dire(f"  {x.photons:6.0f} photons, {x.variante:<22s} RMSE {x.rmse_mM:5.2f} mM (document {doc.get((x.photons, j), np.nan):.1f})"
             f" | η {x.eta_moyen:.3f} ± {x.eta_sd_pct:.1f} % | {x.ms_par_roi_trame:.2f} ms/ROI/trame")
    ecrire('temps_reel_resume', D.to_dict('records'))
    return D


# ================================================================================================ 5. amincissement des mesures
def amincissement(S, g, p, cal, rng):
    q, f = p['amincissement'], p['filtre']; T = q['trames']
    mes = Mesures(os.path.join(ICI, p['donnees']['dossier']), p['donnees']['titration'], cal)
    n_pal = len(q['paliers_mM']); lp = T // n_pal
    C = np.repeat(q['paliers_mM'], lp); C = np.concatenate([C, np.full(T - len(C), C[-1])])
    fb, fy = fond_moyen(cal); n = S.n_fen
    Kc = transition(np.log(g.c), f['s'], f['s_saut'], f['w'])
    reps = 3 if RAPIDE else q['repetitions']; lignes = []; trace = None
    for r in range(reps):
        H = [mes.trame(C[k], q['photons'], rng, q['eta']) for k in range(T)]
        fonds = [(fb * B.sum() / n, fy * Y.sum() / n) for B, Y in H]
        filtres = {'filtre, fond connu': FiltreTempsReel(g, 1, Kc), 'filtre, fond ignoré': FiltreTempsReel(g, 1, Kc)}
        est = {k: np.zeros(T) for k in ('ratio (η = 1)', 'durée de vie seule', *filtres)}; eta = {}
        for k, (B, Y) in enumerate(H):
            est['ratio (η = 1)'][k] = g.lire_ratio(B, Y, fond=fonds[k])
            est['durée de vie seule'][k] = g.lire_flim(B, Y, fond=fonds[k])
            for nom, fl in filtres.items():
                x = fl.mettre_a_jour(B, Y, fonds[k] if 'connu' in nom else None)
                est[nom][k] = x['c'][0]; eta[nom] = (x['eta'], x['eta_sd'])
        if r == 0:
            trace = dict(trame=np.arange(1, T + 1), c_vrai=C, **est)
        for nom, x in est.items():
            for i, c in enumerate(q['paliers_mM']):
                sl = slice(i * lp + lp // 2, (i + 1) * lp)                 # seconde moitié de chaque palier
                lignes.append(dict(repetition=r + 1, methode=nom, palier_mM=c, biais_mM=np.mean(x[sl] - c),
                                   sd_mM=np.std(x[sl] - c), eta_final=eta.get(nom, (np.nan,))[0],
                                   eta_sd_pct=100 * eta.get(nom, (np.nan, np.nan))[1]))
    ecrire('amincissement', lignes); ecrire('amincissement_trace', pd.DataFrame(trace).to_dict('records'))
    D = pd.DataFrame(lignes).groupby(['methode', 'palier_mM'], sort=False).agg(
        biais_mM=('biais_mM', 'mean'), sd_mM=('sd_mM', 'mean'), eta_final=('eta_final', 'mean'),
        eta_final_sd_entre_rep=('eta_final', 'std')).reset_index()
    dire(f"  trames tirées des histogrammes mesurés : {q['photons']} photons, η = {q['eta']} inconnu, fond mesuré "
         f"{100 * fb:.1f} % du bleu, {reps} répétitions ; erreur sur la seconde moitié de chaque palier")
    for _, x in D.iterrows():
        dire(f"  {x.methode:<22s} {x.palier_mM:5.0f} mM : biais {x.biais_mM:+6.2f} mM, sd {x.sd_mM:5.2f} mM"
             + (f" | η final {x.eta_final:.3f} (entre répétitions ± {100 * x.eta_final_sd_entre_rep / x.eta_final:.1f} %)"
                if np.isfinite(x.eta_final) else ''))
    ecrire('amincissement_resume', D.to_dict('records'))
    return D, trace


def titration_relue(S, p, cal):
    """
    Le plancher de la calibration : chaque histogramme complet de la titration (11 à 40 millions de photons) relu
    par le ratio (η = 1), par le FLIM, et par le maximum (c, η) de la vraisemblance hybride, selon le fond utilisé.
    Un η ≠ 1 ici est l'erreur que le FLIM, seul juge de η, ferait à cette concentration avec un nombre infini de trames.
    """
    mes = Mesures(os.path.join(ICI, p['donnees']['dossier']), p['donnees']['titration'], cal)
    cf = np.geomspace(0.5, 300, 6001); gf = Grille(S, cf, np.geomspace(0.6, 1.8, 1201))
    fb, fy = fond_moyen(cal); n = S.n_fen; lignes = []
    for i, c in enumerate(mes.c):
        B, Y, t = mes.B[i], mes.Y[i], cal['titration'][i]
        for nom, fd in (('fond ajusté', (t['fond_bleu_par_case'], t['fond_jaune_par_case'])),
                        ('part moyenne du fond', (fb * B.sum() / n, fy * Y.sum() / n)), ('sans fond', None)):
            L = gf.loglik(B, Y, fd); j = np.unravel_index(np.argmax(L), L.shape)
            lignes.append(dict(c_mM=c, fond=nom, c_ratio_eta1=gf.lire_ratio(B, Y, fond=fd), c_flim=gf.lire_flim(B, Y, fd),
                               c_hybride=cf[j[0]], eta_hybride=gf.eta[j[1]]))
    ecrire('titration_relue', lignes)
    D = pd.DataFrame(lignes)
    for nom, x in D.groupby('fond', sort=False):
        dire(f"  {nom:<21s} ratio (η = 1) : écart max {np.max(np.abs(x.c_ratio_eta1 - x.c_mM)):5.2f} mM | η du FLIM : "
             + ' '.join(f"{e:.3f}" for e in x.eta_hybride) + f"  (à {', '.join(f'{c:g}' for c in x.c_mM)} mM)")


# ================================================================================================ 8. η et le modèle de transition
def eta_transition(S, g, p, rng):
    """
    Le score ℓ(η) ne recueille pas que l'information FLIM : il absorbe aussi l'écart entre la dynamique supposée de c
    et la vraie. À 1 000 photons, 10 ROI : η final selon le noyau K, pour un c statique et des réponses (section 14).
    """
    q, f = p['temps_reel'], p['filtre']; M, T, N = q['n_roi'], q['trames'], q['photons'][0]
    lc = np.log(g.c); lignes = []
    K = {'document (s = 0,04, w = 0,03)': transition(lc, f['s'], f['s_saut'], f['w']),
         's = 0,04 sans sauts': transition(lc, f['s']), 's = 0,01 sans sauts': transition(lc, 0.01)}
    for r in range(1 if RAPIDE else 4):
        for mode in ('c statique (10 mM)', 'réponses (section 14)'):
            C = np.full((M, T), 10.0) if mode.startswith('c statique') else trajectoires_roi(rng, M, T)
            H = [[trame_modele(S, C[m, k], N, rng, q['eta']) for m in range(M)] for k in range(T)]
            for nom, Kc in K.items():
                fl = FiltreTempsReel(g, M, Kc)
                for k in range(T):
                    x = fl.mettre_a_jour(np.array([h[0] for h in H[k]]), np.array([h[1] for h in H[k]]))
                lignes.append(dict(repetition=r + 1, trajectoire=mode, transition=nom, eta_final=x['eta'], eta_sd_pct=100 * x['eta_sd']))
    ecrire('eta_transition', lignes)
    D = pd.DataFrame(lignes).groupby(['trajectoire', 'transition'], sort=False).agg(
        eta=('eta_final', 'mean'), entre_rep=('eta_final', 'std'), sd=('eta_sd_pct', 'mean')).reset_index()
    for _, x in D.iterrows():
        dire(f"  {x.trajectoire:<22s} {x.transition:<30s} η final {x.eta:.3f} (vrai {q['eta']}) ; sd a posteriori {x.sd:.1f} %"
             + (f", entre répétitions {100 * x.entre_rep / x.eta:.1f} %" if np.isfinite(x.entre_rep) else ''))


# ================================================================================================ 6. robustesse (section 16)
def robustesse(S, p, cal):
    q = p['robustesse']; c0, N = q['c_mM'], q['photons']
    cf = np.geomspace(0.5, 300, 4001); ef = np.geomspace(0.6, 1.8, 2401)       # toute la grille du filtre, pas de 0,16 %
    gf = Grille(S, cf, ef)
    def lire(B, Y):
        L = gf.formes(B, Y)[:, None] + gf.ratio(B.sum(), Y.sum())
        return dict(ratio=gf.lire_ratio(B, Y), flim_bleu=float(cf[np.argmax(gf.formes(B, Y, canaux=('B',)))]),
                    hybride=float(cf[np.unravel_index(np.argmax(L), L.shape)[0]]))
    def attendu(Sv, eta=1.0, af=None):
        b, y = Sv.mu(c0, 1.0, eta); A = N / (b.sum() + y.sum()); b, y = A * b, A * y
        if af is not None:
            pb, py = af                                     # parts de l'autofluorescence dans chaque canal
            b = b + pb / (1 - pb) * b.sum() * S.decroissance(q['autofluo_bleu'], 'bleu')
            y = y + py / (1 - py) * y.sum() * S.decroissance(q['autofluo_jaune'], 'jaune')
        return b, y
    fb, _ = fond_moyen(cal)
    cas = [('aucune', attendu(S), (0, 0, 0)),
           (f"η = {q['eta']} (diffusion, profondeur)", attendu(S, q['eta']), (-8.8, 0, 0)),
           (f"autofluorescence = {100 * q['autofluo_part']:.0f} % du jaune", attendu(S, af=(0, q['autofluo_part'])), (-2.7, 0, None)),
           (f"autofluorescence = {100 * q['autofluo_part']:.0f} % des deux canaux",
            attendu(S, af=(q['autofluo_part'], q['autofluo_part'])), (0, -9.6, None)),
           (f"durées de vie {100 * q['taus_plus_courts']:.0f} % plus courtes",
            attendu(Senseur(cal, echelle_tau=1 - q['taus_plus_courts'])), (0, -8.4, -7.2)),
           (f"{100 * q['accepteurs_inactifs']:.0f} % d'accepteurs inactifs en plus",
            attendu(Senseur(cal, part_active=1 - q['accepteurs_inactifs'])), (8.2, 9.0, 8.8)),
           (f"décalage d'horloge de {q['decalage_ps']:+g} ps non modélisé", attendu(Senseur(cal, decalage_ps=q['decalage_ps'])), (None, None, None))]
    b, y = attendu(S)
    cas.append((f"fond mesuré non modélisé ({100 * fb:.1f} % du bleu)", (b + fb / (1 - fb) * b.sum() / len(b), y), (None, None, None)))
    lignes = []
    for nom, (B, Y), doc in cas:
        r = lire(B, Y)
        lignes.append(dict(perturbation=nom, ratio_calibration_fixe=r['ratio'] - c0, flim_bleu_taus_fixes=r['flim_bleu'] - c0,
                           hybride_eta_libre=r['hybride'] - c0, doc_ratio=doc[0], doc_flim_bleu=doc[1], doc_hybride=doc[2]))
        fmt = lambda x: '   —  ' if x is None else f'{x:+6.1f}'
        dire(f"  {nom:<44s} ratio {r['ratio'] - c0:+6.2f} ({fmt(doc[0])}) | FLIM bleu {r['flim_bleu'] - c0:+6.2f} "
             f"({fmt(doc[1])}) | hybride {r['hybride'] - c0:+6.2f} ({fmt(doc[2])}) mM")
    ecrire('robustesse', lignes)


# ================================================================================================ 7. échantillon inconnu
def echantillon_inconnu(S, g, p, cal):
    B, Y = inconnu(os.path.join(ICI, p['donnees']['dossier']), p['donnees']['inconnu'], cal)
    fb, fy = fond_moyen(cal); n = S.n_fen
    fond = (fb * B.sum() / n, fy * Y.sum() / n)
    L = g.loglik(B, Y, fond); P = np.exp(L - L.max()); P /= P.sum()
    m, lo, hi = resume_c(P.sum(1), g.c); e, esd = resume_eta(P.sum(0), g.eta)
    r = dict(N_B=B.sum(), N_Y=Y.sum(), R=Y.sum() / B.sum(), c_ratio_eta1=g.lire_ratio(B, Y, fond=fond),
             c_flim_bleu=g.lire_flim(B, Y, fond, ('B',)), c_flim_2_canaux=g.lire_flim(B, Y, fond),
             c_hybride=m[0], hybride_bas=lo[0], hybride_haut=hi[0], eta_hybride=e, eta_sd_pct=100 * esd)
    ecrire('inconnu', [r])
    dire(f"  inconnu mesuré : N_B = {r['N_B']:.3g}, N_Y = {r['N_Y']:.3g}, R = {r['R']:.3f} | ratio (η = 1) {r['c_ratio_eta1']:.1f} mM |"
         f" FLIM bleu {r['c_flim_bleu']:.1f}, 2 canaux {r['c_flim_2_canaux']:.1f} mM | hybride {r['c_hybride']:.1f} "
         f"[{r['hybride_bas']:.1f} – {r['hybride_haut']:.1f}] mM, η = {r['eta_hybride']:.3f} ± {r['eta_sd_pct']:.1f} %")
    return r


def main():
    t0 = time.time(); os.makedirs(SORTIES, exist_ok=True)
    p = P; cal = charger_calibration(p); f = p['filtre']
    rng = lambda i: np.random.default_rng(p['simulation']['graine'] + i)
    S = Senseur(cal)
    g = Grille(S, grille_c(f), grille_eta(f))
    gp = Grille(S, grille_c(f), grille_eta(f), portes_log(S.n_fen, f['portes']))
    dire(f"{p['donnees']['nom']} — calibration : {cal['source']} ; déviance/ddl {cal['deviance'] / cal['ddl']:.3f}" + (' ; --rapide' if RAPIDE else ''))
    dire("\nLes valeurs « document » sont celles du PDF, calibré sur SCLM_mCer3 : pour un autre senseur, elles ne servent"
         "\nque de repère.")
    dire("\n1. Vérifications (exemples chiffrés du document)"); verifications(S, g)
    dire("\n2. Information de Fisher et borne de Cramér-Rao (section 10)"); cr = information(S, g, p, rng(2))
    dire("\n3. Exemple complet, trame par trame (section 15)"); ex = exemple(S, g, p, rng(3))
    dire("\n4. Version temps réel (section 14)"); tr = temps_reel(S, gp, g, p, rng(4))
    dire("\n5. Trames tirées des histogrammes mesurés (amincissement)"); am, trace = amincissement(S, g, p, cal, rng(5))
    dire("   titration complète relue (plancher de la calibration)"); titration_relue(S, p, cal)
    dire("\n6. Robustesse à 10 mM (section 16) : biais calculé (document)"); robustesse(S, p, cal)
    dire("\n7. Échantillon inconnu mesuré le 2026-10-16"); echantillon_inconnu(S, g, p, cal)
    dire("\n8. η selon le modèle de transition (1 000 photons, 10 ROI)"); eta_transition(S, g, p, rng(8))
    mes = Mesures(os.path.join(ICI, p['donnees']['dossier']), p['donnees']['titration'], cal)
    figures.donnees(S, mes, cal, p, rng(7), os.path.join(SORTIES, 'fig_donnees.png'))
    figures.exemple(ex, p, os.path.join(SORTIES, 'fig_exemple.png'))
    figures.cramer_rao(cr, p['donnees']['nom'], os.path.join(SORTIES, 'fig_cramer_rao.png'))
    figures.temps_reel(tr, os.path.join(SORTIES, 'fig_temps_reel.png'))
    figures.amincissement(trace, am, p['donnees']['nom'], os.path.join(SORTIES, 'fig_amincissement.png'))
    dire(f"\nterminé en {time.time() - t0:.0f} s")
    with open(os.path.join(SORTIES, 'journal.txt'), 'w') as fj:
        fj.write('\n'.join(JOURNAL) + '\n')


if __name__ == '__main__':
    main()
