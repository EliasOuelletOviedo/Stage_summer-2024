#!/usr/bin/env python3
"""
exemples.py — l'exemple de la section 15 (un soma, filtre hybride temps réel et lissé) sous d'autres conditions :
5 nombres de photons par trame × 3 scénarios de η, sur 10 fois plus de trames avec la même forme de courbe (8 mM,
montée vers 22 mM, retour vers 10 mM ; A −30 % au fil de la course). η vrai : 1, 1,1, ou une dérive linéaire de 1 à 1,1
de la première à la dernière trame ; le filtre ne le connaît jamais. C'est celui de run.py, qui suppose η constant
(la marche aléatoire de 0,2 % par trame de s_eta_lisseur est plus fine que le pas de la grille de η, 1,83 % : elle ne
déplace pas la probabilité) : le scénario de dérive montre ce que ce filtre fait d'un η qui bouge.

Sorties dans sorties_exemples/ : fig_exemple_<scénario>_<photons>.png (les deux panneaux de fig_exemple.png),
exemple_<scénario>_<photons>.csv.gz (la table par trame), resume.csv, fig_resume.png, journal.txt.

Lancer depuis ce dossier :  python3 exemples.py  [--rapide]      (~10 min ; --rapide : 300 trames, 2 nombres de photons)
"""
import os, sys, time
import numpy as np
import pandas as pd

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ICI)
from src.modele import charger_params, charger_calibration, Senseur
from src.filtre import Grille, grille_c, grille_eta
from src.exemple import METHODES_EXEMPLE, course_exemple, bilan_exemple, etirer
from src import figures

RAPIDE = '--rapide' in sys.argv
JOURNAL = []


def dire(*a):
    s = ' '.join(str(x) for x in a); print(s); JOURNAL.append(s)


def main():
    t0 = time.time()
    p = charger_params(); q, f = p['exemples'], p['filtre']
    sortie = os.path.join(ICI, q['sorties']); os.makedirs(sortie, exist_ok=True)
    S = Senseur(charger_calibration(p))
    g = Grille(S, grille_c(f), grille_eta(f))
    e0 = etirer(p['exemple'], 1 if RAPIDE else q['facteur_trames'])
    T = e0['trames']
    dire(f"{T} trames : {e0['c_base_mM']:g} mM, montée vers {e0['c_haut_mM']:g} mM à la trame {e0['trame_montee']} "
         f"(τ {e0['tau_montee_trames']:g} trames), retour vers {e0['c_fin_mM']:g} mM à la trame {e0['trame_descente']} "
         f"(τ {e0['tau_descente_trames']:g} trames) ; A {100 * e0['derive_A']:+.0f} %")
    lignes = []
    for i, sc in enumerate(q['scenarios']):
        eta = np.linspace(sc['eta_debut'], sc['eta_fin'], T)
        dire(f"\n{sc['titre']}")
        for j, N in enumerate(q['photons'][:2] if RAPIDE else q['photons']):
            e = dict(e0, photons=N)
            D, _, duree = course_exemple(S, g, f, e, np.random.default_rng(q['graine'] + 100 * i + j), eta)
            nom = f"{sc['nom']}_{N}"
            D.to_csv(os.path.join(sortie, f'exemple_{nom}.csv.gz'), index=False, float_format='%.6g')
            figures.exemple(D, f"{N} photons par trame, {T} trames — {sc['titre']}, inconnu du filtre — A −30 % (photoblanchiment)",
                            os.path.join(sortie, f'fig_exemple_{nom}.png'), eta_lim=(0.8, 1.3))
            b = bilan_exemple(D)
            fin = slice(T - T // 10, T)                          # le dernier dixième : η appris, c sur son plateau final
            b.update(scenario=sc['nom'], photons=N, trames=T,
                     eta_vrai_final=eta[-1], eta_filtre_final_sd_pct=100 * D.eta_filtre_sd.iloc[-1],
                     erreur_eta_lisse_rms_pct=100 * np.sqrt(np.mean((D.eta_lisse / D.eta_vrai - 1) ** 2)),
                     **{f'{col}_erreur_dernier_dixieme': np.sqrt(np.mean((D[col][fin] - D.c_vrai[fin]) ** 2))
                        for _, col in METHODES_EXEMPLE},
                     secondes_filtre_et_lisseur=duree)
            lignes.append(b)
            dire(f"  {N:6d} photons : erreur typique (mM) ratio η = 1 {b['c_ratio_eta1_erreur']:6.2f} | durée de vie seule "
                 f"{b['c_flim_erreur']:6.2f} | filtre {b['c_filtre_erreur']:5.2f} | lissé {b['c_lisse_erreur']:5.2f} ; biais lissé "
                 f"{b['c_lisse_biais']:+.2f} ; η final {b['eta_filtre_final']:.3f} (vrai {eta[-1]:.3f}) ± {b['eta_filtre_final_sd_pct']:.1f} %,"
                 f" η lissé à {b['erreur_eta_lisse_rms_pct']:.1f} % près ; couverture 95 % filtre {100 * b['couverture_filtre']:.0f} %, "
                 f"lissé {100 * b['couverture_lisse']:.0f} %")
    R = pd.DataFrame(lignes)
    R = R[['scenario', 'photons', 'trames'] + [c for c in R.columns if c not in ('scenario', 'photons', 'trames')]]
    R.to_csv(os.path.join(sortie, 'resume.csv'), index=False, float_format='%.6g')
    figures.resume_exemples(R, q['scenarios'], METHODES_EXEMPLE, os.path.join(sortie, 'fig_resume.png'))
    dire(f"\nterminé en {time.time() - t0:.0f} s")
    with open(os.path.join(sortie, 'journal.txt'), 'w') as fj:
        fj.write('\n'.join(JOURNAL) + '\n')


if __name__ == '__main__':
    main()
