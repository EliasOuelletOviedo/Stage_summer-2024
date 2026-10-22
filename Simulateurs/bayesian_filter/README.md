# Filtre bayésien hybride (ratio + FLIM) — SuperClomeleon

Construit le 2026-10-22 d'après « Mesurer le chlore avec SuperClomeleon : FLIM, ratio et filtre bayésien expliqués.pdf »
(dans ce dossier). Le filtre des sections 9 à 14 tourne sur des données simulées d'après la titration SCLM_mCer3 du
2026-10-16 (8 [Cl⁻] de 0 à 131 mM, deux canaux, 80 MHz) : le ratio fait la mesure, la durée de vie apprend η, la
transmission du jaune relative au bleu, propre à chaque tranche.

Écrit en Python, comme le module de calibration qui accompagne ces données (`sclm_hybride.py`, dans
`../FLIM-FRET_2-channels/sorties/declins_purs_2026-10-16_SCLM_mCer3/`).

## Lancer

```bash
pip install -r requirements.txt            # numpy, scipy, pandas, matplotlib
# une seule fois : ajuster la calibration sur les histogrammes mCer3 → calibration.json (~1 min)
python3 step0_calibration.py
# les vérifications (13)
python3 test/run_tests.py
# la course complète (~2,5 min)
python3 run.py                             # --rapide : moins de répétitions (~25 s)
```

## Données : ce qui vient des mesures

Les 900 `.sdt` de `../FLIM-FRET_2-channels/Data/2026-10-16/SCLM_mCer3/` ne sont pas versionnés. Ce dossier part de
leurs histogrammes sommés, versionnés par `step17_bruts_mCer3.jl` (`irf_bruts_mCer3.csv`, `inconnu_bruts_mCer3.csv`) :

- **calibration** (`step0`) : IRF paramétrique, ajustement conjoint de Poisson des 16 histogrammes, isothermes. Elle
  retrouve l'annexe du document : τ₁ 0,0868, τ₂ 0,854, τ_D 3,755, τ_A 3,483 ns ; λ 0,764, G 2,34, g₂ 0,758, κ 0,658 ;
  x_max 0,548, x_min 0,017, K 62,5 mM ; déviance/ddl 2,341. Elle donne aussi ce que le document ne dit pas : **le fond
  vaut 2,3 % des photons bleus de la fenêtre** (jaune ≈ 0), et le décalage d'horloge de chaque acquisition.
- **générateur « modèle »** (`src/simulation.py`) : Poisson autour du modèle calibré, pour toute trajectoire de [Cl⁻] ;
  η, l'amplitude A et sa dérive, le fond et les perturbations de la section 16 s'y règlent sans que le filtre les connaisse.
- **générateur « amincissement »** : Poisson autour des histogrammes **mesurés**, mis à l'échelle de la trame (jaune × η).
  Aucune hypothèse de modèle : fond, IRF réelle, décalages d'horloge et défauts du modèle sont ceux des données. Seules
  les 8 concentrations de la titration sont disponibles (trajectoires en paliers).

## Fichiers

| Fichier | Rôle | lignes |
| --- | --- | --- |
| `params.toml` | le seul endroit pour les constantes (grilles, transition, scénarios, graine) | 92 |
| `step0_calibration.py` | le pont unique avec les données : appelle `sclm_hybride.py`, écrit `calibration.json` | 77 |
| `calibration.json` | la calibration ajustée (constantes, isothermes, IRF, fonds et décalages par acquisition) ; écrit par step0 | — |
| `src/modele.py` | modèle direct μ_B, μ_Y (sections 3 à 7) : IRF ⊗ déclins périodiques, isothermes, η, A ; perturbations | 130 |
| `src/filtre.py` | vraisemblance découpée ratio / FLIM (section 11), grille (c, η), transition, filtre et lisseur (13), temps réel à M ROI, portes et bascule (14), lectures de référence | 190 |
| `src/information.py` | Fisher et Cramér-Rao, nuisances par complément de Schur (section 10) | 50 |
| `src/simulation.py` | les deux générateurs, les trajectoires des sections 14 et 15 | 91 |
| `src/figures.py` | rendu | 108 |
| `run.py` | la course en 8 parties ; un CSV par partie, `journal.txt`, 5 figures | 394 |
| `test/run_tests.py` | 13 vérifications : calibration = annexe, découpage exact, portes (93,7 % de l'information FLIM gardée), η retrouvé, bascule, couverture du lisseur, amincissement | 165 |

## Sorties (`sorties/`)

`journal.txt` (tout ce qui suit, chiffres compris), `verifications.csv`, `cramer_rao.csv`, `cramer_rao_verification.csv`,
`exemple.csv` (la trace de la section 15, trame par trame), `exemple_repetitions.csv`, `temps_reel.csv`,
`temps_reel_resume.csv`, `amincissement.csv`, `amincissement_resume.csv`, `amincissement_trace.csv`,
`titration_relue.csv`, `robustesse.csv`, `inconnu.csv`, `eta_transition.csv` ; figures `fig_donnees`, `fig_exemple`,
`fig_cramer_rao`, `fig_temps_reel`, `fig_amincissement`.

## Résultats, comparés au document

**Ce qui se reproduit à l'identique.** Les exemples chiffrés des sections 6, 7, 9, 11 et 12 (R₁(c), ln L = −5 697,8 /
−5 688,3 / −5 710,7, la « ligne » (c, η) du ratio, l'a posteriori 10,1 [7,6 – 12,7] mM) ; tout le tableau de
Cramér-Rao à 10 mM et 10 000 photons (forme bleue 10,19, jaune 39,94, deux canaux 9,88, ratio 1,28, ratio + formes
1,27 mM) ; sd(ln η) = 1,78 % pour 10⁶ photons (1,8 %) ; 400 tirages : ratio ± 1,26 mM (1,25), durée de vie seule
12,9 ± 11,5 mM (12,8, ≈ 11) ; les biais de la section 16 (η = 1,2 : −8,73 contre −8,8 ; accepteurs inactifs : +8,2 /
+9,0 / +8,9 contre +8,2 / +9,0 / +8,8 ; autofluorescence du jaune sur le ratio −2,74 contre −2,7) ; la référence
« η connu » de la section 14 (1,72 / 1,17 / 0,78 mM contre 1,7 / 1,1 / 0,8).

**Section 15 (un soma, η = 1,25 inconnu, 5 000 photons/trame, A −30 % pendant la course).** Moyenne sur 20 répétitions
(document entre parenthèses) :

| méthode | biais (mM) | erreur typique (mM) |
| --- | --- | --- |
| ratio avec la calibration (η = 1) | −10,9 (−11,5) | 11,2 (11,8) |
| durée de vie seule, trame par trame | +13,6 | 44,6 (41) |
| filtre hybride, temps réel | −0,1 ± 1,4 (+0,3) | 2,1 ± 0,6 (1,4) |
| filtre hybride, lissé après coup | −0,3 ± 1,0 (−0,5) | 1,1 ± 0,4 (1,0) |

η final 1,244 ± 0,021 ; couverture des intervalles à 95 % : 96 % (filtre), 98 % (lissé). Le document montre un seul
tirage ; ici, celui de `fig_exemple.png` (répétition 1) finit à η = 1,222, d'où un biais de −1,1 mM sur ce tirage-là.

**Section 14 (300 trames, η = 1,25 inconnu, 3 répétitions).** Erreur typique après la trame 100, en mM :

| photons par trame | 1 ROI | 10 ROI partageant η | η connu |
| --- | --- | --- | --- |
| 1 000 | 3,76 (2,6) | 2,35 (2,5) | 1,72 (1,7) |
| 3 000 | 2,26 (1,5) | 1,59 (1,4) | 1,17 (1,1) |
| 10 000 | 1,38 (1,5) | **1,47 (0,9)** | 0,78 (0,8) |

Coût : 0,6 ms par ROI et par trame sur la grille 260 × 61, 0,2 à 0,3 ms avec les portes et la bascule.

## Ce que la simulation apprend en plus du document

1. **η n'est pas appris que du FLIM.** Le score ℓ_m(η) de la section 14 absorbe aussi l'écart entre la dynamique supposée
   de c et la vraie : un η statique peut « expliquer » une partie d'un changement de ratio. Partie 8, 1 000 photons,
   10 ROI, η vrai 1,25 : c statique et transition du document → η = 1,18 (−6 %, ≈ 5 sd a posteriori) ; c qui répond et
   transition trop raide (s = 0,01, sans sauts) → η = 1,45. Dans un essai hors course, avec c figé et K = identité,
   η revenait à 1,24–1,26.
   Le module de référence (`sclm_hybride.py`, `demo_temps_reel`) fait la même chose : 1,21 à 1 000 photons pour trois
   graines, et 1,231 à 1,250 à 3 000 photons (le 1,250 ± 1,0 % du document est la graine 0).
2. **La bascule à 1 % fige η trop tôt.** À 10 000 photons, sd(ln η) passe sous 1 % vers la trame 40, pendant la ligne de
   base, alors que l'effet ci-dessus tire encore η vers le bas : η est figé à 1,21–1,23 dans les trois répétitions, et
   10 ROI font moins bien qu'une seule (1,47 contre 1,38 mM). La sd a posteriori ne voit pas ce biais.
3. **Le fond doit être dans le modèle.** Les données en ont 2,3 % dans le bleu. Ignoré, il fausse la forme du déclin :
   η +6 % et +3 à +8 mM sur les trames amincies (partie 5), +4,9 mM pour l'hybride à 10 mM (partie 6). Le filtre accepte un
   fond connu (photons par case et par trame), soustrait des totaux et ajouté aux formes.
4. **Le plancher de la calibration.** Relus en entier (11 à 40 millions de photons), les histogrammes de la titration
   donnent la bonne [Cl⁻] par le ratio à 0,5 mM près, mais le FLIM, seul juge de η, y lit η entre 0,991 et 1,046 selon la
   concentration (1,046 à 20 mM). Sur les trames amincies avec fond connu, le filtre garde ainsi −1,6 / −0,3 / −2,7 mM
   de biais à 10 / 20 / 40 mM (η final 1,238 ± 1,7 %), contre −9,7 à −21 mM pour le ratio lu à η = 1.
5. **Le décalage d'horloge.** +10 ps non modélisé : +3,5 mM pour l'hybride, +7,1 mM pour le FLIM bleu, rien pour le ratio.
   Les acquisitions de la titration ont dérivé de −5 à +15 ps : il faut mesurer le décalage de chaque session.
6. **L'échantillon « inconnu »** du même jour (16,5 + 33,2 millions de photons) : ratio (η = 1) 37,7 mM, FLIM 38,6 mM,
   hybride 36,3 [35,9 – 37,7] mM avec η = 0,988 ± 0,8 %, cohérent avec η = 1 dans les conditions de la calibration.

## Pas fait

- L'autofluorescence n'est pas modélisée dans le filtre (le document la donne corrigée « si on la modélise ») : la
  colonne « hybride » de la partie 6 en montre l'effet sans modèle, avec des formes d'autofluorescence illustratives.
- Ni la ré-estimation des durées de vie en tranche, ni l'estimation du décalage d'horloge par le filtre, ni le temps mort.
- Les remèdes aux points 1 et 2 (apprendre η sur une ligne de base avec une transition ajustée, bascule sur un critère
  qui tienne compte du biais) restent à choisir.
