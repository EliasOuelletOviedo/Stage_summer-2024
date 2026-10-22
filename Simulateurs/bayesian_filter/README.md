# Filtre bayésien hybride (ratio + FLIM) — SuperClomeleon

Construit le 2026-10-22 d'après « Mesurer le chlore avec SuperClomeleon : FLIM, ratio et filtre bayésien expliqués.pdf »
(dans ce dossier). Le filtre des sections 9 à 14 tourne sur des données simulées d'après la titration SCLM_mCer3 du
2026-10-16 (8 [Cl⁻] de 0 à 131 mM, deux canaux, 80 MHz) : le ratio fait la mesure, la durée de vie apprend η, la
transmission du jaune relative au bleu, propre à chaque tranche.

Écrit en Python, comme le module de calibration qui accompagne ces données (`sclm_hybride.py`, dans
`../FLIM-FRET_2-channels/sorties/declins_purs_2026-10-16_SCLM_mCer3/`).

## Lancer

```bash
pip install -r requirements.txt            # numpy, scipy, pandas, matplotlib, sdtfile
# une seule fois : ajuster la calibration sur les histogrammes mCer3 → calibration.json (~1 min)
python3 step0_calibration.py
# les vérifications (13)
python3 test/run_tests.py
# la course complète (~2,5 min)
python3 run.py                             # --rapide : moins de répétitions (~25 s)

# le second jeu : SCLM régulier (800 nm) → calibration_SCLM.json, sorties_SCLM/ (~2 + 2,5 min)
python3 step0_calibration.py --params params_SCLM.toml
python3 run.py --params params_SCLM.toml

# comparaison de 4 techniques sur les .sdt réels de SCLM_mCer3 → sorties_boxplots/ (~4 min, --rapide ~1 min)
python3 boxplots.py
```

Un jeu de données = un fichier de paramètres : seul le bloc `[donnees]` diffère entre `params.toml` (SCLM_mCer3) et
`params_SCLM.toml` (SCLM régulier) ; filtre, scénarios et graines sont les mêmes, pour comparer les deux senseurs.

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
| `params.toml` | le seul endroit pour les constantes (données, grilles, transition, scénarios, graine, boxplots) ; SCLM_mCer3 | 118 |
| `params_SCLM.toml` | les mêmes, pour la titration SCLM régulière à 800 nm | 94 |
| `step0_calibration.py` | le pont unique avec les données : appelle le `sclm_hybride.py` du jeu, écrit sa calibration (`calibrer` sert aussi à boxplots.py) | 86 |
| `calibration.json`, `calibration_SCLM.json` | la calibration ajustée (constantes, isothermes, IRF, fonds et décalages par acquisition) ; écrites par step0 | — |
| `src/modele.py` | modèle direct μ_B, μ_Y (sections 3 à 7) : IRF ⊗ déclins périodiques, isothermes, η, A ; perturbations | 137 |
| `src/filtre.py` | vraisemblance découpée ratio / FLIM (section 11), grille (c, η), transition, filtre et lisseur (13), temps réel à M ROI, portes et bascule (14), lectures de référence | 190 |
| `src/information.py` | Fisher et Cramér-Rao, nuisances par complément de Schur (section 10) | 50 |
| `src/simulation.py` | les deux générateurs (lit les colonnes voie1/voie2 de step17 et bleu/jaune de step16), les trajectoires des sections 14 et 15 | 97 |
| `src/figures.py` | rendu | 151 |
| `src/estimateurs.py` | les 4 lectures d'un échantillon (temps de vie unique, FLIM du système, ratio, hybride à η du contexte), en bloc | 98 |
| `boxplots.py` | partage des .sdt, calibration sur une moitié, échantillons et contextes tirés de l'autre, boxplots | 158 |
| `run.py` | la course en 8 parties ; un CSV par partie, `journal.txt`, 5 figures | 398 |
| `test/run_tests.py` | 13 vérifications : calibration = annexe, découpage exact, portes (93,7 % de l'information FLIM gardée), η retrouvé, bascule, couverture du lisseur, amincissement | 165 |

## Sorties (`sorties/` pour SCLM_mCer3, `sorties_SCLM/` pour SCLM régulier)

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
12,9 ± 11,5 mM (12,8, ≈ 11) ; les biais de la section 16 (η = 1,2 : −8,76 contre −8,8 ; accepteurs inactifs : +8,2 /
+9,0 / +8,8 contre +8,2 / +9,0 / +8,8 ; autofluorescence du jaune sur le ratio −2,74 contre −2,7 ; autofluorescence des
deux canaux sur le FLIM bleu −9,50 contre −9,6 ; durées de vie 5 % plus courtes, FLIM bleu −8,39 contre −8,4) ; la référence
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
5. **Le décalage d'horloge.** +10 ps non modélisé : +3,6 mM pour l'hybride, +7,1 mM pour le FLIM bleu, rien pour le ratio.
   Les acquisitions de la titration ont dérivé de −5 à +15 ps : il faut mesurer le décalage de chaque session.
6. **L'échantillon « inconnu »** du même jour (16,5 + 33,2 millions de photons) : ratio (η = 1) 37,7 mM, FLIM 38,6 mM,
   hybride 36,3 [35,9 – 37,7] mM avec η = 0,988 ± 0,8 %, cohérent avec η = 1 dans les conditions de la calibration.

## Second jeu : SCLM régulier (800 nm), `sorties_SCLM/`

Mêmes scénarios, mêmes graines ; données `../FLIM-FRET_2-channels/sorties/declins_purs_2026-10-16_SCLM/`
(`irf_bruts_purs_800nm.csv`, `inconnu_irf_bruts_800nm.csv`, de `step16_declins_purs.jl` ; 64 ps par canal, une période
de 195 canaux, fenêtre de 182 canaux). Le fichier à 850 nm s'utiliserait en changeant deux noms dans `params_SCLM.toml`.

**La calibration ajuste bien, mais ne se lit pas comme celle de mCer3.** Déviance/ddl 2,00, résidus à ± 4 (`fig_donnees`).
τ₁ 0,063, τ₂ 0,723, τ_D 2,715, τ_A 3,260 ns ; λ 0,720, G 1,59, g₂ 0,254, κ 0,154 ; fond 3,8 % du bleu et 1,8 % du jaune.
Le donneur n'est pas mCerulean3 : τ_D = 2,72 ns (le τ₃ = 2,72 ns de l'analyse FLIM-FRET de ces données), et la
« population 2 » **croît** avec le chlore (f₂ de 0,33 à 0,55 entre 0 et 131 mM, isotherme ajustée de 0,335 vers 0,643) :
c'est sans doute une composante propre au donneur, pas une population en FRET. Le modèle reste une calibration empirique
valable de la forme et du ratio en fonction de [Cl⁻] sur 0 à 131 mM — c'est tout ce dont le filtre a besoin — mais
x = f₁ + f₂ n'y est plus la fraction en FRET (x de 0,864 à 0,641, x_min ajusté 0,475). Ce qui suppose cette lecture
physique est donc sans objet pour ce senseur : la vérification § 6, et la perturbation « 10 % d'accepteurs inactifs »
(qui multiplie f₁ et f₂ par 0,9 : +14,9 / +59,7 / +28,6 mM). Les valeurs « document » du journal sont celles de mCer3.

**Ce que change le senseur.** La forme du déclin renseigne davantage : à 10 mM et 10 000 photons, σ_c = 5,6 mM par les
formes (9,9 pour mCer3), 1,19 mM par le ratio (1,28) ; sd(ln η) = 1,26 % pour 10⁶ photons (1,78 %). η s'apprend donc
plus vite et le filtre fait mieux :

| | SCLM régulier | SCLM_mCer3 |
| --- | --- | --- |
| section 15, temps réel / lissé (mM, 20 répétitions) | 1,45 / 0,87 | 2,13 / 1,14 |
| section 15, η final | 1,244 ± 0,016 | 1,244 ± 0,021 |
| section 14, 1 000 photons : 1 ROI / 10 ROI / η connu (mM) | 2,30 / 1,84 / 1,75 | 3,76 / 2,35 / 1,72 |
| section 14, 10 000 photons : 1 ROI / 10 ROI / η connu (mM) | 0,92 / 0,75 / 0,75 | 1,38 / 1,47 / 0,78 |
| η final, 1 000 photons, 10 ROI, c statique, transition du document | 1,228 | 1,176 |
| amincissement, filtre avec fond connu, biais à 10 / 20 / 40 mM | +0,0 / −1,4 / −1,3 | −1,6 / −0,3 / −2,7 |
| amincissement, filtre sans fond, biais à 10 / 20 / 40 mM | +2,4 / +0,7 / +1,9 | +3,3 / +6,7 / +7,9 |
| titration relue : η du FLIM (fond ajusté) | 0,976 à 1,026 | 0,991 à 1,046 |

Les points 1 et 2 ci-dessus valent encore, atténués : η reste tiré vers le bas (1,228 à 1 000 photons) ou vers le haut
avec une transition trop raide (1,361), et la bascule fige η à 1,224 à 3 000 photons (trame 68). À 10 000 photons, 10 ROI
rejoignent la référence η connu (0,75 mM).

**Ce qui est propre à ce jeu.** L'autofluorescence (formes illustratives) coûte bien plus à l'hybride : +10,0 mM
(5 % du jaune) et +9,8 mM (5 % des deux canaux), contre +0,7 et −9,0 pour mCer3. L'échantillon « inconnu » (1,9 + 5,2
millions de photons) se lit 43,3 mM par le ratio (η = 1), 31,7 / 35,9 mM par le FLIM (bleu / deux canaux), 36,7 mM par
l'hybride, avec η = 0,931 ± 0,2 % : un écart de 7 % à η = 1 que la titration relue (η du FLIM à ± 3 %) n'explique pas
entièrement ; l'intervalle à 95 % de l'hybride ([36,8 – 36,8]) est plus étroit que le pas de la grille (2,5 %).

## Comparaison de quatre techniques sur les .sdt réels (`sorties_boxplots/`)

`boxplots.py`, SCLM_mCer3 seulement, à partir des 798 acquisitions de 1 s de `Data/2026-10-16/SCLM_mCer3/` (100 par
concentration, 98 à 2,5 mM ; leur somme redonne exactement `irf_bruts_mCer3.csv`, vérifié à chaque course).

- **Partage** : numéros impairs (c001, c003, …) pour calibrer, pairs pour tester — entrelacés dans le temps, aucun
  fichier des deux côtés (`repartition_fichiers.csv`).
- **Calibration** sur la somme des fichiers impairs : le modèle du document (step0 ; `calibration_moitie.json`,
  déviance/ddl 1,65, τ₁ 0,0865, τ₂ 0,852, τ_D 3,752, τ_A 3,483 ns, K 62,7 mM) pour trois techniques ; pour le temps de
  vie unique, une exponentielle (⊗ IRF bleue, fond compris) ajustée au canal bleu de chaque concentration
  (`calibration_temps_de_vie.csv` : τ de 3,404 ns à 0 mM à 3,653 ns à 131 mM), inversée par interpolation monotone.
- **Techniques**, toutes par maximum de vraisemblance de Poisson avec le même fond connu (2,37 % du bleu, 0,11 % du
  jaune), c cherché entre 0 et 300 mM :
  - *temps de vie unique* : une exponentielle sur le bleu, puis τ(c) ;
  - *FLIM du système* : les formes des deux canaux (morceau « quand » de la section 11) ;
  - *ratio* : N_Y/N_B inversé sur R₁(c), η = 1 ;
  - *hybride FRET-FLIM* : η n'est pas supposé, il est fixé par le **contexte de l'expérience** — un échantillon de
    **100 × N** photons bleus au même état (l'expérience mise en commun : beaucoup de trames, tout le champ), sur lequel
    le FLIM du système donne ĉ et η̂ = R / R₁(ĉ), le maximum de la vraisemblance hybride à η libre ; l'échantillon de
    N photons est ensuite lu par les deux morceaux à η = η̂. Le contexte est tiré par Poisson autour de la forme de
    l'ensemble des fichiers de test de la concentration (100 × 100 000 photons dépasse ce qu'ils contiennent à 0 mM).
    Distribution de η̂ : `eta_contexte.csv`.
- **Échantillons** : 20 000 par concentration et par nombre de photons (800 000 en tout, 3,2 millions de lectures,
  autant de contextes), chacun tiré d'un ou deux fichiers de test pris au hasard et aminci (binomial) à N photons
  **bleus** attendus dans la fenêtre, N = 1 000, 3 000, 10 000, 30 000, 100 000 ; le jaune suit avec la même part
  (≈ 2,5 N photons jaunes). Amincir un comptage de Poisson donne un comptage de Poisson : chaque échantillon a la
  statistique exacte d'une acquisition plus courte, avec les écarts réels d'un fichier à l'autre (horloge, dérive). Les
  échantillons d'un même fichier ne sont pas indépendants entre eux à fort N (50 fichiers de test par concentration).
- **Sorties** : `boxplots_<N>.png` (un panneau par concentration, les 4 techniques ; boîte 25–75 %, trait médian,
  moustaches 2,5–97,5 %, pas de points isolés), `fig_resume.png` (demi-largeur 2,5–97,5 % et |médiane| contre N),
  `ecarts_resume.csv` (moyenne, écart-type, RMSE, quantiles, part des lectures au bord de la grille), `eta_contexte.csv`,
  `journal.txt`.

Écart c mesuré − c réel, médiane ± demi-largeur de l'intervalle 2,5–97,5 % (mM) :

| photons bleus | temps de vie unique | FLIM du système | ratio (η = 1) | hybride (η du contexte) |
| --- | --- | --- | --- | --- |
| **10 mM** — 1 000 | −0,6 ± 150 | −0,7 ± 47 | −0,1 ± 4,0 | −0,9 ± 5,0 |
| 10 000 | −0,4 ± 17 | −0,9 ± 9,7 | −0,1 ± 1,4 | −0,9 ± 1,6 |
| 100 000 | −0,4 ± 6,3 | −0,9 ± 2,9 | −0,1 ± 0,5 | −0,9 ± 0,6 |
| **40 mM** — 1 000 | −2,3 ± 150 | −2,0 ± 150 | −0,6 ± 8,8 | −2,5 ± 12 |
| 10 000 | −2,7 ± 66 | −2,5 ± 29 | −0,6 ± 2,9 | −2,5 ± 3,8 |
| 100 000 | −2,3 ± 17 | −2,5 ± 8,6 | −0,6 ± 1,0 | −2,5 ± 1,3 |
| **131 mM** — 1 000 | −1,0 ± 150 | +15 ± 145 | −1,1 ± 33 | +15 ± 67 |
| 10 000 | +2,2 ± 133 | +15 ± 120 | −1,0 ± 10 | +15 ± 20 |
| 100 000 | +2,7 ± 51 | +15 ± 48 | −0,9 ± 3,5 | +15 ± 6,4 |

Ce qu'on y lit :

1. **L'hybride a la précision du ratio, et le biais du FLIM.** Sa dispersion vaut 1,1 à 1,3 fois celle du ratio (deux
   fois à 131 mM, où R₁(c) est plat) : l'incertitude de η̂ s'y ajoute, et elle suit la loi du document — à 10 mM et
   N = 10 000, le contexte (10⁶ photons bleus, ≈ 3,5 × 10⁶ au total) fixe η à 0,89 % près (1,8 % × √(1/3,5) ≈ 0,96 %).
   Mais η̂ hérite du plancher de la calibration FLIM : sur ces données, η̂ vaut 0,983 à 1,039 selon la concentration
   (1,037 à 20 mM, 1,039 à 131 mM) au lieu de 1, et la lecture garde le biais de la partie FLIM (+3,4 mM à 20 mM, −2,5 à
   40, +4,4 à 80, +15 à 131 mM), constant en N. À η = 1 connu, l'hybride se confondait avec le ratio.
2. **Le ratio**, lu ici à η = 1 (exact dans ces conditions), atteint la borne de Cramér-Rao : à 10 mM et 10 000 photons
   bleus, ± 0,68 mM d'écart-type, soit 1,28 mM × √(10 000/35 000). Biais sous 1 mM partout. En tranche (η ≠ 1), il
   serait biaisé sans apprentissage de η (partie 3 de run.py) ; ce banc ne le montre pas.
3. **Le FLIM du système est 5 à 10 fois plus dispersé que le ratio**, mais 1,5 à 3 fois moins que le temps de vie
   unique : le modèle à trois populations et les deux canaux tirent davantage des mêmes photons qu'une exponentielle
   unique. Son biais (+3,4 mM à 20 mM, +15 à 131 mM) est le plancher de la calibration (partie 5 de run.py).
4. **Le temps de vie unique est le moins précis mais le moins biaisé** sur ces 8 points (|médiane| ≤ 3 mM) : sa courbe
   τ(c) est mesurée aux concentrations mêmes du test, au lieu de passer par les isothermes. Ailleurs qu'en ces 8 points,
   il dépendrait de l'interpolation.
5. **Bords** : une lecture sous 0 mM est ramenée à 0, d'où ~ 50 % des lectures à 0 à 0 mM pour toutes les techniques,
   et des boîtes tronquées à −c aux faibles concentrations. À 1 000 et 3 000 photons, le temps de vie unique et le FLIM
   touchent aussi le haut de la grille (300 mM) : jusqu'à 50 % des lectures (colonne `part_au_bord_pct`). Leurs
   moustaches à ces N mesurent surtout la grille.

## Pas fait

- L'autofluorescence n'est pas modélisée dans le filtre (le document la donne corrigée « si on la modélise ») : la
  colonne « hybride » de la partie 6 en montre l'effet sans modèle, avec des formes d'autofluorescence illustratives.
- Ni la ré-estimation des durées de vie en tranche, ni l'estimation du décalage d'horloge par le filtre, ni le temps mort.
- Les remèdes aux points 1 et 2 (apprendre η sur une ligne de base avec une transition ajustée, bascule sur un critère
  qui tienne compte du biais) restent à choisir.
