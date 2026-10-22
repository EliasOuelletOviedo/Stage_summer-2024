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

# comparaison de 4 techniques sur les .sdt réels de SCLM_mCer3 → sorties_boxplots/ (~4,5 min, --rapide ~1 min)
python3 boxplots.py
# la même, avec η de chaque échantillon uniforme de 0,9 à 1,1 → sorties_boxplots_eta/ (~4,5 min)
python3 boxplots.py --eta

# l'exemple de la section 15 : 5 nombres de photons × 3 scénarios de η, 3 000 trames → sorties_exemples/ (~3 min)
python3 exemples.py
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
| `params.toml` | le seul endroit pour les constantes (données, grilles, transition, scénarios, graine, boxplots) ; SCLM_mCer3 | 133 |
| `params_SCLM.toml` | les mêmes, pour la titration SCLM régulière à 800 nm | 95 |
| `step0_calibration.py` | le pont unique avec les données : appelle le `sclm_hybride.py` du jeu, écrit sa calibration (`calibrer` sert aussi à boxplots.py) | 86 |
| `calibration.json`, `calibration_SCLM.json` | la calibration ajustée (constantes, isothermes, IRF, fonds et décalages par acquisition) ; écrites par step0 | — |
| `src/modele.py` | modèle direct μ_B, μ_Y (sections 3 à 7) : IRF ⊗ déclins périodiques, isothermes, η, A ; perturbations | 137 |
| `src/filtre.py` | vraisemblance découpée ratio / FLIM (section 11), grille (c, η), transition, filtre et lisseur (13), temps réel à M ROI et portes (14 ; η jamais figé), lectures de référence | 179 |
| `src/information.py` | Fisher et Cramér-Rao, nuisances par complément de Schur (section 10) | 50 |
| `src/simulation.py` | les deux générateurs (lit les colonnes voie1/voie2 de step17 et bleu/jaune de step16), les trajectoires des sections 14 et 15 | 97 |
| `src/figures.py` | rendu | 174 |
| `src/estimateurs.py` | les 4 lectures d'un échantillon (temps de vie unique, FLIM du système, ratio, hybride à η du contexte), en bloc | 98 |
| `boxplots.py` | partage des .sdt, calibration sur une moitié, échantillons et contextes tirés de l'autre (η = 1 ou variable), boxplots | 173 |
| `run.py` | la course en 8 parties ; un CSV par partie, `journal.txt`, 5 figures | 357 |
| `src/exemple.py` | une course de la section 15 (trajectoire étirée, η constant ou par trame), partagée par run.py et exemples.py | 66 |
| `exemples.py` | l'exemple de la section 15 sous 15 conditions (photons × η) | 80 |
| `test/run_tests.py` | 13 vérifications : calibration = annexe, découpage exact, portes (93,7 % de l'information FLIM gardée), η retrouvé, η jamais figé, couverture du lisseur, amincissement | 167 |

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

**Section 14 (300 trames, η = 1,25 inconnu, 3 répétitions).** Erreur typique après la trame 100, en mM. η est mis à
jour à chaque trame jusqu'à la fin : la bascule du document (η figé et filtre sur le ratio seul une fois
sd(ln η) < 1 %) est retirée (voir le point 2 plus bas) ; la référence « η connu » fait la même mise à jour, formes
comprises, avec la grille de η réduite à 1,25.

| photons par trame | 1 ROI | 10 ROI partageant η | η connu |
| --- | --- | --- | --- |
| 1 000 | 3,76 (2,6) | 2,35 (2,5) | 1,70 (1,7) |
| 3 000 | 2,26 (1,5) | 1,51 (1,4) | 1,16 (1,1) |
| 10 000 | 1,38 (1,5) | 0,82 (0,9) | 0,77 (0,8) |

η final pour 10 ROI : 1,215 ± 1,1 %, 1,226 ± 0,5 %, 1,247 ± 0,4 % à 1 000, 3 000, 10 000 photons. Coût : 0,6 à 1 ms par
ROI et par trame avec les portes (la bascule le ramenait à 0,2 ms, au prix d'un η figé).

## Ce que la simulation apprend en plus du document

1. **η n'est pas appris que du FLIM.** Le score ℓ_m(η) de la section 14 absorbe aussi l'écart entre la dynamique supposée
   de c et la vraie : un η statique peut « expliquer » une partie d'un changement de ratio. Partie 8, 1 000 photons,
   10 ROI, η vrai 1,25 : c statique et transition du document → η = 1,18 (−6 %, ≈ 5 sd a posteriori) ; c qui répond et
   transition trop raide (s = 0,01, sans sauts) → η = 1,45. Dans un essai hors course, avec c figé et K = identité,
   η revenait à 1,24–1,26.
   Le module de référence (`sclm_hybride.py`, `demo_temps_reel`) fait la même chose : 1,21 à 1 000 photons pour trois
   graines, et 1,231 à 1,250 à 3 000 photons (le 1,250 ± 1,0 % du document est la graine 0).
2. **La bascule à 1 % figeait η trop tôt ; elle est retirée.** À 10 000 photons, sd(ln η) passait sous 1 % vers la
   trame 40, pendant la ligne de base, alors que l'effet ci-dessus tirait encore η vers le bas : η restait figé à
   1,21–1,23, et 10 ROI faisaient moins bien qu'une seule (1,47 contre 1,38 mM). Sans bascule, η continue d'apprendre
   pendant les réponses et finit à 1,247 ; 10 ROI rejoignent la référence η connu (0,82 contre 0,77 mM). L'effet du
   point 1 reste visible à 3 000 photons : η = 1,226 ± 0,5 % (≈ 4 sd sous 1,25).
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
| section 14, 1 000 photons : 1 ROI / 10 ROI / η connu (mM) | 2,30 / 1,81 / 1,71 | 3,76 / 2,35 / 1,70 |
| section 14, 10 000 photons : 1 ROI / 10 ROI / η connu (mM) | 0,89 / 0,73 / 0,73 | 1,38 / 0,82 / 0,77 |
| η final, 1 000 photons, 10 ROI, c statique, transition du document | 1,228 | 1,176 |
| amincissement, filtre avec fond connu, biais à 10 / 20 / 40 mM | +0,0 / −1,4 / −1,3 | −1,6 / −0,3 / −2,7 |
| amincissement, filtre sans fond, biais à 10 / 20 / 40 mM | +2,4 / +0,7 / +1,9 | +3,3 / +6,7 / +7,9 |
| titration relue : η du FLIM (fond ajusté) | 0,976 à 1,026 | 0,991 à 1,046 |

Le point 1 ci-dessus vaut encore, atténué : η reste tiré vers le bas (1,228 à 1 000 photons) ou vers le haut avec une
transition trop raide (1,361). Sans bascule, 10 ROI font 1,81 / 1,20 / 0,73 mM à 1 000 / 3 000 / 10 000 photons
(η connu : 1,71 / 1,12 / 0,73), contre 1,84 / 1,45 / 0,75 quand la bascule figeait η à 1,224 à 3 000 photons.

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
- **Sorties** : deux séries de graphiques par nombre de photons — `boxplots_<N>.png` (écart signé) et
  `boxplots_abs_<N>.png` (|écart|, échelle log commune aux 8 panneaux ; un |écart| sous 10⁻³ mM, c'est-à-dire une
  lecture exacte au bord c = 0, est ramené à 10⁻³) ; un panneau par concentration, les 4 techniques ; boîte 25–75 %,
  trait médian, moustaches 2,5–97,5 %, pas de points isolés — `fig_resume.png` (demi-largeur 2,5–97,5 % et |médiane| contre N),
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

### Avec η variable (`sorties_boxplots_eta/`, `python3 boxplots.py --eta`)

Mêmes fichiers, même calibration, mêmes échantillons, mais chaque échantillon reçoit son propre η, tiré uniformément de
0,9 à 1,1 : ses photons jaunes (amincis avec la part × η) et ceux de son contexte (la même expérience) en sont
multipliés. Aucune lecture ne connaît η : le ratio le prend égal à 1, l'hybride l'apprend de son contexte de 100 × N
photons ; les deux lectures par la forme (temps de vie unique, FLIM du système) n'en dépendent pas et ne changent pas.

Demi-largeur de l'intervalle 2,5–97,5 % de l'écart (mM), ratio / hybride :

| photons bleus | 0 | 2,5 | 5 | 10 | 20 | 40 | 80 | 131 mM |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 000 | 2,8 / 1,6 | 4,5 / 3,1 | 5,9 / 4,3 | 7,2 / 4,8 | 9,9 / 7,8 | 15 / 12 | 31 / 28 | 58 / **69** |
| 3 000 | 2,4 / 0,9 | 4,1 / 2,2 | 5,2 / 2,5 | 6,2 / 2,8 | 8,4 / 4,5 | 13 / 6,8 | 26 / 16 | 47 / 37 |
| 10 000 | 2,2 / 0,5 | 3,9 / 1,2 | 4,8 / 1,4 | 5,7 / 1,6 | 7,7 / 2,5 | 12 / 3,8 | 23 / 8,7 | 42 / 21 |
| 30 000 | 2,2 / 0,3 | 3,8 / 0,7 | 4,6 / 0,8 | 5,5 / 1,0 | 7,5 / 1,5 | 11 / 2,3 | 22 / 5,3 | 40 / 13 |
| 100 000 | 2,1 / 0,2 | 3,7 / 0,4 | 4,6 / 0,5 | 5,4 / 0,6 | 7,3 / 0,9 | 11 / 1,3 | 22 / 2,9 | 40 / 8,2 |

Médiane de |écart| (mM), ratio / hybride :

| photons bleus | 2,5 | 5 | 10 | 20 | 40 | 80 | 131 mM |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 000 | 2,3 / 1,3 | 2,4 / 1,5 | 2,9 / 1,8 | 3,9 / 3,6 | 6,0 / 4,5 | 12 / 9,3 | 21 / **22** |
| 10 000 | 2,2 / 0,4 | 2,4 / 0,8 | 2,8 / 1,0 | 3,7 / 3,4 | 5,8 / 2,6 | 11 / 4,7 | 20 / 15 |
| 100 000 | 2,2 / 0,2 | 2,4 / 0,7 | 2,8 / 0,9 | 3,6 / 3,4 | 5,8 / 2,5 | 11 / 4,4 | 20 / 15 |

1. **L'hybride est plus précis que le ratio seul**, à toutes les concentrations dès 3 000 photons bleus, et partout à
   1 000 sauf à 131 mM. L'erreur du ratio ne baisse plus avec les photons : elle vient de η (± 10 % de η ≈ ± 5,5 mM à
   10 mM) et plafonne vers 2 à 40 mM de demi-largeur selon c ; celle de l'hybride continue de suivre 1/√N.
2. **Ce qui reste à l'hybride est le plancher de la calibration FLIM**, le même qu'à η = 1 : son η̂ est décalé de
   −1,5 à +3,9 % selon la concentration (`eta_contexte.csv`, colonne `erreur_eta_pct_mediane` : +3,7 % à 20 mM, +3,9 %
   à 131 mM), d'où des biais de +3,4 mM à 20 mM et +15 mM à 131 mM. À 20 mM, la médiane de |écart| de l'hybride
   (3,4 mM) n'est donc qu'à peine sous celle du ratio (3,6 à 3,9 mM), et à 131 mM et 1 000 photons, l'hybride fait un peu
   moins bien que le ratio (69 contre 58 mM de demi-largeur, 22 contre 21 mM de médiane).
3. Les deux lectures par la forme donnent, aux fluctuations de tirage près, les mêmes boîtes qu'à η = 1 (à 10 000
   photons et 10 mM : 9,8 contre 9,7 mM de demi-largeur pour le FLIM, 17,4 contre 17,1 pour le temps de vie unique) :
   elles ne voient pas η, mais restent 5 à 10 fois plus dispersées que l'hybride.

## L'exemple de la section 15 sous d'autres conditions (`sorties_exemples/`)

`python3 exemples.py` (~3 min). La course de la section 15 (un soma ; 8 mM, montée vers 22 mM, retour vers 10 mM ;
A −30 %) sur **3 000 trames** au lieu de 300, même forme de courbe (montée à la trame 601, τ 40 trames ; descente à
1 401, τ 250), pour 1 000, 3 000, 10 000, 30 000 et 100 000 photons par trame et trois η vrais, inconnus du filtre :
1, 1,1, et une dérive linéaire de 1 à 1,1 de la première à la dernière trame. Une figure par cas
(`fig_exemple_<scénario>_<photons>.png`, les deux panneaux de `fig_exemple.png`), la table par trame
(`exemple_<scénario>_<photons>.csv.gz`), `resume.csv`, `fig_resume.png`, `journal.txt`.

Erreur typique sur toute la course (mM), ratio lu à η = 1 / filtre temps réel / lissé :

| photons par trame | η = 1 | η = 1,1 | η dérive de 1 à 1,1 |
| --- | --- | --- | --- |
| 1 000 | 4,8 / 3,6 / 3,1 | 6,7 / 2,4 / 2,4 | 5,3 / 5,0 / 4,0 |
| 3 000 | 2,9 / 1,8 / 2,0 | 6,0 / 1,2 / 0,8 | 4,0 / 2,5 / 2,1 |
| 10 000 | 1,5 / 0,7 / 0,5 | 5,6 / 1,0 / 1,0 | 3,5 / 2,0 / 1,7 |
| 30 000 | 0,9 / 0,5 / 0,4 | 5,6 / 0,6 / 0,4 | 3,2 / 2,0 / 1,5 |
| 100 000 | 0,5 / 0,4 / 0,3 | 5,6 / 0,4 / 0,3 | 3,2 / 1,3 / 0,5 |

1. **η constant et inconnu (1,1)** : le filtre retrouve η (1,084 à 1,098 dès 3 000 photons) et lit c à 0,3–1 mM près, là
   où le ratio calibré garde ~ 5,6 mM d'erreur quel que soit le nombre de photons. Couverture des intervalles lissés à
   95 % : 95 à 98 %.
2. **Peu de photons et longue ligne de base** : à 1 000 et 3 000 photons, η finit trop bas (0,951 et 0,970 pour η = 1 ;
   1,062 pour η = 1,1) — l'effet du point 1 de « Ce que la simulation apprend », plus fort ici avec 600 trames de ligne
   de base statique — et les intervalles lissés ne contiennent la vraie valeur que 28 à 65 % du temps : ils sont trop
   étroits, parce qu'ils supposent η connu à la précision de sa distribution a posteriori, qui ne voit pas ce biais.
3. **η qui dérive : ce filtre ne le suit pas.** Il suppose η constant (section 14 du document) : la marche aléatoire de
   `s_eta_lisseur` (0,2 % par trame) est plus fine que le pas de la grille de η (1,83 %), le poids vers le point voisin
   vaut ~ 10⁻¹⁹, et η reste en pratique constant. Le filtre choisit un η de compromis (0,98 à 1,04 dans le lisseur) :
   trop haut au début, trop bas à la fin ; à 10 000 et 30 000 photons, l'erreur du filtre temps réel plafonne à 2,0 mM
   (2,8 mM sur le dernier dixième), et la couverture lissée tombe à 28–77 %. Sa sd a posteriori « ± 0,0 % » veut dire
   que toute la probabilité est sur un seul point de la grille, pas que η est connu à mieux que le pas. L'hybride
   reste plus juste que le ratio calibré (4,7 à 5,7 mM sur le dernier dixième), mais pour suivre une dérive il faudrait
   une vraie marche aléatoire sur η (noyau qui garde la variance sur la grille, ou grille plus fine), au prix d'un η
   moins précis quand il ne bouge pas.

## Pas fait

- L'autofluorescence n'est pas modélisée dans le filtre (le document la donne corrigée « si on la modélise ») : la
  colonne « hybride » de la partie 6 en montre l'effet sans modèle, avec des formes d'autofluorescence illustratives.
- Ni la ré-estimation des durées de vie en tranche, ni l'estimation du décalage d'horloge par le filtre, ni le temps mort.
- Le remède au point 1 (apprendre η sur une ligne de base avec une transition ajustée) reste à choisir. Sans la
  bascule, le filtre temps réel coûte 2 à 5 fois plus de calcul à 10 ROI ; un filtre léger ne reviendrait qu'avec un
  critère qui tienne compte du biais de η, pas seulement de sa sd a posteriori.
