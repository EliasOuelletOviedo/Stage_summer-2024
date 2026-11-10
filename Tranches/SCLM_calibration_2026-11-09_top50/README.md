# SCLM — calibration en tranche, 2026-11-09 — somas les plus brillants (50 % par FOV)

Même tranche, mêmes données et même analyse que `../SCLM_calibration_2026-11-09/` (voir son README pour le
contexte, les données et les contrôles), **à une différence près : dans chaque acquisition, seuls les ROIs dont le
total de photons (CFP + YFP, cycles gardés) est dans les 50 % les plus élevés sont gardés** (`step3_select.py
--top-photons 50`). Les ROIs Cellpose, les déclins et les contrôles des étapes 1 à 4 sont repris tels quels du
dossier d'origine (copiés : ils ne dépendent pas de la sélection) ; les étapes 5 à 9 sont relancées.

```bash
cd ../sclm_calibration
python3 step3_select.py ../SCLM_calibration_2026-11-09_top50 --top-photons 50
python3 step5_irf.py ../SCLM_calibration_2026-11-09_top50
python3 step6_somas.py ../SCLM_calibration_2026-11-09_top50
python3 step7_calibration.py ../SCLM_calibration_2026-11-09_top50
python3 step8_cinetiques.py ../SCLM_calibration_2026-11-09_top50
python3 step9_flim_deux_canaux.py ../SCLM_calibration_2026-11-09_top50
python3 step9_flim_deux_canaux.py ../SCLM_calibration_2026-11-09_top50 --transfert 1
```

Pour refaire les étapes 1 à 4 dans ce dossier plutôt que de les copier : `step1_decode.py`, `step2_cellpose.py`,
`step4_declins.py` comme dans le README du pipeline (≈ 15 min).

## Sélection

Somas gardés : basal 355 (sur 706), 0 mM 265 (527), 20 mM 286 (570), 131 mM 220 (439) ; time-lapses : basal→0
80/159, 0→20 57/113, 20→131 56/112. Le seuil (médiane des photons des ROIs de l'acquisition) va de 16 000
photons (basal/fov2) à 1,8 million (basal→0, 65 cycles). Taux YFP enregistré par impulsion, médiane des somas : 0,008
(basal), 0,051 (0 mM), 0,029 (20 mM), 0,007 (131 mM) ; la correction du temps mort change le rapport de +2,1 % à
0 mM (`results/temps_mort.png`).

## Contrôles

Ceux des étapes 1 à 4 sont ceux du dossier d'origine (`qc/inventaire.csv`, `validation_sdt.csv`,
`recalage_bliq.csv`, `derive.png`). **IRF** (`qc/irf.png`, `qc/irf_somas.png`) : décalage +55,6 ± 1,1 ps d'un FOV à
l'autre, aucun élargissement voulu ; trois populations à durées de vie communes 2,08 / 0,51 / 0,04 ns ; déviance/ddl
médiane d'un soma 1,12 (774 somas ≥ 20 000 photons), 3,7 en mono-exponentielle.

## Ce que disent les résultats

### Calibration en tissu (`results/calibration.png`, `results/calibration.csv`)

| Lecture | 0 mM | 20 mM | 131 mM | X∞ | K (mM) | [Cl⁻]i basale (mM) |
|---|---|---|---|---|---|---|
| **rapport YFP/CFP brut** | 5,56 | 3,10 | 1,64 | 1,16 [0,41 ; 1,55] | **16 [3 ; 29]** | **4,0 [0,9 ; 7,5]** |
| rapport corrigé du temps mort | 5,71 | 3,14 | 1,64 | 1,16 | 15 [3 ; 28] | 4,4 [1,0 ; 7,8] |
| τ CFP mono-exp. (ns) | 1,426 | 1,617 | 1,766 | 1,82 | 21 [9 ; 39] | 3,5 [0,1 ; 13,4] |
| τm CFP, temps moyen (ns) | 1,440 | 1,583 | 1,706 | 1,75 | 24 [12 ; 45] | 7,0 [3,3 ; 19,5] |
| τ CFP, 3 populations (ns) | 1,702 | 1,817 | 1,895 | 1,92 | 18 [7 ; 35] | 3,3 [0,5 ; 12,9] |

Comparé à tous les somas : rapport un peu plus haut à 0 mM (5,56 contre 5,38) et dynamique un peu plus grande
(X₀/X∞ 4,8 contre 4,2) ; K et la [Cl⁻]i basale changent peu (16 contre 14 mM ; 4,0 contre 4,4 mM).

### [Cl⁻]i basale : 3 à 7 mM

- Médiane des somas du basal : **4,0 mM [0,9 ; 7,5]** par le rapport, 3,3 à 7,0 mM par les durées de vie.
- Cellule par cellule (rapport) : médiane 4,0 mM, IQR 2,2–6,4 mM ; 6 % des somas sous 0 mM (`results/basal.png`).
- **Robustesse** : en retirant un FOV à la fois, la [Cl⁻]i basale va de 2,7 à 5,1 mM (K de 10 à 21 mM).

### Précision d'une cellule seule (auto-cohérence, `results/autocoherence.csv`)

| Imposé | Relu par le rapport : médiane [IQR] | Relu par τ mono |
|---|---|---|
| 0 mM | 0,0 [−1,7 ; 2,7] mM | −0,1 [−2,4 ; 3,2] |
| 20 mM | 19,9 [12,5 ; 53,5] mM | 18,8 [12,3 ; 32,8] |
| 131 mM | 124 [96 ; 176] mM (13 % hors courbe) | 54 [25 ; 143] (30 % hors courbe) |

Même conclusion qu'avec tous les somas : une cellule seule n'est lisible que sous ~10 mM.

### Cinétiques des changements de solution (`results/transition_*.png`)

- **0→20** : 3,67 au premier cycle, 3,65 au dernier (τ ≈ 1,5 ± 0,8 min) ; **20→131** : 1,42 → 1,25, traîne de
  τ = 5,3 ± 0,4 min (médiane par soma 4,3 min). Les FOV de la condition d'arrivée ont été pris à ≥ 90 % de ces
  traînes. Le gros du changement précède toujours le premier cycle.
- **basal→0** : pas de plateau en 49 min (rapport 4,98 → 6,10, τ CFP 1,44 → 1,22 ns) ; deux exponentielles de 6,1
  et 16,7 min pour le rapport, avec la même descente initiale qu'avec tous les somas.

### FLIM à deux canaux, soma par soma (étape 9)

Calibration du modèle de `Simulateurs/bayesian_filter` sur les somas réunis de 0, 20 et 131 mM
(`results/deux_canaux_calibration.png`) : déviance/ddl **35** (48 avec tous les somas) ; τ1 8,5 ps, τ2 0,61, τD 2,31,
τA 2,96 ns ; λ 0,67, G 0,19, g2 0, κ 0,04, k0 0 ; f1 de 0,91 à 0,67, f2 de 0,05 à 0,15 entre 0 et 131 mM.
Lectures, médiane (±1,48·MAD entre somas ; écart-type des médianes de FOV ; % hors courbe), en mM — les [Cl⁻]
négatives sont gardées, les lectures hors courbe comptent comme +∞ :

| Lecture | basal | 0 mM | 20 mM | 131 mM |
|---|---|---|---|---|
| τ mono-exp. → isotherme de l'étape 7 | 3,5 (±7,5 ; 6,6 ; 5 %) | 0,0 (±3,7 ; 0,9) | 20,0 (±12,9 ; 15,4 ; 4 %) | 131 (±59 ; 87 ; 30 %) |
| temps de vie unique | 0,0 (±11,9 ; 6,3) | −3,4 (±7,0 ; 2,3) | 19,6 (±19,7 ; 24,1) | 115 (±107 ; 62) |
| ratio, η = 1 | 2,9 (±3,3 ; 1,3) | −1,5 (±3,1 ; 1,4) | 19,9 (±15,7 ; 27,6) | 118 (±41 ; — ; 8 %) |
| FLIM du système | 1,8 (±4,3 ; 1,9) | 0,4 (±3,2 ; 0,6) | 17,7 (±13,8 ; 18,0) | 201 (±97 ; 130 ; 25 %) |
| **hybride, η partagé par le FOV** | **1,9 (±3,3 ; 2,1)** | 0,6 (±3,0 ; 1,1) | 18,2 (±13,5 ; 21,2) | 171 (±101 ; 191 ; 3 %) |

(« — » : à 131 mM, le FOV 4, très pâle, a un rapport médian de 1,08, sous X∞ : sa médiane relue est démesurée.)

- **Ce qui change avec les somas brillants** : la lecture hybride suit toujours le rapport (ρ de Spearman −0,95 à
  −1,00), mais aussi, maintenant, le τ mono (ρ 0,58 à 0,77 sous 20 mM, contre 0,28–0,36 avec tous les somas). Une
  bonne part du désaccord cellule à cellule entre rapport et durée de vie venait donc du **bruit de photons des
  somas pâles** sur le τ ; chez les somas brillants, les deux lectures du FRET voient en grande partie la même
  variation d'une cellule à l'autre.
- Les dispersions entre somas baissent un peu sous 20 mM (hybride ±3,0 à 0 mM, ±13,5 à 20 mM contre ±3,6 et
  ±14,9) ; l'intervalle à 95 % d'un soma (hybride) fait 0,5 à 1,5 mM sous 20 mM. La dispersion entre somas reste
  bien plus grande que ce que les photons expliquent.
- Basal : 1,8 à 3,5 mM selon la lecture (hybride 1,9 mM). À 0 mM, l'hybride lit 0,6 mM (η des FOV 1,04–1,14).
- **Un seul temps de transfert** (`--transfert 1`, sorties `*_1t`) : déviance/ddl 130 ; τT 0,29, τD 1,93, τA 2,94 ns,
  G = 0, κ 2,6, f de 0,60 à 0,23 (K 36 mM). Mêmes constats qu'avec tous les somas : plus resserré jusqu'à 20 mM
  (FLIM du système ±10,2 ; FOV 12,2 à 20 mM) mais **biaisé en haut** (hybride 65 mM, FLIM du système 73 mM à
  131 mM), parce que la forme des déclins à 131 mM n'y est pas reproduite. Basal 1,4 à 2,8 mM.

### Réserves

Celles du dossier d'origine s'appliquent toutes (champs qui diffèrent autant que les conditions — à 20 mM, les FOV
3 et 5, pâles, ont un rapport médian de 2,0–2,1 contre 3,4–3,7 pour les autres, ρ rapport–intensité 0,69 ; K
déterminé par le seul point à 20 mM ; équilibre à 0 mM incertain ; pH ; IRF paramétrique et bords de la fenêtre).
Propre à cette sélection : garder les somas les plus brillants **de chaque acquisition** ne compare pas des
cellules de même expression d'une condition à l'autre (la puissance et l'extinction de la YFP changent l'échelle),
et écarte en partie les cellules pâles, qui sont peut-être aussi les moins saines.

## Sorties

Mêmes fichiers que le dossier d'origine : `qc/` (étapes 1 à 5), `masks/` (copiés), `selection/` (top 50 %),
`results/` (`somas.csv`, `conditions.csv`, `calibration.csv`, `basal_cellules.csv`, `autocoherence.csv`,
`transitions_*`, `deux_canaux_*` et `deux_canaux_*_1t`, figures).
