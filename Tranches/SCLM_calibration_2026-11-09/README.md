# SCLM — calibration en tranche, 2026-11-09 ([Cl⁻] imposée par ionophores)

Une tranche (`Slice 1`), neurones exprimant SCLM seul. Images en ACSF (**basal**), puis
[Cl⁻] imposée par ionophores à **0, 20 et 131 mM** ; à chaque changement de solution, un
time-lapse. But : l'isotherme du capteur **dans le tissu** (rapport YFP/CFP et durée de vie
CFP selon la [Cl⁻]) et la [Cl⁻]i basale des neurones.

Analysé avec le pipeline `../sclm_calibration/` (voir son README pour la méthode) :

```bash
cd ../sclm_calibration
python3 step1_decode.py ../SCLM_calibration_2026-11-09
~/.venvs/cellpose/bin/python step2_cellpose.py ../SCLM_calibration_2026-11-09
python3 step3_select.py ../SCLM_calibration_2026-11-09 --tout
python3 step4_declins.py ../SCLM_calibration_2026-11-09
python3 step5_irf.py ../SCLM_calibration_2026-11-09
python3 step6_somas.py ../SCLM_calibration_2026-11-09
python3 step7_calibration.py ../SCLM_calibration_2026-11-09
python3 step8_cinetiques.py ../SCLM_calibration_2026-11-09
python3 step9_flim_deux_canaux.py ../SCLM_calibration_2026-11-09
python3 step9_flim_deux_canaux.py ../SCLM_calibration_2026-11-09 --transfert 1   # un seul temps de transfert
```

**Les résultats portent sur tous les ROIs de Cellpose (`--tout`), sans tri visuel.** Après un
tri (`step3_select.py` sans `--tout`), relancer les étapes 5 à 8.

## Ce que contient `data/`

Deux liens (non versionnés) : `data/FLIM` → `~/Documents/Maîtrise/Données/labo/2026-11-09_FLIM`
(SPCM, SPC-QC-104, 9 Go) et `data/Bliq` → `…/labo/2026-11-09` (images Bliq simultanées).

| Acquisition | Heure | Contenu |
|---|---|---|
| basal, fov1–5 | 14:41–14:58 | 6 cycles de 10 s par FOV (fov1 : 7) |
| basal→0 | 15:23:53 (+49 min) | 71 cycles, un toutes les ~41 s (dossiers `_1`, `_2`, `_3` bout à bout) |
| 0 mM, fov1–5 | 16:15–16:26 | 6 cycles par FOV |
| 0→20 | 16:35:51 (+7 min) | 12 cycles |
| 20 mM, fov1–5 | 16:44–16:53 | 6 cycles par FOV |
| 20→131 | 17:00:00 (+10 min) | 16 cycles |
| 131 mM, fov1–5 | 17:12–17:22 | 6 cycles par FOV |
| IRF | 2026-11-10 10:58 | `irf_850nm_16X.sdt`, même carte, mêmes décalages du TDC |

Laser à 80,25 MHz ; images 512 × 1024 à 0,509 µm/px (16X), 128 × 128 dans les `.sdt` d'SPCM.

- **Canaux** : IN1 = YFP (= Bliq C1), IN2 = CFP (= Bliq C2), établi par les données :
  IN1/IN2 vaut 5,1 à 0 mM et 1,7 à 131 mM, et le déclin de IN2 s'allonge avec le Cl⁻.
  **C'est l'inverse des sessions FLIMApp du 2026-10-08/09.** Bliq C3 ne montre que les
  rayures du polygone : ignoré.
- **Les 20 FOV fixes sont 20 champs différents** (aucune paire ne se ressemble d'une condition à
  l'autre, NCC < 0,3), et les trois time-lapses sont trois autres champs. La calibration
  compare donc des **populations de cellules**, pas les mêmes cellules d'une condition à l'autre.
- **Puissance** : en basal, les deux canaux ont ~4 fois moins de photons qu'ensuite (YFP
  0,2–0,8 MHz par champ, contre 1,0–1,4 MHz à 0 mM) : réglage changé entre les deux.

## Contrôles

- **Géométrie** : pour les 23 acquisitions, les photons de IN1 ressemblent le plus à l'image
  Bliq C1 simultanée au décalage nominal (14 lignes, 28 pixels) (NCC 0,71–0,97, `qc/recalage_bliq.csv`).
- **Microtemps** : le déclin du dernier cycle de chaque FOV colle à celui du `.sdt` écrit par
  SPCM, à moins de 0,8 % du pic pour 18 FOV sur 20 ; 2 à 4 % pour 131 mM fov1 et fov2, sur le
  front de montée (l'image `.sdt` de fov1 n'a pas les lignes vides des autres : SPCM a rangé
  ces photons autrement) (`qc/validation_sdt.csv`).
- **Dérive** : ≤ 0,08 px dans les FOV, 0,19 px pour 0→20, 3,4 px sur les 49 min de basal→0
  (masques déplacés cycle par cycle) (`qc/derive.png`).
- **Plan focal** : cycles écartés car pris à un autre plan (retouches de la mise au point) :
  basal/fov2 c2–c3 (et c4, interrompu après 1,8 s), et 6 cycles de basal→0 (vers 10–14 min et
  28 min). basal/fov2 garde 3 cycles (`qc/inventaire.csv`).
- **IRF** (`qc/irf.png`, `qc/irf_somas.png`) : décalage ajusté de +55 ± 1 ps d'un FOV à l'autre
  (même repère de temps : attendu ≈ 0), aucun élargissement voulu par les données. Avec trois
  populations à durées de vie communes, la déviance/ddl médiane d'un soma est **1,10** (992
  somas) ; une mono-exponentielle donne 3,5 : elle ne décrit pas un soma, et τ mono n'est
  qu'une lecture.
- **Temps mort** : au plus 0,03 photon YFP par impulsion (médiane à 0 mM) ; la correction change
  le rapport de +1,3 % et la [Cl⁻]i basale de +0,2 mM (`results/temps_mort.png`).

## Ce que disent les résultats

### Calibration en tissu

Somas gardés : basal 706, 0 mM 527, 20 mM 570, 131 mM 439. Médianes par condition et
isotherme X(c) = X∞ + (X₀ − X∞)·K/(K + c) (exacte sur 3 points ; IC 95 % par bootstrap
hiérarchique, FOV puis somas) (`results/calibration.png`, `results/calibration.csv`) :

| Lecture | 0 mM | 20 mM | 131 mM | X∞ | K (mM) | [Cl⁻]i basale (mM) |
|---|---|---|---|---|---|---|
| **rapport YFP/CFP brut** | 5,38 | 2,99 | 1,67 | 1,27 [0,57 ; 1,65] | **14 [3 ; 28]** | **4,4 [1,0 ; 7,6]** |
| rapport corrigé du temps mort | 5,45 | 3,01 | 1,68 | 1,26 | 14 [3 ; 27] | 4,6 [1,1 ; 7,9] |
| τ CFP mono-exp. (ns) | 1,414 | 1,597 | 1,724 | 1,77 | 19 [8 ; 34] | 2,5 [−0,8 ; 12,5] |
| τm CFP, temps moyen (ns) | 1,432 | 1,572 | 1,689 | 1,73 | 23 [12 ; 41] | 6,8 [2,4 ; 20] |
| τ CFP, 3 populations (ns) | 1,695 | 1,803 | 1,875 | 1,90 | 18 [7 ; 29] | 2,8 [−0,6 ; 12,6] |

Le rapport a une dynamique de ×4,2 (X₀/X∞) ; τ CFP s'allonge de 22 % entre 0 et 131 mM. Les deux
lectures du FRET vont dans le même sens, soma par soma (`results/rapport_vs_tau.png`).

### [Cl⁻]i basale : 3 à 7 mM

- Médiane des somas du basal : **4,4 mM [1,0 ; 7,6]** par le rapport ; les trois lectures de la
  durée de vie donnent 2,5 à 6,8 mM, avec des IC plus larges. Toutes s'accordent sur une
  [Cl⁻]i basse, de l'ordre de quelques mM.
- Cellule par cellule (rapport) : médiane 4,4 mM, IQR 2,2–7,5 mM ; 6 % des somas se lisent sous
  0 mM (`results/basal.png`, `results/basal_cellules.csv`).
- **Robustesse** : en retirant un FOV à la fois, la [Cl⁻]i basale reste entre 3,2 et 5,6 mM ;
  sans les deux FOV pâles de 20 mM, elle vaut 6,7 mM (rapport) et 3,1 mM (τ).

### Précision d'une cellule seule (auto-cohérence)

Les somas des conditions imposées, relus un par un par la courbe (`results/autocoherence.csv`) :

| Imposé | Relu par le rapport : médiane [IQR] | Relu par τ mono |
|---|---|---|
| 0 mM | 0,0 [−1,6 ; 3,1] mM | −0,1 [−2,8 ; 2,8] |
| 20 mM | 19,9 [12,2 ; 48,8] mM | 18,9 [8,9 ; 36,2] |
| 131 mM | 114 [83 ; 165] mM (16 % hors courbe) | 36 [10 ; 94] (42 % hors courbe) |

Une lecture sur une seule cellule n'a de sens que sous ~10 mM, là où la courbe est raide :
c'est justement le domaine du basal. Au-delà de 20 mM, une cellule seule ne dit presque rien.

### Cinétiques des changements de solution (`results/transition_*.png`)

- **0→20** : le rapport est déjà à 3,63 dès les 10 premières secondes, plat seconde par seconde,
  puis ne bouge presque plus (3,59 ; τ ≈ 1,7 min). **20→131** : 1,44 au premier cycle (les FOV à
  20 mM sont à ~3,0), puis une traîne de τ = 5,1 ± 0,4 min jusqu'à 1,26 (médiane par soma 3,8
  min). Les FOV de la condition d'arrivée ont été pris à ≥ 91 % de ces traînes.
- **Le gros du changement précède donc le premier cycle** des deux time-lapses, alors que t = 0
  est censé être le changement de solution : soit l'échange prend moins de 10 s (ionophores
  déjà dans les membranes), soit la solution a été changée avant le lancement. À vérifier.
- **basal→0** : pas de plateau en 49 min. Pendant les ~8 premières minutes, le rapport monte puis
  redescend sous sa valeur de départ, et τ CFP monte (1,42 → 1,50 ns) : c'est l'inverse d'une
  sortie de Cl⁻ (entrée transitoire de Cl⁻, ou acidification à l'arrivée des ionophores). Ensuite
  le rapport monte et τ baisse régulièrement jusqu'à la fin (rapport 6,1, τ 1,22 ns), **au-delà**
  des FOV à 0 mM pris plus tard (5,4 et 1,41 ns).

### FLIM à deux canaux, soma par soma (étape 9)

Le modèle de `Simulateurs/bayesian_filter` (populations du donneur, fuite, émission sensibilisée, excitation
directe, η), calibré sur les déclins réunis des somas de 0, 20 et 131 mM, puis lu soma par soma par maximum de
vraisemblance, sans filtre temporel (fond libre par canal, η partagé par FOV ; choix du 2026-11-10)
(`results/deux_canaux_calibration.png`, `deux_canaux_lectures.png`, `deux_canaux_tau_mono.png`).

**Calibration.** Déviance/ddl **48** (le simulateur obtenait ~2 en solution) : les résidus se concentrent sur le
front de montée (±25 à 50), le reste du déclin est à ±5–10. τ1 8,5 ps, τ2 0,58, τD 2,29, τA 2,98 ns ; λ 0,71,
G 0,17 (dix fois moins qu'en solution), g2 0, κ 0,04, k0 0 : le jaune est surtout la fuite du donneur. La
population 1 (τ1 ≈ 0, en FRET, presque sans photons bleus) passe de 0,91 à 0,67 entre 0 et 131 mM, la population
2 de 0,05 à 0,16. Comme en solution pour SCLM régulier, c'est une calibration **empirique** de la forme et du ratio
en fonction de [Cl⁻], pas une lecture physique des populations. R₁(c) va de 5,2 à 1,4, demi-chemin vers 15–20 mM.

**Lectures** (médiane (IQR des somas), mM ; `results/deux_canaux_resume.csv`). **Les [Cl⁻] négatives sont gardées
partout** : la grille de c descend jusqu'à la limite du modèle (−15,7 mM, où f1 + f2 = 1 ; aucun soma ne l'atteint),
et aucune lecture n'est bornée à 0 — les tronquer biaiserait vers le haut les lectures basses. En haut, la grille va
à 10⁵ mM : une lecture dont le maximum reste au bord n'a pas de maximum fini (forme au-delà de l'asymptote du
modèle) et est « hors courbe », comptée comme +∞ dans les médianes et les quantiles (l'exclure les biaiserait vers
le bas ; « — » : quartile hors courbe). Ainsi, les médianes du τ mono retombent exactement sur la calibration de
l'étape 7. Une grille arrêtée à 300 mM plafonnait l'hybride à 131 mM (135 au lieu de 166 mM).

| Lecture | basal | 0 mM | 20 mM | 131 mM |
|---|---|---|---|---|
| τ mono-exp. (étape 6 → isotherme de l'étape 7) | 2,5 (15,2) ; 5 % hors courbe | 0,0 (6,0) ; 1 % hors courbe | 20,0 (30,9) ; 4 % hors courbe | 131,1 (—) ; 37 % hors courbe |
| temps de vie unique (IRF paramétrique) | −0,9 (18,3) | −4,2 (11,3) | 17,3 (39,3) | 100,5 (169,6) |
| ratio, η = 1 | 3,9 (5,8) | −0,9 (5,7) | 22,0 (39,4) | 112,9 (95,4) ; 9 % hors courbe |
| FLIM du système (η libre par soma) | 1,8 (7,3) ; 2 % hors courbe | 0,2 (4,8) | 18,3 (32,5) ; 5 % hors courbe | 236,3 (—) ; 33 % hors courbe |
| **hybride, η partagé par le FOV** | **2,9 (5,7)** | 1,4 (5,7) | 20,0 (31,2) | 166,1 (220,0) ; 3 % hors courbe |

- **Basal** : toutes les lectures donnent une [Cl⁻]i de −1 à 4 mM (hybride 2,9 mM), comme l'étape 7 (4,4 mM).
- **Le τ mono-exponentiel** et le temps de vie unique du modèle sont la même lecture (τ presque identiques soma
  par soma) : bonnes à 0 et 20 mM, ils perdent pied au-dessus (la courbe τ(c) s'aplatit ; 37 % des somas à
  131 mM hors de la courbe). Le FLIM du système aussi (33 % hors courbe à 131 mM) ; l'hybride n'en perd que 3 %
  (médiane 166 mM à 131 mM, avec une forte dispersion vers le haut).
- **À 0 mM, l'hybride lit 1,4 mM** : η du FOV y vaut 1,05–1,16 au lieu de 1. Quand c était bornée à 0, les somas
  dont le ratio dépasse R₁(0) tiraient η vers le haut et l'hybride lisait 3,6 mM ; garder les négatifs réduit cet
  effet de bord sans l'annuler.
- **Figure** (`deux_canaux_lectures.png`, sans le τ unique) : chaque boîte sans ses lectures hors courbe ni celles à
  plus de 5 σ de sa moyenne (retrait répété) ; axe y commun, le plus petit intervalle qui contient la moyenne ± 5 σ
  de chaque boîte.
- **η du FOV** varie de 0,74 à 1,21 ; à 20 mM, les deux FOV pâles ont η̂ 0,91–0,93 contre ~1,0 pour les autres :
  η absorbe en partie l'écart entre champs relevé plus haut.
- **Ce qui disperse les somas n'est pas le bruit de photons.** L'intervalle à 95 % d'un soma (hybride) fait 0,7
  à 2 mM sous 20 mM, alors que l'IQR entre somas d'une même [Cl⁻] imposée est de 6 à 31 mM. Dans une condition,
  la lecture hybride suit le rapport brut presque exactement (ρ de Spearman −0,95 à −0,99) mais à peine le τ
  mono (ρ 0,28 à 0,36) : la variation d'une cellule à l'autre est une variation du **rapport** que la durée de
  vie ne voit pas — sans doute pas du Cl⁻ (expression, maturation de la YFP, pH, autofluorescence).

**Un seul temps de transfert** (`--transfert 1`, sorties `deux_canaux_*_1t`) : une population en FRET au lieu de
deux (donneur f·e^{−t/τT} + (1 − f)·e^{−t/τD}, une isotherme). La calibration ajuste nettement moins bien
(déviance/ddl **182** contre 48) : τT 0,27, τD 1,91, τA 2,96 ns, f de 0,61 à 0,23 (K 36 mM), λ 0,62, et le jaune
devient tout entier de l'excitation directe (G = 0, κ 2,4). Lectures, médiane (±1,48·MAD entre somas ; écart-type
des médianes de FOV ; % hors courbe) :

| | 2 populations | 1 population |
|---|---|---|
| FLIM du système, basal / 0 / 20 / 131 mM | 1,8 / 0,2 / 18,3 (±13 ; FOV 25) / 236 (33 % hc) | 1,2 / −0,1 / 14,2 (±10 ; FOV 11) / 67 (7 % hc) |
| hybride, basal / 0 / 20 / 131 mM | 2,9 / 1,4 / 20,0 (±15 ; FOV 20) / 166 (3 % hc) | 2,4 / 1,2 / 17,3 (±12 ; FOV 16) / 57 (0 % hc) |
| ratio (η = 1), mêmes conditions | 3,9 / −0,9 / 22,0 / 113 | 3,8 / −1,0 / 22,2 / 116 |

Avec un seul temps de transfert, les lectures FLIM sont **plus resserrées** jusqu'à 20 mM (dispersion et écart
entre FOV plus petits, moins de somas hors courbe) mais **biaisées vers le bas en haut de gamme** (60 mM au lieu
de 131) : le modèle ne décrit pas la forme des déclins à 131 mM. Le basal ne bouge presque pas (1–4 mM).

### Réserves

1. **Les champs diffèrent autant que les conditions.** À 20 mM, les médianes des FOV vont de 2,0
   à 3,5 ; les FOV pâles ont un rapport plus bas **et** un τ CFP plus long (ρ rapport–intensité
   = 0,66 entre somas à 20 mM). Ce n'est donc pas un simple fond sur le rapport : soit une vraie
   différence entre champs (état des cellules, pH, clamp incomplet), soit de l'autofluorescence
   qui touche les deux lectures. Les FOV du basal, imagés à plus faible puissance, sont les plus
   pâles : si la pâleur abaisse le rapport, leur [Cl⁻] est surestimée. Une tranche sans capteur
   (autofluorescence) trancherait.
2. **K est mal déterminé** (3–28 mM pour le rapport) : il repose sur le seul point à 20 mM (K
   passe de 10 à 23 mM selon les FOV gardés). Une concentration de plus vers 5–10 mM est ce qui
   resserrerait le plus la calibration ; la [Cl⁻]i basale, proche de X₀, en dépend moins.
3. **Équilibre à 0 mM incertain** : basal→0 n'a pas de plateau en 49 min. La dérive tardive va
   au-delà des FOV à 0 mM, donc elle n'est sans doute pas que du Cl⁻ (pH ? 70 cycles
   d'illumination ?), mais si les FOV à 0 mM n'étaient pas à l'équilibre, X₀ est sous-estimé et
   la [Cl⁻]i basale aussi.
4. **pH** : les ionophores imposent aussi le pH ; en basal, le pH n'est pas celui des solutions
   de calibration, et la YFP y est sensible.
5. **IRF** : utilisée telle quelle. Sa courbe CFP (ligne 2 du `.sdt`) porte 19 % de son aire
   dans une traîne lente : gardée, elle raccourcit les durées de vie absolues (τ longue
   ~2,0 ns ; ~2,9 ns avec une IRF sans traîne). Les lectures de calibration étant relatives,
   K et la [Cl⁻]i basale en dépendent peu ; les τ absolus, oui. Attribution des deux courbes
   (ligne 0 = IN1, ligne 2 = IN2) faite dans l'ordre des entrées, à confirmer. Les résidus des
   déclins réunis montrent un écart aux deux bords de la fenêtre (avant la montée, et après
   11,5 ns) : effet du TDC près de ses limites, faible au niveau d'un soma.

## Sorties

| | |
|---|---|
| `qc/inventaire.csv` | chaque cycle : heure, durée, photons, dérive, plan focal, gardé |
| `qc/validation_sdt.csv`, `qc/recalage_bliq.csv` | contrôles du décodage |
| `qc/derive.png`, `qc/apercu_acquisitions.png` | dérive des transitions ; somme de chaque acquisition |
| `qc/irf.json`, `qc/irf_fovs.csv`, `qc/irf.png`, `qc/irf_somas.png` | IRF, τ communes, amplitudes par FOV |
| `masks/<acq>_masks.npy`, `masks/<acq>_overlay.png` | ROIs Cellpose |
| `results/somas.csv` | une ligne par ROI × FOV : photons, rapport, μ, τ mono, τm, populations |
| `results/conditions.csv`, `results/calibration.csv` | résumé par condition ; isothermes et IC |
| `results/basal_cellules.csv`, `results/autocoherence.csv` | lectures cellule par cellule |
| `results/calibration.png`, `basal.png`, `rapport_vs_tau.png`, `temps_mort.png` | figures |
| `results/transitions_*.csv`, `results/transition_<nom>.png` | cinétiques |
| `results/deux_canaux_*` | FLIM à deux canaux : calibration (json), lectures par soma, η par FOV, figures |
