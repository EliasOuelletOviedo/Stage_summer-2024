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
python3 step10_validation_croisee.py ../SCLM_calibration_2026-11-09              # tirage aléatoire des ROIs
python3 step10_validation_croisee.py ../SCLM_calibration_2026-11-09 --par-fov    # un FOV par pli
python3 step11_fond.py ../SCLM_calibration_2026-11-09
python3 step12_banc_essai.py ../SCLM_calibration_2026-11-09      # banc d'essai (~15 min ; --sans-calcul ensuite)
python3 step13_additif.py ../SCLM_calibration_2026-11-09          # signal du tissu (après l'étape 12)
python3 step14_ancrage.py ../SCLM_calibration_2026-11-09          # rapport seul contre rapport ancré (après l'étape 13)
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

### En bref (état après le banc d'essai des étapes 12 et 13)

- **Lecture à retenir : le rapport YFP/CFP**, cellule par cellule et dans le temps. À photons égaux il est 6 à 9 fois
  plus précis que toute lecture de la forme du déclin (2,4 mM-éq. à 1 000 photons contre 15 à 27), son bruit est le
  bruit de photons, et il voit un changement dans une même cellule 4 à 5 fois mieux. La forme du déclin CFP sert de
  **vérification à l'échelle d'un champ** (photons réunis), lue par gabarits mesurés plutôt que par ajustement.
- **Les champs à 20 mM diffèrent réellement** : le rapport et la forme du déclin CFP, deux mesures indépendantes, leur
  donnent les mêmes écarts (r = 0,91). L'« erreur par champ » de l'étape 10 n'était pas une erreur de lecture : le
  clamp n'est pas uniforme, et c'est lui qui limite K (3 à 23 mM selon les champs retenus) et la [Cl⁻]i basale.
- **Le tissu ajoute une composante prompte dans chaque pixel** (elle arrive avec l'impulsion, dans les deux canaux) :
  3 à 24 % des photons CFP proches du pic pour un soma médian. Elle fausse les lectures de forme des cellules pâles
  et faisait paraître ces lectures plus robustes d'un champ à l'autre qu'elles ne le sont. Elle se mesure sans
  tranche témoin (régression entre cellules) et se retire.
- **[Cl⁻]i basale : quelques mM** — 3,5 à 4,5 mM par le rapport, 6 à 9 mM par la forme CFP une fois le signal du tissu
  retiré. Les sections qui suivent gardent l'ordre dans lequel l'analyse a avancé ; celles des étapes 9 à 11 sont à
  lire avec les révisions des étapes 12 et 13.


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

### Validation croisée : quelle lecture est la plus juste ? (étape 10)

Cinq plis : le pli k retire le k-ième FOV de chaque [Cl⁻] imposée, recalibre (isothermes de l'étape 7, modèle de
l'étape 9 en 2 et en 1 population) sur les 4 autres FOV de chaque condition, et lit les 3 FOV retirés ; le basal est
lu par chaque pli. Critère fixé d'avance : erreur de la médiane de chaque FOV lu par rapport à la [Cl⁻] imposée
(erreur quadratique moyenne sur les 5 FOV ; `results/validation_croisee_resume.csv`, `validation_croisee.png`).

| Lecture | 0 mM | 20 mM | 131 mM | basal (médiane des plis [min ; max]) |
|---|---|---|---|---|
| rapport brut (isotherme) | 2,9 | 30,6 | ∞ | 3,9 [3,1 ; 5,9] |
| τ mono (isotherme) | 1,6 | 10,7 | ∞ | 4,6 [3,4 ; 5,2] |
| τm (isotherme) | 2,2 | **8,4** | ∞ | 9,6 [7,5 ; 10,3] |
| τ 3 populations (isotherme) | 1,7 | 10,7 | ∞ | 4,6 [3,4 ; 5,3] |
| temps de vie unique | 6,1 | 13,8 | 76 | −1,1 [−1,5 ; 0,0] |
| ratio, η = 1 | 3,2 | 30,0 | ∞ | 3,8 [2,7 ; 4,8] |
| FLIM du système, 2 pop. | 0,8 | 26,3 | 717 | 2,2 [1,1 ; 2,4] |
| **FLIM du système, 1 pop.** | **0,7** | 11,6 | 82 | 1,3 [0,6 ; 1,6] |
| hybride, 2 pop. | 3,6 | 25,0 | 296 | 3,0 [1,9 ; 3,6] |
| hybride, 1 pop. | 2,6 | 17,0 | 95 | 2,6 [1,5 ; 2,7] |

(erreur quadratique moyenne des médianes de FOV, en mM ; « ∞ » : au moins un FOV dont la médiane est hors courbe.)

- **À 0 mM**, tout marche (erreur ≤ 3 mM sauf le temps de vie unique) ; le FLIM du système fait le mieux (0,7 mM).
- **À 20 mM**, l'erreur vient de deux FOV : fov3 et fov5, les deux FOV pâles, sont lus de 29 à 75 mM par toutes les
  lectures quand ils sont hors calibration. Les lectures fondées sur le **rapport** en souffrent le plus (rapport,
  ratio : 30 mM d'erreur ; ces deux FOV à 58–75 mM) ; celles fondées sur la **durée de vie** et le FLIM du système à
  une population, le moins (8 à 12 mM ; ces FOV à 29–37 mM). L'écart entre champs est donc surtout dans le rapport.
- **À 131 mM, aucune lecture n'est juste** (erreur ≥ 64 mM, ou FOV hors courbe) : le haut de gamme n'est pas
  lisible en tranche avec cette calibration.
- **Meilleure lecture pour le domaine utile (0–20 mM) : le FLIM du système à une population** (0,7 et 11,6 mM,
  biais −0,1 et −1,7 mM), suivie des lectures de durée de vie par isotherme (τ mono, τ 3 populations, τm). Le
  modèle à une population, moins bon pour décrire les déclins (déviance/ddl 182 contre 48), lit mieux sous 20 mM ;
  il échoue au-dessus (biais −42 mM à 131 mM).
- **[Cl⁻]i basale** : 1 à 5 mM selon la lecture (τm à part, 9,6 mM ; temps de vie unique ~0) ; le FLIM du système à
  une population donne 1,3 mM [0,6 ; 1,6], le τ mono 4,6 mM [3,4 ; 5,2].

Réserves : 5 FOV par condition seulement, et l'erreur à 20 mM tient à deux d'entre eux ; les τ communes du modèle
« 3 populations » (étape 5) et l'IRF paramétrique sont ajustées sur tous les FOV (fuite faible vers les plis).

**Révisé aux étapes 12 et 13** : l'erreur par champ à 20 mM ne mesure pas la justesse d'une lecture. Les deux champs
« mal lus » ont réellement un autre état du capteur (le rapport et la forme CFP le voient pareil), et les lectures
de forme n'y paraissaient meilleures que parce que le signal prompt du tissu, plus fort dans ces champs pâles, tirait
leur lecture vers le bas.

### Validation croisée par tirage aléatoire des ROIs (étape 10, mode par défaut)

Dans chaque [Cl⁻] imposée, les ROIs sont mélangés (graine fixe) et répartis en 5 plis ; chaque pli calibre sur les
4/5 des ROIs et lit le 1/5 restant (η de l'hybride : par FOV, sur les ROIs lus). Critère : médiane des ROIs lus de
chaque pli contre la [Cl⁻] imposée (erreur quadratique moyenne sur les 5 plis ; `results/validation_aleatoire_*`).

| Lecture | 0 mM | 20 mM | 131 mM | basal (médiane des plis [min ; max]) |
|---|---|---|---|---|
| rapport brut (isotherme) | 0,4 | 1,8 | 22 | 4,4 [4,3 ; 4,6] |
| τ mono (isotherme) | 0,5 | 2,2 | ∞ | 4,3 [3,9 ; 4,4] |
| τm (isotherme) | 0,7 | 1,5 | 123 | 8,9 [8,6 ; 9,3] |
| ratio, η = 1 (2 pop.) | 0,9 | 3,0 | 21 | 3,8 [3,7 ; 4,0] |
| FLIM du système, 2 pop. | 0,3 | 2,7 | 126 | 1,9 [1,5 ; 2,3] |
| FLIM du système, 1 pop. | 0,4 | 6,3 (biais −6,2) | 64 | 1,2 [0,9 ; 1,4] |
| hybride, 2 pop. | 1,5 | 1,1 | 139 | 2,9 [2,7 ; 3,0] |
| hybride, 1 pop. | 1,2 | 2,6 | 75 | 2,4 [2,2 ; 2,6] |

- **Les erreurs sont bien plus petites qu'en retirant des FOV** (rapport brut 1,8 mM contre 30,6 mM par FOV) : les champs servent à la fois à calibrer et à
  lire, et la médiane d'un cinquième des somas pris au hasard ressemble à celle des quatre autres cinquièmes. L'écart
  entre les deux validations mesure ce que coûte un **nouveau champ** : c'est l'écart entre champs, pas
  l'échantillonnage des cellules, qui domine l'erreur.
- Ici, les lectures par **isotherme sur les médianes** (rapport, τ mono, τm) sont les plus justes à 0 et 20 mM : par
  construction, leur courbe passe par les médianes, et c'est la médiane qui est jugée. Les lectures du modèle à deux
  canaux sont calibrées sur les déclins réunis, dominés par les somas brillants : le FLIM du système à 1 population
  y garde un biais à 20 mM, alors qu'il était le meilleur en retirant des FOV.
- À 131 mM, le rapport et le ratio sont les moins mauvais (erreur ~17–22 mM) ; le modèle à deux canaux reste faux.
- [Cl⁻]i basale stable d'un pli à l'autre : 1,2 à 4,4 mM selon la lecture (τm à part, 7–9 mM).

### D'où vient l'écart entre champs ? (étape 11)

Hypothèse testée : un fond non lié au capteur, plus lourd dans les champs pâles, abaisse leur rapport
(`step11_fond.py`, `results/fond.png`, `fond_fovs.csv`, `fond_validation.csv`).

- **Le fond uniforme ajusté à l'étape 9 est petit** : 1,5 à 6 % des photons, plus grand en part dans les FOV pâles
  (ρ −0,79 pour la CFP), mais pas constant en photons par pixel (0,04 à 0,39 par cycle dans le bleu) : ce n'est pas
  un simple bruit de détecteur. Le retirer change le rapport de 1 à 2 % ; la validation par FOV ne bouge pas (à
  20 mM, 29,8 → 30,5 mM d'erreur). **Hypothèse rejetée.**
- **L'écart n'est pas un effet de luminosité.** À intensité égale (somas de 10 à 40 photons/px/cycle), fov3 et fov5
  à 20 mM restent à un rapport de 2,0–2,1 et un τ CFP de 1,65–1,67 ns, contre 3,1–3,5 et 1,52–1,60 ns dans les
  trois autres FOV. Dans un même FOV, les somas brillants ont un rapport un peu plus haut (ρ +0,1 à +0,55 sous
  20 mM), mais cet effet est bien plus petit que l'écart entre champs.
- **Les deux lectures du FRET s'accordent** : dans ces champs, le rapport est plus bas *et* le donneur vit plus
  longtemps, c'est-à-dire moins de FRET. Un artefact optique (fond, diffusion selon la longueur d'onde) ne changerait
  pas la durée de vie. Ces deux champs ont donc vraiment moins de FRET : [Cl⁻]i plus haute que la consigne (clamp
  incomplet ? cellules en moins bon état ?) ou un autre facteur qui éteint le FRET (pH). L'analyse ne peut pas
  trancher ; c'est une propriété des champs, à contrôler à l'acquisition (suivre chaque champ jusqu'au plateau,
  noter la profondeur).

### Banc d'essai : quelle lecture pour ces tranches ? (étape 12)

Regard neuf sur les mêmes photons : 39 lectures (`../sclm_calibration/banc.py`) — rapports de photons ; forme du
déclin CFP sans ajustement (temps moyen, portes, phaseur) ou avec (une, deux exponentielles, trois populations) ;
forme du déclin YFP ; **gabarits** (un déclin lu comme mélange des motifs mesurés à 0 et à 131 mM, sans IRF ni
modèle) ; apprentissage sur 20 portes temporelles ; modèle de l'étape 9 — notées sur des critères fixés d'avance :
séparation des conditions soma par soma, écart entre champs, répétabilité d'un cycle à l'autre, précision selon le
nombre de photons (somme des cycles amincie par tirage binomial), exactitude en validation croisée, changement dans
une même cellule (time-lapse 20→131), dérive pendant la minute d'illumination. « mM-éq. » : écart de la lecture
divisé par sa pente moyenne entre 0 et 20 mM. Tout est dans `results/banc_lectures.csv` ; `banc_classement.png`.

| Lecture | AUC 0 / 20 mM | σ à 1 000 photons | biais à 1 000 ph. | bruit d'un cycle | z, 20→131 | dérive par min | dispersion vraie entre somas | ρ avec la brillance | erreur par champ, 0 / 20 mM |
|---|---|---|---|---|---|---|---|---|---|
| **R** — rapport YFP/CFP | 0,96 | **2,4** | −0,1 | 0,45 | **6,5** | −0,3 | **5,4** | −0,28 | 2,9 / 31 |
| R_tard — photons > 1 ns après le pic | 0,96 | 2,7 | −0,2 | 0,48 | 6,1 | −0,2 | 5,3 | −0,24 | 2,6 / 28 |
| gab_joint — gabarits, formes + rapport | 0,96 | 2,3 | −0,1 | 0,43 | 6,0 | −0,3 | 5,5 | −0,32 | 2,9 / 30 |
| lda_3 — (ln R, temps moyens) | 0,96 | 4,1 | −0,1 | 0,72 | 4,5 | −0,1 | 5,5 | −0,10 | 2,2 / 26 |
| gab_B — gabarits, forme CFP | 0,87 | 14,8 | −0,1 | 2,9 | 1,3 | +0,3 | 10,4 | +0,39 | 1,8 / 12 |
| gab_formes — gabarits, formes des 2 canaux | 0,85 | 10,7 | −1,0 | 1,6 | 2,1 | −0,1 | 8,4 | −0,56 | 4,7 / 22 |
| tm_B — temps moyen CFP | 0,87 | 21,8 | +0,1 | 3,9 | 1,4 | +0,9 | 10,1 | +0,33 | 2,1 / 8,7 |
| tphi_B — phaseur CFP | 0,87 | 22,3 | +2,2 | 4,2 | 1,3 | +0,7 | 10,0 | +0,38 | 1,7 / 12 |
| tau_mono — décalage fixé | 0,86 | 20,7 | +0,5 | 4,3 | 1,3 | +0,6 | 10,1 | +0,40 | 1,5 / 13 |
| tau_mono_s — décalage libre (FLIMApp) | 0,86 | 23,5 | −0,2 | 4,9 | 1,1 | +0,7 | 10,0 | +0,36 | 1,5 / 12 |
| tau_bi — deux exponentielles | 0,86 | 26,9 | +5,6 | 6,4 | 0,8 | +1,2 | 10,6 | +0,35 | 1,6 / 11 |
| tau_3pop — trois populations | 0,87 | 20,5 | −2,8 | 4,0 | 1,2 | +0,6 | 10,1 | +0,37 | 1,5 / 13 |
| mod1_flim — modèle, 1 population | 0,89 | 20,9 | +6,2 | 3,0 | 0,7 | −0,1 | 10,9 | +0,40 | 0,7 / 12 |
| mod2_flim — modèle, 2 populations | 0,92 | 19,7 | +5,8 | 2,4 | 1,0 | −0,1 | — | +0,27 | 0,8 / 26 |
| tm_queue_B — queue du déclin CFP | 0,85 | 52,3 | −1,5 | 9,3 | 1,1 | +1,2 | 10,9 | +0,04 | 2,7 / 7,5 |
| tm_Y — temps moyen YFP | 0,73 | 20,8 | −0,4 | 3,9 | 0,6 | +0,3 | 16,8 | −0,45 | 14 / 31 |
| ridge_portes — régression sur 20 portes | 0,95 | 10,5 | −1,4 | 1,8 | 0,8 | −0,2 | 5,1 | −0,29 | 2,2 / 23 |
| gbr_portes — arbres boostés | 0,96 | 3,6 | +0,1 | 1,1 | 0,3 | −0,1 | 2,3 | −0,14 | 0,2 / 10 |

(mM-éq. sauf AUC, z et ρ ; dérive positive : la lecture monte vers plus de Cl⁻ ; dispersion vraie : covariance des
deux moitiés de cycles d'un soma, donc sans bruit de photons.)

1. **Le rapport est 6 à 9 fois plus précis à photons égaux** (40 à 85 fois moins de photons pour la même précision) :
   2,4 mM-éq. à 1 000 photons contre 14,8 pour la meilleure lecture de forme et 21 à 27 pour les ajustements
   (`banc_precision.png`). Lui ajouter les formes ne gagne que 3 % (gab_joint).
2. **Son bruit est le bruit de photons** : d'un cycle de 10 s au suivant, le bruit observé vaut 1,05 fois celui que
   prédit l'amincissement ; d'un cycle au suivant dans les time-lapses (41 s), 0,9 à 1,0 fois le bruit de Poisson
   (`banc_time_lapses.csv`). Rien d'autre que les photons ne limite le suivi d'une cellule.
3. **Dans une même cellule**, le changement du time-lapse 20→131 est vu à z = 6,5 par le rapport, 0,7 à 1,6 par les
   lectures de forme. Les deux familles suivent pourtant le même processus : d'une cellule à l'autre, leurs
   changements sont corrélés (ρ = −0,44 à −0,55 pour 20→131, −0,67 à −0,77 pour basal→0).
4. **Parmi les lectures de forme, les gabarits font le mieux** : 1,4 à 1,8 fois plus précis que temps moyen,
   phaseur ou ajustements, sans biais à bas photons, sans IRF. Les ajustements n'apportent rien et se déplacent
   quand les photons manquent : deux exponentielles +5,6 mM-éq. à 1 000 photons, modèle de l'étape 9 +6 mM (+25 à
   300 photons, 8 à 11 % d'échecs). Décalage de l'IRF libre ou fixé : pareil.
5. **Sans effet** : la correction du temps mort, des ROIs rognés de 3 px ou réduits à leurs pixels brillants
   (mêmes notes que R) ; le neuropile soustrait dégrade (AUC 0,89). **La forme du déclin YFP ne porte presque rien**
   (AUC 0,73).
6. **L'apprentissage ne se transfère pas.** Le discriminant sur les 20 portes échoue hors de ses champs ; les arbres
   boostés sortent premiers du classement mais **classent au lieu de mesurer** : avec trois concentrations ils
   n'interpolent pas (z = 0,3 sur le vrai changement 20→131) et, les conditions se suivant dans le temps, ils
   peuvent reconnaître la séance plutôt que le chlore.
7. **Les champs à 20 mM diffèrent réellement** (`banc_champs_rapport_forme.csv`, `banc_deux_etats.png`). Photons
   réunis par champ, en part f du motif « 131 mM » : le rapport donne 0,43 / 0,53 / 0,86 / 0,53 / 0,83 et la forme
   du déclin CFP 0,44 / 0,44 / 0,83 / 0,63 / 0,73 — r = 0,91, pente 0,9, écart-type commun 0,18 contre 0,04 pour la
   différence. Deux mesures indépendantes (combien de photons dans chaque canal ; quand ils arrivent dans un seul)
   s'accordent : par l'isotherme de l'étape 7, ces champs se lisent 11 / 16 / 61 / 16 / 52 mM. À 0 mM, rien de tel
   (écarts de 0,08, non corrélés). Selon qu'on croit les trois champs « bas » ou les deux « hauts », K vaut 23 ou
   3 mM et la [Cl⁻]i basale 6,7 ou 1,0 mM : **la calibration est limitée par le clamp à 20 mM, pas par la lecture**.
8. **D'une cellule à l'autre dans un champ**, la dispersion vraie est de ±5 mM-éq. pour le rapport et ±10 pour la
   forme CFP, et les deux ne sont que faiblement liées (corrélation vraie 0,3 ; `banc_correlations.csv`) : ce qui
   distingue les cellules d'un champ n'est pas surtout le chlore. Rapport et forme dépendent de la brillance en
   sens opposés (ρ = −0,28 et +0,33 à +0,40) : voir l'étape 13.
9. **Dérive pendant l'illumination** : en une minute, le rapport glisse de −0,3 mM-éq., les lectures de forme de
   +0,3 à +1,2 (les photons baissent de 12 à 15 %).
10. **Classement** (moyenne des rangs sur les huit critères fixés d'avance) : arbres boostés (réserve du point 6),
    puis R_tard, R, gab_joint, lda_3, gab_B. Sans les deux critères qui supposent que chaque champ est à sa [Cl⁻]
    nominale (invalidés par le point 7) : R_tard, R et gab_joint à égalité avec les arbres, puis lda_3.
    Sur la moitié la plus brillante des somas de chaque champ (`banc_lectures_brillants.csv`), mêmes conclusions.

### Ce que le tissu ajoute dans chaque pixel (étape 13)

Dans un même champ, les somas pâles ont plus de photons très précoces dans le canal CFP et un temps moyen plus
court de 0,07 à 0,13 ns — autant que l'effet de 0 → 20 mM. Si chaque pixel reçoit, en plus du capteur, un signal
qui ne dépend pas de la quantité de capteur, alors pour chaque canal de temps, entre les somas d'un champ :
photons = σ(t)·(photons du soma) + ι(t)·(pixels × cycles). La régression donne σ, le motif du capteur seul, et ι,
ce qui s'ajoute (`results/additif_champs.csv`, `additif.png`).

- **Ce qui s'ajoute est une composante prompte** : elle arrive avec l'impulsion, dans les deux canaux, et dure
  moins de 0,5 ns. Autour du pic (± 0,25 ns), elle vaut 3 à 24 % des photons CFP d'un soma médian (2 à 21 % des
  YFP), le plus en basal, pris à plus faible puissance — ce qu'on attend d'une lumière qui croît moins vite avec
  la puissance que la fluorescence à deux photons (diffusion ou fuite de l'excitation ? autofluorescence brève ?).
- **Vérification hors échantillon** (ι mesuré sur une moitié des somas de chaque champ, retiré de l'autre) :

  | Lecture | ρ avec la brillance, avant → après | erreur par champ à 20 mM | basal (mM) |
  |---|---|---|---|
  | rapport R | −0,28 → −0,08 | 31 → 27 | 4,4 → 3,5 |
  | gabarit CFP | +0,40 → −0,01 | 12 → 27 | 2,4 → 5,9 |
  | temps moyen CFP | +0,33 → +0,03 | 8,7 → 22 | 5,4 → 9,1 |
  | phaseur CFP | +0,38 → −0,02 | 12 → 33 | 1,6 → 6,5 |

  La dépendance à la brillance disparaît, et **l'avantage des lectures de forme d'un champ à l'autre disparaît
  avec elle** : les champs pâles (3 et 5 à 20 mM) portent plus de signal prompt, qui tirait leur lecture de forme
  vers « moins de Cl⁻ » et compensait en partie leur état réellement différent. La dispersion vraie entre somas,
  elle, ne baisse pas (±5 et ±10 mM-éq.) : le signal du tissu n'explique pas ce qui distingue les cellules.
- **Cellule par cellule, sans régression** (`additif_prompt.csv`) : la forme de la composante prompte étant connue,
  on peut la laisser libre dans le mélange de gabarits (gabp_B). La dépendance à la brillance disparaît aussi
  (ρ = −0,07) et la corrélation vraie avec le rapport monte de 0,33 à 0,50, au prix d'une précision deux fois
  moindre (35 mM-éq. à 1 000 photons contre 17). Cette variante enlève aussi la part du signal du capteur qui se
  trouve près du pic : les champs 3 et 5 à 20 mM y sont lus f ≈ 0,5 à 0,6 au lieu de 0,7 à 0,9. Le début et la
  suite du déclin ne disent donc pas exactement la même chose dans ces champs : l'état du capteur ne tient pas en
  un seul nombre (deux populations en FRET, ou [Cl⁻] et pH).
- **Basal** : une fois le signal du tissu retiré, le rapport donne 3,5 mM et les formes 6 à 9 mM (f = 0,25 contre
  0,38, photons réunis par champ). Les réglages des cartes sont les mêmes qu'en calibration (seul le
  déclenchement diffère) ; l'écart reste inexpliqué.

### Rapport seul contre rapport ancré par la forme (étape 14)

Ancrage : pour chaque champ, la forme du déclin CFP (motif du capteur seul, composante prompte retirée) donne la part
f du motif « 131 mM », d'où le rapport attendu R(f) = R₀·(1 − f) + R₁·f ; η̂ = rapport mesuré / R(f), et chaque soma
est lu à R / η̂. Gabarits et droite R(f) calibrés sur les champs d'entraînement à 0 et 131 mM ; même isotherme pour les
deux lectures ; validation croisée par champ. Pour simuler une autre tranche ou un autre réglage, le η du champ lu
est changé par amincissement binomial d'un canal (`results/ancrage_resume.csv`, `ancrage.png`).

| Erreur de la médiane du champ lu (mM) | η = 0,7 | 0,85 | 1 (données) | 1,15 | 1,3 |
|---|---|---|---|---|---|
| 0 mM, rapport seul | 11,9 | 5,9 | **2,9** | 2,8 | 3,9 |
| 0 mM, rapport ancré | 4,3 | 4,3 | 4,2 | 4,3 | 4,3 |
| 20 mM, rapport seul | 135 | 56 | 31 | 19 | 15 |
| 20 mM, rapport ancré | 44 | 44 | 44 | 44 | 43 |
| basal lu, rapport seul | 17,1 | 9,0 | 4,5 | 1,6 | −0,5 |
| basal lu, rapport ancré | 7,9 | 7,8 | 7,6 | 7,7 | 7,4 |

- **L'ancrage fait ce qu'on attend** : la lecture ne dépend plus de η (0,7 à 1,3 : moins de 0,4 mM de déplacement en
  basal, contre −5 à +14 mM pour le rapport seul).
- **Il coûte quand η ne change pas** : 4,2 mM contre 2,9 à 0 mM. η̂ d'un seul champ est incertain de 10 % à 0 mM,
  14 % à 20 mM, 6 à 8 % en basal, et inutilisable à 131 mM (64 %) : non par manque de photons, mais parce que la
  forme varie d'un champ à l'autre. L'ancrage par champ ne gagne que si η change de plus de ~10–15 %.
- **Il porte un biais qui dépend de l'état** : η̂ vaut 1,02 en moyenne à 0 mM, 0,96 à 20 mM mais **1,15 en basal**
  (les cinq champs, 1,08 à 1,25). Rapport et forme ne racontent pas la même chose en basal (pH sans ionophores ?
  puissance quatre fois plus faible ?) ; l'ancrage reporte cet écart sur la lecture : basal 7,6 mM au lieu de 4,5.
  Rien dans ces données ne dit lequel des deux est juste.
- **Conséquence** : ancrer au niveau de la tranche ou de la séance (moyenne des η̂ de tous ses champs : ±6 % avec
  5 champs, ±4 % avec 10), seulement quand η peut avoir changé (autre tranche, autre réglage optique, eNpHR-YFP), et
  calibrer l'ancrage dans l'état physiologique des mesures — sans ionophores, l'écart de 15 % du basal le montre.

### Quelle lecture retenir

- **Pour suivre une cellule ou comparer des cellules : le rapport YFP/CFP**, en photons bruts, converti par
  l'isotherme en tissu. Variante sans risque : ne compter que les photons arrivés plus de 1 ns après le pic
  (R_tard), insensible à la composante prompte pour 10 à 15 % de précision en moins.
- **Pour corriger un η qui a pu changer (autre tranche, réglage, eNpHR) : ancrer le rapport par la forme au niveau
  de la tranche**, pas du champ (étape 14).
- **Pour vérifier le niveau d'un champ ou d'une séance : la forme du déclin CFP par gabarits**, sur les photons
  réunis (au moins 10⁵ : ±1,5 mM-éq. ; avec la composante prompte libre, ±3,5), comparée au rapport. Un désaccord
  signale un champ hors de la courbe de calibration.
- **À éviter cellule par cellule** : les ajustements (biais à bas photons), la forme YFP, l'apprentissage sur trois
  concentrations.
- Ce que l'analyse ne peut pas régler : l'état réellement différent des champs à 20 mM, et ce qui distingue les
  cellules d'un même champ (±5 mM-éq. sur le rapport).

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

6. **Une seule tranche.** Le rapport gagne ici, où la transmission relative des deux canaux varie peu d'un champ
   à l'autre (±3 à 5 % à 0 mM et en basal) ; entre tranches ou préparations elle peut varier davantage, et seule
   la forme du déclin en est indépendante. Avec eNpHR-YFP dans les cellules (tranches d'octobre), le rapport absolu
   est dilué par une YFP sans FRET : il ne vaut alors que pour des variations.
7. **Critères du banc.** Les mM-éq. sont linéaires autour de 0–20 mM ; l'exactitude par champ suppose chaque champ à
   sa [Cl⁻] nominale, ce que la comparaison rapport / forme contredit à 20 mM. Les conditions se suivent dans le
   temps (basal, 0, 20, 131 mM) : toute lecture apprise peut confondre le chlore et l'heure.
8. **Régression entre cellules** (étape 13) : elle suppose les somas d'un champ dans le même état, et cet état
   indépendant de leur brillance ; vrai sous clamp, à vérifier ailleurs. Une tranche sans capteur dirait
   directement ce que le tissu ajoute.

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
| `results/validation_croisee_*` | validation croisée par FOV (étape 10, `--par-fov`) |
| `results/validation_aleatoire_*` | validation croisée par tirage aléatoire des ROIs (étape 10) |
| `results/fond*` | le fond explique-t-il l'écart entre champs ? (étape 11) |
| `results/banc_lectures.csv`, `banc_classement.png` | banc d'essai : toutes les notes de chaque lecture (étape 12) |
| `results/banc_precision.*`, `banc_somas.csv`, `banc_lectures_brillants.csv` | précision selon les photons ; valeurs par soma ; somas brillants |
| `results/banc_champs*.csv`, `banc_deux_etats.png`, `banc_champs.png` | champs : rapport contre forme, phaseur, lectures en mM |
| `results/banc_correlations.csv`, `banc_time_lapses.csv` | corrélations vraies entre lectures ; time-lapses cellule par cellule |
| `results/additif*` | ce que le tissu ajoute dans chaque pixel, notes avant / après (étape 13) |
| `results/ancrage*` | rapport seul contre rapport ancré par la forme, η simulé (étape 14) |
| `qc/banc_instrument.json` | décalages de l'IRF et pic utilisés par le banc |
