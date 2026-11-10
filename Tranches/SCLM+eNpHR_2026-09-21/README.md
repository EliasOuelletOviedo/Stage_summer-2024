# Analyse SCLM + eNpHR — tranches 2026-09-21, ratio FRET C2/C3 **dans le temps**

Même démarche que `../sclm/` — ROIs segmentés puis triés à la main, ratio C2/C3 en
valeurs brutes — mais chaque acquisition est ici un **décours temporel** : le ratio est
calculé image par image au lieu d'une fois par condition.

Les données sont dans `Data/` et tout le reste est écrit à côté, dans ce dossier.

**Les étapes se lancent dans l'ordre** : chacune consomme la sortie de la précédente.
`step4` refuse de tourner tant qu'il n'y a pas de sélection de ROIs.

Seule l'étape 2 a besoin de Cellpose, donc de l'environnement dédié
`~/.venvs/cellpose/bin/python` (Cellpose-SAM 4.2.1, torch 2.13, MPS actif).
Les étapes 1, 3, 4, 5 et 6 n'utilisent que numpy / scipy / matplotlib / pandas / tifffile.

Toutes les commandes se lancent depuis ce dossier.

```bash
python3 step1_align.py                          # fenêtre de stim, recalage, projections
~/.venvs/cellpose/bin/python step2_cellpose.py  # segmentation des somas  (Cellpose requis)
python3 step3_select.py                         # tri manuel des ROIs (interactif)
python3 step4_ratio.py                          # ratio C2/C3 dans le temps + figures
python3 step5_decay.py                          # les 50 meilleures réponses
python3 step6_deltaR_vs_Rbasal.py                # réponse en fonction du R basal
```

FOV11 (train de stimulations) est analysé à part, par `train_FOV11.py`, qui écrit dans
`FOV11_train/` — voir la section « FOV11 » plus bas. Le pipeline ci-dessus l'ignore
(`SEPARATE_FOVS` dans `config.py`).

Les étapes 1 et 4 prennent environ une minute chacune, l'étape 2 une minute, l'étape 5
quelques secondes.

## Ce que contient `Data/`

Découvert automatiquement (`config.survey()`), rangé par ordre chronologique. Une
acquisition n'est retenue que si C2 et C3 existent et ont le même nombre d'images.

| FOV | Acquisition | Images | Heure | Stim détectée |
|---|---|---|---|---|
| FOV1 | `950 nm-001` | 1 | 21/09 14:33 | — |
| FOV1 | `800 nm-001` | 1 | 21/09 14:34 | — |
| FOV1 | `stim-001` | 150 | 21/09 14:40 | 30–90 (29,7–90,9 s) |
| FOV1 | `stim-002` | 88 | 21/09 14:46 | 39–60 (39,1–60,5 s) |
| FOV2 | `950-001` | 1 | 21/09 14:56 | — |
| FOV2 | `800-001` | 1 | 21/09 14:57 | — |
| FOV2 | `stim_VIS-001` | 150 | 21/09 15:00 | 30–89 (29,9–90,1 s) |
| FOV2 | `1064-001` | 1 | 21/09 15:04 | — |
| FOV2 | `stim_IR-001` | 150 | 21/09 15:05 | 30–89 (29,9–90,0 s) |
| FOV2 | `stim_max-001` | 150 | 21/09 15:39 | 30–90 (29,9–91,1 s) |
| FOV3 | `800-001` | 1 | 22/09 11:01 | — |
| FOV3 | `950-001` | 1 | 22/09 11:02 | — |
| FOV3 | `1064-001` | 1 | 22/09 11:05 | — |
| FOV3 | `stim_VIS_10-001` | 150 | 22/09 11:08 | 29–89 (28,6–89,8 s) |
| FOV3 | `stim_VIS_50-001` | 150 | 22/09 11:12 | 29–89 (28,9–90,0 s) |
| FOV3 | `stim_IR-001` | 150 | 22/09 11:15 | 30–90 (29,9–91,1 s) |
| FOV4 | `800-001` | 1 | 22/09 11:21 | — |
| FOV4 | `950-001` | 1 | 22/09 11:24 | — |
| FOV4 | `stim_VIS_10-001` | 150 | 22/09 11:26 | 30–89 (29,8–90,0 s) |
| FOV4 | `stim_VIS_50-001` | 150 | 22/09 11:36 | 30–89 (29,9–90,0 s) |
| FOV4 | `stim_IR-001` | 150 | 22/09 11:41 | 29–89 (28,8–90,2 s) |
| FOV5 | `800-001` | 1 | 22/09 13:32 | — |
| FOV5 | `950-001` | 1 | 22/09 13:34 | — |
| FOV5 | `stim_VIS_1-001` | 180 | 22/09 13:36 | 30–89 (29,9–90,1 s) |
| FOV5 | `stim_VIS_10-001` | 180 | 22/09 13:39 | 30–90 (29,9–91,1 s) |
| FOV5 | `stim_VIS_50-001` | 180 | 22/09 13:42 | 30–90 (29,9–91,1 s) |
| FOV6 | `800-001` | 1 | 22/09 13:48 | — |
| FOV6 | `950-001` | 1 | 22/09 13:50 | — |
| FOV6 | `stim_VIS_1-001` | 180 | 22/09 13:53 | 30–89 (29,9–90,0 s) |
| FOV6 | `stim_VIS_10-001` | 180 | 22/09 13:56 | 30–90 (29,9–91,2 s) |
| FOV6 | `stim_VIS_50-001` | 180 | 22/09 13:59 | 30–90 (29,8–91,0 s) |
| FOV7 | `800-001` | 1 | 22/09 14:08 | — |
| FOV7 | `950-001` | 1 | 22/09 14:10 | — |
| FOV7 | `stim_VIS_10-001` | 180 | 22/09 14:22 | 30–90 (29,8–91,0 s) |
| FOV7 | `stim_VIS_50-001` | 180 | 22/09 14:25 | 30–90 (29,9–91,1 s) |
| FOV7 | `stim_VIS_100-001` | 180 | 22/09 14:28 | 30–89 (30,0–90,3 s) |

Dans les fichiers Bliq, **C2 porte les numéros T impairs et C3 les pairs** : l'image
`k` du décours est donc `C2-T(2k-1)` + `C3-T(2k)`. Les deux canaux ont exactement le
même `DeltaT` dans l'en-tête OME — ils sont acquis simultanément, la paire est sûre.
L'axe des temps vient de ce `DeltaT` (≈ 1,02 s entre images), pas d'un pas supposé.

Les acquisitions à une seule image (800 / 950 / 1064 nm) n'ont pas de décours, donc pas
de ratio. **Seul leur C2, à 800 et 950 nm, est analysé** (`FIXED_CHANNEL`,
`FIXED_WAVELENGTHS` dans `config.py`) : C2 moyen par ROI, mesuré par l'étape 4 dans
`results/images_fixes_C2.csv`, et affiché par l'étape 1 dans
`projections/<FOV>_images_fixes.png`. C3 et 1064 nm sont laissés de côté. Toutes les
images fixes, 1064 nm compris, servent malgré tout au recalage et entrent dans l'image
somme qui sert à la segmentation.

## Un seul jeu de ROIs par FOV

Toutes les acquisitions d'un FOV visent le même champ (vérifié : NCC 0,62–0,96 entre
décours, 0,17–0,95 pour les acquisitions à une seule image, bien plus bruitées ;
dérive ≤ 11 px). Les ROIs sont donc segmentés **une seule fois
par FOV**, sur la somme des projections, et le ROI n° 12 désigne la même cellule dans
`stim_VIS`, `stim_IR` et `stim_max`. C'est ce qui permet de comparer les conditions de
stimulation cellule par cellule.

## Étape 1 — fenêtre de stim, recalage, projections

**Détection de la stim.** La lumière de stimulation illumine tout le champ et fuit
dans C2 : la moyenne pleine image y monte en un plateau rectangulaire, d'un facteur
1,1 à 30 selon l'acquisition.

La détection porte sur les **bords** du plateau, pas sur son niveau. Ce qui distingue
une stim d'un photoblanchiment, c'est que l'allumage et l'extinction sont francs : on
exige donc deux bords, montant puis descendant, d'au moins 6 écarts robustes. Un seuil
sur le niveau (médiane + 8 MAD) échouait dès que le plateau occupait une bonne part de
l'acquisition, puisqu'il gonflait lui-même l'estimation du bruit — jusqu'à ×25 sur
FOV3 — et se hissait hors de sa propre portée.

Le bord est mesuré en comparant les moyennes de **5 images de part et d'autre**, et
non deux images voisines : une marche franche y répond de toute sa hauteur, une dérive
lente de sa pente fois 5 seulement, et le bruit est divisé par √(5/2). Sans ce
fenêtrage, `FOV5/stim_VIS_1-001` — dont la fuite n'est que de +2,9 % — reste sous le
seuil (3,3 écarts sur des différences brutes, 9,1 avec la fenêtre).

Les deux bords trouvés, la fenêtre est délimitée par un seuil à mi-chemin entre le
niveau du plateau et celui du reste, ramenée au plus long bloc contigu, puis
**élargie d'une image de chaque côté** : l'obturateur s'ouvre et se ferme en milieu
d'image, qui est alors à moitié éclairée et tombe pile sur le seuil.

Le plateau doit durer au moins 10 images, soit plus que la demi-fenêtre : sans cette
condition, le fenêtrage étale une image aberrante isolée sur 5 images et en fait un
court plateau (3,4 % de faux positifs dans les essais).

Calibré sur 14 000 tirages de contrôles négatifs (bruit seul, photoblanchiment, dérive
linéaire, dérive avec rebond, 1 à 5 images aberrantes) : **aucun faux positif**. La
détection est fiable à partir d'environ +2,5 % de hausse ; la plus faible stim réelle
du jeu, `FOV5/stim_VIS_1-001`, est à +2,9 % et ressort à 9 fois le bruit.

**Dérive.** Estimée par **blocs de 25 images** et non image par image : à ces niveaux
(8 bits, moyenne ~7 counts) une image isolée est trop bruitée pour donner un pic de
corrélation fiable (NCC ~0,5 contre ~0,94 pour un bloc). Les blocs sont formés des
images hors stim uniquement — le fond uniforme de la stim fausserait la corrélation —
puis le décalage est interpolé linéairement sur toute l'acquisition. En pratique la
dérive intra-acquisition est ≤ 1 px partout.

**Recalage entre acquisitions.** Chaque projection moyenne hors stim est recalée par
corrélation de phase sur celle de l'acquisition de référence du FOV (le décours le
plus long, puis le plus ancien). Garde-fou identique à `../sclm/` : si le NCC est
`< 0.15` le pic n'est pas fiable et **aucun décalage n'est appliqué**.

Le plus gros décalage est `FOV2/stim_max-001` → `stim_VIS-001` (dy = +10 px), cohérent
avec les 39 minutes qui les séparent.

Sorties : `projections/<FOV>.npz`, `<FOV>_sum.tif/.png`, `<FOV>_projections.png`,
`<FOV>_images_fixes.png`, `<FOV>_align.json` (temps, fenêtre de stim et décalage de
chaque image) et `qc/alignement.tsv`.

## Étape 2 — Cellpose

Modèle **Cellpose-SAM** (`cpsam`) sur l'image somme de chaque FOV.

| Paramètre | Valeur | Raison |
|---|---|---|
| `diameter` | 60 px | soma ~30 µm à 0,509 µm/px, mesuré sur les images somme |
| `flow_threshold` | 0.4 | défaut, masques bien formés |
| `cellprob_threshold` | 0.0 | défaut ; plus bas ferait fusionner les somas serrés |
| `min_size` | 300 px | écarte débris et taches de neuropile |
| `max_size_fraction` | 0.25 | empêche un ROI d'avaler le champ |

35 ROIs sur FOV1, 61 sur FOV2, 41 sur FOV3, 59 sur FOV4, 75 sur FOV5, 37 sur FOV6,
56 sur FOV7. `model.eval()` normalise son entrée
**en place** : on lui passe toujours une copie.

## Étape 3 — tri manuel

| Touche | Effet |
|---|---|
| clic sur un ROI | garder / rejeter |
| `a` / `n` | tout garder / tout rejeter |
| `b` | change le fond (somme, puis chaque acquisition) |
| `l` | affiche/masque les numéros |
| `←` `→` | FOV précédent / suivant (enregistre au passage) |
| `q` | enregistre et quitte |

Tous les ROIs sont **gardés par défaut** : il suffit de cliquer ceux à rejeter. La
sélection vaut pour toutes les acquisitions du FOV. Elle est enregistrée dans
`selection/<FOV>_selection.json` et relue au lancement suivant, donc on peut reprendre
où on en était ou repasser sur un seul FOV : `step3_select.py FOV2`.

`python3 step3_select.py --tout` écrit la sélection par défaut (tout gardé) sans
ouvrir de fenêtre — pratique pour faire tourner l'étape 4 tout de suite, quitte à
trier ensuite et à relancer l'étape 4. **C'est ce qui a été fait ici : les résultats
présents reposent sur les 364 ROIs bruts de Cellpose, pas sur un tri visuel.**

## Étape 4 — ratio FRET dans le temps

`ratio = moyenne(C2) / moyenne(C3)` sur les pixels du ROI, **valeurs brutes** (aucune
soustraction de fond), pour chaque image du décours — exactement la définition de
`../sclm/`, appliquée image par image.

Le ROI est défini dans le repère du FOV mais l'image brute n'est pas recalée : c'est
le **masque** qui est déplacé de l'inverse du décalage, ce qui évite d'interpoler les
intensités. Un ROI dont moins de 90 % des pixels restent dans le champ est mis à `NaN`
pour cette image (colonne `note`).

`dR/R0` : `R0` est la **médiane du ratio brut avant la stim** (les 2 premières images
sont écartées, transitoire PMT), par ROI et par acquisition.

### Deux mesures de la réponse, et pourquoi

`dR_R0_apres` moyenne **toute** la période après la stim. C'est la mesure la plus
simple, mais elle suit la dérive lente de la ligne de base : sur FOV4 celle-ci monte
déjà de +0,6 % par 10 s **avant** la stim, ce qui suffit à faire passer
`stim_VIS_10-001` de négatif à positif.

`dR_R0_apres_immediat` ne moyenne que les **10 premières secondes** après la stim
(`RESPONSE_WINDOW_S`). La dérive n'a pas le temps d'y peser, et c'est là que la
réponse est la plus franche. C'est la colonne à utiliser pour comparer des
acquisitions entre elles.

`pente_avant_pct_par_10s`, dans `resume_acquisitions.csv`, donne la pente de la ligne
de base avant la stim : c'est le diagnostic qui dit à quel point la première colonne
est à prendre avec des pincettes pour une acquisition donnée.

### Les images de stim ne sont pas mesurables

Pendant la stim, C2 est noyé par la fuite lumineuse. Ces images restent dans le CSV
avec `stim = True` et `phase = pendant`, mais **le trait s'interrompt dans toutes les
figures** et elles n'entrent pas dans `R0`. Les colonnes `fond_C2` / `fond_C3`
(5e centile de l'image) et `ratio_fond_soustrait` sont des diagnostics : hors stim le
fond vaut 0 et le ratio corrigé est identique au ratio brut ; pendant la stim la
correction ne **rétablit pas** la valeur de base (54 → 18 alors que la base est 5,3
pour `stim_max`), la fuite n'étant pas un simple décalage uniforme. Il n'y a pas de
saturation (aucun pixel à 255), mais le signal FRET reste inexploitable là.

Ce qui est mesurable, c'est **l'état avant et après** la stimulation.

### Sorties dans `results/`

| Fichier | Contenu |
|---|---|
| `ratio_temps.csv` | une ligne par ROI × acquisition × image (171 500 lignes) |
| `resume_par_roi.csv` | `R0`, ratio avant / après, `dR_R0_apres` et `dR_R0_apres_immediat` par ROI × acquisition |
| `resume_acquisitions.csv` | n, moyenne, médiane, écart-type et `pente_avant_pct_par_10s` par acquisition × phase |
| `trace_<FOV>_<acq>.png` | ratio brut et `dR/R0` dans le temps : un trait par ROI, médiane et écart interquartile, bande de stim |
| `trace_<FOV>_toutes.png` | les acquisitions du FOV côte à côte, axe y commun |
| `heatmap_<FOV>_<acq>.png` | `dR/R0` en ROI × temps, ROIs triés par réponse après la stim, avec la **moyenne des ROIs en bande au-dessus**, sur la même échelle de couleurs |
| `avant_apres_stim.csv` | `dR/R0` moyen par ROI et par stim, 10 s avant le début de la stim et 10 s après sa fin, avec le type de stim (VIS / IR) |
| `avant_apres_stim.png` | toutes les stims VIS regroupées dans un panneau, toutes les IR dans l'autre : un trait par ROI et par stim, médiane en gras, même axe y. Les stims sans étiquette (FOV1, `FOV2/stim_max`) n'y entrent pas |
| `images_fixes_C2.csv` | C2 moyen par ROI des images fixes, à 800 et 950 nm |
| `images_fixes_C2.png` | ces mêmes valeurs : un trait par ROI de 800 à 950 nm, médiane en gras, un panneau par FOV |
| `retour_base_meilleurs_neurones.png` | (étape 5) les 50 meilleures réponses sur toute l'acquisition, disposition des figures trace |
| `retour_base_neurones.csv` | (étape 5) les réponses retenues (FOV, acquisition, ROI) et leur amplitude |
| `retour_base_ajustements.csv` | (étape 5) `A`, `tau`, `C`, `R2` de l'ajustement après la stim, réponse par réponse |
| `deltaR_vs_Rbasal.png` | (étape 6) ΔR absolu après toutes les stims VIS en fonction du R basal, un point par ROI et par stim, avec la droite de tendance intra-acquisition ; les statistiques sont imprimées par le script |
| `deltaR_vs_Rbasal.csv` | (étape 6) R basal, R avant, R après, ΔR et ΔR/R0 par ROI et par stim étiquetée |

Les figures lissent sur 5 images la médiane, l'écart interquartile et la carte
`heatmap` ; les traits fins par ROI et **tout le CSV restent en valeurs brutes**.

## Étape 5 — les meilleures réponses sur toute l'acquisition

**Sélection.** Les `N_BEST_TRACES` (50) réponses les plus fortes de tout le jeu, jugées
sur `dR_R0_apres_immediat` — une réponse étant un couple ROI × stim, un même neurone peut
apparaître plusieurs fois. `FOV1/stim-002` est écartée (`MIN_POST_S`) : 88 images, 29 s
seulement après la stim, un calendrier différent des autres.

**Figure.** Exactement la disposition des figures `trace_<FOV>_<acq>.png` (même
fonction `trace_figure`) : ratio C2/C3 brut en haut, `dR/R0` en bas, sur toute la
durée, un trait par réponse, médiane et écart interquartile lissés. La stim est
masquée sur la réunion des fenêtres de toutes les réponses retenues — elles diffèrent
d'une image selon l'acquisition, et sans masque commun la médiane sauterait aux bords.

Les acquisitions durent 150 ou 180 images : au-delà de 152 s, la médiane ne porte plus
que sur les traits des acquisitions longues (23 sur 50). L'ajustement
`A·exp(-t/tau) + C` de la médiane après la stim, limité à la fenêtre où au moins 80 % des
traits contribuent, est imprimé par le script mais n'est pas tracé.

**Biais de sélection.** Les réponses sont choisies sur leur amplitude juste après la
stim puis tracées sur cette même fenêtre : cette amplitude est surestimée (régression
vers la moyenne), et une partie des réponses retenues sont des ROIs bruités plutôt que
de vrais répondeurs — trois viennent de `FOV5/stim_VIS_1-001`, une condition qui ne
répond pas en moyenne. La forme du retour, et donc `tau`, en souffre beaucoup moins.

## Étape 6 — la réponse dépend-elle du R basal ?

Question posée : les cellules à R élevé (peu de Cl⁻) répondent-elles plus fort ?

- **R basal** : ratio moyen au début de la ligne de base (de l'image 3 jusqu'à 10 s
  avant la stim), une fenêtre distincte de celle du ΔR — sinon le bruit de mesure
  créerait à lui seul une corrélation négative (régression vers la moyenne).
- **ΔR** : ratio des 10 s après la stim moins celui des 10 s avant ; **ΔR/R0** : le
  même rapporté au ratio des 10 s avant.
- Tendances **intra-FOV** (valeurs centrées par tranche, Theil-Sen, Spearman sur les
  rangs) : les tranches diffèrent par leur R basal et leur amplitude de réponse.

| | rho intra-FOV | pente |
|---|---|---|
| A · ΔR absolu, VIS 50 % | −0,44 (p < 0,001) | −0,07 par unité de R |
| B · ΔR/R0, VIS 50 % | −0,28 (p < 0,001) | −1,7 % par unité de R |
| B · ΔR/R0, témoins IR + VIS 1 % | −0,29 (p < 0,001) | −1,5 % par unité de R |
| C · VIS 50 % − témoin, même cellule (n = 210) | **−0,04 (p = 0,6)** | −0,7 % [IC95 −2,3 ; +0,9] |

En ΔR absolu, les cellules à R élevé baissent bien davantage, mais c'est en grande
partie mécanique : une même baisse relative de 3 % donne −0,045 à R = 1,5 et −0,09 à
R = 3. En relatif, la tendance est plus faible **et la même existe sans réponse**, dans
les conditions témoins. Une fois retranchée la réponse de chaque cellule à son témoin,
il ne reste aucune dépendance détectable au R basal — sans pouvoir exclure un effet
modeste, vu la largeur de l'intervalle.

Le R brut n'est pas non plus une mesure calibrée du Cl⁻ : sans soustraction de fond, il
est corrélé à la luminosité de la cellule (rho ≈ 0,7 intra-FOV), et il change d'une
acquisition à l'autre pour les mêmes cellules (FOV2 : ≈ 2,7, puis ≈ 5,3 sur
`stim_max-001`).

## FOV11 — trains de stimulations, analysés à part

FOV11 (acquis le 24/09) contient deux trains de quatre stimulations, sur les mêmes
cellules, à 26 minutes d'écart. Chacun dure 270 images (275 s) :

| | pulses (images) | durée des pulses | début des pulses |
|---|---|---|---|
| `train_stim-002` | 33–42, 90–101, 150–161, 208–218 | ≈ 10 s chacun | 33, 91, 152, 211 s |
| `train_stim-003` | 31–38, 90–101, 150–166, 208–229 | ≈ 6, 10, 15 et 20 s | 31, 91, 152, 212 s |

Les étapes 1 à 6 ne traitent qu'une stimulation par acquisition. `train_FOV11.py`
reprend leurs fonctions, ajoute ce qu'il faut pour un train et écrit tout dans
`FOV11_train/`, une figure trace et une figure ΔR / R basal par train :

```bash
python3 train_FOV11.py                 # tout, de la détection aux figures
python3 train_FOV11.py --resegmenter   # idem en relançant Cellpose (renumérote les ROIs)
SCLM_OUT=FOV11_train SCLM_FOVS=FOV11 python3 step3_select.py FOV11   # tri manuel
```

Cellpose n'est lancé que si les masques n'existent pas encore, et la sélection par
défaut n'est écrite que s'il n'y en a pas : relancer le script ne défait pas un tri
manuel. `SCLM_OUT` et `SCLM_FOVS` (lus par `config.py`) redirigent les sorties et
restreignent les FOVs.

**Détection des pulses.** Même principe que la détection principale (bords fenêtrés de
la moyenne pleine image de C2), mais tous les bords sont gardés, et chaque bord montant
est apparié au premier bord descendant après lequel C2 est revenu près de la base. Sur
les 14 décours à stim unique des autres FOVs, elle retrouve les mêmes fenêtres
(±1 image) ; aucun faux positif sur 8 000 témoins simulés ; un train de 4 pulses est
retrouvé en entier dès +3 %. Les pulses de FOV11 font +33 % sur C2.

**Mesures.** R0 est la médiane avant le premier pulse ; les images entre deux pulses
sont en phase `apres`, la colonne `pulse` donne le numéro du pulse qui les précède.
ΔR = 10 s après chaque pulse − 10 s avant ; R basal est mesuré avant le premier pulse,
le même pour les quatre. Chaque cellule donne quatre points par train : les tests se
font sur la moyenne des quatre pulses par cellule (n = 37).

| ΔR/R0 médian | pulse 1 | pulse 2 | pulse 3 | pulse 4 |
|---|---|---|---|---|
| `train_stim-002` | −2,8 % | −2,1 % | −1,8 % | −2,0 % |
| `train_stim-003` | −2,5 % | −2,9 % | −3,5 % | −3,4 % |

Dans le premier train, chaque réponse revient presque à la base avant le pulse suivant.
Dans le second, aux pulses plus longs, les réponses s'additionnent et la ligne de base
descend d'un pulse à l'autre. D'un train à l'autre, les mêmes cellules gardent leur R
basal (ρ = 0,90) et leur réponse (ρ = 0,74).

ΔR contre R basal, par cellule : ρ = −0,68 puis −0,66 en ΔR absolu, **ρ = −0,53 puis
−0,45 en ΔR/R0**. FOV11 n'a aucune condition témoin (ni IR ni VIS 1 %) pour retrancher
la part de la tendance qui existe sans réponse, comme au panneau C de l'étape 6.

### Les neurones qui répondent le mieux

Pour chaque train, les 5 neurones à la **réponse typique** la plus forte — la médiane,
sur les 4 pulses, du ΔR/R0 (10 s après − 10 s avant), pour ne pas retenir un neurone sur
un seul pulse exceptionnel. `meilleurs_neurones_<train>.png` les montre avec la
disposition des figures trace, un neurone par couleur (la même d'un train à l'autre) ;
`meilleurs_neurones.csv` donne leurs réponses pulse par pulse.

| | ROIs (réponse typique) |
|---|---|
| `train_stim-002` | 23 (−8,4 %), 27 (−6,7 %), 11 (−5,9 %), 14 (−4,7 %), 3 (−3,6 %) |
| `train_stim-003` | 27 (−10,5 %), 23 (−7,2 %), 14 (−7,0 %), 3 (−6,3 %), 13 (−5,5 %) |

Tous baissent aux 4 pulses, et quatre neurones (3, 14, 23, 27) sont dans les deux
trains. Leur C3 n'est pas faible (7 à 10 counts, médiane 8) : ce ne sont pas des
cellules sombres dont le ratio serait instable.

### Avant / après chaque stimulation

`avant_apres_stim.png` reprend la présentation de la figure du même nom du pipeline
principal, avec un panneau par train : un trait par ROI et par stimulation, ΔR/R0 moyen
sur les 10 s avant le début de chaque pulse et les 10 s après sa fin (R0 = base avant le
premier pulse) ; valeurs dans `avant_apres_stim.csv`.

| | 10 s avant → 10 s après (médianes) | en baisse |
|---|---|---|
| `train_stim-002` | −0,2 % → **−2,3 %** | 140/148 |
| `train_stim-003` | −1,0 % → **−4,1 %** | 144/148 |

Pour le second train, dont les pulses ont quatre durées différentes,
`avant_apres_stim_train_stim-003_par_stimulation.png` montre chaque stimulation dans
son panneau (un trait par ROI, axe y commun ; les trains concernés sont listés dans
`PAR_PULSE`) :

| stimulation (durée) | 10 s avant → 10 s après (médianes) | écart médian par ROI | en baisse |
|---|---|---|---|
| 1 (≈ 6 s) | 0,0 % → −2,3 % | −2,5 pts | 37/37 |
| 2 (≈ 10 s) | −1,0 % → −4,0 % | −2,9 pts | 35/37 |
| 3 (≈ 15 s) | −1,9 % → −5,7 % | −3,4 pts | 36/37 |
| 4 (≈ 20 s) | −2,9 % → −6,8 % | −3,2 pts | 36/37 |

Le « avant » descend d'un pulse à l'autre : les réponses s'additionnent. L'écart propre
à chaque pulse grandit de 6 à 15 s puis plafonne à 20 s.

### Sauts instrumentaux sur C3 (second train)

Le second train montre cinq sauts brusques de C3 hors de toute stim (à 48, 66, 243, 251
et 266 s), **de même valeur en counts** (0,1 à 0,24) dans les pixels sombres, dans le
neuropile et dans les cellules : un décalage additif sur tout le canal, donc
instrumental. Ils tombent hors des fenêtres de mesure de la réponse, ou à leur dernière
image (pulses 1 et 4) ; la réponse typique (médiane des 4 pulses) n'en dépend pas.

**Correction d'une analyse antérieure.** J'avais d'abord attribué à un artefact la
hausse de C3 qui suit chaque pulse, en la comparant à des pixels de fond sombres. Ces
pixels, proches du zéro du détecteur 8 bits, ont leurs propres sauts — à l'image 26 du
premier train, +0,6 count dans le fond sombre, rien dans les cellules — et ne sont pas
un bon témoin. Sur le neuropile (pixels hors cellules, d'intensité intermédiaire), les
pulses produisent la même signature que dans les cellules, en plus petit : C2 (YFP)
baisse et C3 (CFP) monte. C'est la réponse attendue de SuperClomeleon (le Cl⁻ éteint la
YFP, le transfert d'énergie diminue, la CFP remonte), pas un artefact. Les conclusions
tirées de cette « correction du fond » — notamment que le premier pulse ne répondait
plus — sont donc caduques ; les chiffres de FOV11 sont ceux du ratio brut.

Deux autres réserves : le nom `train_stim` n'indique pas la longueur d'onde ; et l'image
somme montre une discontinuité horizontale vers le milieu du champ, présente dans les
images brutes, que traversent les ROIs 13 à 17.

## Slide de résultats préliminaires

```bash
python3 slide_resultats.py     # présentation/slide_resultats_FR.png et _EN.png, + _FR_dark / _EN_dark sur fond noir
```

Une figure en quatre panneaux, toute sur FOV11 (à lancer après `train_FOV11.py`), sans
autre texte que le titre et les titres de panneaux :

- **A** (en haut à gauche) — images à 800 nm et à 950 nm côte à côte, canal C2 (YFP),
  contraste propre à chaque image : à 800 nm SCLm remplit les corps cellulaires (noyau
  sombre), à 950 nm eNpHR marque les membranes (cellules en anneaux). En bas à droite de
  chaque image, un encart agrandit un carré de 150 px centré sur le ROI 6 (`ZOOM_ROI`, `ZOOM_COTE`,
  cadre jaune sur l'image) ; barre
  d'échelle de 50 µm en bas à gauche de l'image à 950 nm ;
- **B** (en haut à droite) — **tracé recollé** à partir du second train
  (`train_stim-003`) : la ligne de base (avant la première stimulation, 0–30 s), puis la
  stimulation 2 et son retour, décalés dans le temps pour suivre directement la ligne de
  base, comme si cette stimulation commençait à ~31 s (`B_TRAIN`, `B_PULSE`). Ce montage
  saute la stimulation 1 et ses 60 s de retour ; avant la stimulation 2, le ratio était
  encore à −1,1 %, et non à 0 comme le suggère le tracé : la réponse propre à la
  stimulation 2 est de −3,0 points (de −1,1 % à −4,3 %), pas ~−4,5 %. Le sous-titre
  mentionne « essai 2 : ligne de base, puis stimulation 2 ». Le retour est en plus
  décalé de +1,96 point, pour que le maximum de sa médiane lissée soit à 0 : la courbe
  finit sur la ligne de base, alors que le retour réel est incomplet dans la fenêtre. Le
  creux affiché (≈ −3 %) correspond alors à peu près à la réponse propre à la
  stimulation 2 (−3,0 points) ;
  La remontée est ajustée (courbe orange) par base + A·exp(−t/τ) sur la médiane non lissée,
  base fixée à la ligne de base affichée (0), comme pour les τ de D : τ = 23,8 ± 2,9 s. Sur
  cette fenêtre de 49 s, une base libre absorbe le retour (τ = 46 ± 33 s). La remontée
  affichée, lente au début puis plus rapide, n'est pas tout à fait exponentielle ;
- **C** (en bas à gauche) — premier train complet : quatre stimulations, quatre réponses ;
- **D** (en bas à droite) — temps de retour τ des 10 cellules qui répondent le mieux
  (réponse médiane la plus forte sur les 8 stimulations) : 4 barres, une par rang de
  stimulation, chacune réunissant les deux essais (stim 1 de `train_stim-002` et stim 1
  de `train_stim-003`, etc.) ; un point par cellule et par essai, les τ > 60 s rejetés
  (retour trop lent pour la fenêtre de 42 à 53 s : 12 sur 80), soit n = 15, 18, 17 et 18 ;
  barre = médiane ± écart interquartile. Valeurs dans `FOV11_train/results/tau_par_cellule.csv`.

**Cellules retenues** pour B et C : au moins 10 % plus brillantes que le neuropile qui les
entoure (`MIN_CONTRASTE`), 24/37. Le critère est validé par les données : les cellules
écartées répondent moins (VIS 50 %, FOV3–7 : −1,8 % contre −3,1 %, p = 0,03).

**τ.** Ajustement base + A·exp(−t/τ) du retour, de la fin d'une stimulation au début de
la suivante, la base étant fixée au niveau des 10 s qui précèdent la stimulation. Avec un
plateau libre, celui-ci absorbe le retour et τ part à la borne (300 s) pour un tiers des
cellules ; base fixée, seuls 2 τ sur 80 y restent. Aucun autre critère : chaque cellule a
un τ à chaque stimulation, même imparfait.

| | stim 1 | stim 2 | stim 3 | stim 4 |
|---|---|---|---|---|
| τ médian, deux essais réunis, τ ≤ 60 s | 35 s | 26 s | 22 s | 24 s |
| n | 15 | 18 | 17 | 18 |

Les deux essais n'ont pas les mêmes durées de stimulation : les barres réunissent des
stimulations de même rang, pas de même durée.
Le retour de la stim 1 de l'essai 2 est traversé par un saut instrumental de CFP (66 s),
qui fausse probablement son τ ; la stim 4 de l'essai 2, la plus longue, revient plus
lentement.

C2 = YFP et C3 = CFP : une baisse du ratio est une entrée de Cl⁻. Courbes lissées sur
5 images (moyenne glissante) ; ratio brut. Laissés de côté : le contrôle IR et la dose,
la dépendance au R basal (non concluante), la reproductibilité entre essais.

## Ce que disent les résultats

`dR/R0` sur les 10 premières secondes après la stim, médiane sur les ROIs. La fuite
est l'excédent de C2 pendant la stim, **en counts** et non en facteur : le gain du
détecteur et le zoom sont identiques partout (0,95 et 1,0, vérifiés dans les
métadonnées OME), donc les counts se comparent d'une acquisition à l'autre, alors
qu'un facteur se trouve gonflé quand la ligne de base est sombre.

| FOV | Acquisition | Type | Fuite C2 | Pente avant | **Réponse** |
|---|---|---|---|---|---|
| FOV1 | `stim-001` | ? | 11,1 | −0,30 | −0,6 % |
| FOV1 | `stim-002` | ? | 127 | −0,08 | −3,1 % |
| FOV2 | `stim_IR-001` | IR | 5,1 | −0,01 | **+0,3 %** |
| FOV2 | `stim_VIS-001` | VIS | 9,0 | −0,10 | −1,4 % |
| FOV2 | `stim_max-001` | ? | 155 | −0,06 | −5,0 % |
| FOV3 | `stim_IR-001` | IR | 23,6 | +0,58 | **+0,4 %** |
| FOV3 | `stim_VIS_10-001` | VIS 10 % | 0,46 | −0,02 | −3,5 % |
| FOV3 | `stim_VIS_50-001` | VIS 50 % | 2,2 | −0,04 | −5,9 % |
| FOV4 | `stim_IR-001` | IR | 23,3 | +0,69 | **+2,1 %** |
| FOV4 | `stim_VIS_10-001` | VIS 10 % | 1,3 | +0,56 | −1,3 % |
| FOV4 | `stim_VIS_50-001` | VIS 50 % | 6,3 | −0,18 | −2,9 % |
| FOV5 | `stim_VIS_1-001` | VIS 1 % | 0,08 | +0,13 | **+0,4 %** |
| FOV5 | `stim_VIS_10-001` | VIS 10 % | 0,48 | +0,17 | −1,8 % |
| FOV5 | `stim_VIS_50-001` | VIS 50 % | 2,5 | +0,01 | −4,7 % |
| FOV6 | `stim_VIS_1-001` | VIS 1 % | 0,28 | −0,67 | **−0,5 %** |
| FOV6 | `stim_VIS_10-001` | VIS 10 % | 1,2 | +0,65 | −1,5 % |
| FOV6 | `stim_VIS_50-001` | VIS 50 % | 6,1 | −0,28 | −1,7 % |
| FOV7 | `stim_VIS_10-001` | VIS 10 % | 0,40 | +0,24 | −1,1 % |
| FOV7 | `stim_VIS_50-001` | VIS 50 % | 1,9 | −0,03 | −1,6 % |
| FOV7 | `stim_VIS_100-001` | VIS 100 % | 3,9 | −0,04 | −1,8 % |

### La réponse est spectrale, pas photothermique

Dans une même tranche, au même gain, l'IR déverse sur le détecteur **18 à 51 fois plus
de lumière** que le VIS à 10 % — et ne produit rien :

| | IR | VIS 10 % |
|---|---|---|
| FOV3 | 23,6 counts → **+0,4 %** | 0,46 count → **−3,5 %** |
| FOV4 | 23,3 counts → **+2,1 %** | 1,3 count → **−1,3 %** |

Un photoblanchiment ou un échauffement croîtraient avec la lumière délivrée. Ici c'est
l'inverse, et d'un facteur 20 à 50. La dépendance est donc bien **spectrale**.

### La dose-réponse est réelle mais sature

Les ROIs sont communs aux acquisitions d'un même FOV : la comparaison se fait **par
paire, cellule par cellule**. Test des signes sur la différence par ROI, calculable
depuis `resume_par_roi.csv` :

| FOV | 1 % → 10 % | 10 % → 50 % | 50 % → 100 % |
|---|---|---|---|
| FOV3 | — | −4,0 pts, 33/41, p = 1·10⁻⁴ | — |
| FOV4 | — | −1,7 pts, 45/59, p = 7·10⁻⁵ | — |
| FOV5 | −2,4 pts, 51/74, p = 2·10⁻³ | −3,7 pts, 50/73, p = 2·10⁻³ | — |
| FOV6 | −0,8 pts, 30/37, p = 2·10⁻⁴ | −0,5 pts, 22/37, **n.s.** | — |
| FOV7 | — | −1,1 pts, 33/56, **n.s.** | −0,05 pt, 29/56, **n.s.** |

Le passage de 1 % à 10 % est net dans les deux tranches qui l'ont. Celui de 10 % à
50 % ne l'est que dans trois tranches sur cinq. Et sur FOV7, la seule qui monte à
100 %, **doubler la puissance au-delà de 50 % n'ajoute rien** (−0,05 pt).

La réponse sature donc, à une puissance qui dépend de la tranche : dès 10 % sur FOV6
et FOV7, plus tard sur FOV3 et FOV5. FOV6 et FOV7 sont aussi les tranches les moins
réactives (−1,7 et −1,8 % au maximum, contre −5,9 % sur FOV3).

Attention, le **réglage nominal n'est pas une mesure de la lumière délivrée** : à 50 %
nominal la fuite va de 1,9 à 6,3 counts selon la tranche, soit un facteur 3. Sur les
13 acquisitions VIS, le réglage nominal prédit d'ailleurs mieux la réponse
(r = −0,61 sur le log) que la fuite mesurée (r = −0,50) — les deux restent médiocres,
l'amplitude dépendant surtout de la tranche.

### La relaxation est exponentielle, mais incomplète

Les 50 meilleures réponses (45 neurones, 5 FOVs : surtout `FOV3/stim_VIS_50` et
`FOV5/stim_VIS_50`) partent d'environ −11 % juste après la stim et remontent comme une
exponentielle, **tau ≈ 12 s** (R² = 0,63). Mais la médiane plafonne à **−6,9 %** : plus
de la moitié de la réponse subsiste 60 s après la stim, et encore 90 s après sur les
acquisitions de 180 images. Réponse par réponse, seuls 13 ajustements sur 50 sont
exploitables (R² > 0,2), tau médian 24 s — un décours isolé est trop bruité.

Ce plateau n'est pas un artefact de dérive : les pentes des lignes de base avant stim de
ces acquisitions vont de −0,18 à +0,24 % par 10 s, de signes opposés, soit au plus
±1,4 % sur une minute — loin des −6,9 %. **La fenêtre
d'observation est probablement trop courte** pour voir la fin de la relaxation.

En haut, le ratio brut mélange des tranches de niveaux de base différents : les six
traits autour de 5–6 sont ceux de `FOV2/stim_max-001`.

### Le reste

Sur les conditions qui répondent, la bande de moyenne des heatmaps passe d'un blanc
uniforme avant la stim à un bleu franc juste après, qui s'éclaircit sans revenir au
blanc en 60 s. `trace_<FOV>_toutes.png` met les acquisitions d'un même FOV côte à côte
sur un axe y commun.

Images fixes : à 950 et 1064 nm, C3 est quasiment vide — un ratio C2/C3 y reviendrait
à diviser par du bruit, d'où l'analyse de C2 seul. Dans cinq FOVs sur sept (FOV3 à FOV7),
le C2 des ROIs est nettement plus fort à 950 qu'à 800 nm (rapport des médianes ×1,9 à
×2,8) ; il est à peu près égal sur FOV2 (×1,05) et un peu plus faible sur FOV1 (×0,85) — voir `images_fixes_C2.png`.

Trois réserves. Ces chiffres portent sur les ROIs bruts de Cellpose (voir étape 3),
dont certains tombent sur des cellules non marquées. Les acquisitions de FOV1 ne
portent pas d'étiquette de longueur d'onde, donc elles n'entrent pas dans la lecture
spectrale. Et il manque une acquisition témoin sans aucune stimulation, qui mesurerait
directement la dérive de base — les pentes avant-stim du tableau en sont un substitut,
pas un remplacement.
