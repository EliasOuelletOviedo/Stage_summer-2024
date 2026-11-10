# Analyse SCLM — tranches 2026-08-17, ratio C2/C3

Pipeline en 4 étapes sur `/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-08-17_tranches/SCLM`.
Le canal **C1 est ignoré** ; seuls C2 et C3 sont utilisés.

**Les étapes doivent être lancées dans l'ordre** : chacune consomme la sortie de la
précédente. `step4` refuse de tourner tant que le tri manuel de `step3` n'est pas fait.

Seule l'étape 2 a besoin de Cellpose, donc de l'environnement dédié
`~/.venvs/cellpose/bin/python` (Cellpose-SAM 4.2.1, torch 2.13, MPS actif).
Les étapes 1, 3 et 4 n'utilisent que numpy / scipy / matplotlib / pandas / tifffile et
tournent aussi avec le `python3` du système.

Toutes les commandes se lancent depuis ce dossier (`Tranches/sclm/`).

```bash
python3 step1_sum.py                            # recalage + 15 images somme
~/.venvs/cellpose/bin/python step2_cellpose.py  # segmentation des somas  (Cellpose requis)
python3 step3_select.py                         # tri manuel des ROIs (interactif)
python3 step4_ratio.py                          # ratios C2/C3 + boxplots
```

Les sorties vont dans `../SCLM_2026-08-17/`.

## Appariement des FOVs

Vérifié par corrélation croisée invariante en translation (`config.py: FOV_MAP`).

| Tranche | Basal / High_K | Return |
|---|---|---|
| Slice 1 | `FOV_1-001` … `FOV_5-001` | `FOV-001` … `FOV-005` (diagonale) |
| Slice 2 | `FOV-001` … `FOV-005` | `FOV-001` … `FOV-005` (diagonale) |
| Slice 3 | `FOV-001` … `FOV-005` | 1→001, 2→002, 3→003, **4→aucun**, 5→**004** |

- `Slice 2/High_K/FOV-006` est un timelapse de 1764 images, pas un FOV fixe : **exclu**.
- Les dossiers `Basal_to_high_K*` et `High_K_to_basal*` sont des enregistrements de
  transition (~1764 images chacun) : **non utilisés** ici.
- Slice 3 / Basal `FOV-004` n'a **pas** de Return : son image somme n'additionne que
  Basal + High_K (4 images au lieu de 6) et la condition Return y vaut `NaN`.

## Étape 1 — recalage et somme

Chaque condition est recalée sur **Basal** par corrélation de phase sur `C2+C3`,
translation entière uniquement (pas d'interpolation, les intensités sont préservées).
Garde-fou : si le score NCC est `< 0.15`, le pic de corrélation n'est pas fiable et
**aucun décalage n'est appliqué** (statut `NON RECALE` dans `qc/registration_log.tsv`) —
un décalage arbitraire de 100 px serait bien pire que pas de recalage.

L'image somme est l'addition des 6 images recalées (2 canaux × 3 conditions), les
bords entrants comptant pour 0.

## Étape 2 — Cellpose

Modèle **Cellpose-SAM** (`cpsam`, le modèle par défaut de Cellpose 4), réglé pour des
somas de neurones en tranche :

| Paramètre | Valeur | Raison |
|---|---|---|
| `diameter` | 30 px | soma ~15 µm à 0,509 µm/px (métadonnée OME) |
| `flow_threshold` | 0.4 | défaut, masques bien formés |
| `cellprob_threshold` | 0.0 | défaut ; plus bas ferait fusionner les somas serrés |
| `min_size` | 120 px | écarte débris et taches de neuropile |
| `max_size_fraction` | 0.25 | empêche un ROI d'avaler le champ |

`model.eval()` normalise son entrée **en place** : on lui passe toujours une copie.

## Étape 3 — tri manuel

| Touche | Effet |
|---|---|
| clic sur un ROI | garder / rejeter |
| `a` / `n` | tout garder / tout rejeter |
| `b` | change le fond (somme, Basal C2, High_K C2, Return C2) |
| `l` | affiche/masque les numéros |
| `←` `→` | FOV précédent / suivant (enregistre au passage) |
| `q` | enregistre et quitte |

Tous les ROIs sont **gardés par défaut** : il suffit de cliquer ceux à rejeter.
Vert plein = gardé, rouge pointillé = rejeté. La sélection est enregistrée dans
`selection/<FOV>_selection.json` et le script peut être relancé pour reprendre où on
en était, ou sur un seul FOV : `step3_select.py Slice2_FOV4`.

## Étape 4 — ratio et figures

`ratio = moyenne(C2) / moyenne(C3)` sur les pixels du ROI, **valeurs brutes**
(aucune soustraction de fond), calculé indépendamment dans chaque condition.

Un ROI dont moins de 90 % des pixels restent dans le champ après recalage est mis à
`NaN` pour cette condition (colonne `note` du CSV).

Sorties dans `../SCLM_2026-08-17/results/` :

| Fichier | Contenu |
|---|---|
| `ratios_C2_C3.csv` | une ligne par ROI × condition (`mean_C2`, `mean_C3`, `ratio`, `note`) |
| `resume_stats.csv` | n, moyenne, médiane, écart-type, Q1, Q3 par tranche × condition |
| `boxplot_Slice{1,2,3}.png` | distribution du ratio par condition, une figure par tranche |
| `boxplot_combine.png` | les 3 tranches groupées sur un même axe |
| `trace_Slice{1,2,3}.png` | **trajectoire de chaque ROI** : un trait mince relie les 3 conditions |
| `trace_combine.png` | les 3 tranches en petits multiples, axe y commun |

Les deux familles de figures répondent à des questions différentes : le boxplot montre
la distribution à chaque condition, le tracé apparié montre si les ROIs répondent tous
dans le même sens ou si la médiane masque des comportements opposés. Un ROI dont une
condition manque voit son trait s'interrompre au lieu d'être écarté.
