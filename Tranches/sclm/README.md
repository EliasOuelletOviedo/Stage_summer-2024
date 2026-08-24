# Analyse SCLM — tranches 2026-08-17, ratio C3/C2

Pipeline en 4 étapes sur `/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-08-17_tranches/SCLM`.
Le canal **C1 est ignoré** ; seuls C2 et C3 sont utilisés.

Environnement : `~/.venvs/cellpose/bin/python` (Cellpose-SAM 4.2.1, torch 2.13, MPS actif).
Toutes les commandes se lancent depuis ce dossier (`Tranches/sclm/`).

```bash
~/.venvs/cellpose/bin/python step1_sum.py       # recalage + 15 images somme
~/.venvs/cellpose/bin/python step2_cellpose.py  # segmentation des somas
~/.venvs/cellpose/bin/python step3_select.py    # tri manuel des ROIs (interactif)
~/.venvs/cellpose/bin/python step4_ratio.py     # ratios C3/C2 + boxplots
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

`ratio = moyenne(C3) / moyenne(C2)` sur les pixels du ROI, **valeurs brutes**
(aucune soustraction de fond), calculé indépendamment dans chaque condition.

Un ROI dont moins de 90 % des pixels restent dans le champ après recalage est mis à
`NaN` pour cette condition (colonne `note` du CSV).

Sorties dans `../SCLM_2026-08-17/results/` :
`ratios_C3_C2.csv` (une ligne par ROI × condition), `resume_stats.csv`,
`boxplot_Slice1.png`, `boxplot_Slice2.png`, `boxplot_Slice3.png`, `boxplot_combine.png`.
