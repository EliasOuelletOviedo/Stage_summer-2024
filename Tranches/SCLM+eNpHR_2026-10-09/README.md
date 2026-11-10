# SCLM + eNpHR — tranche du 2026-10-09, cinq sessions FLIMApp

Les sessions gardées sont listées dans `tri.csv` (`../sclm_flim/step0_tri.py`) : une par
sous-dossier, chacune analysée par le pipeline commun `../sclm_flim/` (étapes 1 à 7 ; résultats dans
`<session>/results/`). `FOV2_1158_4x10s_100pct` n'a pas de ligne de base (1re stim à 0,1 s) :
l'étape 5 n'y écrit que les valeurs brutes.

Ce fichier ne porte que sur l'**étape 8**, le filtre bayésien hybride (ratio + FLIM) de
`Simulateurs/bayesian_filter`, appliqué aux somas de chaque session :

```bash
cd ../sclm_flim
for s in 1119_3x10s_100pct FOV2_1143_2x10s_25pct FOV2_1150_2x10s_50pct FOV2_1158_4x10s_100pct FOV3_1334_4x20s_50pct; do
    python3 step8_filtre_bayesien.py ../SCLM+eNpHR_2026-10-09/$s      # 1 à 20 min chacune (FOV2 : 62 somas)
done
cd ../SCLM+eNpHR_2026-10-09                                           # étape 9 : toutes les sessions d'un coup
python3 ../sclm_flim/step9_redescente.py 1119_3x10s_100pct FOV2_1143_2x10s_25pct FOV2_1150_2x10s_50pct FOV2_1158_4x10s_100pct FOV3_1334_4x20s_50pct
```

Méthode et réserve générale dans le README de `../sclm_flim/` : trames de ~3 s (blocs FLIMApp
d'1 s regroupés par 3, comme à l'étape 5), calibration de la titration SCLM à 800 nm, transférée
à la session par un ajustement conjoint des somas (IRF de la QC-104, décalages d'horloge, fond et
YFP de eNpHR par soma, durées de vie réajustées).

## Le transfert de la calibration

| session | somas | photons bleus / soma / 3 s | déviance/ddl, τ de la solution | déviance/ddl, τ réajustés | τD réajusté (ns) | η | c médian (mM), τ solution → réajustés |
|---|---|---|---|---|---|---|---|
| 1119, 3 × 10 s, 100 % | 28 | 875 | 41,5 | 7,6 | 1,75 | 0,83 | 11 → 15 |
| FOV2 1143, 2 × 10 s, 25 % | 62 | 424 | 27,2 | 4,9 | 1,50 | 1,33 | 1,4 → 63 |
| FOV2 1150, 2 × 10 s, 50 % | 62 | 438 | 23,1 | 4,6 | 1,55 | 1,37 | 3,4 → 71 |
| FOV2 1158, 4 × 10 s, 100 % | 62 | 459 | 2,5 | 1,3 | 1,61 | 1,35 | 5,4 → 65 |
| FOV3 1334, 4 × 20 s, 50 % | 19 | 445 | 25,5 | 3,5 | 1,50 | 0,95 | 1,3 → 14 |

(La déviance/ddl baisse avec les photons : 78 s seulement pour FOV2 1158. En solution, τD = 2,72 ns.)

1. **Les durées de vie de la solution ne décrivent pas les somas**, et celles qu'on réajuste
   s'en éloignent (τD 1,5 à 1,75 ns contre 2,72) : elles absorbent ce que le modèle ne contient
   pas. Le fond retenu vaut 11 à 16 % des photons bleus et 7 à 9 % des jaunes, la YFP directe
   (eNpHR) est retenue pour la moitié des somas environ.
2. **La [Cl⁻] absolue dépend du modèle** : 1,4 ou 63 mM pour la même session FOV2 1143 ; une
   partie des somas touche le bord de la grille (300 mM). Elle n'est pas à lire.
3. **Le transfert est reproductible sur un même champ** : les trois sessions de FOV2 redonnent
   les mêmes τ (τD 1,50–1,61 ns, τA 2,79–2,84 ns), le même η (1,33–1,37) et le même c médian
   (63–71 mM). Le modèle se trompe de la même façon d'une session à l'autre.

## Les réponses

Δc lissé = c moyen 10 s après la stim − 10 s avant, médiane des somas (étape 5 pour les
fenêtres et les témoins ; `<session>/results/filtre_resume.csv`, `filtre_reponses.png`) :

| session | stims : Δc (mM), somas en hausse | témoins : Δc (mM) | par soma, stims contre témoins | ρ avec ΔR/R0 (étape 5) |
|---|---|---|---|---|
| 1119 | +0,05, +0,08, +0,16 (19 à 24 / 28) | +0,01 à +0,03 | +0,14 contre +0,04 (p = 0,15) | −0,56 |
| FOV2 1143 (25 %) | +0,14, +0,00 (46 et 32 / 62) | +0,15, −0,09 | +0,09 contre +0,03 (p = 0,05) | −0,46 |
| FOV2 1150 (50 %) | +0,14, +0,08 (44 et 39 / 62) | +0,01, +0,03 | +0,11 contre −0,00 (p = 3·10⁻⁶) | −0,34 |
| FOV2 1158 (100 %) | —, +0,01, −0,02, −0,17 | (aucun : stims à 20 s d'intervalle) | — | — |
| FOV3 1334 (50 %, 20 s) | +0,37, +0,12, +0,21, +0,06 (17, 13, 16, 10 / 19) | +0,01 à +0,10 | +0,24 contre +0,06 (p = 0,03) | −0,59 |

- Le sens est celui du ratio : c monte après les stims, un peu plus que pendant les témoins, et
  Δc suit ΔR/R0 soma par soma (ρ de −0,34 à −0,59). Les réponses sont 2 à 4 fois plus petites
  qu'au 2026-10-08 (+0,3 à +0,6 mM), avec 7 à 15 fois moins de photons par soma.
- En temps réel, le filtre ne sépare plus les stims des témoins (p de 0,12 à 0,71) : à
  ~450 photons bleus par soma et par tranche, c'est le lissage qui fait la mesure.
- FOV2 1158 : sans ligne de base ni intervalle, ses fenêtres « avant » tombent dans la réponse
  à la stim précédente ; ses Δc ne se comparent pas aux autres.
- η est appris sans surprise (0,83 à 1,37 selon le champ, stable sur la session), mais son
  incertitude (± 0,15 à 0,7 %) est trop étroite pour un modèle qui ajuste mal.

## La redescente de la [Cl⁻] (étape 9) : pas lisible à ces photons

`filtre_redescente_sessions.png` (ici) met côte à côte la moyenne des stims de chaque
session ; `<session>/results/filtre_redescente.png` et `filtre_redescente_somas.png`
détaillent chaque stim et chaque soma. FOV2 1158 n'y est pas : 10 s entre deux stims,
aucune redescente à ajuster.

| session | dérive de c lissé (%/100 s) | forme après la stim (médiane des somas) | τ des somas mesurables (s) | τ du ratio (étape 7, s) |
|---|---|---|---|---|
| 1119 | +2,5 | monte jusqu'à +0,3 mM vers 80 s, revient vers 140 s | 25 (11 / 28 somas) | 8 (3 à 41) |
| FOV2 1143 (25 %) | +0,3 | ~+0,1 mM pendant 90 s, puis sous la base | 29 (26 / 62) | 63 (non mesurable à la stim 1) |
| FOV2 1150 (50 %) | +0,2 | monte jusqu'à +0,5 mM vers 85 s, redescend à moitié ensuite | 23 (16 / 62) | non mesurable |
| FOV3 1334 (50 %, 20 s) | +2,3 | saut de +0,1 mM, plateau sur 150 s | 21 (9 / 19) | 3 (1 à 190) |

- **Aucune redescente n'est mesurable sur la médiane des somas**, ni par stim ni sur la
  moyenne des stims : après la stim, c reste au-dessus de sa base, ou continue à monter
  pendant une minute ou plus, au lieu de redescendre en ~15 s comme au 2026-10-08.
- Le ratio de l'étape 7 n'y voit pas mieux (τ de 1 à 190 s selon la stim, non mesurable
  sur FOV2) : à ~450 photons bleus par soma et par tranche, la réponse est trop petite
  pour en lire la forme, et le lissage du filtre peut l'étaler sur plusieurs dizaines de
  secondes. Les τ des somas « mesurables » (21 à 29 s) ne portent que sur une partie
  d'entre eux, choisie par l'ajustement : ils ne disent pas le temps de redescente du champ.

## Sorties (dans chaque `<session>/`)

| Fichier | Contenu |
|---|---|
| `qc/filtre_transfert.png`, `.json`, `_somas.csv` | la calibration en solution sur la somme des somas, avec ses durées de vie et avec celles réajustées ; paramètres partagés ; c, YFP directe, fonds et déviance par soma |
| `results/filtre_bayesien.png` | [Cl⁻] de chaque soma, lissée et en temps réel, et η du champ |
| `results/filtre_reponses.png`, `filtre_reponses.csv`, `filtre_resume.csv` | Δc par soma × stim et témoin, résumé par stim, Δc contre ΔR/R0 |
| `results/filtre_par_trame.csv` | c (et intervalle à 95 %) par soma et par tranche, η (non versionné) |
| `results/filtre_redescente.png`, `filtre_redescente_somas.png`, `filtre_redescente.csv` | étape 9 : Δc aligné sur la fin de chaque stim, C + A·exp(−t/τ), carte des somas, τ contre celui du ratio |
| `../filtre_redescente_sessions.png` | étape 9 : la moyenne des stims de chaque session, côte à côte |
