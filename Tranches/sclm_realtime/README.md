# Pipeline SCLM — sessions FLIMApp en mode Realtime par passes (sans images)

Pour les sessions FLIMApp où la carte SPC-QC-104 range les photons **par passe** (un scan de la
NI, borné par les marqueurs M0/M3) et non par pixel : mode ROI (spirales des galvos de
FLIMApp, une ROI par passe) ou mode ROI éteint (une seule série). Le laser opto (1064 nm) suit
le protocole, en PWM : ses phases portent le code de routage réservé (0) et FLIMApp en jette
les photons. Ici, on les garde pour mesurer ce que le 1064 nm ajoute.

Pour le Realtime **en images** (`[clamp] images = true`, M1/M2 enregistrés, ROIs Cellpose), voir
`../sclm_flim/`.

Premier jeu : `../SCLM_realtime_2026-11-11/` (README de résultats).

## Format attendu

```
<jeu>/data/FLIM/        session FLIMApp, telle qu'écrite par l'appli
    run.toml            protocole, PWM, réglages SPC, mode ROI
    frames.csv          ce que FLIMApp a mesuré en direct (une ligne par passe et par ROI)
    visits.csv          commande de chaque créneau de la NI
    readback.bin/.txt   signaux relus par la NI (10 kHz) : galvos, commandes, lignes
    irf.csv             IRF de FLIMApp par canal, 256 canaux de 12,5/256 ns
    spc/<série>.spc     flux FIFO de chaque entrée (3T0089_IN1, _IN2) + _acquisition.ini
```

## Lancer

Depuis ce dossier, le jeu de données en argument :

```bash
python3 step1_decode.py ../<jeu>       # ~40 s : photons par passe, genre et canal de temps ; contrôle frames.csv
python3 step2_fuite.py ../<jeu>        # ce que le 1064 nm ajoute, ce qui en passe dans la mesure
python3 step3_series.py ../<jeu>       # ~1 min : lectures passe par passe, déclins de la session, dérive, bruit
python3 step4_reponses.py ../<jeu>     # réponses aux stims, témoins, modèles cinétiques, ordre de grandeur en mM
python3 step5_absolu.py ../<jeu>       # [Cl⁻] absolue avec la calibration en tranche du 2026-11-09 ?
python3 step6_cinetique.py ../<jeu>    # par ROI : latence, τ d'activation et τ de retour, par stim
python3 step7_hybride.py ../<jeu>      # [Cl⁻] dans le temps : niveau ancré sur τ, variations du rapport
python3 step8_images.py ../<jeu>       # ~1 min : images reconstruites sans M1/M2 (polygone)
~/.venvs/cellpose/bin/python step9_cellpose.py ../<jeu>   # ROIs Cellpose sur l'image reconstruite
python3 step10_rois.py ../<jeu>        # ~1 min : par ROI, [Cl⁻] ancrée, avant / après, cinétique
```

Les étapes 3 à 5 traitent la session comme **une seule série** (mode ROI éteint) ; l'étape 6
sépare les ROIs par code de routage. Pour une session en mode ROI, étendre les étapes 3 à 5 de
la même façon.

## Ce que fait chaque étape

**1 — décodage.** `spc.py` relit les `.spc` (format FIFO_150, traduit de la QC-104 par
FLIMApp) en deux lectures : les marqueurs (bornes des passes, M0 → M3 suivant), puis les
photons, rangés par blocs de mots. Chaque photon d'une passe va dans son genre (**ACQ** : code
≠ 0, acquisition ; **OPTO** : code 0, phase opto du PWM ; **PAUSE** : code 0 entre la fin de la
passe et la suivante) et son canal de temps (256 canaux de 12,5/256 ns, ceux de FLIMApp).
Égalités comme FLIMApp : un photon pile sur M0 est dans la passe, pile sur M3 il n'y est plus.
Les stims viennent de la commande 1 **relue** (`readback.bin`), créneau par créneau (créneau s
= passe s + 1). La colonne `command1` de `frames.csv` est décalée d'une passe au bord des blocs.
**Contrôles intégrés** : photons d'acquisition de chaque passe identiques à `frames.csv` ;
attribution des canaux (`config.CANAL`, vérifiée sur les temps moyens) ; galvos relus.

**2 — fuite du 1064 nm.** Photons des phases opto rapportés à ceux des pauses (850 nm seul, la
même passe) : stim contre les 40 passes de part et d'autre. Profil des photons dans la période
du PWM (1 µs) : ce qui passe dans les phases d'acquisition, retard d'allumage. Forme des photons
ajoutés : fond plat (non synchrone du 850 nm) contre forme du déclin.

**3 — séries.** Par passe, sur les photons d'acquisition : rapport R = YFP/CFP, temps moyen tm
(fenêtre [17, 255)), τ mono de chaque canal (`flim.py`), et la durée de vie de FLIMApp en direct.
Le modèle de `flim.py` est calculé sur une grille fine d'une **période exacte du laser**
(12,461 ns, pas les 12,5 ns de la fenêtre), puis intégré dans les canaux de FLIMApp ; décalage
de l'IRF ajusté une fois sur toute la session, puis fixé. Dérive sur les moments calmes (un
niveau libre par intervalle). Bruit autour d'une moyenne glissante de 41 passes, comparé au
bruit de photons (Monte-Carlo pour τ) ; écart d'Allan.

**4 — réponses.** Modèles sur toute la série, base quadratique + réponse (M1 pompe et
relaxation, M2 avec retard, M3 avec part permanente ; AIC), puis Δ = 10 s après − 10 s avant
sur les lectures redressées, contre des témoins (fausses stims dans la base et en fin
d'intervalle). Ordre de grandeur en mM avec les pentes de l'isotherme du 2026-11-09.

**5 — absolu.** Isothermes du 2026-11-09 sur le rapport et τ CFP, et forme du déclin CFP contre
les motifs du capteur seul à 0 et 131 mM de ce jour-là (étape 13 de `../sclm_calibration/`).

**6 — cinétique par ROI.** Par code de routage : ln R redressé, puis modèle marche-retour avec
latence (δ, τ d'activation, τ de retour) sur chaque stim et sur toutes ensemble (τ communs, IC
par bootstrap des stims). Avec une stim de 10 s, τ_act ne se mesure que s'il est nettement plus
court que la stim.

**7 — lecture hybride.** Ancres dans les fenêtres calmes : τ CFP moyen → c_τ (isotherme de τ
mono du 2026-11-09), échelle du rapport η = R / R_iso(c_τ) ; chaque passe lue au rapport,
R_iso⁻¹(R / η), η interpolé entre les ancres (ancrage continu) ou fixé à la première
(ancrage unique). Diagnostic : Δτ contre Δ ln R entre ancres, comparés à la direction du Cl⁻
(isothermes) et à celle d'une perte d'accepteur (−τ_D·E). Réponses aux stims en mM sur ln R
redressé (base de M1).

**8 — images sans M1/M2.** Quand le polygone balaie le champ pendant une session par passes,
la phase de ses lignes se retrouve dans les photons : sur une trame (540 lignes), leur
repliement modulo la période de ligne donne la projection de l'image sur x ; chaque trame est
recalée sur un gabarit de cette projection (corrélation croisée circulaire, sous-canal), les
décalages déroulés donnent la phase, donc la ligne et la position de chaque photon. Image
roulée pour mettre le retour du galvanomètre en bas ; blocs de 10 s recalés sur l'image de la
session (dérive de la tranche ou erreur lente du suivi).

**9 — ROIs** (venv Cellpose) : Cellpose-SAM sur √photons, bandes du retour écartées.

**10 — par ROI.** Photons d'acquisition rangés par passe et par ROI (position de l'étape 8,
corrigée du décalage du bloc) ; contrôle : la somme redonne l'étape 1. [Cl⁻] de chaque ROI par
ancrage unique sur son τ CFP ; avant / après (Wilcoxon entre ROIs) ; modèle de l'étape 6 par ROI ;
`hybride_ancrage_unique_rois.png`, la figure de l'ancrage unique de l'étape 7 refaite par soma.
Chaque fenêtre de temps est lue sur ses photons sommés, puis inversée par l'isotherme. Ne pas
faire la moyenne des lectures passe par passe : à 100–400 photons CFP par passe et par soma,
l'inversion convexe de l'isotherme la biaise vers le haut (+0,06 à +0,47 mM).

## Canaux

`config.CANAL = {"yfp": 1, "cfp": 2}` : le 2026-11-11, IN1 = YFP et IN2 = CFP (comme le
2026-11-09, l'inverse d'octobre). FLIMApp suit donc la durée de vie du **YFP** sur ch1. Vérifier
cette attribution à chaque nouvelle session (l'étape 1 s'arrête si le CFP n'a pas le temps moyen
le plus court).

## Pièges de ce mode

- **M1 et M2 ne sont pas enregistrés** en Realtime par passes : si le polygone balaie pendant la
  session, les photons viennent de tout le champ. L'étape 8 reconstruit les images en suivant
  la phase des lignes trame par trame (elle dérive de dizaines de µs/s : un modèle à période
  fixe ne suffit pas, le fondamental à 18 kHz non plus). Le Realtime en images reste plus sûr.
- Le laser à 850 nm ne s'arrête pas pendant les pauses (pas d'obturateur) : leurs photons
  servent de référence sans 1064 nm.
- La dernière passe n'a pas de pause complète.
- `concentration_ch*` de `frames.csv` est une conversion fictive de FLIMApp.

## Fichiers

| | |
|---|---|
| `config.py` | paramètres, lecture d'une session (run.toml, readback, frames, visits, IRF) |
| `spc.py` | décodage des `.spc` en passes (genres, canaux de temps, profils du PWM, codes) |
| `flim.py` | durées de vie sur la période exacte du laser, temps moyen |
| `figures.py` | style commun des figures |
| `step1` … `step10` | les étapes (la 9 avec le venv de Cellpose) |
