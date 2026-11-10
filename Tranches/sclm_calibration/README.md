# Pipeline SCLM — calibration en tranche, [Cl⁻] imposée par ionophores (SPCM + Bliq)

Pour les expériences où la [Cl⁻] intracellulaire est **imposée** dans la tranche
(ionophores, solutions à 0, 20, 131 mM de Cl⁻) après des images en conditions basales :
on mesure, soma par soma, le rapport YFP/CFP et la durée de vie de la CFP à chaque
[Cl⁻] connue, on en tire l'isotherme du capteur **dans le tissu**, et on relit les somas du
basal pour estimer leur [Cl⁻]i. Les time-lapses pris aux changements de solution donnent la
vitesse de mise à l'équilibre.

Pourquoi un pipeline à part : les données viennent du logiciel **SPCM** de B&H (et non de
FLIMApp), une acquisition par FOV, en cycles de 10 s ; et la question est une calibration,
pas une réponse à une stimulation. Même démarche que `../sclm_flim/` pour le reste : ROIs
Cellpose triés, rapport en valeurs brutes, durées de vie par vraisemblance de Poisson.
Le code est ici, une seule fois ; chaque jeu garde ses données et ses résultats dans son
propre dossier (premier jeu : `../SCLM_calibration_2026-11-09/`).

## Format attendu

```
<jeu>/data/FLIM/<condition>/fov<k>_c<n>.spc        un cycle de 10 s par fichier (SPCM, FIFO_IMAGE)
<jeu>/data/FLIM/<condition>/fov<k>_c<n>.sdt|.set   images du dernier cycle, réglages (contrôle)
<jeu>/data/FLIM/<a>_to_<b>[_<partie>]/…_c<nn>.spc  time-lapse d'un changement de solution
<jeu>/data/Bliq/<condition>/fov<k>-001/Bliq VMS/C1|C2|C3/*.tif   images Bliq simultanées (10 s)
```

Conditions : `basal`, puis `<x> mM`. Une transition coupée en plusieurs dossiers
(`basal_to_0_1`, `_2`, `_3`) est mise bout à bout. `data/` peut contenir des liens.

## Lancer

Depuis ce dossier, le jeu de données en argument :

```bash
python3 step1_decode.py ../<jeu>                          # ~6 min (9 Go) : photons, dérive, plan focal, contrôles
~/.venvs/cellpose/bin/python step2_cellpose.py ../<jeu>   # ~1 min : somas de chaque acquisition
python3 step3_select.py ../<jeu>                          # tri manuel (fenêtre, ← → entre acquisitions) ; --tout = tout garder
                                                          # --top-photons 50 = les 50 % de ROIs les plus brillants de chaque acquisition
python3 step4_declins.py ../<jeu>                         # ~4 min : déclins par ROI, cycle et canal
python3 step5_irf.py ../<jeu>                             # ~30 s : contrôles de l'IRF, τ communes du déclin CFP
python3 step6_somas.py ../<jeu>                           # ~30 s : mesures par soma dans les FOV fixes
python3 step7_calibration.py ../<jeu>                     # ~1 min : isothermes, [Cl⁻]i basale, auto-cohérence
python3 step8_cinetiques.py ../<jeu>                      # ~30 s : cinétiques des changements de solution
python3 step9_flim_deux_canaux.py ../<jeu>                # ~2 min : modèle FLIM à deux canaux de bayesian_filter, soma par soma
python3 step9_flim_deux_canaux.py ../<jeu> --transfert 1  # la même avec un seul temps de transfert (sorties *_1t)
python3 step10_validation_croisee.py ../<jeu>             # ~10 min : validation croisée, ROIs tirés au hasard
python3 step10_validation_croisee.py ../<jeu> --par-fov   # ~10 min : validation croisée, un FOV de chaque condition par pli
```

Seule l'étape 2 a besoin de Cellpose. Après un tri manuel (étape 3), relancer les étapes 5
à 8 : l'étape 4 mesure tous les ROIs une fois pour toutes. `step1_decode.py --sans-decodage`
refait dérive, plan focal et figures sans relire les `.spc`.

## Ce que fait chaque étape

**1 — décodage** (`spcm.py`). Format FIFO de la QC-104 établi par FLIMApp
(`format_fifo_qc104.jl`) : débordement = bit 31 seul, marqueur = bit 30 (M1 ligne, M2
trame), photon = bits 30-31 nuls, entrée aux bits 28-29 (0 = IN1). Ligne au M1, trame au
M2 ; pixel de 51,82 ns (celui de Bliq) ; l'image Bliq est le rectangle des lignes 14–525 et
des pixels 28–1051, vérifié acquisition par acquisition (meilleur décalage 0, 0). Le canal
de 4 ps du TDC est tiré uniformément dans sa largeur puis rangé modulo la période en 256
canaux (sans ce tirage : peigne de ±4 %). **Contrôles** : le déclin du dernier cycle de
chaque FOV colle à celui du `.sdt` écrit par SPCM (`qc/validation_sdt.csv`) ; ressemblance
des photons à l'image Bliq simultanée (`qc/recalage_bliq.csv`). **Dérive** : décalage de
chaque cycle sur la somme (corrélation de phase), appliqué aux masques à l'étape 4.
**Plan focal** : un cycle dont l'image lissée ressemble à moins de 0,95 (NCC) aux autres
cycles (FOV) ou à ses voisins (transition) est d'un autre plan — retouche de la mise au
point — et écarté. Un cycle interrompu de moins de 3 s est écarté.

**2 — Cellpose** (`cpsam`) sur la racine de la somme des photons CFP + YFP de chaque
acquisition, diamètre 30 px (~15 µm), réglages des tranches 2026-08-17.

**3 — tri manuel**, acquisition par acquisition (clic = garder/rejeter, `b` = fond, `←`/`→`) ; ou sans fenêtre :
`--tout`, ou `--top-photons P` (les P % de ROIs au plus grand total de photons CFP + YFP, dans chaque acquisition).

**4 — déclins.** Second passage sur les photons : un déclin de 256 canaux par ROI, par
cycle et par canal, le neuropile (hors somas dilatés de 10 px, au-dessus du 25e centile),
et les photons de chaque ROI par seconde de cycle.

**5 — IRF** (`flim.py`). L'IRF est celle mesurée par SPCM pour le jeu (`data/FLIM/irf*.sdt`),
utilisée telle quelle : fond plat retiré (médiane avant le pic), traîne gardée ; ligne du
`.sdt` de chaque canal dans `config.IRF_LIGNE`. Contrôles par FOV (déclin CFP des somas réunis,
3 exponentielles) : décalage (≈ 0 attendu, l'IRF étant dans le repère des données),
élargissement voulu par les données, déviance des somas. Puis trois durées de vie communes à
tous les FOV (populations du déclin CFP, comme le modèle « 1-2-3 » de
`Simulateurs/bayesian_filter`), amplitudes, fond et décalage par FOV.

**6 — mesures par soma** (FOV fixes, cycles gardés réunis) : photons, rapport YFP/CFP
**brut** et corrigé du **temps mort** (`mesures.py` : au plus un photon par impulsion et par
entrée sur la QC-104 ; pixel par pixel, n → −E·ln(1 − n/E)), taux par impulsion μ, τ CFP
mono-exponentielle (décalage libre, le modèle de FLIMApp), temps moyen τm, amplitudes des
trois populations et leur τ moyen ; pente du rapport sur les cycles.

**7 — calibration.** Isotherme à un site X(c) = X∞ + (X₀ − X∞)·K/(K + c) pour chaque
lecture, sur les médianes des somas des trois [Cl⁻] imposées (exacte : 3 points, 3
paramètres) ; incertitudes par **bootstrap hiérarchique** (FOV tirés avec remise, puis somas
dans chaque FOV) ; inversion c = K·(X₀ − X)/(X − X∞) pour le basal (médiane et cellule par
cellule). **Auto-cohérence** : les somas des conditions imposées relus un par un par la
courbe — la dispersion de ces [Cl⁻] relues est la précision d'une lecture sur une cellule.

**8 — cinétiques** des transitions : par cycle (somas réunis, chaque soma) et par seconde,
X(t) = X∞ + (X₀ − X∞)·exp(−t/τ_éq) avec t = 0 au changement de solution ; deux
exponentielles pour basal→0. Fraction de la transition accomplie quand les FOV de la
condition d'arrivée ont commencé.

**9 — FLIM à deux canaux** (`deux_canaux.py`, modèle de `Simulateurs/bayesian_filter`, fonctions de son
`sclm_hybride.py`). IRF paramétrique (gaussienne ⊗ traînes exponentielles + fond) ajustée sur l'IRF mesurée ;
calibration par ajustement de Poisson conjoint des déclins bleus et jaunes réunis de tous les somas de chaque
[Cl⁻] imposée (τ1, τ2, τD, τA, λ, G, g2, κ, k0 partagés ; A, f1, f2, fonds, décalages par condition ; η = 1),
puis isothermes de f1 et f2. Chaque soma est lu par maximum de vraisemblance sur une grille de c, **sans filtre
temporel**, avec un fond libre par canal : pour chaque c, amplitude et fond de chaque canal (Newton, Poisson)
donnent la vraisemblance des **formes** ; le rapport des amplitudes, comparé à η·R₁(c), donne la **quantité**
(le découpage formes + totaux de la section 11). Quatre lectures, comme `estimateurs.py` : temps de vie unique
(→ τ(c) des déclins de calibration), FLIM du système (η libre par soma), ratio (η = 1), **hybride** (η partagé
par le FOV, profilé sur ses somas, comme θ dans le filtre). Comparées au τ mono-exponentiel de l'étape 6.

**10 — validation croisée.** Cinq plis, de deux façons : ROIs tirés au hasard dans chaque condition (par
défaut : nouvelles cellules des mêmes champs), ou `--par-fov` : le pli k retire le k-ième FOV de chaque [Cl⁻] imposée, recalibre les
isothermes de l'étape 7 et le modèle de l'étape 9 (2 et 1 population) sur les autres FOV, et lit les FOV retirés ;
le basal est lu par chaque pli. Critère : erreur de la médiane de chaque FOV lu hors calibration par rapport à la
[Cl⁻] imposée (erreur quadratique moyenne par condition), dispersion entre somas, part hors courbe.

## Canaux (2026-11-09)

| Entrée QC-104 | Image Bliq | Fluorophore | Établi par |
|---|---|---|---|
| IN1 (champ 0) | C1 (NCC 0,96) | **YFP** | 4–5 fois plus de photons ; déclin plus lent ; IN1/IN2 baisse avec le Cl⁻ |
| IN2 (champ 1) | C2 (NCC 0,89) | **CFP** | déclin qui s'allonge avec le Cl⁻ (le donneur perd son FRET) |
| — | C3 (gain 0,2) | — | rayures du polygone, pas de cellules : ignoré |

C'est l'inverse des sessions FLIMApp du 2026-10-08/09 (IN1 = CFP) : vérifier les canaux à
chaque nouveau jeu (`config.ENTREE`).

## Fichiers

| | |
|---|---|
| `config.py` | constantes, lecture d'un jeu (`acquisitions`), IRF du jeu (`charger_irf`) |
| `spcm.py` | décodage des `.spc` de SPCM (QC-104) et rangement des photons |
| `recalage.py` | corrélation de phase, NCC |
| `flim.py` | modèles de déclin (mono, n exponentielles, amplitudes à τ fixées), temps moyen |
| `mesures.py` | comptes corrigés du temps mort, lectures de la durée de vie CFP |
| `deux_canaux.py` | modèle FLIM à deux canaux de bayesian_filter : IRF paramétrique, calibration conjointe, lecture d'un soma |
| `figures.py` | style commun des figures |
| `step1` … `step10` | les étapes |
