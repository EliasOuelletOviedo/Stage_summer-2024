# Suivre une dérive de η à très bas photons

Le test : **10 ROI, 1 000 photons par ROI et par trame, 3 000 trames.**
- η, le facteur de transmission jaune/bleu, dérive linéairement de 1 à 1,1 ; le filtre ne le connaît pas.
- L'amplitude perd 30 % par photoblanchiment.
- Chaque ROI a sa propre réponse en [Cl⁻] (montée, plateau, retour).

Calibration : SCLM_mCer3 (`params.toml`). Les données sont les mêmes pour toutes les variantes (même graine).

## Les variantes (`essai.py`, `moteur.py`)

| | dynamique de η | η partagé par les ROI | dynamique de c |
|---|---|---|---|
| référence | marche aléatoire, moyenne de modèles sur s_η | non | s = 4 % par trame |
| 1 | dérive paramétrique ln η_k = ln η0 + D·k/(T−1), θ = (η0, D) statique | non | s = 4 % |
| 2 | marche aléatoire, moyenne sur s_η | oui (projection p(η) Π p(c_m \| η)) | s = 4 % |
| 1-2 | dérive paramétrique | oui (factorisation exacte) | s = 4 % |
| 1-2-3 | dérive paramétrique | oui | moyenne de modèles sur s_c ∈ {0,4 ; 1,3 ; 4} % |

Grille de la dérive : η0 de 0,7 à 1,5 (81 valeurs) et D de −0,15 à +0,15 (31 valeurs).

## Résultats (`sorties/resultats_essai_3000.json`)

| | rms de c (mM) | biais (mM) | couverture à 95 % | rms de η | η final (vrai 1,1) |
|---|---|---|---|---|---|
| référence | 4,06 | −3,37 | 64 % | 6,7 % | 1,012 |
| 1 | 3,58 | −2,88 | 75 % | 5,8 % | 1,034 |
| 2 | 3,67 | −3,37 | 52 % | 5,9 % | 1,015 |
| 1-2 | 3,01 | −2,67 | 70 % | 4,7 % | 1,037 |
| **1-2-3** | **1,67** | **−1,13** | **94 %** | **2,1 %** | **1,075 ± 0,8 %** |

Pour la référence et la variante 1, l'erreur rms est la moyenne des rms par ROI. Les figures, elles, donnent la rms sur les 10 ROI ensemble : 4,09 et 3,62 mM.

**L'ingrédient décisif est 3.** La moyenne de modèles choisit s_c = 0,4 % par trame (poids 1,00).

- **Pourquoi :** avec s = 4 %, le c de chaque trame est presque libre. Un η partagé, estimé à travers des milliers de c presque libres, reste biaisé quel que soit le nombre de trames (problème des paramètres incidents).
- **Ce qui reste :** même avec 1-2-3, η a environ 2 % de retard sur la vérité, d'où un biais de −1,1 mM sur c.

Figures (`sorties/figures/`, style de `fig_exemple_eta`) :
- **`fig_derive_<v>_roi1.png` :** la ROI 1. Le panneau a montre la durée de vie trame par trame, le ratio à η = 1, le filtre temps réel avec son intervalle à 95 % et la vraie valeur. Le panneau b montre η ± 2 sd et les poids des modèles.
- **`fig_derive_<v>_10roi.png` :** les 10 ROI, avec η en bas.
- **`fig_derive_comparaison.png` :** η et l'erreur sur c des 5 variantes.

## Toutes les modifications ensemble (`tout.py`)

« tout » ajoute à 1-2-3 :
- **[1] ancre sur la forme bleue seule :** la forme jaune, sensible à un changement de l'accepteur seul, est retirée de la vraisemblance ;
- **[2] décalage d'horloge appris ;**
- **[3] fond additif appris par ROI ;**
- **[7] grille de θ adaptative :** recentrée à ±6 sd, avec τ choisi pour décorréler ;
- **[8] composante de saut uniforme (p = 0,005).**

Les paramètres lents [2][3] sont estimés par EM : c est pris aux quantiles 10/30/50/70/90 % de sa loi filtrée, avec 4 itérations à la trame 100 (qui retraitent les trames 0-99), puis toutes les 100 trames sur les 600 dernières.

Deux scénarios :
- **A :** les données ci-dessus.
- **B :** les mêmes, plus un fond additif (4 % du bleu, 2 % du jaune, × U(0,5 ; 1,5) par ROI) et le photoblanchiment de l'accepteur seul, de 1 à 0,95.

| | A : rms / biais / couverture | B : rms / biais / couverture |
|---|---|---|
| référence | 4,06 / −3,37 / 64 % | 16,2 / +15,2 / 14 % |
| 1-2-3 | **1,67 / −1,13 / 94 %** | 17,9 / +17,3 / 1 % |
| tout sans [2] | 2,24 / −1,81 / 78 % | **2,48 / −1,20 / 77 %** |
| tout | 3,31 / −2,94 / 49 % | 3,43 / −2,31 / 56 % |

- **Scénario B :** la référence et 1-2-3 s'effondrent (η part à 1,31-1,37). Les modifications [1][3] gardent une erreur 7 fois plus petite.
- **Scénario A, sans défaut :** elles coûtent 0,6 mM.

Ce que la mise au point a montré :
- **ρ_A, un facteur propre à l'accepteur appris en ligne, n'est pas identifiable.** Il partait à 1,2-1,26 pour 0,98 vrai, parce que dans la forme jaune ρ_A et c agissent presque tous deux par la part de fuite. D'où l'ancre bleue [1].
- **Le décalage d'horloge ne s'apprend pas en ligne.** La forme bleue de mCer3 y est extrêmement sensible : un décalage non corrigé de 3 ps fait lire 6,4 / 11,8 / 22,4 mM à 5 / 10 / 20 mM, et 25 ps fait lire 23 / 34 / 60 mM. Appris, le décalage reste biaisé de +1 à +4 ps, parce qu'il se confond avec c dans la forme bleue. Il faut le mesurer par une référence en début de séance, à mieux que 1 ps.

## Relancer

Depuis ce dossier :

    python3 essai.py 3000 derive            # les 5 variantes, bilans → sorties/resultats_essai_3000.json (~1 h 30)
    python3 traces.py ref 1 2 123           # les mêmes avec les traces → sorties/traces/ (mêmes bilans exactement)
    python3 figures_traces.py               # figures → sorties/figures/
    python3 tout.py 3000 A tout_sans_delta tout       # « tout » → sorties/resultats_tout_A_3000.json
    python3 tout.py 3000 B tout_sans_delta tout 123 ref
    python3 tout.py 3000 A abl1 abl3 abl7 abl8        # 1-2-3 + une seule modification
