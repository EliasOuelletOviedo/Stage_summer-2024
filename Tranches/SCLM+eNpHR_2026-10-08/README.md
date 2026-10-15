# SCLM + eNpHR — tranche du 2026-10-08, FLIM en temps réel (une session)

Analysé avec le pipeline commun `../sclm_flim/` (voir son README pour la méthode) :

```bash
cd ../sclm_flim
python3 step1_decode.py ../SCLM+eNpHR_2026-10-08
~/.venvs/cellpose/bin/python step2_cellpose.py ../SCLM+eNpHR_2026-10-08
python3 step3_select.py ../SCLM+eNpHR_2026-10-08 --tout
python3 step4_declins.py ../SCLM+eNpHR_2026-10-08
python3 step5_analyse.py ../SCLM+eNpHR_2026-10-08
python3 step6_ratio_corrige.py ../SCLM+eNpHR_2026-10-08
python3 step7_retour.py ../SCLM+eNpHR_2026-10-08
```

Les figures de présentation (fond noir) ont leur propre script, dans leur dossier :

```bash
python3 présentation/presentation.py      # depuis ce dossier, après l'étape 7
python3 présentation/slide_resultats.py   # slide de résultats FR / EN, fond clair et noir
python3 présentation/presentation.py retour overlay   # seulement les figures dont le nom contient ces mots
```

**Les résultats ci-dessous portent sur les 34 ROIs bruts de Cellpose (`--tout`), pas
sur un tri visuel** — comme au 2026-09-21. Après un tri (`step3_select.py` sans
`--tout`), relancer l'étape 5 suffit.

## Ce que contient `data/`

| | |
|---|---|
| Session FLIM | 15:38:47 → 15:47:26, 8 647 images à 16,7 Hz, 172 blocs de 50 images (3,0 s) |
| Photons | ch1 (CFP) 94 M, ch2 (YFP) 354 M dans l'image ; aucune perte de FIFO |
| Protocole | 4 stims de 10,2 s, commande 1 à 100 % : 30, 150, 270, 390 s (horloge des cartes : +90 ms) |
| Stim | commande 1 = puissance du 1064 nm, avec sa porte (commande 2 à 5 V) d'après `bench.toml` |
| ROIs FLIMApp | 8 polygones dessinés pendant l'expérience (`manual-10` … `manual-17`, ici F1…F8) |
| Images fixes | 800 et 950 nm à 15:15–15:17 (« avant ») et 15:49–15:51 (« après »), canal C2 (C1 quasi vide) |

Le `debug_report.txt` montre des pertes de SYNC **avant** la session (14:58–15:34,
laser coupé ou câble), aucune pendant ; pendant la session, seules des entrées du
journal de FLIMApp ont été perdues (file pleine), pas des photons. Le champ ne dérive pas : ≤ 0,4 px
sur toute la session (`qc/derive.png`), les masques restent fixes.

## Contrôles

- **Décodage** : les photons des 8 ROIs FLIMApp, reconstruits depuis les `.spc`, sont
  identiques à `frames.csv` pour les 172 blocs et les deux canaux (0 écart sur 2 752,
  `qc/validation_flimapp.csv`) ; les débuts de bloc tombent à 0,05 ms près.
- **Durées de vie** : l'ajustement hors ligne suit celui du direct, bloc par bloc,
  r = 0,988 (CFP) et 0,989 (YFP), avec un écart constant de +8 et +11 ps
  (`qc/validation_duree_de_vie.png`) : le décalage de l'IRF n'est pas paramétré pareil.
- **Cohérence des deux lectures du FRET** : d'une cellule à l'autre, R0 et τ0 CFP sont
  anticorrélés (ρ = −0,72, p = 2·10⁻⁶) — les cellules à ratio élevé ont bien le donneur
  le plus éteint.
- Bruit par bloc de 3 s (somas, médiane) : 44 ps sur τ CFP, 1,2 % sur le ratio.

## Images fixes : les « avant » ne sont pas au plan de la session

Ramenées au repère FLIM (échelle x 1,015, translation nulle — même scanner), les images
« après » se superposent au FLIM (NCC 0,84 à 800 nm) ; les « avant », non (0,40 et 0,34),
sans décalage à trouver, ni global ni local. Entre elles, les deux « avant » se
ressemblent (0,60, même plan à deux longueurs d'onde), mais 800 avant / 800 après
ne donne que 0,38 : **la mise au point a changé entre 15:17 et la session**
(`qc/recalage_fovs.png`). On ne peut donc pas comparer une même cellule avant / après
sur les images Bliq ; seules les images « après » sont mesurées par ROI.

## Ce que disent les résultats

Réponse = moyenne des 10 s après la stim − moyenne des 10 s avant (blocs de stim
exclus, la stim ajoutant +3,7 % de photons dans ch2 et +1,1 % dans ch1). Somas Cellpose :
médiane sur 34 ROIs, p de Wilcoxon sur l'écart.

| | ΔR/R0 (points) | Δτ CFP (ps) | Δ photons CFP | Δ photons YFP |
|---|---|---|---|---|
| stim 1 | **−1,56** (27/34 ↓, p = 3·10⁻⁵) | +13,4 (p = 0,1) | +1,8 % | +0,3 % |
| stim 2 | **−1,24** (27/34, p = 9·10⁻⁴) | +11,6 (p = 0,02) | +1,4 % | +0,2 % |
| stim 3 | **−1,32** (27/34, p = 1·10⁻⁶) | +11,1 (p = 0,2) | +1,3 % | 0,0 % |
| stim 4 | **−1,00** (26/34, p = 3·10⁻⁴) | +9,8 (p = 0,04) | +0,6 % | −0,7 % |
| témoins 1–4 | −1,00, −0,33, −0,38, −0,33 | −6,6 à +4,1 | −0,3 à +2,1 % | −0,5 à +1,6 % |

Par ROI, la moyenne des 4 stims contre celle des 4 témoins : **ΔR −1,27 contre −0,42
point (p = 0,004)**, **Δτ CFP +11,0 contre −3,7 ps (p = 0,016)**. Les deux lectures vont
dans le même sens : moins de FRET après chaque stim.

Dans le neuropile et le champ entier, la marche est plus nette encore, faute de dérive
entre les stims : ΔR −1,8 à −3,1 points contre −0,5 à +0,2 pour les témoins ; Δτ CFP
+4 à +13 ps contre −7 à +4 (`results/decomposition.png`).

### Une descente graduelle, et des somas qui reviennent à leur base

Sans correction, les stims semblent s'empiler en escalier : en fin de session, ratio
−8,3 %, τ CFP +79 ps, photons CFP +9,7 % (somas). L'étape 6 sépare les deux choses :
la **descente graduelle** du ratio est estimée hors des stims (la base, puis de 30 s
après chaque stim jusqu'à la suivante), avec une pente commune et un niveau libre par
intervalle, pour que les marches n'entrent pas dans la pente (méthode dans
`../sclm_flim/step6_ratio_corrige.py`). Ratio image par image, en % des 10 s d'avant :

| | descente graduelle | 10 s après la stim | 110 s après (10 s avant la suivante) |
|---|---|---|---|
| 34 somas réunis | **−1,96 %/100 s** (−9,7 % sur la session) | −1,1 à −1,7 % | **−0,3 à +0,1 %** |
| neuropile | −0,20 %/100 s | −2,3 à −3,0 % | −1,3 à −1,7 % |
| champ entier | −1,00 %/100 s | −1,8 à −2,3 % | −0,6 à −1,1 % |

- **Dans les somas, la réponse revient entièrement à la base** avant la stim suivante :
  l'escalier y était fait de la descente graduelle, pas des stims. Cette descente existe
  avant toute stim (−2,01 %/100 s sur les 30 premières secondes) : elle vient de
  l'imagerie, pas de la stimulation.
- **La part qui ne revient pas est dans le neuropile**, qui, lui, ne descend presque pas
  entre les stims.
- Ces chiffres tiennent quelle que soit la durée écartée après chaque stim (20 à 60 s :
  reste des somas entre −0,4 et +0,3 %).

### Temps de remontée des somas : ~13 s, le même pour toutes les cellules

Étape 7 : mono-exponentielle C + A·exp(−t/τ) sur le ratio corrigé des somas, de la fin
de chaque stim jusqu'à 108 s après (tranches de 1 s, en % des 10 s d'avant), équilibre C
libre. Imposer C = base d'avant allonge τ quand le retour n'est pas complet (stim 1 :
12 s libre, 31 s imposé).

| somas réunis | stim 1 | stim 2 | stim 3 | stim 4 | moyenne des 4 retours |
|---|---|---|---|---|---|
| τ (± erreur type) | 12,8 ± 3,0 s | 17,9 ± 3,5 s | 11,4 ± 2,5 s | 11,8 ± 2,9 s | **13,2 ± 1,6 s** |
| A | −1,7 % | −2,1 % | −1,8 % | −1,7 % | −1,8 % |
| équilibre | −0,3 % | −0,1 % | +0,2 % | +0,1 % | 0,0 % |

- **Par soma**, sur la moyenne de ses 4 retours : 22 mesurables sur 34 (les autres sans
  baisse nette), τ médian 14,4 s [10,8–18,5].
- **Cet étalement est celui du bruit de mesure.** En simulant chaque soma avec son
  amplitude et son propre bruit mais le même τ pour tous (13,2 s), on obtient 12,8 s
  [10,0–16,5] ; Kolmogorov-Smirnov p = 0,20. Rien n'indique que les cellules reviennent
  à des vitesses différentes ; ce jeu ne permet pas non plus de l'exclure pour des
  écarts plus fins que ce bruit.
- Un retour isolé (un soma, une stim) est trop bruité pour être lu seul : erreur type
  médiane ≈ 100 % de τ.
- τ des somas réunis selon la durée écartée pour estimer la dérive (20 à 60 s) : 12,3 à
  14,6 s.

`results/ratio_corrige_champ.png` et `ratio_corrige_cellpose.png` montrent le ratio brut
avec la dérive ajustée, puis corrigé. τ CFP n'est pas corrigé de la même façon pour
l'instant ; `results/recuperation.csv` (étape 5) porte sur le ratio non corrigé.

### Réserves sur l'interprétation

Une entrée de Cl⁻ par eNpHR donnerait bien moins de FRET. Plusieurs observations
invitent pourtant à ne pas conclure seul sur ce jeu :

1. **La stim excite l'accepteur.** Pendant la stim, ch2 (YFP) gagne 3,7 % de photons :
   le 1064 nm excite la YFP à deux photons, et peut donc la photoblanchir, ce qui
   donnerait une baisse **durable** du FRET. Dans les somas, une fois la descente
   graduelle retirée, la réponse revient à la base : ce n'est pas ce qu'on y voit. Dans
   le neuropile, en revanche, une marche de −1,3 à −1,7 % reste après chaque stim.
2. **Le neuropile répond plus fort que les somas** (−2,3 à −3,0 % contre −1,1 à −1,7 %)
   et ne revient qu'à moitié.
3. **La lumière CFP gagnée ne colle pas à la durée de vie gagnée.** Un donneur moins
   éteint gagne en intensité ce qu'il gagne en durée de vie. Ici, sur le champ entier,
   les photons CFP montent de +1,1 à +2,4 % par stim quand τm ne monte que de +5 à +8 ps
   (+0,4 à +0,5 % de 1,48 ns) : la plus grande part du gain de CFP vient d'autre chose
   que le FRET.
4. **Les photons YFP baissent peu ou pas** (−1,4 à +0,3 %, la plus forte baisse à la
   stim 4), alors qu'une entrée de Cl⁻ (YFP éteinte) comme un photoblanchiment devraient
   les faire baisser — en partie masqué par la YFP de eNpHR dans ch2.
5. **Les cellules les plus riches en eNpHR répondent moins.** Indice = C2 à 950 nm / C2 à
   800 nm (images « après ») : ρ = +0,59 avec ΔR (p = 3·10⁻⁴), ρ = −0,34 avec Δτ CFP
   (p = 0,05) (`results/expression_eNpHR.png`). Pour ΔR, une partie s'explique par la
   YFP de eNpHR qui dilue le ratio ; pour τ CFP, non. L'indice mêle toutefois l'excitation
   de la YFP de SCLm et de eNpHR à 950 nm : c'est une piste, pas une mesure d'expression.

Les stims « IR » du 2026-09-21 ne donnaient aucune réponse ; si c'est la même source,
cette différence est à éclaircir. Le témoin qui trancherait : le même protocole sur une
tranche SCLm sans eNpHR (ou une stim qui n'excite pas la YFP), pour mesurer ce que la
lumière seule fait au FRET.

## Sorties

| Fichier | Contenu |
|---|---|
| `projections/images.png` | photons ch1, ch2, ratio, temps moyen CFP et YFP, images Bliq recalées |
| `projections/session.json` | blocs, fenêtres de stim, fuite, dérive, recalage des images fixes |
| `masks/cellpose_overlay.png` | les 34 ROIs Cellpose et les 8 ROIs FLIMApp (cyan) |
| `qc/par_image.png` | photons du champ entier image par image (60 ms) et ratio |
| `results/mesures_blocs.csv` | une ligne par ROI × bloc : photons, ratio, τ, τm, χ², écarts à la base |
| `results/avant_apres_stim.csv` | une ligne par ROI × stim (ou témoin) : avant, après, écart, pour chaque mesure |
| `results/resume_stims.csv` | médianes, nombre en baisse, p de Wilcoxon par famille × stim |
| `results/resume_par_roi.csv` | R0, τ0, réponse moyenne aux stims et aux témoins, indice eNpHR |
| `results/recuperation.csv` | avant / après / 10 s avant la stim suivante, par ROI et par stim |
| `results/images_fixes_C2.csv` | C2 moyen par ROI sur chaque image fixe (`meme_plan` dit si elle est utilisable) |
| `results/trace_ratio.png`, `trace_duree_de_vie.png` | ratio et τ CFP de chaque soma dans le temps |
| `results/decomposition.png` | ratio, τ CFP, photons CFP et YFP : somas, neuropile, champ entier |
| `results/reponse_alignee.png` | les 4 stims et les 4 témoins alignés sur leur début |
| `results/avant_apres_stim.png` | un trait par ROI, avant → après, pour chaque stim et les témoins |
| `results/heatmap.png` | ΔR/R0 et Δτ CFP par ROI × temps |
| `results/trace_flimapp.png` | les 8 ROIs dessinées pendant l'expérience |
| `results/expression_eNpHR.png` | indice 950/800 nm contre la réponse |
| `results/ratio_corrige_champ.png`, `ratio_corrige_cellpose.png` | ratio image par image, brut avec la dérive ajustée, puis corrigé |
| `results/derive_ratio.csv` | pente commune et pente propre de chaque intervalle, par famille |
| `results/reponse_ratio_corrige.csv` | ratio corrigé 10 s avant, 10 s après et avant la stim suivante |
| `results/ratio_corrige_par_image.csv` | ratio brut et corrigé, image par image, champ et somas |
| `results/retour_equilibre_cellpose.png` | somas réunis : les 4 retours ajustés, et chacun en détail |
| `results/temps_de_remontee.png` | distribution des τ : par soma (IC 95 %), contre le bruit seul, stim par stim |
| `présentation/` | figures sur fond noir pour les diapos et le script qui les produit (`presentation.py`, réglages en tête : couleurs, polices, titres, titres d'axes) : images nues (image, contours et barre de 50 µm, sans titre ni marge) : FLIM sans ROI, avec les ROIs Cellpose, avec les ROIs FLIMApp, avec les deux (`cellpose_overlay_*`), et les mêmes sur le 950 nm « avant » (`cellpose_overlay_950nm_avant_*` : autre plan focal que la session, les contours n'y tombent pas exactement sur les mêmes cellules) ; somas réunis : photons et ratio image par image, ratio corrigé, retour à l'équilibre, temps de remontée |
| `présentation/slide_resultats_{FR,EN}[_dark].png` | slide de résultats sur le modèle du 2026-09-21, le ratio en haut et la durée de vie CFP en bas : A ratio et C durée de vie sur la session, avec le retour ajusté après chaque stim (C_k + A_k·exp(−t/τ) sur la médiane non lissée, un τ commun aux 4 stims : 8,2 ± 2,4 s pour le ratio, 19 ± 12 s pour la durée de vie) ; B et D τ de retour par stim, un point par cellule, retours mesurables et τ ≤ 40 s (B : ratio, τ de l'étape 7 ; D : durée de vie, même ajustement sur chaque cellule, en blocs de 3 s). A et C corrigés de la descente graduelle (pente de la cellule médiane : ratio −1,68 %/100 s, τ CFP +16,7 ps/100 s, méthode de l'étape 6). Un τ isolé de la durée de vie est deux fois plus bruité que celui du ratio (erreur type médiane 2 τ, contre 1 τ) : D montre surtout la dispersion du bruit (`slide_resultats.py`) |
| `présentation/duree_de_vie_redescente.png` | durée de vie CFP des somas sur la session (médiane ± écart interquartile, corrigée de la descente graduelle comme dans les slides), blocs de stim compris, et les redescentes après chaque stim ajustées ensemble hors stim : C_k + A_k·exp(−t/τ), un τ commun aux 4 stims (19 ± 12 s), une amplitude et un équilibre par stim (A = +3, +11, +18, +9 ps) ; une à une, les redescentes sont trop bruitées pour être ajustées (`presentation.py`) |
| `présentation/duree_de_vie_redescente_2exp.png` | la même chose avec deux durées de vie, le modèle des calibrations (a1·exp(−t/τ1) + (1 − a1)·exp(−t/τ2), τ moyen pondéré par l'intensité), ajusté bloc par bloc sur la somme des déclins des 34 somas (~260 000 photons ; un soma seul, ~6 000, ne suffit pas à trois paramètres libres) : τ1, τ2 et τ moyen, chacun rapporté à sa droite de dérive, et leurs redescentes (τ commun). Avant la 1re stim : τ1 = 2,73 ns, τ2 = 0,50 ns, τ moyen = 1,89 ns ; sur toute la session a1 = 0,25, contre 0,49 à 0,71 pour SCLM purifiée : une composante courte étrangère s'ajoute dans le canal CFP. Seul τ moyen répond dans le même sens à presque chaque stim (A = +25, +15, +12, −5 ps, τ = 11 ± 25 s), plus bruité que la mono-exponentielle ; τ1 (±77 ps d'un bloc à l'autre) et τ2 (±15 ps) changent de signe d'une stim à l'autre : en blocs de 3 s, les deux composantes ne se séparent pas (`presentation.py`) |
| `présentation/duree_de_vie_redescente_2exp_canal2.png` | la même chose pour le canal 2 (YFP, qui contient aussi la YFP de eNpHR). Avant la 1re stim : τ1 = 3,30 ns, τ2 = 1,11 ns, τ moyen = 2,85 ns ; a1 = 0,55 sur toute la session. Pendant les stims, la fuite du 1064 nm allonge τ1 de ~300 ps et τ moyen de ~150 ps : ce n'est pas une réponse. Après chaque stim, la YFP descend légèrement puis remonte, à l'inverse de la CFP : τ moyen −13, −9, −1, −6 ps, τ1 −33, −15, −18, −16 ps (τ commun ≈ 14 ± 15 s) ; τ2 change de signe. Bosse de ~+50 ps sur τ moyen vers 235 s, hors stim, d'origine inconnue (`presentation.py`) |
| `présentation/correlation_ratio_duree_de_vie.png` | ratio et durée de vie CFP, par soma : à la ligne de base, R₀ contre τ₀ (Spearman ρ = −0,72, p = 2·10⁻⁶) ; en réponse aux stims (moyenne des 4, étape 5), ΔR/R₀ contre Δτ (ρ = −0,44, p = 0,01), contre ρ = 0,00 (p = 0,98) pour les témoins des mêmes cellules ; droites de Theil-Sen (`presentation.py`) |
| `results/temps_de_remontee.csv` | A, τ, erreur type, équilibre : somas réunis, chaque soma, chaque soma × stim |
