# Simulateur Kalman temps réel — étape 1

Construit le 2026-10-16 d'après « Simulateur Kalman temps réel — étape 1.pdf » (dans ce dossier).
C'est un simulateur neuf, indépendant du simulateur FLIM-FRET à deux canaux. Il génère des histogrammes SCLM à
10 Hz le long d'une montée de Cl⁻ (5 → 30 mM), puis les lit comme des données réelles avec trois lectures
(J/B, 5a et leur fusion par filtre de Kalman étendu), avec et sans dérive du gain du canal jaune.

Révisé le même jour d'après « Simulateur Kalman — corrections de l'étape 1.pdf » : réglage de référence
σ_w,g = 3·10⁻³, balayage étendu aux trois budgets, modèle de dérive (g, ġ), P₁₁/err² comme critère de réglage,
origine du biais de 5a et de l'écart sur les amplitudes, critères de sortie vers l'étape 2. La section
« Corrections » ci-dessous en rend compte ; les résultats de la première course y sont rappelés quand ils changent.

## Lancer

```bash
# une seule fois : extraire la titration dans params.toml (environnement du simulateur précédent, lecteur .sdt)
julia --project=../FLIM-FRET_2-channels step0_params.jl
# les tests (7 groupes, 32 vérifications)
julia --project=. -t auto test/runtests.jl
# la course complète (~4 min avec 8 fils d'exécution, vidéos comprises)
julia --project=. -t auto run.jl            # --rapide : 2 répétitions, sans vidéo
```

## Fichiers

| Fichier | Rôle | lignes |
| --- | --- | --- |
| `params.toml` | le seul endroit pour les constantes, critères de sortie compris ; le bloc « step0 » est écrit par `step0_params.jl` | — |
| `step0_params.jl` | le pont unique avec la titration : ajustement global → Langmuir, IRF, fond | 90 |
| `src/params.jl` | lecture de `params.toml`, Langmuir, courbe monotone (valeur, dérivée, inverse) | 107 |
| `src/generate.jl` | trajectoire, dérive, histogrammes (six étapes, mode sans bruit, graine par pas) | 94 |
| `src/estimate.jl` | J/B, 5a, calibrations, covariance de Fisher, R attendu | 158 |
| `src/kalman.jl` | prédiction, mise à jour (forme de Joseph), diagnostics ; état (C, g) ou (C, g, ġ) ; aucune photophysique | 74 |
| `run.jl` | les courses (tous les réglages du filtre sur les mêmes lectures), les CSV, les vérifications, les critères | 276 |
| `src/figures.jl` | rendu après simulation : vidéos, résumé, réglage, Fisher, amplitudes | 186 |
| `test/runtests.jl` | les quatre tests du document, deux tests du filtre seul, le test d'amplitudes | 117 |

Les cinq fichiers du document font 709 lignes (plafond fixé à 800), 895 avec le rendu.

Sorties dans `sorties/` :

- `courses.csv` : le fichier plat, une ligne par pas (répétition 1 de chaque scénario, aux trois budgets, N = 300
  compris comme borne basse), pour le filtre de référence et le modèle (g, ġ) à σ_ġ = 10⁻⁴ ;
- `resume.csv` : un résumé par budget, scénario et lecture ou réglage du filtre, moyenne sur 8 répétitions ;
- `criteres.csv` : les critères de sortie, pour chaque réglage du filtre ;
- `verification_R.csv`, `biais_lectures.csv`, `amplitudes.csv`, `coherence.csv` : les vérifications ;
- `fig_resume.png`, `fig_reglage.png`, `fig_fisher.png`, `fig_amplitudes.png`, et trois vidéos au réglage de
  référence (`video_sans_derive`, `video_rampe`, `video_saut`) avec leur dernière image en PNG. Les `.mp4` ne sont
  pas versionnés (`.gitignore`) ; `run.jl` les régénère à l'identique.

## Les paramètres

Tout vient de l'ajustement global « libre » de la titration SCLM du 2026-10-16 à 800 nm (étape 15 du simulateur
précédent : `global_parametres.csv`, `global_cellules.csv`). Les IRF sont celles de `IRF_800_16X.sdt`, et le fond
est le niveau plat avant le pic de ces IRF. Le simulateur ne relit jamais les histogrammes de la titration.

| Quantité | Valeur | Langmuir sur les 8 points |
| --- | --- | --- |
| τ₂, τ₃ (bleu) | 0,721 ; 2,722 ns | — |
| τ_J, τ_A (jaune) | 1,059 ; 3,253 ns | — |
| J/B (photons, fond exclu) | 5,95 à 0 mM → 1,65 à 131 mM | K = 28,9 mM, écart max 0,03 |
| p (fraction d'amplitude de τ₂, bleu) | 0,707 → 0,576 | K = 56 mM |
| q (fraction d'amplitude de τ_J, jaune) | 0,077 → 0,278 | K = 219 mM |
| d (fraction d'intensité de la diffusion) | bleu 7,9 → 3,4 % ; jaune < 0,5 % | bleu K = 16 mM ; jaune : une droite |
| fond par acquisition de 0,1 s | bleu 187, jaune 327 photons | — |

Trois remarques sur ces valeurs.

- **Le J/B du document (7,03 → 1,60) ne vient pas de la titration.** Ces valeurs sont celles du banc test_0/test_131
  (étape 14). La titration à 8 points donne 5,95 → 1,65. J'ai pris celles de la titration pour que J/B et les
  amplitudes viennent du même ajustement.
- **La diffusion entre comme une fraction d'intensité d(C) par canal**, soit une copie de l'IRF (τ = 1 ps) comme dans
  l'ajustement global. Elle n'entre pas comme un débit fixe : à 10 kcps, les 8 400 photons/s de diffusion mesurés
  pendant la titration feraient presque tout le signal. La Langmuir du jaune dégénère en droite (K sur sa borne de
  10⁴ mM), parce que les 8 points ne sont pas monotones. Elle reste sous 0,5 % du canal.
- **Le fond pèse lourd** : 187 photons de fond pour 1 000 photons bleus à 10 kcps. C'est le régime réel, et c'est
  pourquoi la lecture J/B soustrait le fond.

## Écarts au document, et pourquoi

1. **Convolution par FFT sur 2n = 512 points, pas 256.** Le TAC couvre 16,4 ns et la période laser 12,5 ns. Une
   convolution circulaire sur les 256 canaux replierait donc à 16,4 ns, pas à la période. Le déclin est périodisé
   analytiquement (série géométrique), y compris pour les retards négatifs (la fin de l'impulsion précédente). La
   convolution linéaire sur 512 points est alors exacte, sans bord, ce qui est le but de l'étape 4.
2. **5a n'a que ses 4 paramètres et ne modélise pas la diffusion**, que le générateur contient. La calibration,
   construite avec le même générateur, absorbe l'écart qui en résulte sur p et q (voir correctif 3), et
   l'aller-retour est exact. L'invariance à g est exacte au bruit numérique près sur un histogramme que 5a décrit
   exactement ; sur le générateur complet, il reste un déplacement d'environ 0,002 mM-équivalent pour g = 1,2, soit
   10⁻⁴ du bruit d'un pas (test 2, en deux volets).
3. **L'observable de 5a est une combinaison de p et q.** Les deux canaux ne partagent aucun paramètre, donc 5a se
   sépare en un ajustement par canal. Règle écrite d'avance : le jaune entre dans l'observable s'il réduit σ_C d'au
   moins 10 % à 20 mM. Il la réduit de 16 % (p seul : 15,8 mM, q seul : 24,5 mM, combinaison : 13,3 mM à 1 000
   photons), soit 29 % de variance en moins. L'observable est z_S = w_p·p + w_q·q avec w ∝ (p′/σ_p², q′/σ_q²),
   normalisée pour que z_S′ = 1 mM⁻¹ à 20 mM. Les signes opposés des poids (−335, +138) ne demandent pas
   d'explication de plus : p baisse avec le Cl⁻ et q monte, et chaque poids a le signe de sa pente, ce qui fait
   s'additionner les deux sensibilités.
4. **Le filtre reçoit le R de Fisher attendu à son estimation (Ĉ⁻, ĝ⁻), pas celui du tirage.** Avec la variance
   calculée sur chaque histogramme, ĝ se stabilisait à 0,986 au lieu de 1 : var(J/B) ∝ (J/B)², donc les tirages
   hauts pèsent moins que les tirages bas, et des poids corrélés aux mesures biaisent la moyenne que fait le filtre.
5. **Calibration par budget.** Comme 5a est mal spécifié (point 2), son ajustement sans bruit dépend du rapport
   fond/signal. Chaque budget a donc sa courbe, toujours à g = 1.
6. **R diagonal retiré des courses et des figures** (corrections, « ce qui sort »). Le code garde R complet. Pour
   mémoire, la première course donnait un écart R complet / R diagonal de 4 % au plus sur la RMSE et de 8 % sur
   P₁₁/err², avec une corrélation J/B–5a de −0,10 à 1 000 photons : les deux lectures sont quasi orthogonales,
   parce que J/B lit des totaux et 5a la forme du déclin. C'est une ligne de méthodes, pas un résultat.

## Vérifications

**Les quatre tests du document passent**, plus trois autres (32 vérifications) :
1. Aller-retour sans bruit : J/B et 5a retrouvent 20 mM (nœud de la grille) et 23,7 mM (hors nœud) à 0,01 mM près.
2. Invariance de 5a : exacte (< 10⁻⁸) sur un histogramme que 5a décrit exactement, < 10⁻³ σ sur le générateur complet.
3. Sensibilité de J/B : jaune × 1,2 multiplie J/B par 1,2 à 10⁻⁹ près.
4. Filtre à g constant : sur 8 courses, le Kalman est à moins de 5 % de J/B (sur toute la course et sur les 30
   dernières secondes) et ĝ à moins de 2 % de 1, pour la référence comme pour l'ancien défaut. Une course seule
   varie de ±5 % : ĝ erre autour de 1, d'autant plus que σ_w,g est grand.
5. Filtre seul, mesures fabriquées, g = 1,05 inconnu : ĝ et Ĉ tombent dans 3√P de la vérité. Sans 5a, l'erreur sur Ĉ
   est plus de 3 fois plus grande et P₁₁ ne se resserre pas.
6. Modèle (g, ġ), mesures fabriquées, rampe de 10 % : ĝ et la pente tombent dans 3√P de la vérité, et la marche
   aléatoire raide traîne au moins 5 fois plus.
7. Amplitudes : sans diffusion dans le générateur, 5a retrouve p et q injectés à 10⁻⁶ près, à 0, 20, 80 et 131 mM.

![Fisher contre Monte-Carlo](sorties/fig_fisher.png)

Fig. 1 : **le R de Fisher est juste.** Variances de J/B et de 5a, covariance et corrélation, prédites par Fisher
(trait) et mesurées sur 4 000 tirages à 20 mM (points, ±1σ), aux trois budgets. À 1 000 photons : var J/B 0,02139
contre 0,02126 (0,6 %), var 5a 180,2 contre 174,7 (3 %), ρ −0,097 contre −0,098. C'est la validation du socle
statistique (dérivées, convolution, vraisemblance de Poisson).

**Cohérence des calibrations** : à 20 mM, sans bruit, J/B et 5a retrouvent tous deux 20 mM (`coherence.csv`).

## Corrections du 2026-10-16

### Correctif 1 — le réglage de σ_w,g

La référence passe de σ_w,g = 10⁻⁴ à 3·10⁻³ par pas. Le balayage va maintenant de 10⁻⁴ à 3·10⁻², aux trois
budgets, et tous les réglages tournent sur les mêmes histogrammes : les comparaisons sont appariées.

![Réglage](sorties/fig_reglage.png)

Fig. 2 : **le critère de réglage.** En haut, P₁₁ annoncé ÷ erreur quadratique réelle contre σ_w,g, une courbe par
dérive, avec la bande [0,8 ; 2] des critères de sortie (gris). En bas, le biais des 30 dernières secondes et la RMSE
sur toute la course, à 1 000 (trait plein) et 3 000 photons (tirets).

| σ_w,g (1 000 photons) | RMSE sans dérive | biais, rampe 10 % | biais, rampe 20 % | P₁₁/err² (min → max) | rattrapage du saut (moitié) |
| --- | --- | --- | --- | --- | --- |
| 10⁻⁴ (ancien défaut) | 1,19 | −1,45 | −2,89 | 0,33 → 1,52 | 23 s |
| 3·10⁻⁴ | 1,19 | −1,31 | −2,65 | 0,37 → 1,53 | 18 s |
| 10⁻³ | 1,22 | −0,64 | −1,45 | 0,72 → 1,64 | 8,7 s |
| **3·10⁻³ (référence)** | 1,38 | **−0,03** | −0,31 | 1,61 → 1,80 | 1,6 s |
| 10⁻² | 1,88 | +0,29 | +0,21 | 1,79 → 1,85 | 0,2 s |
| 3·10⁻² | 2,42 | +0,46 | +0,44 | 1,71 → 1,74 | 0,4 s |

- **L'optimum est bien à 3·10⁻³ à 1 000 et 3 000 photons, à 10⁻³ à 300.** Au-delà de 3·10⁻³, le sur-desserrage se
  paie cher : la RMSE sans dérive monte à 1,88 puis 2,42 mM (celle de J/B est 2,00), ρ₁₂ passe de 0,68 à 0,88, et le
  biais change de signe. Il n'y a rien à gagner à 10⁻² : P₁₁/err² est déjà dans la bande à 3·10⁻³.
- **Le critère du correctif 4 désigne le même réglage.** σ_w,g tel que P₁₁/err² sous la dérive la plus forte (rampe
  de 20 %) soit ≥ 1 et le plus proche de 1 : 10⁻³ à 300 photons, 3·10⁻³ à 1 000 et à 3 000. Le réglage dépend bien
  du budget, mais la grille (un facteur 3 entre deux valeurs) ne sépare pas 1 000 de 3 000.
- **Le prix de la référence** : sans dérive, la RMSE passe de 1,19 à 1,38 mM à 1 000 photons (0,82 à 0,97 à 3 000),
  et ĝ erre davantage d'une course à l'autre (±0,04 dans la vidéo du saut). Le Kalman reste sous J/B partout.

### Correctif 2 — modéliser la dérive comme une rampe

`kalman.jl` accepte un état (C, g, ġ) : g_{k+1} = g_k + Δt·ġ_k + w_g, ġ_{k+1} = ρ·ġ_k + w_ġ, avec Q sur g laissé à
10⁻⁴. Trois valeurs de σ_ġ (10⁻⁵, 10⁻⁴, 10⁻³ s⁻¹ par pas), plus la variante Ornstein-Uhlenbeck (ρ = 1 − Δt/τ,
τ = 30 s, σ_ġ = 10⁻⁴). La règle pour retenir le modèle à vitesse était écrite avant la course (`params.toml`) : il
doit battre la référence en RMSE sur les rampes, sans perdre plus de 10 % de RMSE sur le saut et sans dérive, ni
plus de 0,1 mM de biais sur le saut.

| Modèle de g (1 000 photons) | rampe 10 % : biais, RMSE | rampe 20 % : biais, RMSE | saut : biais, RMSE, rattrapage | sans dérive : RMSE | bat la référence |
| --- | --- | --- | --- | --- | --- |
| marche, σ_w,g = 10⁻⁴ | −1,45 ; 1,55 | −2,89 ; 2,49 | −0,99 ; 1,47 ; 23 s | 1,19 | — |
| **marche, σ_w,g = 3·10⁻³** | −0,03 ; 1,36 | −0,31 ; 1,44 | +0,21 ; 1,37 ; 1,6 s | 1,38 | — |
| (g, ġ), σ_ġ = 10⁻⁵ | +0,14 ; 1,34 | +0,12 ; 1,38 | +0,26 ; 1,43 ; 6,3 s | 1,36 | oui à 1 000, non à 3 000 |
| (g, ġ), σ_ġ = 10⁻⁴ | +0,30 ; 1,40 | +0,34 ; 1,42 | +0,59 ; 1,46 ; 3,1 s | 1,41 | non |
| (g, ġ), σ_ġ = 10⁻³ | +0,35 ; 1,85 | +0,36 ; 1,80 | +0,39 ; 1,84 ; 0,7 s | 1,86 | non |
| OU, σ_ġ = 10⁻⁴, τ = 30 s | −0,06 ; 1,30 | −0,35 ; 1,41 | +0,24 ; 1,37 ; 4,8 s | 1,31 | oui à 1 000, non à 3 000 |

- **Le modèle à vitesse supprime bien le traînage** : à σ_ġ = 10⁻⁵, son biais ne dépend plus de la pente (+0,14 à
  10 %, +0,12 à 20 %), là où la marche aléatoire passe de −0,03 à −0,31. Le biais qui reste (+0,1 à +0,2 mM) est
  celui de l'observable de 5a, présent aussi sans dérive.
- **Mais il ne bat pas la marche desserrée partout.** À 1 000 photons, (g, ġ) à σ_ġ = 10⁻⁵ et la variante OU passent
  la règle, de peu (RMSE 1,34 contre 1,36 sur la rampe de 10 %). À 3 000 photons, aucun ne la passe : sur le saut, le
  modèle à vitesse dépasse (biais +0,33 contre +0,13 mM, au-delà de la tolérance de 0,1) et rattrape plus lentement
  (5,2 s contre 1,9 s). C'est le dépassement que le document craignait.
- **Décision, selon la règle écrite d'avance : garder la marche aléatoire à 3·10⁻³.** Elle est plus simple, elle
  passe tous les critères aux deux budgets, et le gain du modèle à vitesse est dans le bruit là où il existe. La
  variante OU est la meilleure des variantes à vitesse ; elle redeviendra intéressante si les dérives réelles
  (expérience de profondeur) s'avèrent plus rapides que 20 % en 90 s.

### Correctif 3 — le biais de 5a

**Le +3,31 mM n'est pas un biais d'estimateur.** `biais_lectures.csv` sépare les deux quantités sur 4 000 tirages à
C fixé : le biais de l'observable brute z_S, que reçoit le filtre, et celui de la concentration inversée, calculé sur
les tirages qui tombent dans la plage.

| 1 000 photons | 5 mM | 10 mM | 20 mM | 30 mM |
| --- | --- | --- | --- | --- |
| biais de l'observable z_S (mM-équivalent) | +0,08 | +0,30 | +0,26 | +0,34 |
| biais de la concentration inversée (mM) | +3,9 | +3,1 | +2,8 | +3,6 |
| tirages hors plage (échecs) | 23 % | 8,6 % | 1,3 % | 0,6 % |

L'observable n'est biaisée que de +0,1 à +0,3 mM-équivalent : c'est le vrai biais d'estimateur à faible comptage,
et c'est lui que le filtre absorbe (le Kalman sans dérive est biaisé de +0,28 mM à 1 000 photons, +0,12 à 3 000).
Les +3 mM viennent de l'inversion. Avec σ ≈ 13 mM par pas, la courbe de calibration s'aplatit vers les fortes
concentrations et étire la queue haute : même à 20 mM, où 1,3 % des tirages seulement échouent, la moyenne des
concentrations inversées est décalée de +2,8 mM. À basse concentration, la troncature (« hors plage = échec »)
s'y ajoute. Le biais de 5a seul ne transfère donc pas 3 mM sur ĝ ; corriger p̂ par bootstrap n'est pas utile
avant l'étape 2.

**L'écart sur les amplitudes vient entièrement de la diffusion.** Le test qui tranche : générer sans bruit sur toute
la plage, avec la diffusion du générateur et sans elle.

![Amplitudes](sorties/fig_amplitudes.png)

Fig. 3 : p et q retrouvés par 5a contre injectés (Langmuir de la titration), générateur avec et sans diffusion ; en
bas, l'écart relatif.

- **Sans diffusion, 5a retrouve exactement p et q injectés** (écart max 10⁻¹⁵, test 7). La normalisation de p est
  donc la même des deux côtés : fractions d'amplitude des deux populations, diffusion à part.
- **Avec diffusion, l'écart est de +9,5 à +10,4 % sur p et de +14 à +6 % sur q.** 5a, qui n'a pas de composante
  pour la copie de l'IRF, la range dans la composante rapide.
- **Le critère proposé par le document se serait trompé.** Il disait : écart constant → normalisation, écart qui varie
  avec la concentration → diffusion. Sur p, l'écart est presque constant (9,5 à 10,4 %) alors que d_bleu passe de
  7,8 à 3,5 % : la composante rapide perd elle aussi du poids quand le Cl⁻ monte, et le rapport des deux bouge peu.
  C'est le test sans diffusion qui tranche, pas la forme de l'écart.
- **Conséquence** : p̂ et q̂ de 5a ne sont pas des fractions de population et ne doivent pas être cités comme telles.
  Pour citer une fraction, il faudrait une cinquième composante, la diffusion, à amplitude libre dans 5a. La
  concentration, elle, n'est pas touchée : la calibration porte la même diffusion.

### Correctif 4 — P₁₁/err² comme critère de réglage

La figure 2 en est la mise en œuvre. Au mauvais réglage, le filtre est faux et confiant (P₁₁/err² = 0,33 à 1 000
photons et 0,15 à 3 000 sous une rampe de 20 %). À la référence, il est entre 1,6 et 1,8 partout : légèrement
sous-confiant, donc sûr. Le critère ne demande ni vérité terrain ni connaissance de la dérive en expérience, à
condition de pouvoir estimer l'erreur réelle, ce que la simulation fait et que l'expérience devra approcher (par
exemple sur un étalon).

### Diagnostics d'innovation, après la relance

- **L'excursion pendant la montée est toujours là** : innovation moyenne du rapport −0,26σ pendant la montée,
  +0,01σ ailleurs, à tous les réglages. C'est la signature du modèle de prédiction trop rigide sur C, laissée
  telle quelle pour l'étape 2 (Q dépendant de u).
- **La variance des innovations du rapport ne s'est pas resserrée** : 0,72 à 1 000 photons et 0,51 à 3 000, pour
  tous les σ_w,g de 10⁻⁴ à 3·10⁻³, contrairement à ce que le document attendait. Ce n'est pas σ_w,g qui la fixe
  mais σ_w,C = 1,2 mM par pas : trop large hors de la montée, il gonfle la covariance prédite du rapport. Même
  remède qu'au point précédent.

## Critères de sortie vers l'étape 2

Chiffrés d'avance (`params.toml`, section `[criteres]`), vérifiés automatiquement (`criteres.csv`) aux deux budgets
de fonctionnement. N = 300 n'est pas un régime de fonctionnement (5a y échoue 25 % du temps, J/B jusqu'à 10 %) : il
reste dans les CSV comme borne basse et sort des figures de conclusion.

| Critère | 1 000 photons | 3 000 photons |
| --- | --- | --- |
| 1. biais du Kalman sous rampe de 10 % sous 0,3 mM | **oui** (−0,03) | **oui** (−0,04) |
| 2. P₁₁/err² entre 0,8 et 2 sous toutes les dérives | **oui** (1,61 à 1,80) | **oui** (1,63 à 1,77) |
| 3. RMSE du Kalman sous celle de J/B partout | **oui** (au plus 0,69 × J/B) | **oui** (au plus 0,90 × J/B) |
| 4. écart sur les amplitudes de 5a expliqué | **oui** : la diffusion, que 5a ne modélise pas (test 7) | idem |

**Les quatre critères tiennent au réglage de référence.** Parmi les autres marches aléatoires, seul σ_w,g = 10⁻²
passe aussi les trois premiers à 1 000 photons, avec une RMSE 36 % plus forte sans dérive ; aucun réglage à 10⁻³ ou
moins ne passe le critère 2. Les modèles à vitesse (σ_ġ = 10⁻⁵ et 10⁻⁴) et la variante OU les passent aussi aux deux
budgets (`criteres.csv`) : le choix entre eux et la référence se fait sur la règle du correctif 2.

## Résultats

![Résumé](sorties/fig_resume.png)

Fig. 4 : (a) biais des 30 dernières secondes contre la dérive, référence et ancien défaut ; (b) le même biais pour
chaque modèle de g ; (c) RMSE contre le budget sous rampe de 10 % ; (d) ĝ autour du saut de 5 %, moyenne de 8
courses, pour chaque modèle de g. 1 000 photons sauf en (c).

![Rampe](sorties/video_rampe.png)

Fig. 5 : dernière image de `video_rampe.mp4` (rampe de 10 %, référence, 1 000 photons). Les deux autres vidéos :
`video_sans_derive.mp4` et `video_saut.mp4`.

### La prédiction du document, au réglage de référence (1 000 photons, 8 répétitions)

| | J/B | 5a | Kalman, σ_w,g = 3·10⁻³ | Kalman, ancien défaut 10⁻⁴ |
| --- | --- | --- | --- | --- |
| **Sans dérive** : RMSE | 2,00 mM | 12,5 mM | 1,38 mM | **1,19 mM** |
| biais, 30 dernières s | −0,01 | +3,3 (inversion) | +0,28 ± 0,30 | +0,22 ± 0,32 |
| **Rampe 10 %** : RMSE | 3,07 | 12,5 | **1,36** | 1,55 |
| biais, 30 dernières s | −3,23 | +3,2 | **−0,03 ± 0,29** | −1,45 ± 0,28 |
| **Rampe 20 %** : biais | −5,62 | +3,2 | **−0,31** | −2,89 |
| **Saut 5 %** : biais | −2,03 | +3,3 | **+0,21**, rattrapé à moitié en 1,6 s | −0,99, en 23 s |
| erreur max pendant la montée | 7,7 | 71 | 4,1 | 3,6 |

### Ce que disent les résultats

**1. Les deux volets de la prédiction tiennent au réglage de référence.** Sans dérive, la fusion fait mieux que J/B
(RMSE 1,38 contre 2,00 mM). Sous dérive, elle reste non biaisée : −0,03 mM sous une rampe de 10 %, contre −3,23 pour
J/B, et elle rattrape un saut de 5 % en 1,6 s. À 3 000 photons, même tableau : −0,04 contre −3,23 mM.

**2. La seule différence avec la première course est Q.** Le mécanisme était juste, le réglage de σ_w,g non : à 10⁻⁴,
une marche aléatoire ne peut pas suivre une rampe de 1,1·10⁻⁴ par pas et la dérive passe dans Ĉ. À 3·10⁻³ elle suit,
au prix de 16 % de RMSE sans dérive.

**3. Le filtre ne ment pas sur sa précision à la référence.** P₁₁/err² est entre 1,6 et 1,8 dans toutes les
conditions, aux deux budgets : légèrement sous-confiant. À l'ancien défaut, il tombait à 0,15.

**4. Il n'y a pas de retard de poursuite.** L'erreur maximale pendant la montée vaut 4,1 mM pour le Kalman contre
7,7 pour J/B à 1 000 photons : le filtre est meilleur pendant la transitoire, pas pire. Desserrer σ_w,g la fait
passer de 3,6 à 4,1 mM, sans changer la conclusion.

**5. 5a seul reste inutilisable à ces budgets, mais c'est une bonne ancre.** Il est 5,6 à 5,9 fois plus bruité que
J/B, et sa concentration inversée est biaisée de +3 mM par l'inversion elle-même. Son observable brute, celle que
reçoit le filtre, n'est biaisée que de +0,1 à +0,3 mM-équivalent, et elle est rigoureusement insensible à g (biais
de 3,31 / 3,30 / 3,23 / 3,16 mM pour 0 / 5 / 10 / 20 % de dérive).

### Réserves

- Les valeurs absolues dépendent des Langmuir, ajustées sur 8 points dont deux (5 et 10 mM) sortaient du mélange à
  deux états à l'étape 15. Les rapports entre lectures et le rôle de Q sont plus robustes que les mM.
- La dérive simulée démontre le mécanisme, pas l'amplitude réelle en tranche, qui viendra de l'expérience de
  profondeur. Le choix entre marche aléatoire et modèle à vitesse en dépend : il a été fait ici sur des rampes de
  5 à 20 % en 90 s et un saut de 5 %.
- La référence 3·10⁻³ vient du document, et le critère du correctif 4 la confirme sur la même course. Le balayage
  ayant été vu avant ce choix, il reste à la confirmer sur une autre trajectoire (une autre montée, une autre
  dérive).
- 8 répétitions par scénario : les biais de ±0,3 mM sont à environ 1σ de 0 (écarts-types dans `resume.csv`).

## Pour l'étape 2

Le générateur et les lectures ne bougent pas, et le filtre garde la marche aléatoire à σ_w,g = 3·10⁻³ (10⁻³ si le
budget descend vers 300 photons). Deux points restent ouverts, tous deux du côté de σ_w,C : l'excursion des
innovations pendant la montée et leur variance trop faible au repos (0,72). Un Q sur C qui dépend de la commande u,
large pendant la stimulation et étroit au repos, répond aux deux. Le modèle (g, ġ) est prêt dans `kalman.jl` si
les dérives mesurées en profondeur s'avèrent plus rapides que celles simulées ici.
