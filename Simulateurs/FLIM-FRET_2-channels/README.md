# Simulateur FLIM-FRET à deux canaux

Ce dossier implémente le plan `Simulateur FLIM-FRET à deux canaux — plan de conception.md` :
l'erreur sur [Cl⁻] (en mM) de chaque façon de lire le capteur, selon le budget de photons
et la concentration vraie, avec la borne de Cramér-Rao comme référence absolue.

- **Données et IRF** : celles de `tests/analyse_SCLM_2.jl`, soit la titration du 2026-07-20
  (`Maîtrise/Données/2026-07-20_SCLM`, 3 réplicats × 9 concentrations, canal bleu) et
  `irf_20260710_16x_sclm.sdt`, médiane soustraite comme le fait `load_irf_from_sdt`.
- **Estimateur 2** : **ton** `vec_to_lifetime` (la copie `lifetime_analysis.jl` de ce dossier),
  appelé tel quel. Le fichier n'est pas modifié ; il est chargé dans un module isolé par fil
  d'exécution.

```bash
cd Simulateurs/FLIM-FRET_2-channels
julia --project=. step0_spectres.jl             # spectres FPbase et matrice η (télécharge une fois)
julia --project=. step1_titration.jl            # ta titration → paramètres, et le test de M2 lent (20 s)
julia --project=. -t 8 step2_validation.jl      # 28 tests internes (40 s)
julia --project=. step3_calibration.jl          # calibrations, K apparents, identifiabilité, bornes (1 min)
julia --project=. -t 8 step4_monte_carlo.jl     # grille principale : 351 000 réalisations (≈ 10 min sur le M1)
julia --project=. -t 8 step5_nuisances.jl       # balayages de nuisance (≈ 45 min)
julia --project=. step6_figures.jl              # figures et résumé chiffré (2 min)
julia --project=. -t 8 step7_cible_mobile.jl    # cible mobile et boucle fermée (optionnel, ≈ 10 min)
julia --project=. step8_synthese.jl             # deux diapositives de synthèse : FR / EN, fond clair / noir (1 min)
# révision du 2026-10-15 (voir la section du même nom)
julia --project=. step9_batterie.jl             # batterie H1–H5 sur l'ingestion, IRF et t₀ sur tes données (1 min)
julia --project=. -t 8 step10_facteur8.jl       # le facteur 8 sur tes données amincies (1 min)
julia --project=. -t 8 step11_incertitude.jl    # incertitude des paramètres : bandes, variance paramétrique, classement (≈ 12 min)
julia --project=. step12_canal_jaune.jl         # canal jaune : η_YD depuis tes TIFF, ingestion testée, prédictions (1 min)
julia --project=. -t 8 step13_fusion.jl         # fusion ratio + forme sous dérive de gain (2 min)
```

L'étape 5 compte désormais aussi Φ_A effectif, τ_A ±20 %, le retard du jaune, le blanchiment et un
balayage « reference » (≈ 80 min en tout) ; les étapes 9 à 13 lisent les sorties des étapes 1 à 5.

## Figures de synthèse

Deux diapositives, 16:9, 2560 × 1440 px, comme tes slides, en français et en anglais, fond clair et
fond noir : `sorties/figures/synthese_{FR,EN}{,_dark}.png` (ce que dit le simulateur) et
`sorties/figures/synthese_revision_{FR,EN}{,_dark}.png` (ce que la révision du 2026-10-15 a
vérifié et changé). Les textes ci-dessous sont écrits pour la présentation.

### Diapositive 1 — ce que dit le simulateur

![Synthèse](sorties/figures/synthese_FR.png)

> Nous avons construit un simulateur de l'imagerie FLIM à deux canaux du capteur mCer3–Topaz.
> Il est calé sur la titration en solution du 2026-07-20, et il compare les façons de lire
> [Cl⁻] à budget de photons égal, avec la borne de Cramér-Rao comme référence.
>
> **A.** Un modèle à deux états décrit la titration. Une fraction x des capteurs transfère
> presque toute son énergie (τ_DA = 0,15 ns, E = 0,96) ; dans les autres, le donneur n'est pas
> éteint (τ_D0 = 3,69 ns). Le Cl⁻ ne change que x. Ces deux durées de vie, communes aux 27
> histogrammes, suffisent à les décrire, avec un K vrai de 30 mM.
>
> **B.** Sans ajustement supplémentaire, le modèle retrouve les K apparents mesurés :
> 18,7 mM pour le τ mono-exponentiel (mesuré : 18,5) et 40,6 mM pour le ratio des images
> (mesuré : 44,2).
>
> **C.** Parce que le transfert est presque total, l'état actif est quasi éteint dans le bleu,
> et sa signature temporelle tient sous la réponse instrumentale. La forme des déclins ne porte
> qu'environ 3 % de l'information sur [Cl⁻] ; le rapport jaune/bleu en porte presque tout.
>
> **D.** Le rapport jaune/bleu atteint la borne : σ ≈ 1 mM à 20 mM pour 10⁴ photons bleus. Les
> lectures par la forme sont 6 à 8 fois moins précises, la lecture actuelle (τ
> mono-exponentiel + Hill) environ 60 fois.
>
> **E.** Le ratio dépend toutefois de l'échelle relative des deux canaux : 5 % de dérive de gain
> le biaisent de 3 à 4 mM, alors que la forme y est insensible. À l'inverse, un τ décalé de 5 %
> entre la solution et la tranche biaise la forme de 7 à 15 mM et laisse le ratio intact.
>
> **F.** En boucle fermée, l'intégrateur moyenne le bruit et il ne reste que le biais. Il est
> presque nul avec le ratio, de −0,75 mM avec la forme et de −1,6 mM avec la lecture actuelle
> (−7,5 mM à 10³ photons).
>
> **À retenir.** Pour ce capteur, l'information est dans le rapport des deux canaux FLIM. Mais
> ce rapport n'est exact que si son échelle est étalonnée et suivie : la diapositive suivante
> montre comment la forme des déclins peut la surveiller en continu. Les barres et bandes grises
> donnent l'incertitude due à l'IRF et au décalage t₀.

### Diapositive 2 — la révision : ce qui tient, ce qui change

![Révision](sorties/figures/synthese_revision_FR.png)

> Nous avons ensuite cherché à mettre ce simulateur en défaut, par ordre de ce qui pouvait
> invalider ses conclusions.
>
> **A.** Le paramètre dont tout dépend est la durée de vie courte, τ_DA = 0,15 ns. Les données
> réelles la préfèrent nettement à toute valeur supérieure à 0,3 ns. Une titration simulée à 1 ns
> est retrouvée à 1,00 ns, et à 0,94–1,03 ns avec une IRF fausse de 10 % : une erreur d'IRF ne
> fabrique pas un transfert quasi total.
>
> **B.** Les acquisitions sont à temps fixé, donc l'intensité totale est une prédiction
> indépendante. Le modèle prévoit une hausse de 62 % entre 0 mM et saturation ; elle est de 61 %.
>
> **C.** Sur les histogrammes réels amincis, remplacer la lecture actuelle (τ mono-exponentiel +
> Hill) par un bi-exponentiel à durées de vie figées divise l'erreur par 7,3 à 10⁴ photons, comme
> prévu (6,9).
>
> **D.** En réajustant la titration avec une IRF plus large ou plus étroite de 10 %, ou un t₀
> décalé, les précisions changent, mais jamais l'ordre des lectures. L'incertitude des paramètres
> ajustés ajoute environ 1 mM à 20 mM, pour toutes les lectures.
>
> **E.** Si le gain relatif des détecteurs dérive, le ratio seul dérive avec lui. Un filtre de
> Kalman qui estime ce gain à partir de la forme des déclins garde l'erreur sous 0,8 mM à 10⁴
> photons par mesure, que la dérive se fasse en une minute ou en une heure.
>
> **F.** Les sensibilités sont complémentaires. Le gain et une fuite mal connue biaisent le ratio ;
> un décalage de quelques picosecondes entre détecteurs biaise la forme à deux canaux ; un τ décalé
> entre la solution et la tranche biaise les lectures de forme. Selon son produit, le blanchiment
> de l'accepteur touche toutes les lectures ou seulement le ratio.
>
> **À retenir.** Le modèle tient et ses conclusions relatives aussi. Les valeurs absolues
> demandent maintenant une mesure : celle du canal jaune à trois concentrations, qui donnera la
> fuite, le rendement de l'accepteur et le retard entre détecteurs. C'est aussi le seul test qui
> distingue un accepteur émetteur d'un accepteur sombre, que le canal bleu ne sépare pas.

Tous les paramètres sont dans `parametres.toml`, chacun avec sa provenance. Pour simuler
SCLM plutôt que mCer3_Topaz, il suffit d'y écrire `proteine = "SCLM"`, `donneur = "cerulean"`,
`Phi_D = 0.62`, puis de relancer à partir de l'étape 2.

## En bref

1. **M2 lent décrit ta titration mCer3_Topaz au premier ordre.** τ_D0 = 3,69 ns et τ_DA = 0,15 ns,
   liés sur les 27 histogrammes : χ²ᵣ = 1,96, contre 1,58 si on libère les τ. Le χ²ᵣ de ton
   ajustement mono, point par point, est celui que M2 prédit (21 contre 19 à 0 mM,
   1,35 contre 1,37 à « F »). L'amplitude du donneur à t = 0 reste constante (±5 %) ; les acquisitions
   étant à temps fixé, les comptages montent comme M2 le prédit (×1,61 contre ×1,62). Les phasors
   restent sur le segment. SCLM et M2-SCLM passent moins bien (§ 1).
2. **Le FRET actif est presque total : τ_DA = 0,15 ns, E = 0,96.** Le plan supposait τ_DA ≈ 2,1 ns.
   L'état actif est donc quasi éteint dans le bleu (3 % des photons bleus à 0 mM, dans un
   pic sous l'IRF), et l'émission sensibilisée (τ_A = 3,0 ns) ressemble au donneur qui fuit
   (τ_D0 = 3,69 ns). **La forme des déclins ne contient presque rien sur x ; le rapport jaune/bleu,
   qui varie d'un facteur 2,3 entre 0 et 150 mM, contient presque tout.**
3. **Variance : le ratio est déjà à la borne.** À 20 mM et 10⁴ photons (normalisation A), le ratio
   donne σ = 0,98 mM pour une borne de 0,93 mM. L'ajustement conjoint à échelle commune (5d) fait
   0,97 mM : la forme n'ajoute que 1 à 3 %. L'ajustement conjoint du plan (5a, échelles bleue et
   jaune libres) donne 5,6 mM, le bi-exponentiel à τ figés 8,0 mM et ton mono 59 mM.
4. **Ton pipeline actuel (mono + Hill) est le maillon faible.** τ mono ne varie que de 0,30 ns
   sur 0–150 mM, et ton ajustement libère le fond et le décalage d'IRF, ce qui triple σ_τ
   (0,157 contre 0,053 ns à 10⁴ photons, d'après la variante 2F, où ils sont figés). À bas photons,
   il surestime en outre τ de 0,2 à 0,3 ns. Jusqu'à 3·10³ photons, 74 à 89 %
   des estimations saturent à 150 mM. À 10⁴ photons, 34 % saturent encore (σ robuste 31 mM) ;
   à 3·10⁴, 3 % (σ robuste 13 mM). Voir § 6. Sur tes propres histogrammes amincis à 10⁴ photons, le
   bi-exponentiel à τ figés fait 7,3 fois mieux (révision, facteur 8).
5. **Biais : la forme achète de la robustesse sur l'intensité, pas sur les durées de vie.**
   - Une fuite η_YD fausse de 25 % biaise le ratio de −7 à +10 mM, et 5a de −4 à +5 mM.
   - Une dérive de gain du jaune ne biaise pas 5a, mais biaise le ratio d'environ 0,7 mM par %.
     5a ne devient préférable qu'au-delà de ±8 % d'incertitude sur ce gain (à 10⁴ photons).
   - Un τ_D0 faux de 5 % (l'indice de réfraction du cytoplasme, ≈ 1,37 contre 1,33 en solution,
     suffit à le produire) laisse le ratio à moins de 0,4 mM de biais. Les estimateurs à durées
     de vie figées (5b, 5c, 4) dérivent en revanche de 7 à 15 mM, et de 14 à 44 mM pour un écart
     de 10 %.
6. **K apparents : le modèle retrouve tes deux ajustements de Hill.** Avec K vrai = 30,3 mM (isotherme
   sur x), il prédit 18,7 mM pour ton τ mono (tu as mesuré 18,5 mM) et 40,6 mM pour B/J
   (ton ratio TIFF C2/C1 : 44,2 mM). Le τ mono concorde ; l'écart sur B/J disparaît si la fuite
   η_YD vaut 0,7 à 0,8 fois sa valeur des spectres (révision, changement 3). C2/C1 est bien un
   rapport bleu/jaune (étiquettes inversées).
7. **En boucle fermée à 20 mM**, à 10⁴ photons par échantillon à 10 Hz, l'erreur statique vaut
   0,0 mM avec le ratio ou 5d, −0,75 mM avec 5a et −1,6 mM avec ton pipeline (−7,5 mM à 10³ photons).

Ces chiffres comparent des estimateurs **sous un modèle contrôlé**. La révision du 2026-10-15
(section suivante) teste ce modèle contre tes données et ajoute le photoblanchiment, le décalage
entre détecteurs et l'incertitude des paramètres. Le classement tient ; les valeurs absolues
attendent la mesure du canal jaune.

## Révision du 2026-10-15 : les changements de la note

La note `Simulateur FLIM-FRET changements à apporter.pdf` ordonne les changements par ce qu'ils
peuvent invalider. Tout ce qui se fait sans banc est fait. Le canal jaune est prêt à être ingéré :
le code est testé sur des bancs simulés et attend les fichiers .sdt des deux canaux.

```bash
julia --project=. step1_titration.jl            # + le test d'intensité et le mode d'acquisition des .sdt
julia --project=. -t 8 step2_validation.jl      # 28 tests (dont 8 nouveaux)
julia --project=. -t 8 step5_nuisances.jl tau_A_relatif Phi_A_relatif decalage_jaune_ps \
      blanchiment_accepteur blanchiment_donneur blanchiment_accepteur_sombre reference   # ≈ 40 min
julia --project=. step9_batterie.jl             # batterie H1–H5, IRF et t₀ sur tes données (1 min)
julia --project=. -t 8 step10_facteur8.jl       # le facteur 8 sur tes données amincies (1 min)
julia --project=. -t 8 step11_incertitude.jl    # mondes IRF/t₀, variance paramétrique, classement (≈ 12 min)
julia --project=. step12_canal_jaune.jl         # η_YD depuis tes TIFF, ingestion du jaune, prédictions (1 min)
julia --project=. -t 8 step13_fusion.jl         # fusion ratio + forme sous dérive de gain (2 min)
julia --project=. step6_figures.jl && julia --project=. step8_synthese.jl
```

### La porte : « intensité constante » voulait dire amplitude constante, et le test passe

- **Les 27 acquisitions sont à temps fixé.** Chaque .sdt porte un temps de collecte de 10 s et
  l'arrêt sur temps (`col_t`, `stopt`, lus par `lire_acquisition`). L'intensité totale garde donc
  son information.
- La mention « intensité constante » de l'ancienne figure portait sur l'**amplitude du donneur à
  t = 0** (le A de l'ajustement), pas sur l'intensité totale. Le libellé est corrigé.
- **Le test :** les comptages nets (fond retiré) montent de **×1,609** entre 0 mM et « F ».
  M2 prédit (1 − E·x_sat)/(1 − E·x₀) = **×1,619**. La courbe complète est dans le panneau en haut
  à droite de `titration_mCer3_Topaz.png`. Le modèle à deux états n'est pas falsifié sur ce point.
- SCLM et M2-SCLM échouent (×2,32 et ×2,30 mesurés, ×1,40 et ×1,45 prédits), ce qui confirme le
  choix de mCer3_Topaz.

### Changement 1 — la batterie H1–H5 (`step9_batterie.jl`, fig. 9)

Chaque jeu est une titration complète simulée avec la même IRF, les mêmes photons et fonds par
histogramme et la même grille de concentrations que le 2026-07-20. Elle passe dans le même
ajustement global que tes données. « Biais du ratio » désigne la lecture par J/B calibrée sur les
paramètres retrouvés, appliquée au monde vrai, à 20 mM. Moyenne de 3 tirages ; bruit de Poisson
sur χ²ᵣ : ± 0,018.

| Jeu | Hypothèse générée | χ²ᵣ | τ_DA retrouvé (vrai) | x₀ retrouvé (vrai) | Biais du ratio | Ce que ça dit |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| H1 | modèle actuel | 0,99 | 0,152 (0,151) | 0,426 (0,429) | −0,6 mM | l'ingestion retrouve ses paramètres |
| H2 | transfert vers un accepteur non émissif | 0,99 | 0,151 (0,151) | 0,430 (0,429) | hors plage | **indiscernable par le canal bleu** |
| H2b | complexe statique non émissif | 0,99 | — | ≈ 0 (0) | — | exclu : pas de déclin court, alors que tes données en ont un |
| H3a | 20 % de x à τ = 1,5 ns | 1,29 | 0,223 | 0,348 (0,429) | −8,6 mM | discernable ; reproduit les signatures de tes données |
| H3b | 10 % de donneur nu (4,1 ns) | 0,99 | 0,151 | 0,384 (0,386) | +0,3 mM | indiscernable, sans conséquence (τ_D0 absorbe : 3,735 ns) |
| H4 | τ_DA = 1,0 ns, x × 1,32 (même x·E) | 0,99 | 1,000 (1,000) | 0,561 (0,564) | −0,2 mM | retrouvé |
| H5 | H4 ajusté avec l'IRF ×0,9 / ×1,1 | 1,13 / 1,81 | 1,03 / 0,94 | 0,559 / 0,568 | −1,5 / +3,1 mM | une IRF fausse ne fabrique pas E = 0,96 |
| — | H1 ajusté avec l'IRF ×0,9 / ×1,1 | 1,02 / 1,32 | 0,165 / 0,115 | 0,393 / 0,517 | −6,1 / +28 mM | IRF trop étroite invisible au χ²ᵣ |
| — | H1 avec la gigue de t₀ observée | 1,15 | 0,143 | 0,426 | +4,4 mM | explique ≈ 0,15 du χ²ᵣ |

**Ce que disent les résultats**

1. **τ_DA = 0,15 ns tient.** Sur tes données, χ²ᵣ vaut 1,96 à 0,15 ns, 2,42 à 0,3 ns et 5,46 à
   1 ns (Δχ² ≈ 21 000). Simulée à 1 ns, la titration est retrouvée à 1,000 ns. Ajustée avec une IRF
   fausse de 10 %, elle donne 0,94 à 1,03 ns : une erreur d'IRF ne fabrique pas un transfert quasi
   total.
2. **Le chiffre clé : x₀ sous perturbation de l'IRF.** Tes données ne distinguent pas une IRF de
   5 à 15 % plus étroite que la mesurée (χ²ᵣ de 1,958 à 1,971). Elles rejettent une IRF plus large
   (×1,02 : 2,10 ; ×1,1 : 2,76). Dans la plage permise, x₀ va de 0,363 à 0,429 (−15 % à 0) et x₀·E
   de 0,344 à 0,411. Décaler t₀ d'un canal entier (48,8 ps) est rejeté sans appel (χ²ᵣ de 20 à 22).
   La barre d'erreur systématique du signal est donc dissymétrique : de −0,07 à 0 sur x₀.
3. **Ce que vaut χ²ᵣ = 1,96.** Ni l'IRF (×0,9 ne change pas le χ²ᵣ de tes données, ×1,1
   l'aggrave) ni un donneur nu (H3b, invisible) ne l'expliquent. Deux contributions sont
   identifiées.
   - **La gigue de t₀ entre échantillons** vaut environ 0,15 à 0,17 (voir ci-dessous).
   - **Un état FRET intermédiaire** (H3a : 20 à 40 % de la population active vers 1,5 à 2,5 ns)
     donne de 1,1 à 1,65. H3a reproduit en plus deux signatures que l'IRF ne reproduit pas.
     D'abord l'écart du τ mono de ton code entre les données et le modèle ajusté : −0,11 ns à 0 mM
     et +0,07 ns à saturation, contre −0,08 et +0,07 ns pour H3a. Le modèle à deux états
     sous-estime d'un tiers la plage du τ mono (0,31 ns contre 0,52 ns mesurés). Ensuite la
     distance des phasors au segment : 0,11 %, contre 0,12 % pour H3a.

   Avec un état intermédiaire et la gigue (« H3a 30 % + gigue »), χ²ᵣ = 1,78.
4. **Le mécanisme à nommer dans toute figure : H2.** Un transfert vers un accepteur qui n'émet pas
   donne exactement les mêmes 27 histogrammes bleus. Seul le canal jaune le tranche : sous H2,
   J/B ne varie pas avec [Cl⁻] (fig. 12).
5. **Conséquence pour les valeurs absolues.** Sous H3a, une lecture par le ratio calibrée sur les
   paramètres ingérés serait biaisée de −5 à −13 mM à 20 mM. Le classement des estimateurs ne
   dépend pas de ces paramètres ; leurs valeurs absolues, si. Une calibration empirique du ratio
   (titration à deux canaux) supprime ce biais.

**t₀ libre par histogramme** (le diagnostic du changement 2, panneau en bas à droite de la fig. 9)

- **Pas de tendance avec [Cl⁻]** : ρ de Spearman = +0,35, p ≈ 0,33 par permutation des
  9 concentrations. τ_DA n'est pas contaminé par un artefact corrélé au signal.
- **Mais une vraie gigue** : 3,6 ps entre échantillons, contre 0,6 ps entre réplicats, ce qui est le
  bruit statistique. Libérer t₀ fait passer χ²ᵣ de 1,96 à 1,79.
- **Effet sur x.** La gigue déplace le x de chaque concentration de ±0,04. Avec t₀ libre, les
  deux non-monotonies de x (5 mM, 146 mM) disparaissent, et l'isotherme passe de K = 30,3 mM,
  n = 1,46 à K = 24,8 mM, n = 1,13.
- **Mais ces paramètres prédisent moins bien tes K apparents** (τ mono : 16,0 mM contre 18,5
  mesuré ; B/J : 35,5 contre 44,2).
- Les deux traitements encadrent donc K vrai. Le t₀ commun reste le nominal, et le monde
  « t₀ par histogramme » entre dans les bandes du changement 2.
- **À l'acquisition :** enregistrer une IRF par échantillon, ou libérer t₀ par échantillon dans
  l'ingestion.

### Le facteur 8 sur tes données (`step10_facteur8.jl`, fig. 10)

Tes 27 histogrammes (≈ 10⁶ photons) sont amincis par tirage binomial, 200 tirages par histogramme et
par budget. Deux lectures sont comparées : 4R, bi-exponentiel à τ_D0, τ_DA et t₀ figés (fond libre,
calibré par l'isotherme sur les histogrammes complets), et 2H, ton pipeline (mono + Hill, calibré
comme le tien). Le jumeau refait tout sur des tirages du modèle ajusté. σ est le bruit de photons
seul, à 20 mM :

| N | 4R réel | 4R jumeau | 4 simulateur | 2H réel | 2H jumeau | 2H simulateur | facteur réel / jumeau / prédit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10³ | 47,2 | 45,8 | 47,5 | 59,4 | 61,9 | 54,5 | 1,3 / 1,4 / 1,1 |
| 3·10³ | 23,6 | 25,4 | 27,3 | 63,3 | 62,7 | 61,5 | 2,7 / 2,5 / 2,3 |
| 10⁴ | 7,3 | 7,8 | 8,7 | 53,4 | 60,2 | 59,7 | **7,3 / 7,8 / 6,9** |
| 3·10⁴ | 3,9 | 4,2 | 4,3 | 11,1 | 22,5 | 26,8 | 2,8 / 5,4 / 6,3 |

**Ce que disent les résultats**

- **Le facteur se matérialise à 10⁴ photons** : 7,3 sur tes données, 7,8 sur le jumeau, 6,9 prédit.
- **Le bi-exponentiel à τ figés se comporte sur tes données exactement comme prévu.** À 10⁴ et
  3·10⁴ photons, son σ réel vaut 0,93 à 0,94 fois celui du jumeau à 20 mM, et 0,92 à 0,99 fois à
  5 mM. Le modèle à deux états décrit donc aussi le bruit de tes données, pas seulement leur forme
  moyenne.
- **Le facteur vient surtout de la saturation du mono + Hill.** Sur tes données à 10⁴ photons, 21 %
  des lectures sont ramenées à 150 mM (36 % pour le jumeau). À 3·10⁴ photons, ton mono fait
  mieux que prévu (11 mM contre 22 à 27) : sur tes données, il voit un τ qui varie de 0,52 ns entre
  0 mM et saturation, contre 0,31 à 0,34 ns prédits. C'est la signature de H3a.
- **Gain immédiat, sans changer d'acquisition** : remplacer mono + Hill par 4R divise σ par
  environ 7 à 10⁴ photons et par environ 3 à 3·10⁴.

### Changement 2 — l'incertitude des paramètres (`step11_incertitude.jl`, fig. 1, 2, 5, 11)

- **L'IRF et t₀ : des bandes sur les figures 1, 2 et 5.** La titration est réajustée avec l'IRF
  ×0,9 et ×1,1, t₀ ± 3,5 ps (la gigue observée, plutôt qu'un canal entier, que les données rejettent)
  et t₀ libre par histogramme. Chaque réajustement définit un monde cohérent : ses paramètres, son
  IRF, et la même IRF donnée à ton code mono. On y refait les cellules des figures 1 et 2 et les
  K apparents de la figure 5.
  - À 20 mM et 10⁴ photons, σ vaut 0,91 à 1,26 mM pour le ratio (0,98 au nominal), 5,6 à 7,4 pour
    5a, 7,7 à 13,5 pour 4 et 54 à 60 pour 2H. **Le classement ne change dans aucun monde.**
  - Les K apparents prédits pour le τ mono vont de 16,0 à 22,6 mM (18,5 mesuré). Ceux de B/J vont de
    35,5 à 41,4 mM : la valeur mesurée (44,2) sort de la plage, ce que le changement 3 explique par
    η_YD.
- **La covariance des ajustements : la variance paramétrique.** Les incertitudes ajustées valent
  σ(τ_D0) = 0,003 ns et σ(τ_DA) = 0,002 ns (covariance multipliée par χ²ᵣ), puis σ(x₀) = 0,008,
  σ(x_sat) = 0,010, σ(K) = 2,4 mM et σ(n) = 0,14. On tire 200 vecteurs, que chaque estimateur,
  calibré sur le point estimé, lit sans bruit.
  - Le σ paramétrique vaut environ 1,1 mM à 20 mM (0,8 à 5–12 mM, 2,3 à 40, 4,3 à 60), **le même
    pour tous les estimateurs**, puisque c'est la calibration x → C qui est incertaine.
  - Pour le ratio à 10⁴ photons, il est comparable au bruit de photons (0,98) : σ total = 1,48 mM.
    Pour les lectures de forme, il est négligeable.
  - Vérification par un vrai Monte-Carlo sur 40 vecteurs : pour le ratio, 1,97 mM² obtenus contre
    2,19 prédits (l'estimation à 40 vecteurs est à ±23 % près).
  - `resultats_principaux.csv` porte désormais des colonnes séparées : `variance_photonique_mM2`,
    `variance_parametrique_mM2`, `biais_parametrique_moyen_mM`, `sigma_total_mM`.
- **Les paramètres supposés** (τ_A ±20 %, Φ_A ±30 %, η_YD ±50 % ; fig. 11, `classement.csv`).
  - Le classement en précision survit : τ de Kendall moyen de 0,99 en calibration exacte, et le
    meilleur estimateur est inchangé dans 72 cas sur 72 (5, 20 et 40 mM ; 10³ et 10⁴ photons).
  - En revanche, si l'on garde la valeur supposée alors qu'elle est fausse, le meilleur change dans
    17 cas sur 72. Un η_YD ou un Φ_A faux de 15 à 50 % porte la RMSE du ratio calibré sur le modèle à
    5 à 30 mM à 20 mM (jusqu'à 110 mM à 40 mM), et les lectures de forme passent devant.
  - L'avantage du ratio suppose donc sa calibration mesurée, c'est le changement 3.

### Changement 3 — le canal jaune (`step12_canal_jaune.jl`, `src/Jaune.jl`, fig. 12)

- **`parametres.toml` à trois statuts.** Chaque paramètre est « ajusté », « supposé » ou
  « mesuré » ; le chargeur refuse un « mesuré » sans date ni incertitude, et `decrire_provenance`
  les liste.
- **η_YD tiré de tes images, sans banc.** Le K apparent de B/J (44,2 mM mesuré contre 40,6 prédit)
  donne η_YD × **0,80** (0,62 à 0,83 selon les mondes IRF/t₀). La plage dynamique de B/J (×2,99
  mesurée contre ×2,48 prédite) donne indépendamment η_YD × **0,71** (0,52 à 0,75). Les deux
  inversions placent la fuite 20 à 30 % sous la valeur des spectres (0,130 et 0,116 contre 0,164),
  et la sortie est écrite dans `sorties/calibration/eta_YD_inverse.toml`.
  - Réserve : tes TIFF Bliq ne sont pas photoniques, et un décalage d'intensité déplacerait les deux
    nombres. La valeur nominale n'est donc pas remplacée.
  - Si la fuite vaut vraiment 0,75 fois la valeur supposée, un ratio calibré sur le modèle est biaisé
    de +10 mM à 20 mM (fig. 11).
- **L'ingestion du jaune est prête : `ajuster_bleu_jaune`.** C'est un ajustement conjoint du bleu
  et du jaune aux trois concentrations, avec τ_D0 et τ_DA de la titration. Il rend
  ρ_D = η_YD/η_BD, ρ_A = η_YA·Φ_A/(η_BD·Φ_D) (d'où Φ_A effectif), τ_A, ε (borné à [0 ; 0,15]), le
  retard Δt du jaune et les fonds.
  - Aller-retour sur 20 bancs simulés dont les valeurs vraies diffèrent des supposées (η_YD × 0,8,
    τ_A = 2,8 ns, Φ_A = 0,50, ε = 0,04, Δt = 30 ps) : sans biais, écarts normés de 0 ± 1.
  - Précision à 10⁶ photons bleus par concentration : ρ_D à ±1,3 %, ρ_A à ±1,6 %, τ_A à
    ±0,011 ns, Δt à ±0,7 ps. ε n'est connu qu'à ±0,016, donc seulement borné.
  - L'ajustement en deux temps (bleu, puis jaune) sous-estime l'erreur sur ρ_D d'un facteur 12 : il
    ne faut pas l'utiliser.
  - `parametres_mesures` écrit le bloc TOML au statut « mesuré » (exemple dans
    `sorties/canal_jaune/parametres_mesures_exemple.toml`).
- **Ce que la mesure confrontera.**
  - **J/B (gains égaux)** vaut 1,36 à 0 mM, 0,85 à 30 mM et 0,55 à saturation, soit une plage de
    ×2,48. L'enveloppe des paramètres supposés va de ×1,72 à ×4,32. Sous H2, J/B reste constant à
    0,49 (×1,00), et sous H3a la plage vaut ×2,30. La plage ne dépend pas du gain relatif des
    détecteurs : c'est le test qui sépare M2 de H2.
  - **Le jaune à saturation** est presque pure fuite : il doit décroître en ⟨τ⟩ = 3,62 ns, proche de
    τ_D0 et non de τ_A, contre 3,30 ns à 0 mM.
  - **La montée** de 0,15 ns tient dans la largeur de l'IRF (130 ps à mi-hauteur), mais elle n'est
    pas statistiquement inobservable. Le jaune seul à 0 mM la mesure à ±0,018 ns avec 10⁵ photons
    jaunes (±0,024 si Δt est inconnu), et à ±0,006 ns avec 10⁶. Une montée ≥ 0,5 ns (H4 : 1 ns,
    ±0,13 ns) serait sans ambiguïté. La limite pratique est la connaissance de l'IRF du jaune, à
    mesurer au diffuseur.

### Changement 4 — fusionner le ratio et la forme (`step13_fusion.jl`, `src/Fusion.jl`, fig. 7b)

**Ce qui a été ajouté**

- L'**estimateur 7** (`DeuxLectures`) rend les deux lectures, ratio et forme 5a, avec leurs
  variances de Fisher.
- La **dérive de gain** du jaune est simulée dans le générateur (`derive_gain`) : une marche
  aléatoire sur log g, d'écart-type 5 % après τ_g.
- Un **filtre de Kalman étendu** suit l'état ([Cl⁻], log g). La mesure z₁ = log(N_J/N_B) a pour
  modèle h₁ = log g + log R(C), et la mesure z₂ = x̂ de 5a a pour modèle h₂ = x(C).

**Deux leçons de mise au point**

- **La forme doit entrer en x, pas en mM.** Inverser l'isotherme à N fini ajoute un biais de Jensen
  (−1,5 mM à 10⁴ photons, −4,8 à 10³), et la fusion en héritait.
- **Le bruit de processus de [Cl⁻] doit dominer celui du gain**, en unités de log R. Sinon le
  filtre attribue au gain toute variation du ratio et jette sa précision. q_C correspond ici à
  1 mM/s.

**Résultats** ([Cl⁻] = 20 + 8·sin(2πt/5 min) mM, 10 Hz pendant 30 min ; RMSE en mM sur la
seconde moitié, moyenne de 6 réalisations) :

| τ_g de la dérive | ratio seul, 10⁴ | forme seule, 10⁴ | **fusion, 10⁴** | gain : dérive / erreur de la fusion | fusion, 10³ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 min | 22,2 | 0,84 | **0,77** | 9,8 % / 1,2 % | 2,07 |
| 5 min | 5,2 | 0,86 | **0,68** | 3,8 % / 0,9 % | 1,69 |
| 20 min | 1,6 | 0,82 | **0,51** | 1,8 % / 0,7 % | 1,07 |
| 1 h | 1,9 | 0,86 | **0,42** | 1,5 % / 0,5 % | 0,93 |

- **La fusion bat les deux lectures seules dans tous les cas.** Après un saut de gain de 5 % à 10⁴
  photons, le ratio seul reste biaisé de −3,1 mM. La fusion revient à 1/e de son pic en **13 s**
  (130 mesures) pour mCer3_Topaz, et en 9 s pour le capteur hypothétique à E = 0,43. L'estimation de
  la note (40 à 60 mesures) est optimiste d'un facteur 2 à 3 pour ce réglage. Toute dérive plus lente
  que la minute est néanmoins observée et corrigée.
- **La forme seule, filtrée, devient compétitive** dès que le gain dérive de plus de 1 à 2 % : la
  moyenne glissante de 5a, sans biais en x, récupère la précision. Pour le capteur à E = 0,43, la
  fusion n'apporte presque rien sur la forme seule (0,56 à 0,69 mM contre 0,70 à 0,74), ce qui
  confirme la note : avec un transfert moins total, la forme suffit.

### Changement 5 — les deux biais exclus (`step5_nuisances.jl`, fig. 13, 14)

- **Les modèles.** Le blanchiment fait décroître exponentiellement, pendant l'acquisition, le
  nombre de donneurs et d'accepteurs intacts, indépendamment, β étant la fraction blanchie en fin
  d'acquisition. L'histogramme accumulé est celui des populations moyennées. Pour l'accepteur
  blanchi, deux produits sont simulés. Soit il n'absorbe plus : le donneur devient nu, à 4,1 ns. Soit
  il reste extincteur sans émettre : le bleu est inchangé. Le retard Δt du détecteur jaune va de
  −200 à +200 ps. Les biais sont **ajoutés**, c'est-à-dire mesurés contre le monde nominal avec les
  mêmes germes (balayage « reference »).

Biais ajouté à 20 mM et 10⁴ photons (mM) :

| Perturbation | 1 ratio | 5d | 5a | 4 bi-exp. | 2H |
| --- | ---: | ---: | ---: | ---: | ---: |
| gain du jaune +5 % | −3,1 | −3,0 | 0,0 | 0 | 0 |
| retard du jaune +20 ps | 0,0 | −0,1 | **−3,8** | 0 | 0 |
| 10 % d'accepteurs blanchis, non absorbants | +2,7 | +2,8 | +5,0 | +5,3 | +9,5 |
| 10 % d'accepteurs blanchis, extincteurs sombres | **+1,8** | +1,7 | +1,0 | 0,0 | ≈ 0 |
| 10 % de donneurs blanchis | −0,1 | −0,1 | +0,3 | +0,2 | +0,6 |
| η_YD réel = 0,75 × supposé | **+10,3** | +9,8 | −3,8 | 0 | 0 |
| τ_D0 en tranche +5 % | −0,2 | +0,1 | +13,0 | **+15,2** | +60 |

- **Δt.** Le ratio, 5d et les estimateurs à un canal n'y sont pas sensibles : moins de 0,07 mM pour
  le ratio à 200 ps, et 1 mM pour 5d au-delà de 130 ps. Les ajustements de forme à deux canaux le
  lisent comme de l'amplitude : 2 à 3 mM par 10 ps pour 5a et 5c.
  - **Tolérance à exiger du banc : ±4 à 5 ps si la fusion utilise 5a** (±3 ps pour 5c).
    L'ajustement conjoint du jaune mesure Δt à ±0,7 ps, et un diffuseur vu par les deux canaux
    devrait faire aussi bien.
- **Blanchiment : tout dépend du produit.**
  - *S'il n'absorbe plus*, le donneur se déséteint et x baisse pour tous. La forme est alors **plus**
    biaisée que le ratio (5,3 contre 2,7 mM par 10 %) : elle n'est pas protégée.
  - *S'il reste extincteur sombre*, le bleu ne change pas : 4, 3 et le phasor sont immunisés, et seul
    le ratio (1,8 mM par 10 %) et, à moitié, 5a sont touchés. C'est le cas où le blanchiment joue
    contre la recommandation actuelle, comme le craignait la note, et c'est aussi ce que suggèrent
    tes tranches (ratio instable, durée de vie stable).
  - Le blanchiment du donneur est négligeable pour tous.
- **Mesure décisive :** la durée de vie du donneur sous illumination répétée à 594 nm. Si elle monte,
  l'accepteur blanchi n'absorbe plus ; si elle reste stable, c'est un extincteur sombre, et la baisse
  de J/B se comporte comme une dérive de gain que la fusion peut suivre (à vérifier en simulation).

### Ce que la révision change dans la conclusion

1. **Le cœur du modèle tient.** τ_DA = 0,15 ns survit à la batterie, et le test d'intensité passe.
   L'excès de χ²ᵣ s'explique surtout par la gigue de t₀ entre échantillons et par une hétérogénéité
   du FRET (un état intermédiaire), pas par l'IRF.
2. **Le facteur 8 est réel** à 10⁴ photons (7,3 sur tes données). Recommandation immédiate :
   remplacer mono + Hill par 4R dans le pipeline actuel.
3. **Les valeurs absolues d'un ratio calibré sur le modèle ne sont pas défendables.** η_YD est
   probablement 20 à 30 % sous les spectres, et H3a ainsi que l'IRF pèsent aussi. Seule une
   calibration du ratio mesurée au banc l'est. Le classement des estimateurs, lui, est robuste.
4. **La recommandation devient la fusion** (ratio et forme, gain dans l'état du filtre) plutôt que
   « ratio, en surveillant le gain ».
5. **Le banc, par priorité :**
   1. le canal jaune à trois concentrations : la plage de J/B tranche entre M2 et H2, et la mesure
      donne η_YD, Φ_A effectif, τ_A et Δt ;
   2. un diffuseur vu par les deux canaux : Δt à ±3–5 ps, et l'IRF du jaune ;
   3. la durée de vie du donneur sous blanchiment de l'accepteur ;
   4. une IRF par échantillon, ou t₀ libre par échantillon.

### Écarts à la note

| Note | Ici | Pourquoi |
| --- | --- | --- |
| les changements dans les modules existants | nouveaux fichiers `Batterie.jl`, `Reel.jl`, `Incertitude.jl`, `Jaune.jl`, `Fusion.jl` et étapes 9 à 13 ; les modules existants reçoivent les briques (Φ_A effectif, Δt du jaune, IRF dilatée, états du donneur, blanchiment, dérive de gain, t₀ par histogramme) | lisibilité ; aucune refonte |
| t₀ décalé d'un canal (48,8 ps) | calculé, mais exclu des bandes ; remplacé par ± la gigue observée (3,5 ps) | χ²ᵣ de 20 à 22 : un décalage d'un canal n'est pas un monde plausible |
| inversion de η_YD dans `Ingest.jl` | dans `Jaune.jl`, sortie dans `sorties/calibration/` | elle utilise l'interface des estimateurs, définie après `Ingest.jl` |
| blanchiment par fluorophore | deux produits possibles pour l'accepteur blanchi | ils mènent à des conclusions opposées |
| la forme comme seconde mesure du filtre | lue en x, avec h₂ = x(C) | biais de Jensen sinon |
| facteur 8 : bi-exponentiel à τ figés | 4R, avec un fond libre | tes histogrammes ont ≈ 1 % de fond |
| — | une deuxième inversion de η_YD, par la plage dynamique du TIFF | indépendante de la première, elle la confirme |
| — | des signatures en plus du χ²ᵣ dans la batterie (τ mono, phasors) | le χ²ᵣ seul ne distingue pas l'IRF d'un état intermédiaire |


## 1. Ta titration et le test de M2 lent (`step1_titration.jl`)

Ajustement global du canal bleu, sur les canaux 14–242 (ton gating ; les canaux 1–13 et 245–256 de tes .sdt
sont vides, limites du TAC). Le modèle de chaque histogramme est
μ = A·[x·B(τ_DA) + (1 − x)·B(τ_D0)] + fond, où B est l'IRF convoluée avec le noyau du § 4.
τ_D0, τ_DA et le décalage t₀ sont liés sur les 27 histogrammes ; x, A et le fond sont libres pour
chacun.

| | mCer3_Topaz | SCLM | M2-SCLM |
| --- | ---: | ---: | ---: |
| τ_D0 (ns) | 3,689 ± 0,002 | 2,824 ± 0,003 | 2,749 ± 0,003 |
| τ_DA (ns) | 0,151 ± 0,001 | 0,426 ± 0,002 | 0,375 ± 0,002 |
| χ²ᵣ, τ liés (27 histogrammes) | **1,96** | 2,96 | 2,71 |
| χ²ᵣ, τ libres par histogramme | 1,58 | 2,63 | 2,38 |
| χ²ᵣ de ton mono (moyenne) | 7,5 | 19,3 | 18,4 |
| intensité A du donneur selon [Cl⁻] | constante (CV 5 %, sans tendance) | monte de 0,70 à 1,2 | monte de 0,74 à 1,4 |
| isotherme sur x : x₀ / x_sat | 0,43 / 0,05 | 0,52 / 0,26 | 0,54 / 0,26 |
| K_D vrai, n (isotherme sur x) | **30,3 mM, 1,46** | 15,1 mM, 1,18 | 19,5 mM, 1,23 |

Les incertitudes sont statistiques ; l'erreur systématique due à l'IRF est plus grande (τ_DA libre
varie de 0,08 à 0,35 ns d'un histogramme à l'autre). Tout est dans
`sorties/titration/` ; la figure `sorties/figures/titration_mCer3_Topaz.png` réunit
les déclins avec le modèle global, les résidus, x selon [Cl⁻], ton χ²ᵣ mono observé et prédit,
l'intensité et les phasors.

**Pourquoi mCer3_Topaz est la protéine par défaut.** C'est la seule des trois que M2 lent décrit
proprement. C'est aussi celle dont parle le plan : donneur mCer3, Φ_D = 0,87, et ton χ²ᵣ mono
moyen de 7,5. Le K_D de 14,4 mM et le n de 1,21 cités par le plan sont en revanche ceux de SCLM. Le
donneur de SCLM est vraisemblablement Cerulean, dont le déclin est multi-exponentiel par lui-même,
ce qui expliquerait ses χ² plus élevés et son intensité qui dérive.

**Les diagnostics du plan, sur mCer3_Topaz :**

- *Profil de χ²ᵣ mono.* Le plan attendait un profil en cloche, pire au mélange 50/50. Il supposait
  x₀ ≈ 1 ; or x₀ = 0,43, donc le mélange le plus prononcé est à 0 mM, et M2 lent prédit un profil
  **monotone**. C'est ce que tu observes, et quantitativement : prédit 21,4, 13,3, 14,9, 9,8, 7,7,
  2,8, 2,0, 1,3, 1,35 ; observé 18,7, 12,6, 11,9, 9,0, 6,4, 3,1, 2,1, 2,1, 1,37 (0 mM → F, moyennes
  des réplicats).
- *Validation croisée par les photons.* M2 prédit que l'amplitude du donneur à t = 0 ne dépend pas
  de [Cl⁻] : seule la répartition entre τ_DA et τ_D0 change. Ajustée librement, elle reste constante
  à 5 % près, sans tendance. Les acquisitions étant à temps fixé (10 s, arrêt sur temps), les
  comptages nets doivent alors monter comme (1 − E·x)/(1 − E·x₀) : ×1,619 prédit entre 0 mM et
  « F », ×1,609 mesuré. Le modèle explique donc à lui seul la hausse des comptages.
- *Phasors.* Les 27 points restent à 0,1 % (médiane), et au plus 0,32 %, de la longueur du segment
  τ_D0–τ_DA. Ils ne le parcourent toutefois que sur 4,5 %. Une vraie dérive de τ_D0 de 0,3 ns, comme
  celle des ajustements à τ libres (3,53 ns à 0 mM, 3,83 ns à 146 mM), écarterait les points de
  1,4 à 1,6 %, soit 5 à 15 fois plus. Cette dérive est donc plus probablement un artefact de
  corrélation entre paramètres libres (pic court sous l'IRF) qu'un vrai glissement.
- *τ_D0 contre mCer3 nu.* Pas de mCer3 seul dans ce jeu. FPbase donne 4,1 ns pour mCerulean3, contre
  3,69 ns ici. L'écart de 0,4 ns peut venir d'un transfert résiduel vers la forme protonée, ou de
  conditions de mesure différentes. Dans M2, il est absorbé par τ_D0 ; il n'est pas mis dans x_sat.
- *« F » compté comme C = ∞*, comme dans tes scripts. Sans ce point, x_sat n'est plus déterminé
  (l'isotherme donne x_sat < 0) : si « F » est du fluorure et non une saturation en Cl⁻, il faudra
  une titration qui monte plus haut.

## 2. Le modèle et ses paramètres

Sous M2 lent, chaque canal est une somme de trois exponentielles constantes (τ_DA, τ_D0, τ_A) dont
les amplitudes sont affines en x ; le Cl⁻ n'agit que sur x. Les formules du plan sont reprises telles
quelles (`src/Photophysics.jl`), plus le terme de risque ε_dark·Φ_dark·(1 − x), émis comme
l'accepteur. Le commutateur `mecanisme = "M2_rapide"` (un seul τ qui glisse) est en place : il est
testé et il fait sortir les phasors de la droite (test 8).

| Paramètre | Valeur | Provenance |
| --- | --- | --- |
| τ_D0, τ_DA | 3,689 ns ; 0,1515 ns | ta titration, ajustement global (§ 1) |
| x₀, x_sat, K_D, n | 0,429 ; 0,049 ; 30,3 mM ; 1,46 | isotherme à deux plateaux sur les 27 x |
| τ_A | 3,0 ns | inféré des YFP T203Y, comme le plan ; à mesurer sur mTopaz seul |
| Φ_D, Φ_A | 0,87 ; 0,60 | FPbase (mCerulean3, Topaz) |
| η_BD, η_BA, η_YD, η_YA | 0,336 ; 0,0019 ; 0,164 ; 0,590 | spectres FPbase × transmissions FPbase de tes filtres, détection plate |
| ε | 0,03 | milieu de la fourchette du plan (0,02–0,05) |
| ε_dark·Φ_dark | 0 | forme protonée non émissive ; balayée jusqu'à 0,05 |
| IRF, fenêtre | `irf_20260710_16x_sclm.sdt` ; canaux 14–242 | ton IRF, médiane soustraite ; ton gating |
| période, canaux | 12,5 ns ; 256 canaux de 48,8 ps | laser à 80 MHz |

**Les filtres ne sont pas ceux que supposait le plan** (480/40 et 540/40). Le bleu est un
**Chroma D470/40m** (« 473-em » dans tes noms de fichiers) et le jaune un **Chroma ET530/30m**,
d'après `tests/plot.jl` et `tests/plot_GECO_2_channels.jl` (canal 2 = 530/30). C'est une
inférence, à confirmer. Avec les rectangles idéaux du plan, η_BD vaudrait 0,529 au lieu de 0,336.

**Le contrôle « fuite ≈ 24 % du canal jaune ».** La fuite du donneur représente 36 % du jaune à
0 mM, 50 % à 20 mM et 89 % à saturation. L'ordre de grandeur concorde à bas Cl⁻, mais je ne sais
pas dans quelles conditions tu avais obtenu 24 %.

## 3. Les estimateurs (`src/Estimators.jl`)

Interface commune `observable(estimateur, h_bleu, h_jaune, contexte) → (observable, statut)`. Chaque
estimateur a sa propre courbe de calibration sans bruit, et [Cl⁻] s'obtient en l'inversant.

| # | Estimateur | Canaux | Paramètres libres | Observable |
| --- | --- | --- | --- | --- |
| 1 | Ratio d'intensités | B + J | 0 | N_jaune / N_bleu |
| 2 | Mono-exponentiel MLE, **ton code** | B | τ (+ décalage d'IRF et fond, comme dans ton code) | τ apparent |
| 2F | Ton code, décalage et fond figés | B | τ | τ apparent |
| 2H | Ton code + Hill, **ton pipeline** | B | idem 2 | τ, inversé par une Hill ajustée sur ta grille de titration |
| 3 | Bi-exponentiel MLE, τ libres | B | τ₁, τ₂, A₁, A₂ | fraction d'amplitude courte |
| 4 | Bi-exponentiel MLE, τ figés | B | A₁, A₂ | fraction d'amplitude courte |
| 5a | Conjoint deux canaux, τ figés | B + J | x, échelle B, échelle J | x |
| 5b | 5a avec des τ figés faux | B + J | idem | x (balayages des τ) |
| 5c | Conjoint, τ_A libre | B + J | x, τ_A, échelles B et J | x |
| 5d | Conjoint, **une seule échelle** | B + J | x, échelle | x (utilise aussi le rapport J/B) |
| 6 | Phasor | B | 0 | position sur l'axe de la titration |

J'ai ajouté **5d, 2F et 2H** au plan, pour les raisons suivantes :

- **5d** sépare ce que l'ajustement conjoint gagne grâce à la forme (5a : échelles libres, donc
  insensible au rapport J/B) de ce qu'il gagne grâce au rapport d'intensités. Sans 5d, on ne
  verrait pas que la forme n'apporte presque rien ici.
- **2H** est ton pipeline complet : la calibration y est ta Hill, pas la courbe exacte.
- **2F** mesure ce que te coûtent tes deux paramètres de nuisance. Il suffit de passer
  `fixed_parameters` à ton propre appel.

Toutes les vraisemblances sont de Poisson. Les modèles linéaires en x (4, 5a, 5b, 5d) profilent
les échelles analytiquement : il ne reste qu'une recherche 1-D, rapide et sans minimum parasite.
Les modèles non linéaires (3, 5c) passent par un Levenberg-Marquardt avec bornes (`src/Ajustement.jl`).

## 4. Choix de modélisation, et pourquoi

- **Le noyau instrumental est exact** (`src/Instrument.jl`). Il combine trois éléments :
  - le déclin périodisé par la série géométrique des impulsions précédentes, sans troncature ;
  - l'IRF supposée uniforme à l'intérieur de chacun de ses canaux (c'est un histogramme) ;
  - l'intégration sur le canal de détection.

  Le résultat est le déclin lissé par un triangle de demi-largeur Δ : il est dérivable partout, y
  compris à l'impulsion, et conserve exactement les photons. La convolution est circulaire, par FFT
  planifiée.
- **La fenêtre TAC** (canaux 14–242) est ajoutée au plan : tes histogrammes réels sont vides hors
  de cette fenêtre, et les photons qui y tomberaient sont perdus.
- **Trois normalisations du budget, comme le plan.**
  - A : temps d'acquisition fixé, avec N = photons bleus attendus à 20 mM. Les estimateurs à un
    canal y voient tes photons bleus actuels, et le canal jaune vient en plus.
  - B : N_bleu + N_jaune = N.
  - C : N_bleu = N.
- **Un échec** est un ajustement non convergé, ou convergé contre une borne d'un paramètre non
  linéaire (τ à sa borne, composantes confondues). Les échecs sont exclus des statistiques et
  comptés dans `taux_echec`. Une estimation hors de la plage de calibration est ramenée à 0 ou à
  150 mM et comptée dans `taux_hors_plage`.
- **Plafond à 150 mM** (≈ ta concentration la plus haute). Sans plafond, une estimation au plateau
  de l'isotherme part à l'infini. Ce plafond domine les statistiques des estimateurs qui saturent :
  pour eux, `biais_median_mM` et `sigma_robuste_mM` (1,4826 × MAD) sont plus parlants.
- **Bornes de Cramér-Rao en mM**, avec dérivées analytiques (le test 7 les compare à l'AD).
  `CRLB_mM` est la borne de **toute l'information** : deux canaux, τ connus, rapport J/B connu.
  `CRLB_propre_mM` est la borne des données que consomme chaque estimateur (totaux seuls pour le
  ratio, bleu seul, deux canaux à échelles séparées…). L'efficacité vaut σ²_borne/σ²_MC, comme
  dans le texte du plan (1 = au bout de ce qui est possible).
- **Ton code dans des modules isolés.** `lifetime_analysis.jl` garde son état dans des globales,
  ce qui le rend non sûr en parallèle. Il est donc chargé une fois par fil. Seuls deux noms
  de FLIMApp sont fournis (`FLIMCore.Reglages`, qui n'apparaît que dans des signatures de fonctions
  d'import jamais appelées, et `DEFAULT_HISTOGRAM_RESOLUTION`). L'appel reprend `fit_decay`
  d'`analyse_SCLM_2.jl` (`first_fit = true`, jusqu'à 5 essais), avec `hist_size_threshold = 0` :
  sans cela, ton code refuse tout histogramme de moins de 500 photons, ce qui donnerait 100 %
  d'échecs à 10² et 3·10².
- **Germes.** Chaque réalisation a son germe (`StableRNG`, stable d'une version de Julia à
  l'autre), le même pour ses deux canaux. Les résultats ne dépendent ni du nombre de fils ni de
  l'ordre d'exécution. Dans les balayages, les mêmes germes sont réutilisés d'une valeur à l'autre
  (nombres aléatoires communs) pour lisser les courbes.

## 5. Validation interne (`step2_validation.jl`)

**28 tests sur 28 réussis** (`sorties/validation/rapport_validation.md`) : les huit du plan, plus le
noyau contre l'intégrale double exacte, la FFT contre la convolution directe, l'adaptateur contre
tes résultats, et huit tests des modules de la révision (tableau du bas).

| Test du plan | Résultat |
| --- | --- |
| 1. Aller-retour sans bruit (7 paramètres libres) | écart relatif 9·10⁻¹⁶ ; 4, 5a, 5d rendent x à 3·10⁻¹⁵ |
| 2. Conservation des photons (donneur, émission sensibilisée) | 2·10⁻¹⁶ |
| 3. Amplitudes = fractions (bleu idéal) | exact |
| 4. Limite saturante : bleu mono-exponentiel, jaune = (η_YD/η_BD)·bleu | exact ; ton code y rend τ_D0 à 0,06 % près |
| 5. Limite dégénérée : jaune = −A à τ_DA, +A à τ_A | rapport −1,0000000000 |
| 6. Efficacité à 10⁶ photons (5a, 5d, 4, 3 ; 20 000 réalisations) | 1,005 ; 0,999 ; 0,996 ; 0,995 (IC 95 % ± 0,02) |
| 7. Dérivées analytiques contre ForwardDiff | 2·10⁻¹⁵ |
| 8. Linéarité phasor sous M2 lent | 7·10⁻¹⁴ (M2 rapide : 0,14) |
| Adaptateur : ton τ mono de SCLM 0 mM c1 | 2,0707331135068343, identique à `raw_data.csv` |

| Test de la révision | Résultat |
| --- | --- |
| 12. Modèle à N états du donneur, réduit à deux états, contre M2 lent | 5·10⁻¹⁶ |
| 13. IRF dilatée : identité à s = 1, aire conservée, largeur croissante | 4·10⁻¹⁶ ; 2·10⁻¹⁶ ; 2,40 < 2,63 < 2,90 canaux |
| 14. Retard Δt du jaune : forme décalée, ∂/∂t₀ contre différences finies | exact ; 4·10⁻⁶ |
| 15. Blanchiment : β = 0 ≡ M2, β → 0 continu, populations conservées | 3·10⁻¹⁶ ; 4·10⁻¹⁰ ; 1·10⁻¹⁶ |
| 16. Blanchiment du donneur : survie moyenne (1 − e^{−κ})/κ | 1·10⁻¹⁶ |
| 17. Ingestion d'une titration synthétique sans bruit (t₀ commun ou par histogramme) | τ à 2·10⁻⁹, x à 5·10⁻¹⁰, t₀ à 4·10⁻¹⁰ ns |
| 18. Ajustement conjoint bleu + jaune d'un banc simulé sans bruit | 8·10⁻⁸ |
| 19. Dérive de gain : écart-type de log g après τ_g | 0,049 contre 0,050 |

Le 0,06 % du test 4 n'est pas une erreur du simulateur. Ton `vec_to_lifetime` place les canaux au pas
du TAC (0,048861 ns), alors que la période et l'IRF sont au pas de 12,5/256 (0,048828 ns), d'où
τ × 1,0007. C'est sans conséquence ici.

## 6. Ce que disent les résultats

Figures dans `sorties/figures/`, chiffres dans `sorties/resultats/resume.md`, fichier plat du plan
dans `sorties/resultats/resultats_principaux.csv` (une ligne par cellule et par estimateur) et
`nuisances.csv`.

### Variance — question 1 (fig. 1, 2, 6)

σ[Cl⁻] à 20 mM, normalisation A, en mM :

| Estimateur | 10² | 10³ | 10⁴ | 10⁵ |
| --- | ---: | ---: | ---: | ---: |
| 1 ratio J/B | 11,2 | 3,1 | 0,98 | 0,30 |
| 5d conjoint, échelle commune | 11,3 | 3,0 | 0,97 | 0,30 |
| 5a conjoint, échelles séparées | 63 | 31 | 5,6 | 1,8 |
| 5c conjoint, τ_A libre | 64 | 38 | 6,8 | 2,1 |
| 4 biexp τ figés | 69 | 47 | 8,0 | 2,3 |
| 3 biexp τ libres | 19 * | 29 * | 12,5 | 3,6 |
| 6 phasor | 72 | 63 | 33 | 5,2 |
| 2F ton code, décalage et fond figés | 73 | 66 | 39 | 6,5 |
| 2 ton code | 71 | 56 | 59 | 9,8 |
| *borne de Cramér-Rao (toute l'information)* | 9,3 | 2,9 | 0,93 | 0,29 |

\* 49 % et 20 % d'ajustements en échec : σ ne porte que sur les réussis, d'où une valeur
trompeusement basse.

- **Le ratio et 5d atteignent la borne** dès 300 photons, avec une efficacité de 0,86 à 0,99.
  Elle tombe à 0,7 à 100 photons.
- **La forme des déclins ne porte qu'environ 3 % de l'information de Fisher sur x** : la borne à
  échelles séparées vaut 5,55 mM à 10⁴ photons, contre 0,93 mM pour la borne commune, et
  (0,93 / 5,55)² ≈ 0,03. Dès 10⁴ photons, 5a exploite cette information à 99 % ; il lui faut
  pourtant environ 36 fois plus de photons que le ratio pour la même précision.
- **Les trois normalisations du plan**, à 10⁴ photons et 20 mM :
  - **A, gain pratique** : ratio 0,98 mM contre 59 mM pour ton pipeline (2H), un facteur 60 ;
  - **B, gain méthodologique** (photons égalisés) : ratio 1,34 mM contre 15,3 mM pour le meilleur
    estimateur bleu seul (4), un facteur 11 ;
  - **C, gain pur du deuxième canal** (canal bleu identique) : 5d 0,92 mM contre 8,7 mM pour 4,
    un facteur 9.
- **Le long de l'isotherme** (fig. 1, 10⁴ photons), σ du ratio est minimal entre 8 et 16 mM
  (≈ 0,75 mM), puis passe à 1,2 mM à 30 mM et 9 mM à 100 mM. À 0 mM, la borne est infinie :
  l'isotherme y a une pente nulle, puisque n > 1.
- **Efficacité propre** (fig. 6), c'est-à-dire la part de l'information de leurs propres données
  qu'ils extraient :
  - dès 10⁴ photons, 5a, 5c, 4 et 3 en tirent 80 à 100 % ;
  - sous 3·10³ photons, ils tombent à 0,2–0,6, car leurs estimations s'empilent contre les
    plateaux de la calibration ;
  - à 10⁴ photons, ton code (2) n'en tire que 1,5 %, et le phasor 5 %.

### Stabilité à bas photons — question 3 (fig. 2, 3)

- **Ratio et 5d** : biais de +1,5 mM à 100 photons, de moins de 0,3 mM dès 300 ; aucun échec.
- **3 (bi-exp, τ libres)** est le seul qui cesse vraiment de converger : 49 % d'échecs à 100
  photons, 36 % à 300, 20 % à 10³, 11 % à 3·10³ et 2 % à 10⁴. 5c échoue dans 5 % des cas à
  100 photons.
- **Les estimateurs à forme seule** ne divergent pas, ils saturent. Sous 3·10³ photons, leurs
  estimations s'empilent contre 0 ou 150 mM, d'où un biais de +3 à +45 mM à 20 mM. C'est le
  « biais qui disparaît à 10⁴ » du plan : celui de 5a passe de +7 mM à 10³ photons à −0,1 mM
  à 10⁴.
- **Ton code** surestime τ de 0,2 à 0,3 ns sous 3·10³ photons (fond et décalage libres), d'où la
  saturation décrite plus haut. Son biais vaut +37 mM à 10⁴, +7 mM à 3·10⁴ et +1 mM à 10⁵.
  **Figer le décalage et le fond (2F)** divise σ_τ par 3 (0,053 contre 0,157 ns à 10⁴ photons) et
  supprime le biais de τ à bas photons : 3,529 ns à 10³ photons pour 3,526 ns vrai, contre 3,769 ns
  avec ton ajustement libre. En mM, le σ robuste passe de 31 à 17 mM à 10⁴ photons. C'est mieux,
  mais encore 17 fois le ratio, puisque τ mono ne bouge que de 0,3 ns sur toute la plage. À photons
  infinis, ta Hill ne coûte que ±0,9 mM entre 0 et 60 mM (+2,9 mM à 100 mM) : l'écart entre
  le modèle mono et M2 n'est pas le problème, c'est le bruit.

### Biais sous nuisance — question 2 (fig. 4, 4b)

Il s'agit du biais **ajouté** à 20 mM et 10⁴ photons quand la tranche diffère de la calibration :
estimateurs et calibrations construits sur les valeurs nominales, données tirées du monde perturbé.
Il est mesuré par rapport au monde nominal tiré avec les mêmes germes. En mM, valeurs basse / haute :

| Nuisance | 1 ratio | 5d | 5a (5b pour les τ) | 5c | 4 | 3 | 2 ton code |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| η_YD × 0,75 / × 1,25 | +10 / −7 | +10 / −7 | −3,8 / +4,6 | −3,0 / +3,9 | 0 | 0 | 0 |
| gain du jaune × 0,95 / × 1,05 | +3,6 / −3,1 | +3,5 / −3,0 | +0,2 / 0,0 | +0,2 / 0,0 | 0 | 0 | 0 |
| gain du jaune × 0,90 / × 1,10 | +7,9 / −5,8 | +7,7 / −5,7 | +0,3 / −0,1 | +0,2 / −0,1 | 0 | 0 | 0 |
| ε = 0 / 0,10 (nominal 0,03) | +1,0 / −2,2 | +1,0 / −2,2 | +0,2 / −0,4 | −0,3 / +0,8 | 0 | 0 | ≈ 0 |
| ε_dark·Φ_dark = 0,01 / 0,05 | −1,3 / −5,9 | −1,3 / −5,8 | −0,2 / −1,2 | +0,4 / +2,0 | 0 | 0 | ≈ 0 |
| τ_D0 × 0,95 / × 1,05 | +0,2 / −0,1 | −0,1 / +0,1 | −9 / +13 | −7 / +8,5 | −9 / +15 | +0,6 / +0,9 | −43 / +60 |
| τ_D0 × 0,90 / × 1,10 | +0,4 / −0,3 | −0,3 / +0,3 | −17 / +33 | −14 / +18 | −17 / +44 | +1,0 / +1,6 | −58 / +86 |
| τ_DA × 0,90 / × 1,10 | −0,2 / +0,2 | −0,2 / +0,2 | +1,5 / −1,1 | +2,8 / −2,2 | +2,5 / −2,0 | +0,5 / −0,1 | +5 / −7 |
| τ_A × 0,90 / × 1,10 | −0,1 / +0,1 | −0,4 / +0,4 | −8,7 / +17 | +0,1 / +0,1 | 0 | 0 | ≈ 0 |
| τ_A × 0,85 / × 1,15 | −0,2 / +0,2 | −0,6 / +0,6 | −12 / +35 | 0 / +0,1 | 0 | 0 | ≈ 0 |

« 0 » : l'estimateur n'utilise pas ce paramètre. « ≈ 0 » : dans le bruit Monte-Carlo de ton code
(σ/√1000 ≈ 2 mM).

- **La fuite spectrale (fig. 4, la figure centrale du plan)** confirme l'hypothèse du plan, mais
  de peu. Une erreur de 25 % sur η_YD biaise le ratio 2 à 3 fois plus que 5a. Si η_YD est bien
  étalonné, la fuite fait passer σ du ratio de 0,56 mM (sans fuite) à 0,95 mM, et améliore 5a
  (8,9 → 5,7 mM), comme le plan l'avait prévu. Le ratio reste toutefois 6 fois plus précis :
  avec ±25 % d'erreur sur η_YD, ratio et 5a ont une RMSE comparable (7 à 10 mM contre 7 mM).
  Avec η_YD étalonné à quelques % près, le ratio l'emporte nettement.
- **Une dérive de gain est le seul cas où 5a l'emporte**, et c'est l'argument central du plan.
  - Une erreur de 5 % sur le gain relatif du jaune biaise le ratio et 5d de ±3–4 mM (±6–8 mM à
    10 %, +20 / −11 mM à 20 %), et laisse 5a intact (< 0,3 mM).
  - En RMSE moyenne sur la fenêtre de clamping, à 10⁴ photons, 5a reste entre 5,6 et 6,1 mM. Il
    ne passe devant le ratio qu'au-delà de ±8 % d'erreur de gain environ (ratio : 3,3–3,7 mM à
    ±5 %, 6,0–7,8 mM à ±10 %).
  - À 10³ photons, 5a ne passe jamais devant le ratio jusqu'à ±20 % (26–31 mM contre 10–21 mM).
  - Le photoblanchiment de l'accepteur n'est pas une simple perte de gain : voir les réserves.
- **Les durées de vie sont le talon d'Achille des lectures en forme.** Il suffit d'un τ_D0 décalé
  de 5 % entre la solution de calibration et la tranche pour biaiser 5b, 5c et 4 de 7 à 15 mM.
  C'est un ordre de grandeur réaliste : l'indice de réfraction du cytoplasme (≈ 1,37 contre 1,33)
  raccourcit les durées de vie d'environ 5 %. Le ratio et 5d n'en bougent pas (< 0,5 mM). Une
  erreur de 10 % sur τ_A, la valeur supposée, coûte −9 à +17 mM à 5a ; 5c, qui ajuste τ_A,
  n'en bouge pas. 3, qui ajuste ses τ, est insensible à toutes ces erreurs, au prix d'un
  σ de 12,5 mM.
- **ε et ε_dark·Φ_dark sont des nuisances mineures** (≤ 2 mM sur la plage du plan), sauf
  ε_dark·Φ_dark sur le ratio : −1,3 mM par centième. Le risque à 800 nm que signalait le plan
  pèse donc surtout sur les lectures en intensité.
- **Le cas dégénéré du plan** (η_YD = 0 ; ε laissé à 0,03) : σ du ratio = 0,56 mM et σ de
  5a = 8,9 mM, contre 0,95 et 5,7 mM avec la fuite réelle.

### K apparents — le livrable secondaire (fig. 5)

| Lecture | K apparent (mM) | K app. / K vrai | Mesuré par toi |
| --- | ---: | ---: | ---: |
| x directement (3, 4, 5a, 5d) | 30,3 (n = 1,46) | 1,00 | — |
| ratio J/B, phasor, ⟨τ⟩ pondéré en intensité | 21,8 (n = 1,46) | 0,72 | — |
| ton τ mono (2) | 18,7 (n = 1,34) | 0,62 | **18,5 (n = 1,30)**, Hill sur τ mono |
| ta Hill sur ta grille de titration (2H) | 18,5 (n = 1,33) | 0,61 | |
| ratio inversé B/J | 40,6 (n = 1,46) | 1,34 | **44,2 (n = 1,06)**, ton C2/C1 TIFF |

- **Lectures pondérées en intensité.** K_app = K·(D(0)/D(∞))^{1/n}, où D est l'intensité du
  donneur à chaque plateau. C'est la formule du plan, généralisée à x₀ ≠ 1 et à n ≠ 1 ; elle
  donne 21,8 mM, exactement ce que trouvent le ratio et le phasor.
- **Ton τ mono** n'est pas exactement pondéré en intensité, d'où 18,7 mM et n = 1,34.
- **Le K de 14,4 mM cité par le plan** est celui de SCLM. Pour mCer3_Topaz, le K vrai vaut
  30 mM : il tombe dans la fourchette de 26 à 36 mM du plan, mais par un autre chemin. Le plan
  divisait par τ_DA/τ_D0 ≈ 0,4–0,55 ; ce qui compte ici est le rapport des intensités du donneur
  aux deux plateaux (0,62), puisque τ_DA/τ_D0 ne vaut que 0,04.
- **D'où vient n > 1.** Le n de 1,46 se trouve déjà sur x, avant toute lecture. La coopérativité
  apparente vient donc du capteur, ou du point « F », ou de la force ionique ; elle ne vient pas
  du biais de l'ajustement mono. Le plan supposait n = 1 sur x.

### La carte (carte_meilleur_estimateur.png)

- **Sur toute la grille** (13 concentrations × 7 budgets × 3 normalisations), le ratio et 5d sont
  ex æquo, à 5 % près, pour la plus petite RMSE.
- **Parmi les estimateurs de forme seule**, le meilleur est presque partout 5a, avec une RMSE 2 à
  12 fois plus grande.
- **Donc** : l'ajustement conjoint ne « gagne » ni sous un seuil de photons, ni au-delà d'une
  fuite. Il ne gagne que si l'échelle J/B devient incertaine (balayage du gain) ; c'est la
  réponse au « quand l'utiliser » du plan.

### Si le capteur transférait moins (fig. 8)

σ à 20 mM, 10⁴ photons, calibration exacte, en mM, pour des capteurs identiques sauf τ_DA :

| τ_DA (ns) | 0,15 (ta titration) | 0,5 | 1,0 | 1,5 | 2,1 (hypothèse du plan) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 ratio | 0,98 | 1,05 | 1,29 | 1,62 | 2,31 |
| 5d | 0,97 | 1,02 | 1,20 | 1,45 | 2,02 |
| 5a | 5,6 | 3,9 | 3,6 | 3,7 | 4,7 |
| 4 | 8,0 | 4,7 | 4,0 | 4,1 | 5,4 |
| 6 phasor | 33 | 6,7 | 4,6 | 4,7 | 6,1 |

- **Même avec le capteur du plan** (τ_DA = 2,1 ns), le ratio reste 2 fois plus précis que 5a.
- **La forme apporte le plus vers τ_DA ≈ 1 ns.** 5d n'y gagne que 13 % sur le ratio.
- **Un capteur à FRET moins efficace** rendrait la forme plus informative, mais il dégraderait
  davantage le ratio qu'il n'aiderait la forme.

### Cible mobile et boucle fermée (fig. 7)

- **Boucle ouverte**, à 10⁴ photons par échantillon :
  - RMSE sur la trajectoire de 0,97 mM pour le ratio, 0,96 pour 5d, 6,2 pour 5a, 10,3 pour 4 et
    71 pour ton pipeline ;
  - un Kalman dont la variance de mesure R est la borne évaluée à l'estimation courante
    (usage 2 du plan) ramène le ratio à 0,32 mM ;
  - à 10³ photons : 3,2, 3,2, 36, 52, 115 et 0,82 mM (Kalman).
- **Boucle fermée** (consigne 20 mM, intégrateur PI, relaxation de 20 s) :
  - erreur statique de +0,03 mM pour le ratio et +0,02 mM pour 5d, −0,75 mM pour 5a, −1,3 mM pour
    4 et −1,6 mM pour ton pipeline à 10⁴ photons ;
  - −7,5 mM pour ton pipeline à 10³ photons.

  L'intégrateur moyenne le bruit sur ~20 s : il reste le **biais** à la consigne, comme le
  plan le soupçonnait. Les gains du PI sont arbitraires ; seule l'erreur statique est
  interprétable.

### Ce que ça change pour toi

1. **Lis [Cl⁻] sur le rapport J/B des deux canaux FLIM**, ou par 5d. Les comptages TCSPC sont des
   photons, contrairement à tes TIFF Bliq. C'est la borne à 1–3 % près, et environ 60 fois mieux
   que ton τ mono à 10⁴ photons. Depuis la révision : de préférence en fusion avec la forme (gain
   dans l'état d'un filtre de Kalman, fig. 7b), et, en attendant le canal jaune, remplace au moins
   mono + Hill par le bi-exponentiel à τ figés (facteur 7 vérifié sur tes données).
2. **Ce que cela exige : étalonner l'échelle J/B, puis la surveiller.** En pratique, mesurer la
   fuite η_YD/η_BD sur un échantillon donneur seul (mCer3 seul), puis suivre le gain relatif des
   deux détecteurs. Seuil à retenir : tant que ce gain relatif est connu à mieux que ±8 % environ
   (à 10⁴ photons), le ratio reste meilleur. La forme (5a), insensible à l'échelle, sert de
   garde-fou lent plutôt que d'estimateur principal.
3. **Si tu gardes une lecture en durée de vie**, fige au moins le décalage d'IRF et le fond (2F).
   Sache qu'elle reste bien moins précise que le ratio, et qu'elle est très sensible à un τ
   décalé entre la solution et la tranche.
4. **Mesures qui manquent**, par ordre d'impact sur ces conclusions :
   - le gain relatif des deux détecteurs et sa stabilité ;
   - τ_A sur mTopaz seul ;
   - τ_D0 en tranche contre en solution ;
   - le décalage temporel entre détecteurs. Il n'affecte pas le ratio mais touche 5a et 5c,
     dont la montée du jaune tient sous l'IRF ;
   - la confirmation des filtres.

### Réserves

- **Photoblanchiment** (exclusion n° 1 du plan, simulé depuis la révision, fig. 14). Le résultat
  dépend du produit du blanchiment de l'accepteur. S'il n'absorbe plus, toutes les lectures sont
  biaisées, la forme plus que le ratio. S'il reste extincteur, seul le ratio l'est (1,8 mM par 10 %
  d'accepteurs blanchis). Les taux sont supposés indépendants de l'état (un donneur éteint passe
  moins de temps excité).
- **Décalage entre détecteurs** (exclusion n° 2, simulé depuis la révision, fig. 13). Il ne pèse
  pas sur le ratio, mais 5a et 5c prennent 2 à 3 mM par 10 ps.
- **Valeurs supposées.** τ_A = 3,0 ns et ε = 0,03 sont supposés, les filtres sont déduits de tes
  scripts, et Φ vient de FPbase. La sensibilité à chacune est dans le tableau des nuisances.
- **τ_DA = 0,15 ns** survit à la batterie de la révision : tes données rejettent τ_DA ≥ 0,3 ns,
  et une IRF fausse de 10 % ne fabrique pas un transfert quasi total. Son incertitude systématique
  (IRF plus étroite, gigue de t₀) va de 0,12 à 0,19 ns. La fig. 8 montre que le classement tient
  jusqu'à 2,1 ns.
- **M2 lent reste une description au premier ordre** (χ²ᵣ 1,96 à 10⁶ photons par histogramme).
  La révision attribue l'essentiel de l'excès à une gigue de t₀ entre échantillons (3,6 ps) et à un
  état FRET intermédiaire (H3a), qui rendraient les valeurs absolues du modèle incertaines de
  quelques mM.
  L'IRF, médiane soustraite comme dans ton code, garde un plancher de bruit écrêté (2,3 % de son
  aire hors du pic). Il est utilisé à l'identique pour générer et pour ajuster : c'est cohérent,
  pas tout à fait réaliste.
- **Protéine purifiée en solution.** En tranche, l'autofluorescence, le pH, la diffusion et
  l'indice de réfraction s'ajoutent, et le capteur des tranches est SCLM, que M2 décrit moins bien.
- **Une seule machine** (M1, 8 fils). Le `time_limit` de 60 ms de ton optimiseur n'a jamais
  produit de NaN ici ; sur une machine plus lente ou chargée, le taux d'échec de l'estimateur 2
  pourrait changer.

## 7. Écarts au plan

| Plan | Ici | Pourquoi |
| --- | --- | --- |
| Filtres 480/40 et 540/40, rectangles | D470/40m et ET530/30m, transmissions FPbase | tes vrais filtres (à confirmer) |
| Une fraction x par concentration | une fraction x par histogramme | les réplicats d'un même tube collent à ±0,005, mais d'un tube à l'autre l'écart est de ±0,02 ; l'isotherme en tient compte |
| 256 canaux tous utiles | fenêtre TAC 14–242 | tes histogrammes sont vides hors de cette fenêtre |
| 5 estimateurs (+ 4 et phasor) | + 5d, 2F, 2H | voir § 3 |
| K_D = 14,4 mM, n = 1 au départ | K vrai = 30,3 mM, n = 1,46, tirés de l'isotherme sur x | le plan demandait de ne pas injecter un K apparent ; 14,4 mM était celui de SCLM |
| — | balayages du gain du jaune et du τ_DA du capteur | le gain teste l'argument central (robustesse de la forme) ; τ_DA montre quand la forme deviendrait utile |
| Kalman (usage 2 de la borne) | démontré dans la cible mobile | R = borne de Cramér-Rao à l'estimation courante |

## 8. Sorties

| Fichier | Contenu |
| --- | --- |
| `parametres.toml` | tous les paramètres, leur statut (ajusté, supposé, mesuré) et leur provenance |
| `sorties/spectres/` | spectres et filtres FPbase (cache), `eta.toml` |
| `sorties/titration/` | ajustement global par protéine : paramètres (`parametres_*.toml`), un histogramme par ligne (`titration_*.csv`), courbes (`ajustement_global_*.csv`) |
| `sorties/validation/rapport_validation.md` | les 28 tests |
| `sorties/calibration/` | courbes de calibration, K apparents, identifiabilité, bornes (par C, N et normalisation) |
| `sorties/resultats/resultats_principaux.csv` | grille principale (fichier plat du plan) |
| `sorties/resultats/nuisances.csv` | balayages de nuisance |
| `sorties/resultats/resume.md` | tableaux chiffrés |
| `sorties/figures/` | figures 1 à 14 et 7b, carte, titrations, deux synthèses (`synthese_{FR,EN}{,_dark}.png`, `synthese_revision_*`, étape 8) |
| `sorties/batterie/` | batterie H1–H5 (`batterie.csv`), tes données réajustées sous IRF et t₀ perturbés (`reel_perturbations.csv`), profil de χ²ᵣ en τ_DA, t₀ par histogramme, chiffres clés (`batterie_resume.toml`) |
| `sorties/reel/` | le facteur 8 : résumé (`facteur8.csv`), calibrations, toutes les lectures (`amincis_tirages.csv`, ignoré par git) |
| `sorties/incertitude/` | mondes IRF/t₀ (`mondes_grille.csv`, `mondes_K_apparents.csv`), variance paramétrique, décomposition par Monte-Carlo, survie du classement (`classement.csv`) |
| `sorties/canal_jaune/` | inversions de η_YD, aller-retour de l'ingestion du jaune, prédictions de J/B, observabilité de la montée, exemple de bloc « mesuré » |
| `sorties/calibration/eta_YD_inverse.toml` | η_YD tiré du K apparent et de la plage dynamique de tes images TIFF |
| `sorties/resultats/fusion_*.csv` | fusion ratio + forme : résumé, trajectoire, réponse à un saut de gain |
| `src/` | le simulateur : `Types`, `Photophysics`, `Instrument`, `Spectra`, `Generate`, `Fisher`, `Ajustement`, `Utilisateur`, `Ingest`, `Parametres`, `Estimators`, `Calibrate`, `Experiments`, et pour la révision `Batterie`, `Reel`, `Incertitude`, `Jaune`, `Fusion` ; `Figures`, `Resume` |
