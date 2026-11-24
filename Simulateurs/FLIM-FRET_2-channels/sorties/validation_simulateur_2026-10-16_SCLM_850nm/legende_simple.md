**Figure A.1 — Validation du simulateur de déclins de fluorescence.**
Le simulateur génère les histogrammes de temps d'arrivée des photons (comptage de photons uniques
corrélé en temps, canaux de 64 ps, période laser de 12,5 ns) émis par le donneur du capteur de chlorure
SuperClomeleon (Cerulean, donneur ; Topaz, accepteur ; excitation à deux photons à 850 nm). Le déclin du
donneur est modélisé comme un mélange de deux populations : des donneurs en FRET (durée de vie 0,655 ns)
et des donneurs sans FRET (2,57 ns). On y ajoute la lumière diffusée et le bruit de fond du détecteur,
dans les proportions mesurées. Le tout est convolué par la réponse instrumentale mesurée, puis tiré
selon une loi de Poisson. Les durées de vie proviennent d'un ajustement global de la titration en
solution (8 concentrations de Cl⁻, de 0 à 131 mM). La fraction FRET est la part d'amplitude de la
population en FRET.

**(a) Retrouver une fraction FRET connue.** Pour chaque fraction injectée (0,05 à 0,95) et chaque
nombre de photons par déclin (1 000, 10 000 et 100 000), 500 déclins sont simulés. La fraction est
ensuite estimée par maximum de vraisemblance de Poisson, durées de vie, réponse instrumentale et
origine des temps étant connues. Les points donnent la moyenne, les barres ±1 écart-type ; ils sont
décalés horizontalement selon le nombre de photons, pour la lisibilité. Le tireté marque l'égalité
entre fraction retrouvée et injectée.
- Les estimations suivent l'égalité : l'écart moyen reste sous 0,002 à 10 000 photons et sous 0,0004
  à 100 000.
- La seule exception est à 1 000 photons pour la plus petite fraction (+0,032), où l'interdiction
  d'amplitudes négatives tronque la distribution.
- L'écart-type décroît comme 1/√N : 0,01–0,11 à 1 000 photons, 0,004–0,04 à 10 000, 0,001–0,014 à
  100 000.
- Ce panneau vérifie la cohérence interne : le même modèle génère et ajuste.

**(b) Comparer une simulation à une mesure réelle.** Les points noirs sont le déclin mesuré du donneur
à 20 mM de Cl⁻ (somme de 10 acquisitions de 1 s, 1 387 370 photons). La courbe orange est un déclin
simulé avec les mêmes durées de vie, amplitudes, lumière diffusée, fond et réponse instrumentale, et
exactement le même nombre de photons. L'origine des temps est placée à mi-hauteur du front montant ;
les points avant 0 viennent de la queue de l'impulsion précédente.
- **Le panneau inférieur** donne l'écart de chaque canal au modèle commun, en écarts-types du bruit de
  photons ((n − μ)/√μ). La bande grise (±2) contient 95 % des écarts dus au seul bruit.
- **La simulation est du bruit pur :** χ²ᵣ de 0,8 pour le tirage montré, 1,0 en moyenne sur 200
  tirages.
- **La mesure garde un écart systématique** (χ²ᵣ = 4,0), concentré près du front montant et en fin de
  fenêtre. Il trahit une description imparfaite de la réponse instrumentale.
- **Cet écart ne se voit qu'avec beaucoup de photons.** Ramenée par sous-échantillonnage binomial à
  1 000, 10 000 et 100 000 photons, la mesure donne χ²ᵣ = 0,98, 1,02 et 1,21, contre 0,99, 0,99 et 1,01
  pour la simulation. Jusqu'à ~10⁴ photons par déclin, le régime des simulations, mesure et simulation
  sont indiscernables ; au-delà de ~10⁵, le simulateur est légèrement optimiste pour les lectures
  fondées sur la forme du déclin.
