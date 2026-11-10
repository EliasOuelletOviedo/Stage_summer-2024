# Simulateur FLIM-FRET à deux canaux — plan de conception

Oct 15, 2026 · @Elias

## Ce que le simulateur doit trancher

Une seule sortie, en mM : l'erreur sur \[Cl⁻\] de cinq estimateurs, en fonction du budget de photons et de la concentration vraie. Tout le reste est un moyen d'y arriver.

Trois questions distinctes, qui n'ont pas la même réponse :

1. **Variance.** À budget de photons égal, quel estimateur a le plus petit σ\[Cl⁻\] ? Référence absolue : la borne de Cramér-Rao, qui dit combien d'information existe dans les deux histogrammes indépendamment de l'estimateur.
2. **Biais.** Sous fuite spectrale du donneur dans le canal jaune, excitation directe de l'accepteur, et surtout sous mauvaise spécification du modèle d'ajustement, lequel reste juste ? Un estimateur de faible variance mais biaisé de 8 mM est inutilisable pour un asservissement.
3. **Stabilité à bas photons.** À 300 ou 1 000 photons, lesquels convergent encore, et lesquels développent un biais qui disparaît à 10⁴ ? C'est le régime réel de tes tornades round-robin, et c'est là que le nombre de paramètres libres décide de tout.

**Plage d'intérêt.** La titration simulée couvre 0 à 100 mM, avec une densité de points accrue entre 5 et 35 mM : c'est ta fenêtre de clamping, et c'est aussi la zone où la sensibilité d\[signal\]/d\[Cl⁻\] est la plus faible. Un estimateur qui gagne à 2 mM et perd à 25 mM ne t'aide pas.

**Livrable secondaire, non demandé mais gratuit.** La cartographie K\_vrai → K\_apparent pour chacun des cinq estimateurs. Sous le modèle retenu ci-dessous, chaque estimateur déforme la constante de dissociation d'un facteur différent et connu. Si tes titrations expérimentales donnent déjà des K\_D apparents différents entre ratio et durée de vie, c'est une validation expérimentale directe du modèle.

## Le modèle physique : M2 lent

**Le chlorure n'agit sur aucune durée de vie. Il agit sur une fraction de population.** Les trois durées de vie du système sont des constantes ; toute la dépendance en \[Cl⁻\] passe par une amplitude. C'est la décision structurante du simulateur et elle découle de trois faits.

### Pourquoi ce mécanisme et pas un autre

**L'absorption de l'accepteur change.** R₀ dépend de l'intégrale de recouvrement J(λ) = ∫F\_D(λ)·ε\_A(λ)·λ⁴dλ, donc du coefficient d'extinction de l'accepteur, pas de son rendement quantique. Si le chlorure ne faisait que quencher l'émission de mTopaz, R₀ et E seraient inchangés et τ du donneur ne bougerait pas. Il bouge. Donc c'est l'absorption qui change : la protonation du chromophore T203Y déplace son absorption de 514 vers \~400 nm, le recouvrement avec l'émission de mCer3 s'effondre, et l'accepteur cesse d'être un accepteur.

**L'échange est gelé à l'échelle de la nanoseconde.** Le clignotement par protonation des YFP relaxe sur des dizaines de microsecondes. Devant un déclin de 3 ns, c'est un mélange statique, pas une moyenne temporelle.

**Tes données le confirment.** Déclin court et visiblement non mono-exponentiel à bas Cl⁻, mono-exponentiel propre à haut Cl⁻. Un échange rapide prédirait du mono-exponentiel partout, avec un τ qui glisse. Un mélange statique prédit exactement ce que tu vois.

### Les trois durées de vie et le seul paramètre libre

| Symbole | Signification | Dépend de \[Cl⁻\] ? |
| --- | --- | --- |
| τ\_D0 | mCer3 sans accepteur actif | non |
| τ\_DA | mCer3 avec accepteur actif, = τ\_D0(1−E) | non |
| τ\_A | mTopaz déprotoné | non |
| x | fraction de capteurs à accepteur actif | **oui — c'est le signal** |

Avec la fraction déprotonée donnée par une isotherme de Langmuir à deux plateaux :

```latex
x(C) = x_{\mathrm{sat}} + (x_0 - x_{\mathrm{sat}})\,\frac{K_D^{\,n}}{K_D^{\,n} + C^{\,n}}
```

Les deux plateaux x₀ et x\_sat absorbent à la fois la maturation incomplète de l'accepteur, la saturation incomplète et le pH de repos. Ces trois effets sont **mathématiquement indistinguables** ici : ils n'entrent que par le produit qui fixe les plateaux. Inutile donc de porter un paramètre de maturation séparé — il ne serait pas identifiable.

### Émission du donneur

Par excitation du donneur, en photons par ns :

```latex
I_D^{\mathrm{em}}(t) = \frac{\Phi_D}{\tau_{D0}}\Big[\,x\,e^{-t/\tau_{DA}} + (1-x)\,e^{-t/\tau_{D0}}\Big]
```

Les amplitudes à t = 0 valent exactement x et (1−x). **L'analyse pondérée par les amplitudes lit donc la fraction directement**, sans passer par une durée de vie moyenne. C'est le fondement de l'avantage attendu.

### Émission de l'accepteur

Deux termes. L'émission sensibilisée, issue de la fraction active, avec k\_t = 1/τ\_DA − 1/τ\_D0 :

```latex
I_A^{\mathrm{sens}}(t) = \frac{\Phi_A}{\tau_A}\cdot\frac{x\,k_t}{k_A - k_{DA}}\Big(e^{-t/\tau_{DA}} - e^{-t/\tau_A}\Big)
```

dont l'intégrale vaut Φ\_A·x·E, comme attendu. Puisque τ\_DA < τ\_A, le terme en τ\_DA porte l'amplitude négative : c'est la montée, et sa constante de temps est τ\_DA — constante en \[Cl⁻\]. Et l'excitation directe, qui ne touche que les accepteurs déprotonés, donc proportionnelle à x elle aussi :

```latex
I_A^{\mathrm{dir}}(t) = \varepsilon\,x\,\frac{\Phi_A}{\tau_A}\,e^{-t/\tau_A}
```

### Ce que voient les deux canaux

Chaque canal est une combinaison linéaire des mêmes trois exponentielles, avec des poids fixés par la matrice de recouvrement spectral η\[canal, émetteur\] :

```latex
I_{\mathrm{bleu}}(t) = \eta_{BD}\,I_D^{\mathrm{em}} + \eta_{BA}\big(I_A^{\mathrm{sens}} + I_A^{\mathrm{dir}}\big)
```

```latex
I_{\mathrm{jaune}}(t) = \eta_{YD}\,I_D^{\mathrm{em}} + \eta_{YA}\big(I_A^{\mathrm{sens}} + I_A^{\mathrm{dir}}\big)
```

**Conséquence centrale, et c'est l'argument du projet.** Les deux histogrammes sont entièrement décrits par trois durées de vie constantes, deux facteurs d'échelle d'intensité, et **un seul paramètre physique libre, x**. Le problème d'estimation à durées de vie figées a donc un paramètre d'intérêt unique, contre deux ou trois pour un ajustement bi-exponentiel classique. À 1 000 photons, cette différence de dimensionnalité n'est pas marginale.

**Et la fuite spectrale change de camp.** Le terme η\_YD·I\_D porte la même fraction x que le reste, avec la même forme temporelle. Pour un ratio d'intensités, c'est une catastrophe : il ajoute au canal jaune une copie atténuée du canal bleu, qui tire le rapport vers 1 et comprime ta plage dynamique. Pour un ajustement conjoint, c'est une troisième contrainte sur le même paramètre. Le simulateur doit montrer ces deux effets opposés sur la même figure.

## Paramètres et provenance

Chaque valeur porte son origine. Les valeurs inférées sont des points de départ à remplacer dès que tu auras la mesure correspondante — le code doit rendre ce remplacement trivial.

| Paramètre | Valeur de départ | Provenance |
| --- | --- | --- |
| τ\_D0 | à extraire | **de ta titration** : ton déclin mono-exponentiel à Cl⁻ saturant *est* τ\_D0 |
| τ\_DA | à extraire | **de ta titration** : composante courte, liée globalement sur toute la série |
| τ\_A (mTopaz) | 3,0 ns | inféré des YFP T203Y ; à mesurer sur mTopaz seul |
| x₀, x\_sat | à extraire | **de ta titration** : amplitudes aux deux extrêmes |
| K\_D, n | 14,4 mM, 1,0 | ton ajustement de Hill, mais voir l'avertissement ci-dessous |
| η\_BD, η\_YD | calculé | spectre d'émission mCer3 (FPbase) intégré sur 460–500 et 520–560 nm |
| η\_BA, η\_YA | calculé | spectre d'émission Topaz (FPbase) sur les mêmes bandes |
| ε (exc. directe) | \~0,02–0,05 | spectre 2P mCer3 mesuré (Leben 2022) ÷ spectre 2P inféré de Citrine (Drobizhev 2011), à 800 nm |
| Φ\_D, Φ\_A | 0,87 / 0,60 | FPbase (mCerulean3 / Topaz) |

### Deux avertissements sur les valeurs de départ

**N'injecte pas ton K\_D de 14,4 mM comme vérité.** Sous M2, ⟨τ⟩ est une fonction homographique de x, et x est une fonction homographique de \[Cl⁻\] ; la composée est homographique, donc la courbe τ(\[Cl⁻\]) est **rigoureusement hyperbolique, n = 1 exactement**, quel que soit le mélange. Le biais du fit mono-exponentiel ne peut pas fabriquer ton n = 1,21 — il fabrique autre chose :

```latex
K_D^{\mathrm{apparent}} = K_D^{\mathrm{vrai}}\cdot\frac{\tau_{DA}}{\tau_{D0}}
```

pour une lecture pondérée par l'intensité, et aucun décalage pour une lecture pondérée par les amplitudes ; ton MLE mono tombe entre les deux. Avec un rapport τ\_DA/τ\_D0 de 0,4 à 0,55, **ton K\_D vrai serait plutôt de 26 à 36 mM**. Le simulateur prend K\_vrai en entrée et produit les K\_apparents ; il ne consomme pas un K\_apparent comme s'il était vrai.

Le n = 1,21 reste à expliquer séparément. Candidats, par ordre de probabilité : plateau de saturation non atteint dans ton ajustement, second site de liaison, effet de force ionique sur tes solutions de calibration.

**ε mérite un paramètre de risque.** À 800 nm en 2P, tu excites à l'équivalent 2P de 400 nm — précisément où absorbe le chromophore **protoné** de mTopaz. La forme neutre des YFP T203Y est considérée comme non émissive, donc le terme devrait être nul. Mais s'il ne l'est pas, il produit une composante **qui croît avec \[Cl⁻\]** et qui mime du signal. Garde dans le code un terme ε\_dark·Φ\_dark, par défaut à zéro, pour pouvoir chiffrer ce que coûterait une valeur non nulle.

## Générateur direct : de \[Cl⁻\] aux deux histogrammes

Sept étapes, de la concentration aux comptages entiers. L'idée structurante : **tout est une somme d'exponentielles**, donc chaque canal se représente par un vecteur de couples (amplitude, durée de vie), et une seule routine de convolution sert aux deux canaux et à tous les estimateurs.

1. **Fraction active.** x = x(C) par l'isotherme à deux plateaux.
2. **Amplitudes analytiques par canal.** Développer les expressions de la section précédente en une somme sur les trois durées de vie {τ\_DA, τ\_D0, τ\_A}. Le canal bleu comme le jaune se réduisent à trois couples (α, τ), dont l'un porte une amplitude négative. Aucune intégration numérique ici.
3. **Déclin idéal périodique.** Le laser tire toutes les 12,5 ns et le déclin ne s'éteint pas dans la fenêtre. Sommer analytiquement la série géométrique des impulsions précédentes plutôt que de tronquer :

```latex
I^{\mathrm{per}}(t) = \sum_{\alpha,\tau} \alpha\,\frac{e^{-t/\tau}}{1 - e^{-T_{\mathrm{rep}}/\tau}},\qquad 0 \le t < T_{\mathrm{rep}}
```

4. **Convolution avec l'IRF.** IRF chargée depuis ton fichier, rééchantillonnée sur la grille des 256 canaux, normalisée à somme unité, et **convolution circulaire** sur la période — pas linéaire, sinon tu réintroduis le bord que l'étape 3 vient d'éliminer. Un décalage t₀ reste un paramètre du modèle même s'il vaut zéro par défaut : les estimateurs doivent pouvoir le laisser libre.
5. **Binning.** 256 canaux, Δt = 12,5/256 ≈ 48,8 ps. Intégrer le modèle sur chaque canal plutôt que l'échantillonner au centre — l'échantillonnage ponctuel biaise les premiers canaux, là où la montée vit.
6. **Mise à l'échelle photonique.** Normaliser chaque canal pour que la somme des comptages attendus reproduise le budget demandé, selon la normalisation choisie (voir le plan d'expériences). Le rapport bleu/jaune n'est **pas** libre : il sort du modèle, c'est une prédiction.
7. **Bruit de Poisson.** Tirage `rand(Poisson(μ_k))` canal par canal. Même germe pour les deux canaux d'une même réalisation, germes indépendants entre réalisations.

**Deux garde-fous dans le générateur.** Les amplitudes négatives sont physiques dans le modèle mais les comptages attendus doivent rester positifs partout : lever une erreur si un μ\_k devient négatif, ce qui signalerait une incohérence de paramètres. Et le générateur doit pouvoir produire le déclin **sans bruit** : c'est ce mode qui sert à construire les courbes de calibration de chaque estimateur et à valider le code.

## Inférer le canal jaune depuis ta titration en canal bleu

Ta titration bleue contient quatre des neuf paramètres du modèle, et elle les contient **mieux** que des mesures séparées, parce qu'ils y sont liés sur toute la série. Trois paramètres manquent et devront venir d'ailleurs.

### L'ajustement global à faire sur tes données

Ajuster **simultanément** tous les histogrammes de la titration avec τ\_D0 et τ\_DA **liés entre toutes les concentrations**, et une seule amplitude x\_i libre par concentration. C'est l'ajustement que le modèle M2 impose, et il a moins de paramètres libres que ce que tu fais actuellement concentration par concentration.

Il te rend : τ\_D0, τ\_DA, et la série {x\_i}. Puis un ajustement de Langmuir sur {x\_i, C\_i} donne x₀, x\_sat, K\_D, n.

### Ce test passe ou échoue, et son échec est informatif

Si un modèle à deux exponentielles à durées de vie globalement liées décrit toute ta titration avec un χ²ᵣ proche de 1, M2 lent est confirmé et le simulateur repose sur du solide. Si le χ²ᵣ reste élevé, quelque chose d'autre domine et il faut me le dire avant d'aller plus loin.

Deux diagnostics complémentaires, tous deux gratuits sur tes données actuelles :

- **Le profil de χ²ᵣ de ton ajustement mono-exponentiel, point par point.** Sous M2 lent, il doit être au meilleur à Cl⁻ saturant (vrai mono), le pire à concentration intermédiaire (mélange 50/50), et intermédiaire à Cl⁻ nul. Tu m'as donné 7,5 en moyenne ; sors-le point par point. Un profil monotone ou plat contredit le modèle.
- **τ\_D0 comparé à du mCer3 nu.** Ton mono à Cl⁻ saturant est censé être τ\_D0. Un écart résiduel de 0,2–0,3 ns signifie que la saturation est incomplète et te donne x\_sat directement.

### Les trois paramètres que le canal bleu ne contient pas

| Manquant | Source | Sensibilité du résultat |
| --- | --- | --- |
| τ\_A | mTopaz seul, ou littérature YFP (3,0 ns) | modérée — balayer ±15 % |
| ε (excitation directe) | rapport des spectres 2P à 800 nm | faible si ε < 0,05, à vérifier |
| η (matrice spectrale) | spectres FPbase + bandes de filtres | **forte** — c'est le paramètre de nuisance principal |

Les η se calculent à la construction, une fois, depuis l'API REST de FPbase. Le spectre d'émission de mCerulean3 et celui de Topaz s'obtiennent par `https://www.fpbase.org/api/proteins/<slug>/?include_fields=states__em_spectrum&format=json`, qui renvoie des paires \[longueur d'onde, valeur\]. Chaque η est l'intégrale du spectre normalisé sur la bande du filtre, en supposant des filtres rectangulaires idéaux (460–500 et 520–560 nm) et une détection plate.

**Garde ce chiffre en tête comme contrôle :** ton propre calcul sur tes données avait abouti à une fuite du donneur représentant environ 24 % du canal jaune. Si l'intégration des spectres te donne un ordre de grandeur comparable, les deux approches se confirment. Un écart important signale une erreur de ma part ou une bande de filtre différente de ce que je suppose.

### Validation croisée gratuite

Une fois les paramètres extraits, le générateur prédit le **nombre total de photons du canal bleu** en fonction de \[Cl⁻\], indépendamment de la forme du déclin. Compare cette prédiction à tes comptages mesurés : tu obtiens une vérification du modèle qui ne coûte aucune expérience nouvelle. Si elle échoue, le modèle a un défaut avant même qu'on parle du canal jaune.

**Ce qu'il me faut de toi :** les histogrammes de la titration, l'IRF, un exemple de ton ajustement actuel, et les bandes exactes de tes filtres si elles diffèrent de 480/40 et 540/40 en largeur à mi-hauteur nominale.

## Les cinq estimateurs

Ce qui distingue ces estimateurs n'est pas la formule mais **le nombre de paramètres libres** et **les données qu'ils consomment**. C'est là que se joue le comportement à bas photons.

| # | Estimateur | Canaux | Paramètres libres | Observable |
| --- | --- | --- | --- | --- |
| 1 | Ratio d'intensités | bleu + jaune | 0 | R = N\_jaune / N\_bleu |
| 2 | τ mono-exponentiel MLE | bleu | 2 (τ, A) | τ\_apparent |
| 3 | Bi-exponentiel MLE, τ libres | bleu | 4 (τ₁, τ₂, A₁, A₂) | fraction d'amplitude |
| 4 | Bi-exponentiel MLE, τ figés | bleu | 2 (A₁, A₂) | fraction d'amplitude |
| 5 | **Conjoint deux canaux, τ figés** | bleu + jaune | 3 (x, échelle\_B, échelle\_J) | x directement |
| 6 | Phasor | bleu | 0 | position (g, s) |

L'estimateur 4 est ajouté exprès : il isole ce que gagne le fait de **figer les durées de vie**, séparément de ce que gagne le **deuxième canal**. Sans lui, un avantage de l'estimateur 5 serait inattribuable, et la comparaison ne prouverait rien.

### Détails qui comptent

**Vraisemblance.** Poisson, pas moindres carrés — à 300 photons sur 256 canaux, la plupart des canaux sont à 0, 1 ou 2 comptages et l'approximation gaussienne échoue exactement dans le régime d'intérêt. Maximiser la log-vraisemblance de Poisson sur les deux canaux simultanément pour l'estimateur 5.

**Variantes de l'estimateur 5.** Trois versions, parce que le choix de ce qu'on fige *est* le résultat :

- **5a** — τ\_D0, τ\_DA, τ\_A figés aux vraies valeurs. Borne optimiste.
- **5b** — figés à des valeurs volontairement fausses (±5 %, ±10 %). Mesure ce que coûte une calibration imparfaite, qui est le cas réel.
- **5c** — τ\_A laissé libre. Teste le mauvais conditionnement attendu entre τ\_A (3,0 ns) et τ\_DA (\~2,1 ns), qui ne sont séparés que de 0,9 ns. Surveiller le coefficient de corrélation dans la matrice de covariance ; au-delà de 0,9, l'information est illusoire.

**Chemin vers \[Cl⁻\].** Chaque estimateur reçoit sa **propre** courbe de calibration, construite en passant le générateur en mode sans bruit sur la grille de concentrations et en enregistrant la valeur de son observable. L'estimation se fait ensuite par inversion de cette courbe. C'est exactement ce que tu ferais au banc, et c'est ce qui produit les K\_apparents divergents.

**Phasor.** Inclus parce qu'il est rapide, donc compatible avec ta boucle de rétroaction. Mais il faut en attendre peu : il projette le déclin sur deux coordonnées, ce qui suffit pour un mono-exponentiel et pas pour un mélange à trois composantes. Son intérêt réel est que sous M2, les points d'une titration doivent tomber sur un **segment de droite** entre les positions de τ\_D0 et τ\_DA dans le plan phasor, la position le long du segment donnant x. Si tes données expérimentales font ça, c'est une confirmation visuelle immédiate de M2 — et c'est un test que tu peux faire cette semaine sur ta titration existante, sans rien ajuster.

## Borne de Cramér-Rao

La borne dit combien d'information existe dans les deux histogrammes, indépendamment de l'estimateur. Sans elle, tu ne peux pas distinguer « cet estimateur est mauvais » de « cette information n'existe pas ».

Les comptages par canal sont poissoniens et indépendants, donc l'information de Fisher s'écrit directement, en sommant sur les canaux temporels **des deux détecteurs** :

```latex
\mathcal{I}_{ij}(\theta) = \sum_{c \in \{B,\,J\}} \; \sum_{k=1}^{256} \frac{1}{\mu_k^{(c)}}\, \frac{\partial \mu_k^{(c)}}{\partial \theta_i}\, \frac{\partial \mu_k^{(c)}}{\partial \theta_j}
```

Les dérivées sont analytiques : le modèle est une somme d'exponentielles convoluée par une IRF fixe, donc ∂μ/∂θ s'obtient en dérivant les amplitudes et les taux avant convolution, puis en convoluant. Pas de différences finies — elles sont instables sur les paramètres de durée de vie.

**Propagation vers la concentration.** La variance qui t'intéresse est en mM :

```latex
\sigma^2_{[\mathrm{Cl}^-]} \;=\; \left(\frac{dx}{d[\mathrm{Cl}^-]}\right)^{-2} \left[\mathcal{I}^{-1}\right]_{xx}
```

Le facteur dx/d\[Cl⁻\] est le plus petit aux extrémités de l'isotherme, et c'est lui qui crée la zone de faible sensibilité à tes concentrations de travail. La borne en mM, et non en ns, est la seule qui t'informe utilement.

**Trois usages concrets au-delà de la comparaison d'estimateurs :**

1. **Diagnostic d'identifiabilité.** Calculer la borne avec θ = (x, τ\_D0, τ\_DA, τ\_A, échelles) et inspecter les corrélations. Si la corrélation entre τ\_A et τ\_DA dépasse 0,9, tu sais avant d'écrire un estimateur que les laisser libres tous les deux ne marchera pas.
2. **Le R de ton Kalman.** La borne évaluée aux paramètres estimés donne la variance de mesure à chaque instant, gratuitement, et elle varie avec le nombre de photons effectivement collecté et avec la position sur l'isotherme. C'est strictement meilleur qu'un R constant.
3. **Efficacité des estimateurs.** Le rapport σ²\_Monte-Carlo / σ²\_CRLB. Un estimateur à efficacité 0,3 laisse 70 % de l'information sur la table et vaut la peine d'être amélioré ; un estimateur à 0,95 est au bout de ce qui est possible et il faut chercher ailleurs.

## Plan d'expériences numériques

### Les trois normalisations de budget, et pourquoi il en faut trois

Comparer un estimateur qui n'utilise qu'un canal à un estimateur qui en utilise deux attribue à la méthode un gain qui vient simplement du nombre de photons. Tu m'as dit que ce gain t'intéresse quand même — il reste un gain réel dans ton régime. Mais il faut pouvoir le séparer, sinon la conclusion est inattaquable et inutile.

| Normalisation | Ce qu'elle fixe | Ce qu'elle mesure |
| --- | --- | --- |
| **A — temps d'acquisition** | durée de la tornade | le gain pratique, tout compris. **Référence principale** : c'est ta contrainte réelle |
| **B — photons totaux** | N\_bleu + N\_jaune | le gain méthodologique seul, photons égalisés |
| **C — photons du canal bleu** | N\_bleu | le gain pur du deuxième canal, comme bonus par-dessus l'existant |

Chaque figure doit porter sa normalisation en toutes lettres. C'est la première chose qu'un évaluateur cherchera.

### Grilles de balayage

- **Concentration :** 0, 1, 2, 5, 8, 12, 16, 20, 25, 30, 40, 60, 100 mM — densité accrue dans ta fenêtre de clamping.
- **Budget photonique :** 10², 3·10², 10³, 3·10³, 10⁴, 3·10⁴, 10⁵. Ton point de fonctionnement actuel (10⁴ par 0,1 s) est au milieu ; les deux décades sous lui sont le cœur de la question.
- **Répétitions Monte-Carlo :** 2 000 par cellule en dessous de 10³ photons, 1 000 au-dessus. À 2 000 réalisations, l'incertitude relative sur σ est d'environ 1,6 %, suffisante pour distinguer des estimateurs qui diffèrent de 10 %.

### Balayages de nuisance

Tous à 10³ et 10⁴ photons, normalisation A :

1. **Fuite spectrale η\_YD**, de 0 à 1,5× la valeur nominale. Le paramètre le plus important : il doit dégrader fortement le ratio et peu l'ajustement conjoint. Si ce n'est pas ce qu'on observe, l'hypothèse de départ du projet tombe.
2. **Excitation directe ε**, de 0 à 0,10.
3. **Mauvaise spécification des durées de vie figées**, ±5 % et ±10 % sur τ\_D0, τ\_DA, τ\_A séparément. Mesure la robustesse à une calibration imparfaite.
4. **Mauvaise spécification du modèle**, le cas le plus instructif : générer sous M2 et ajuster avec un mono-exponentiel. Chiffre en mM l'erreur que fait ton pipeline actuel.
5. **Terme ε\_dark non nul**, de 0 à 0,05, pour borner le risque d'excitation 2P du chromophore protoné à 800 nm.

### Deux conditions de contrôle

**Le cas dégénéré.** η\_YD = 0 et ε = 0 : le canal jaune devient de l'émission sensibilisée pure. Sert à vérifier que le code retrouve la situation idéale de la littérature, et à mesurer combien la fuite réelle coûte par rapport à cet idéal.

**La cible mobile.** Les estimateurs sont comparés à concentration fixe, mais ton usage réel est un asservissement. Une passe finale optionnelle : une trajectoire de \[Cl⁻\] (échelon de 5 à 25 mM, constante de temps de 20 s), échantillonnée à 10 Hz, pour voir si un estimateur à faible variance mais à biais dépendant de la concentration produit une erreur statique en boucle fermée.

## Sorties et figures

Six figures. Chacune répond à une question que tu peux poser à voix haute, et chacune porte sa normalisation de budget en légende.

| # | Figure | Question à laquelle elle répond |
| --- | --- | --- |
| 1 | σ\[Cl⁻\] en fonction de \[Cl⁻\], une courbe par estimateur, plus la borne CRLB, à 10⁴ photons | Lequel gagne, et où sur l'isotherme ? |
| 2 | σ\[Cl⁻\] en fonction du budget photonique, à 20 mM, en log-log | Lequel tient à bas photons ? La pente doit valoir −0,5 ; un décrochage signale une rupture |
| 3 | Biais en fonction du budget photonique, mêmes conditions | Lesquels développent un biais quand les photons manquent ? |
| 4 | Biais en fonction de η\_YD | Le ratio s'effondre-t-il là où l'ajustement conjoint tient ? **C'est la figure centrale du projet** |
| 5 | K\_apparent / K\_vrai par estimateur | De combien chaque méthode déforme la constante de dissociation |
| 6 | Efficacité (σ²\_MC / σ²\_CRLB) en fonction du budget | Quel estimateur laisse de l'information sur la table, et où |

### Sortie numérique

Un seul fichier plat, une ligne par cellule de la grille : `[Cl]_vrai, estimateur, N_photons, normalisation, paramètre_nuisance, valeur_nuisance, biais_mM, σ_mM, RMSE_mM, CRLB_mM, efficacité, taux_échec`.

Le **taux d'échec** est une colonne à part entière, pas une note de bas de page : à 300 photons avec quatre paramètres libres, une fraction des ajustements ne convergera pas ou convergera sur une frontière. Un estimateur qui échoue 20 % du temps n'est pas « légèrement moins bon », il est inutilisable en boucle fermée, et cette information disparaît si tu ne rapportes que σ sur les réalisations réussies.

### Le résultat qu'on cherche vraiment

Une carte en deux dimensions — budget photonique en abscisse, concentration en ordonnée — indiquant quel estimateur minimise la RMSE dans chaque région. Si l'analyse conjointe gagne partout, c'est une conclusion. Si elle gagne seulement sous un certain seuil de photons, ou seulement au-delà d'une certaine fuite spectrale, c'en est une meilleure, parce qu'elle te dit **quand** l'utiliser et justifie la complexité supplémentaire auprès de ton directeur.

## Architecture Julia

Le type central est une somme d'exponentielles. Tout le reste s'y ramène, ce qui veut dire une seule routine de convolution, une seule routine de dérivation analytique, et des estimateurs qui diffèrent uniquement par les paramètres qu'ils laissent libres.

```julia
struct ExpSum                      # un canal, entièrement décrit
    α::Vector{Float64}             # amplitudes, certaines négatives
    τ::Vector{Float64}             # ns
end

struct SensorParams
    τ_D0::Float64; τ_DA::Float64; τ_A::Float64   # ns, constantes
    x₀::Float64; x_sat::Float64                  # plateaux de l'isotherme
    K_D::Float64; n::Float64                     # mM, sans unité
    Φ_D::Float64; Φ_A::Float64
end

struct OpticsParams
    η::SMatrix{2,2,Float64}        # [canal, émetteur] : BD BA / JD JA
    ε::Float64                     # excitation directe de l'accepteur
    ε_dark::Float64                # terme de risque, défaut 0.0
end

struct InstrumentParams
    n_bins::Int                    # 256
    T_rep::Float64                 # 12.5 ns
    irf::Vector{Float64}           # chargée, normalisée à somme 1
    t₀::Float64                    # décalage, défaut 0.0
end
```

### Modules et responsabilité unique de chacun

| Module | Rôle |
| --- | --- |
| `Spectra.jl` | récupère les spectres FPbase, calcule η, met en cache en CSV. Exécuté une fois |
| `Photophysics.jl` | `x(C, p)` et `expsums(C, p, o) -> (ExpSum, ExpSum)`. Aucune notion d'instrument |
| `Instrument.jl` | périodisation, convolution circulaire, intégration par canal. `forward(es, inst) -> Vector{Float64}` |
| `Generate.jl` | mise à l'échelle photonique, tirage de Poisson, gestion des germes |
| `Estimators.jl` | les six estimateurs, interface commune `estimate(h_bleu, h_jaune, cfg) -> (x̂, converged)` |
| `Fisher.jl` | dérivées analytiques, matrice d'information, CRLB en mM |
| `Calibrate.jl` | courbes de calibration sans bruit, inversion |
| `Ingest.jl` | ajustement global de ta titration bleue, extraction des paramètres |
| `Experiments.jl` | balayages, parallélisation, écriture du fichier plat |

### Trois décisions d'implémentation

**Convolution circulaire par FFT.** La périodisation à 12,5 ns rend la convolution naturellement circulaire. Une FFT sur 256 points, planifiée une fois et réutilisée, est à la fois plus rapide et plus juste qu'une convolution linéaire tronquée.

**Dérivées analytiques, pas automatiques.** La différentiation automatique fonctionnerait, mais les dérivées d'une somme d'exponentielles sont triviales à écrire à la main et le code devient lisible pour un évaluateur. Réserve l'AD comme test : les deux doivent coïncider à 10⁻¹⁰ près, et c'est un bon test unitaire.

**Parallélisation sur les réalisations Monte-Carlo.** `Threads.@threads` sur la boucle des répétitions, germe dérivé de l'indice de réalisation pour la reproductibilité. Les balayages complets représentent de l'ordre de 10⁶ ajustements ; sans parallélisation, c'est des heures.

### Point d'accroche avec ton code existant

`Estimators.jl` doit exposer l'estimateur 2 (mono-exponentiel MLE) derrière la même interface que les autres, mais appelant **ta fonction existante** plutôt qu'une réimplémentation. Deux raisons : tu valides ainsi ton propre code contre une vérité terrain connue, et la comparaison porte sur ce que tu utilises réellement, pas sur une version idéalisée. Si ta fonction a une signature différente, un adaptateur de dix lignes suffit.

## Validation interne

Huit tests. Tant qu'ils ne passent pas tous, aucun résultat du simulateur ne vaut la peine d'être interprété — un simulateur silencieusement faux produit des figures convaincantes.

1. **Aller-retour sans bruit.** Générer sans Poisson, estimer avec le modèle exact : les paramètres doivent revenir à la précision machine. Échoue si la convolution, le binning ou la périodisation ont un défaut.
2. **Conservation des photons.** L'intégrale de l'émission sensibilisée doit valoir Φ\_A·x·E, celle du donneur Φ\_D·\[x·τ\_DA/τ\_D0 + (1−x)\]. Vérifiable analytiquement, donc testable exactement.
3. **Amplitudes égales aux fractions.** Dans le canal bleu idéal, les amplitudes à t = 0 valent x et (1−x). C'est le fondement du modèle ; s'il ne tient pas dans le code, le code est faux.
4. **Limite à Cl⁻ saturant.** Quand x → x\_sat avec x\_sat = 0, le canal bleu doit devenir rigoureusement mono-exponentiel à τ\_D0, et le canal jaune doit disparaître. Reproduit ton observation expérimentale.
5. **Limite dégénérée.** Avec η\_YD = 0 et ε = 0, un ajustement bi-exponentiel du canal jaune doit retrouver exactement une amplitude négative à τ\_DA et une positive à τ\_A, avec un rapport d'amplitudes de 1 en valeur absolue.
6. **Efficacité asymptotique.** À 10⁶ photons, le MLE correctement spécifié doit atteindre la CRLB à 2 % près. Si l'efficacité plafonne plus bas, c'est l'optimiseur ou le calcul de Fisher qui a un défaut, pas la physique.
7. **Dérivées analytiques contre AD.** Coïncidence à 10⁻¹⁰ sur un jeu de paramètres aléatoires.
8. **Linéarité phasor.** Sous M2, la titration doit tracer un segment de droite dans le plan phasor. Un test géométrique simple qui vérifie d'un coup le modèle et le générateur.

**Ce que je ne propose pas.** Reproduire des chiffres publiés — tu as dit que ce n'était pas nécessaire, et les tests 1 à 8 sont de toute façon plus contraignants qu'une comparaison à une référence externe dont on ne contrôle pas les conditions.

## Exclusions assumées

Ce qui est hors de cette version, et ce que chaque exclusion te coûte. Les deux premières ne sont pas neutres.

| Exclu | Effet sur la conclusion | Priorité de réintégration |
| --- | --- | --- |
| **Photoblanchiment** | **sous-estime l'avantage de l'analyse conjointe** | 1 |
| **Décalage temporel inter-canaux** | **surestime l'avantage de l'analyse conjointe** | 2 |
| Tag YFP d'eNpHR | sous-estime l'avantage | 3 |
| Variation de pH | ignore un mode de défaillance réel du capteur | 4 |
| Autofluorescence cellulaire | optimiste pour tous les estimateurs également | 5 |
| Donneur ou accepteur libres | optimiste pour tous également | 5 |
| Anisotropie, distribution de distances | effet faible, surtout sur les premiers canaux | 6 |
| Pile-up, temps mort, optique non idéale | négligeable à 100 kcps sur 80 MHz | 6 |

### Les deux biais qui ne s'annulent pas

**Le photoblanchiment manque dans le sens qui te défavorise.** C'est l'un des deux arguments les plus forts pour la méthode duale : un accepteur qui blanchit change l'amplitude de l'émission sensibilisée mais pas sa constante de temps, alors qu'il déplace le ratio d'intensités sans recours. Tes propres observations en tranche vont dans ce sens — ratio instable sous illumination répétée à 594 nm, durée de vie stable. Donc tout résultat favorable que produira ce simulateur est un **plancher**, pas une estimation centrale.

**Le décalage inter-canaux manque dans le sens qui te favorise.** La montée dure quelques centaines de ps. Un décalage relatif non calibré entre tes deux HPM-100-40-C — longueurs de câble, tensions, lignes à retard internes — se transfère directement en biais sur tout ce que le canal jaune apporte. Le simulateur à décalage nul décrit donc un banc parfaitement calibré, ce que ton banc n'est pas encore. C'est la première chose à ajouter, et la mesure correspondante au banc coûte vingt minutes avec un diffuseur vu simultanément dans les deux canaux.

**Mis ensemble :** les deux premières lignes du tableau pointent en sens opposés, et aucune des deux n'est petite. Ne présente pas un résultat de ce simulateur comme une estimation de performance réelle tant qu'elles ne sont pas toutes les deux dedans. Présente-le comme ce qu'il est : une comparaison d'estimateurs sous un modèle contrôlé, où le classement est fiable et les valeurs absolues ne le sont pas.

### Le commutateur à garder

Garde dans `Photophysics.jl` une option de mécanisme : `:M2_lent` (défaut, mélange statique) et `:M2_rapide` (échange rapide, un seul τ qui glisse). Elle coûte vingt lignes et te donne un test de discrimination de mécanisme prêt à l'emploi le jour où tu auras le canal jaune réel. Si tes histogrammes jaunes montrent une montée dont la constante de temps varie avec \[Cl⁻\], c'est M2 rapide et toute la section « modèle physique » est à refaire.
