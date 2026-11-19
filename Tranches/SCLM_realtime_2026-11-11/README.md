# SCLM — session FLIMApp en temps réel, 2026-11-11 (4 stims à 1064 nm, PWM)

Une session FLIMApp en mode **Realtime par passes** (`2026-11-11_160751`, 580 s) : la carte
SPC-QC-104 date chaque photon pendant que le polygone du microscope balaie le champ, et la NI
envoie le protocole du laser à 1064 nm, **4 stims de 10 s à 100 %**, séparées de 130 s, en
**PWM 50 % à 5 kHz** (un échantillon de 100 µs de laser opto, un échantillon d'acquisition).
Protocole en volts : boucle ouverte, le PI n'agissait pas.

Analysé avec le pipeline `../sclm_realtime/` (voir son README pour la méthode) :

```bash
cd ../sclm_realtime
python3 step1_decode.py ../SCLM_realtime_2026-11-11     # ~40 s : photons rangés par passe, contrôle frames.csv
python3 step2_fuite.py ../SCLM_realtime_2026-11-11      # ce que le 1064 nm ajoute, ce qui passe dans la mesure
python3 step3_series.py ../SCLM_realtime_2026-11-11     # ~1 min : lectures passe par passe, dérive, bruit
python3 step4_reponses.py ../SCLM_realtime_2026-11-11   # réponses, témoins, modèles cinétiques
python3 step5_absolu.py ../SCLM_realtime_2026-11-11     # [Cl⁻] absolue avec la calibration du 2026-11-09 ?
python3 step6_cinetique.py ../SCLM_realtime_2026-11-11  # latence, τ d'activation et τ de retour, par stim
python3 step7_hybride.py ../SCLM_realtime_2026-11-11    # [Cl⁻] dans le temps : niveau de τ, variations du rapport
python3 step8_images.py ../SCLM_realtime_2026-11-11     # ~1 min : images reconstruites sans M1/M2
~/.venvs/cellpose/bin/python step9_cellpose.py ../SCLM_realtime_2026-11-11   # ROIs (Cellpose)
python3 step10_rois.py ../SCLM_realtime_2026-11-11      # ~1 min : par ROI, avant / après, cinétique
```

## Ce que contient `data/`

`data/FLIM` → `~/Documents/Maîtrise/Données/2026-11-11_160751` (lien non versionné, 1,3 Go) :
`run.toml`, `frames.csv` (ce que FLIMApp a mesuré en direct, une ligne par passe),
`visits.csv` (commande de chaque créneau), `readback.bin` (signaux relus par la NI à 10 kHz),
`irf.csv` (l'IRF `irf_850nm_16X.sdt` du 2026-11-09), `spc/3T0089_IN{1,2}.spc` (flux FIFO).

| | |
|---|---|
| Passes | 2416 : scan de 210 ms + pause de 30 ms (créneau de 240 ms) |
| Stims | début à 60,1, 190,2, 320,3 et 450,1 s ; 42, 42, 41 et 42 passes ; commande 0,986 V (100 %), porte 5 V |
| PWM | 5 kHz, 50 % : pendant chaque scan, 100 µs de 1064 nm (code de routage 0, photons jetés par FLIMApp), puis 100 µs d'acquisition (code 1) |
| Canaux | **IN1 = YFP** (6,1 fois plus de photons, ~2,9 ns), **IN2 = CFP** (~1,5 ns) : comme le 2026-11-09, l'inverse d'octobre |
| Photons | YFP 362 kHz, CFP 59 kHz (toute la session) ; par passe, en acquisition : 37 860 YFP et 6 160 CFP ; R = 6,15 |

## Contrôles

- **Décodage exact** : les photons d'acquisition de chaque passe sont ceux de `frames.csv`, pour
  les 2416 passes et les deux canaux (`qc/validation_flimapp.csv`) ; ni GAP ni débordement ;
  2416 M0 et 2416 M3 par carte. 0,04 % des photons de code 1 tombent hors des passes (décalage
  d'un échantillon entre le code et le signal de passe, déjà signalé par FLIMApp).
- **Stims** : la commande relue par la NI (`readback.bin`) et `visits.csv` concordent. La
  colonne `command1` de `frames.csv` est décalée d'une passe au bord de 3 blocs (5 passes) : ne
  pas s'en servir pour situer les stims.
- **Ce qui a été balayé — tout le champ, au polygone.** Le mode ROI de FLIMApp était éteint
  (`[roi] active = false`) : ses galvos sont restés à 0 V (readback). Mais le débit de photons
  oscille à **17 997,3 Hz** (ligne de 55,564 µs) et à son harmonique, et à **33,25 Hz** (trame
  de 540 lignes) : le polygone balayait le champ entier, ce que l'utilisateur confirme. Les
  débits sont d'ailleurs ceux de l'image complète prise 50 min plus tôt (15:14 : YFP 411 kHz,
  CFP 69 kHz, rapport 5,95).
- **Pas de M1/M2 dans les fichiers, mais des images reconstruites.** En mode Realtime par
  passes, FLIMApp coupe exprès les marqueurs de ligne (M1) et de trame (M2)
  (`parametres_clamp`, `routing_mode = 6400` : M0 et M3 seulement) ; il n'y en a aucun dans les
  deux fichiers. L'étape 8 retrouve les images en suivant la phase des lignes du polygone trame
  par trame, sur la structure de l'image elle-même : 12 somas séparés par Cellpose (étapes 9 et
  10). Les étapes 2 à 7 portent sur le champ entier (somas, neuropile, fond), une seule série.

## Ce que disent les résultats

### En bref

- **Le multiplexage par PWM protège la mesure.** Le 1064 nm excite la YFP (+9,0 % de photons
  YFP pendant ses phases, rien sur le CFP), mais ces photons tombent dans les phases jetées. Il
  n'en passe que **0,07 %** dans l'acquisition, dans la première µs après la coupure.
- **Chaque stim baisse le rapport YFP/CFP de ~1,2 %**, 4 stims sur 4 (z = −2,7 à −4,0 ; témoins
  entre −0,3 et +0,3 %). Le CFP monte, le YFP baisse, τ CFP s'allonge de ~10 ps : moins de FRET,
  le sens d'une entrée de Cl⁻.
- **Cinétique** : latence de **~4 s** (IC 3,2–4,8 s), descente en rampe pendant la stim (τ
  d'activation non mesurable sur 10 s), minimum ~4 s après la fin, retour avec **τ ≈ 16 s**
  (IC 10–33 s). Le rapport revient entièrement : aucune part permanente (α = 1,00 ± 0,35).
- **Le rapport est limité par le bruit de photons** passe par passe (résidus/attendu 0,99) :
  1,35 % par passe de 240 ms, ~0,7 % par seconde. τ CFP et tm CFP le sont aussi, mais ~10 fois
  moins précis en mM-équivalent.
- **Ordre de grandeur** : avec les pentes de la calibration du 2026-11-09, une stim vaut
  +0,2 à +0,4 mM selon le rapport et +0,5 à +1,2 mM selon τ CFP.
- **Pas de [Cl⁻] absolue.** Avec cette calibration, le rapport lit −2,6 mM et τ CFP 5,7 mM au
  même moment, et leurs dérives diffèrent de 13 fois. Une calibration dans ce mode est
  nécessaire.
- **Lecture hybride (étape 7)** : niveau ancré sur τ, variations au rapport. Elle a la précision
  du rapport (0,3 mM par seconde, contre 6 mM pour τ seul). Mais sa montée de 5,7 à 14 mM vient
  de τ, et τ monte 6 fois trop pour du Cl⁻ par rapport à la baisse du rapport : signature d'un
  blanchiment de l'accepteur. Seule la première ancre (5,7 ± 0,6 mM) est un niveau défendable.
- **Par cellule (étapes 8 à 10)** : images reconstruites sans M1/M2, 12 somas. Le rapport baisse
  dans la majorité des somas à chaque stim (médiane −0,8 à −1,9 % ; p = 0,002 à la stim 4), soit
  +0,8 mM en moyenne sur 44 paires ROI × stim (rapport ancré sur le τ de chaque ROI ; 34/44 ↑),
  +0,7 mM dérive lente retirée. Somas à 13,7 ± 2,0 mM en moyenne avant la 1re stim (5,6 à 29 mM
  selon le soma). Les constantes de temps ne se lisent pas cellule par cellule (trop peu de
  photons).
- **À régler dans FLIMApp** : ch1 = IN1 = YFP, donc la durée de vie « ch1 » suivie en direct
  est celle du YFP. L'IRF du 2026-11-09 arrive 0,13 à 0,19 ns après les déclins d'aujourd'hui.
  Le Realtime en images (M1/M2 enregistrés) éviterait la reconstruction.

### Le laser à 1064 nm et le PWM (étape 2)

| | YFP (IN1) | CFP (IN2) |
|---|---|---|
| photons ajoutés dans les phases opto (stims 1–4) | +8,3, +9,0, +9,2, +9,0 % (± 0,2–0,3) | −0,5 à +1,1 % (± 0,5) : rien |
| photons ajoutés par passe | ~3 400 | — |
| forme dans le déclin | plate (103 ± 0,4 % plat) : non synchrone du 850 nm | — |
| évolution pendant les 10 s | aucune tendance nette (0,00 à 0,08 point de % par s) | — |
| passé dans l'acquisition (profil du PWM) | 0,066 % des photons d'acquisition, dont 0,05 % dans la 1re µs | 0,0 % |
| retard d'allumage du 1064 nm | ~10 µs après le changement de code | |

Référence : les photons des pauses de la même passe (850 nm seul), qui suivent la fluorescence
au même moment. Le 1064 nm excite donc la YFP en deux photons (YFP du capteur et, s'il y en a,
YFP de eNpHR), pas le CFP ; plat dans le déclin parce que ses impulsions ne sont pas
synchrones du 850 nm. Sur la mesure pendant une stim, la fuite biaiserait le rapport de
+0,07 % au plus (~0,02 mM-éq.), loin sous la réponse ; les fenêtres de 10 s avant et après
n'ont pas de stim. Figure : `results/fuite.png`.

### Lectures passe par passe, dérive et bruit (étape 3)

**Déclins de toute la session** (`results/declins.png`) : CFP mono-exponentiel 1,55 ns
(déviance/ddl 1 484 sur 15 M de photons) ; à trois exponentielles, une composante prompte
(4,7 % des photons CFP, 0,6 % des YFP), 0,72 et 2,84 ns. **L'IRF de FLIMApp (celle du
2026-11-09) doit être avancée de 128 à 193 ps** selon le modèle (YFP : 138 à 158 ps) :
les déclins d'aujourd'hui montent avant elle. Le décalage est ajusté une fois sur toute la
session, puis fixé passe par passe.

**Bruit** (`results/bruit.png`) : toutes les lectures hors ligne sont au bruit de photons, par
passe de 240 ms.

| lecture | écart-type par passe | bruit de photons | rapport |
|---|---|---|---|
| R (YFP/CFP) | 0,083 (1,35 %) | 0,084 | 0,99 |
| tm CFP | 30 ps | 31 ps | 0,97 |
| τ CFP mono | 50 ps | 51 ps | 0,98 |
| τ YFP mono | 38 ps | 36 ps | 1,05 |
| τ CFP de FLIMApp en direct | 58 ps | 51 ps | 1,14 |

L'écart d'Allan du rapport suit le bruit de photons jusqu'à ~4 s, puis s'aplatit (~0,3 % de
16 à 32 s) : des fluctuations lentes s'ajoutent au-delà de quelques secondes. Les photons YFP
seuls ont un excès de bruit (χ²/ddl 1,44 dans les modèles de l'étape 4) que le rapport annule :
des fluctuations communes aux deux canaux (puissance, focus).

**Dérive** (moments calmes, un niveau libre par intervalle) : rapport −1,2 ± 0,2 %/100 s, YFP
−1,6 %/100 s, CFP −0,4 %/100 s, tm CFP +0,25 %/100 s ; τ CFP passe de 1,50 à 1,57 ns en 8 min.
La composante prompte ne l'explique pas (5,0 → 4,6 % des photons) : ce sont les durées de vie
du donneur qui s'allongent (0,70 → 0,74 ns et 2,82 → 2,88 ns). Le YFP baisse quatre fois plus
vite que le CFP : signature d'un **blanchiment de l'accepteur** par l'imagerie continue (ou
d'une lente montée de Cl⁻ ; les deux se ressemblent).

**Fin de session** (après 495 s, écartée des étapes 4 à 6) : le rapport monte de 4,8 % parce
que le CFP baisse de 5,1 % sans que le YFP (−0,6 %) ni τ CFP (−0,6 %) ne bougent : pas un
changement de FRET. Plutôt un changement optique du champ (focus ?).

**FLIMApp en direct** : sa durée de vie CFP (ch2) est à +8 ps de l'ajustement hors ligne
(écart-type 24 ps, r = 0,92), 14 % plus bruitée que le bruit de photons ; sa durée de vie
« ch1 » est celle du **YFP** (+89 ps). La colonne `concentration_ch*` est la conversion fictive
de FLIMApp.

### Réponses aux stims (étape 4)

Chaque lecture est redressée par une base quadratique ajustée sur toute la session avec la
réponse (modèle M1). Une pente linéaire tirée des moments calmes était trop raide et fabriquait
un faux dépassement de +0,5 % 40 à 100 s après les stims. Δ = 10 s après la fin − 10 s avant le
début. Témoins : la même mesure autour de fausses stims (deux dans la base, une à la fin de
chacun des trois premiers intervalles).

| lecture | stims 1–4 | témoins (5) | erreur type de bruit de photons |
|---|---|---|---|
| rapport YFP/CFP | −1,20, −1,04, −0,82, −1,09 % | −0,29 à +0,27 % | 0,30 % |
| photons CFP | +0,95, +0,46, +0,37, +0,71 % | −0,29 à +0,60 % | 0,28 % |
| photons YFP | −0,25, −0,58, −0,45, −0,38 % | −0,02 à +0,30 % | 0,11 % |
| τ CFP mono | +13, −6, +16, +16 ps | −17 à +8 ps | 11 ps |
| tm CFP | +9, +9, +4, +1 ps | −10 à +10 ps | 7 ps |

Seul le rapport sépare chaque stim des témoins (z de −2,7 à −4,0, témoins entre −1 et +1).
Le partage entre CFP et YFP varie d'une stim à l'autre (stim 1 : surtout le CFP) au gré de
fluctuations communes aux deux canaux, de la taille de celles des témoins. Les durées de vie
vont dans le même sens (+10 ps en moyenne contre −6 ps pour les témoins), à ~2 σ.

**Modèles cinétiques sur le rapport** (`results/reponses.png`, `reponses_modeles.csv`), sur
toute la série jusqu'à 495 s, résidus au bruit de photons (χ²/ddl 0,99) :

| modèle | paramètres | AIC |
|---|---|---|
| M1 : pompe + relaxation | −0,14 %/s pendant la stim, τ = 30 ± 7 s | +12 |
| **M2 : avec retard** | deux constantes égales de 9,0 s, minimum −1,12 % 4,8 s après la fin | **0** |
| M3 : avec part permanente | α = 1,00 ± 0,35 : aucune part permanente | +14 |

Avec un gain par stim (M1, τ commun 36 ± 8 s) : −0,14, −0,14, −0,06, −0,16 %/s ; la stim 3
répond moins. Les photons YFP seuls laissent voir une part permanente (M3 : α = 0,43 ± 0,22),
peu sûre vu leur excès de bruit.

### Latence, τ d'activation et τ de retour (étape 6)

Modèle marche-retour avec latence sur le rapport redressé (`results/cinetique.png`,
`cinetique.csv`) : rien pendant δ, puis une descente A·(1 − e^(−t/τ_act)) pendant la durée de la
stim, puis un retour en e^(−t/τ_ret). Une stim à la fois, puis les 4 ensemble (δ, τ_act et τ_ret
communs, amplitude propre à chaque stim ; IC par bootstrap des stims).

| | amplitude au minimum | latence δ | τ d'activation | τ de retour |
|---|---|---|---|---|
| stim 1 | −1,41 % | 4,1 ± 3,0 s | non mesurable | 17,6 ± 7,0 s |
| stim 2 | −1,38 % | 4,8 ± 2,9 s | non mesurable | 16,0 ± 6,5 s |
| stim 3 | −1,29 % | 3,4 ± 2,6 s | non mesurable | 8,5 ± 4,5 s |
| stim 4 | −0,97 % | 6,9 ± 0,6 s | non mesurable | 50 ± 43 s (35 s de retour avant 495 s) |
| **les 4 ensemble** | **−1,19 %** [−1,41 ; −1,01] | **4,3 s** [3,2 ; 4,8] | 7,5 s [3 ; 300] : non mesurable | **16,4 ± 3,8 s** [10,3 ; 32,8] |

Le τ d'activation n'est pas mesurable : sur 10 s de stim et à ce rapport signal/bruit, une
descente exponentielle rapide après une latence et une rampe lente ne se distinguent pas. Ce
qui se lit, c'est la **latence (~4 s)** et la **pente de la descente (−0,14 à −0,22 %/s)**. Il
faudrait des stims plus longues (30 à 60 s, jusqu'au plateau) pour mesurer τ_act. Les écarts
types viennent du jacobien et les IC des τ par stim sont log-normaux.

### Ordre de grandeur en mM (étape 4) et [Cl⁻] absolue (étape 5)

La seule calibration en tissu est celle du 2026-11-09 (autre tranche, imagerie au polygone et
SPCM). Ses pentes donnent un ordre de grandeur, selon la [Cl⁻] basale supposée
(`results/reponses_mM.csv`) :

| basale supposée | Δ après une stim, rapport | Δ, τ CFP | Δ, tm CFP | précision par seconde, rapport | τ CFP | tm CFP |
|---|---|---|---|---|---|---|
| 0 mM | +0,20 mM | +0,52 mM | +0,45 mM | 0,12 mM | 1,3 mM | 1,1 mM |
| 4,4 mM | +0,27 mM | +0,79 mM | +0,63 mM | 0,17 mM | 2,0 mM | 1,6 mM |
| 10 mM | +0,39 mM | +1,21 mM | +0,92 mM | 0,25 mM | 3,0 mM | 2,3 mM |

**Une [Cl⁻] absolue ne se lit pas ici** (`results/absolu.png`, `absolu.csv`) :
- le rapport lit **−2,6 mM** au début (R = 6,30, au-dessus du X0 = 5,38 du 2026-11-09 ;
  concentration négative gardée), τ CFP **5,7 mM** au même moment ;
- au fil de la session, τ CFP monte de 5,7 à 14,7 mM, le rapport de −2,6 à −1,9 mM : dérives de
  +2,0 et +0,15 mM/100 s. Rapporté aux pentes du 2026-11-09, **τ CFP bouge 3 à 5 fois plus par %
  de rapport** qu'en calibration, pour la dérive comme pour les stims ;
- la forme du déclin CFP ne se décompose pas sur les motifs à 0 et 131 mM du 2026-11-09 :
  la part f du motif « 131 mM » vaut 0,09, 0,50 ou 0,53 selon le début de la fenêtre, avec
  un résidu systématique de ~3 %.

Deux explications possibles de l'écart entre rapport et τ : une YFP sans FRET (eNpHR-YFP ?) qui
dilue les variations du rapport, ou une échelle des déclins (IRF, composante prompte, trajet)
qui n'est pas celle du 2026-11-09. Seule une calibration dans ce mode (ionophores à deux
concentrations au moins, en Realtime) tranchera.

### Lecture hybride : niveau du temps de vie, précision du rapport (étape 7)

À la demande de l'utilisateur : le niveau absolu vient de τ CFP (isotherme du 2026-11-09, qui
donne des valeurs plausibles), les variations du rapport. Dans chaque fenêtre calme (avant la
1re stim, puis de 60 s après chaque stim jusqu'à la suivante), τ moyen donne c_τ et l'échelle
du rapport η = R / R_iso(c_τ) ; chaque passe est lue au rapport, R_iso⁻¹(R / η). Deux variantes :
**ancrage continu** (η interpolé entre les ancres : le niveau lent suit τ) et **ancrage
unique** (η de la première ancre, avant toute stim : le niveau lent suit le rapport).
Figures : `results/hybride.png` (les deux variantes et le diagnostic) et
`results/hybride_ancrage_unique.png` (le rapport seul, ancré sur τ pendant les 60 s avant la
1re stim puis fixé). `results/avant_apres_stim.png` reprend la figure d'octobre sur cette
série unique : 10 s avant contre 10 s après, valeurs brutes (dérive comprise) — [Cl⁻] ancrée
+0,42, +0,38, +0,30, +0,40 mM (z de +3,2 à +4,3 ; témoins +0,07 mM en moyenne, de −0,04 à
+0,17) ; ΔR/R₀ −1,41, −1,21, −0,94, −1,21 % ; Δτ CFP +18, −2, +20, +17 ps (z ≤ 1,9).

| | ancrage continu | ancrage unique | τ CFP seul |
|---|---|---|---|
| niveau, 0–60 s | 5,9 mM | 5,7 mM | 5,7 ± 0,6 mM (IC 95 %) |
| niveau, 420–450 s (avant la stim 4) | 14,4 mM | 6,9 mM | 14,2 ± 1,0 mM |
| réponse aux stims (10 s après − 10 s avant, dérive lente retirée) | +0,36, +0,40, +0,35, +0,52 mM | +0,35, +0,31, +0,24, +0,35 mM | non lisible |
| précision par passe (240 ms) | 0,54 mM | 0,41 mM | 13 mM |
| précision par seconde | 0,27 mM | 0,20 mM | 6 mM |

La lecture hybride a bien la précision du rapport (20 à 30 fois mieux que τ seul) et le
niveau de τ aux ancres. **Mais la montée lente de 5,7 à 14 mM vient entièrement de τ, et elle
n'a pas la signature du Cl⁻** (panneau en bas à gauche de la figure) : entre les ancres, τ CFP
monte de **18 à 22 ps par % de baisse du rapport**, quand un changement de Cl⁻ en donne
**3,2** (pentes des deux isothermes du 2026-11-09 à 5,7 mM). Une perte d'accepteur en donne
**12** (modèle simple : capteurs dont la YFP est blanchie, sans FRET ; Δτ/Δ ln R = −τ_D·E,
τ_D = 2,72 ns en solution, E = 0,45). Avec la YFP qui blanchit quatre fois plus vite que le
CFP (étape 3), la dérive lente ressemble à un **blanchiment de l'accepteur** par l'imagerie
continue, pas à une charge en Cl⁻. L'ancrage continu la convertit en ~8 mM de Cl⁻ qui
n'existent sans doute pas. L'ancrage unique ne garde que la dérive du rapport (+1,2 mM),
elle-même en partie due au blanchiment. Les réponses aux stims (+0,3 à +0,4 mM) ne dépendent
presque pas de la variante. La moyenne des stims tombe entre les deux directions
(Δτ +10 ± 11 ps pour −1,0 % de rapport) : elle ne tranche pas entre Cl⁻ et accepteur.

En pratique : le temps de vie ne fixe le niveau absolu qu'avant que le blanchiment s'installe
(ici la première minute : 5,7 ± 0,6 mM, si l'isotherme de τ du 2026-11-09 vaut pour ce
montage). Ensuite, ni τ ni le rapport ne séparent seuls le Cl⁻ du blanchiment de l'accepteur.
Les deux ensemble le pourraient, puisqu'ils y répondent dans des proportions différentes : il
faudrait mesurer la direction « blanchiment » (blanchir l'accepteur exprès, sous ionophores)
et décomposer les dérives sur les deux directions.

### Images reconstruites sans M1/M2 et ROIs (étapes 8 à 10)

**Reconstruction (étape 8)**, `results/images.png`, `qc/images_suivi.png`. Sur une trame entière
(540 lignes, 30 ms), les photons repliés modulo la période de ligne donnent la projection de
l'image sur x, très structurée (des somas brillants sur un fond sombre). Elle se décale au fil
du temps avec la phase du polygone. Chaque trame (19 391) est recalée sur un gabarit de cette
projection par corrélation croisée, et les décalages, déroulés d'une trame à l'autre, donnent la
phase des lignes. Elle dérive de 79 µs/s en moyenne par rapport à l'horloge de la carte (ligne
suivie : 17 995,93 Hz), avec des écarts de ±1,8 µs d'une trame à l'autre. Chaque photon reçoit
alors sa ligne (y) et sa position dans la ligne (x). Contrôles :
- **pas de ligne perdue** : les blocs de 10 s, recalés sur l'image de la session, ne bougent pas
  en y (−0,5 à +0,4 px) ;
- **pas de motif de facette** : décalage x d'une ligne à l'autre de 0,26 canal seulement ;
- **dérive en x** : le champ bouge en x de −19 à +26 px au fil de la session, avec des sauts vers
  290 s et 535 s (le second avec le changement optique de la fin). Tranche ou erreur lente du
  suivi, elle est corrigée bloc par bloc.

L'image a le flou d'un suivi trame par trame (les somas restent nets, les dendrites un peu
étirées en x), assez pour des ROIs grossiers.

**ROIs (étape 9)** : Cellpose-SAM (diamètre 50 px) sur √photons de la session, bandes du retour
du galvanomètre écartées : **12 somas**, 47,5 % des photons (`results/rois.png`). Les photons
de chaque passe sont rangés par ROI avec la correction de dérive du bloc. Contrôle : la somme
des ROIs et du reste redonne exactement les photons de chaque passe (étape 1).

**Par ROI (étape 10)** : 100 à 400 photons CFP par passe et par soma (contre 6 160 pour le
champ). Ancre sur le τ CFP de chaque soma (60 s avant la 1re stim) : de 5,6 à 29 mM selon les
somas (IC 95 % de ±3 à ±11 mM). Le **ROI 2** (τ 1,70 ns, trop près du plateau de l'isotherme à
1,77 ns : 76 ± 29 mM) n'est pas lisible en mM et reste seulement dans les lectures du rapport.
Les somas ont un τ CFP plus long que le champ entier (1,50 à 1,70 ns contre 1,50) : le niveau
de 5,7 mM du champ venait surtout du neuropile et du fond.

Chaque fenêtre (10 s avant, 10 s après, 2 s ou 10 s dans les figures) est lue sur ses photons
sommés, puis inversée par l'isotherme. Avec 100 à 400 photons par passe, la moyenne des
lectures passe par passe serait biaisée vers le haut, puisque l'inversion est convexe. Sur la
fenêtre de l'ancre elle-même, ce biais va de +0,06 à +0,47 mM selon le soma. Une première
version l'avait ; les écarts aux stims n'en changent presque pas (+0,83 → +0,81 mM réunies).

| stim | ΔR/R₀ (médiane, 12 ROIs) | Δ[Cl⁻] (médiane, 11 ROIs) | Δ[Cl⁻] dérive retirée (moyenne) | Δτ CFP (médiane) |
|---|---|---|---|---|
| 1 | −1,86 % (3/12 ↑, p = 0,03) | +0,61 mM (8/11 ↑, p = 0,03) | +0,59 mM (7/11 ↑) | +11 ps (p = 0,5) |
| 2 | −1,08 % (4/12 ↑, p = 0,05) | +0,78 mM (8/11 ↑, p = 0,02) | +0,71 mM (8/11 ↑) | +6 ps (p = 0,3) |
| 3 | −0,75 % (3/12 ↑, p = 0,15) | +0,53 mM (9/11 ↑, p = 0,03) | +0,66 mM (9/11 ↑) | +5 ps (p = 0,5) |
| 4 | −1,48 % (2/12 ↑, p = 0,002) | +0,87 mM (9/11 ↑, p = 0,005) | +0,83 mM (9/11 ↑) | +23 ps (p = 0,7) |

Les stims réunies donnent **+0,81 mM [+0,50 ; +1,12]**, 34 paires sur 44 en hausse,
p = 6·10⁻⁷ (Wilcoxon) (`results/avant_apres_rois.png` ; dans `avant_apres_rois_reunies.png`, le 0
est la moyenne des ROIs dans les 10 s avant chaque stim, et chaque ROI garde son écart à cette
moyenne). Dérive lente retirée (ln R redressé par la base de M1 de chaque soma, comme l'étape
7) : **+0,70 mM [+0,41 ; +0,99]**, 33/44 ↑, p = 4·10⁻⁶. C'est le double du champ (+0,37 mM
brut, +0,31 mM dérive retirée) : les somas sont ancrés plus haut sur l'isotherme, où chaque %
de rapport vaut plus de mM.

**[Cl⁻] des somas dans le temps** (`results/hybride_ancrage_unique_rois.png`, la figure
`hybride_ancrage_unique.png` refaite par soma). Chaque soma est ancré une fois sur son τ CFP
(60 s avant la 1re stim), puis son échelle reste fixe. La figure trace des fenêtres de 2 s et
des fenêtres glissantes de 10 s.
- **Moyenne des 11 somas lisibles** : ancre à 13,7 ± 2,0 mM (IC 95 %, photons seulement ; sans
  l'erreur de l'isotherme ni celle de l'IRF, communes à tous). Bruit de 0,66 mM par seconde.
- **Bruit par soma** : de 0,9 à 3,3 mM par seconde.
- **Dérive** : la même montée lente que le champ, d'environ +2 mM entre 60 et 460 s, puis une
  baisse après 495 s (le changement optique de la fin). Elle vient sans doute en partie du
  blanchiment de la YFP (voir l'étape 7).
- **Réponses visibles** : une bosse à chaque stim sur la moyenne. Soma par soma, les plus nettes
  sont la marche de +2,8 mM du ROI 8 à la 1re stim et quelques réponses de 1 à 3 mM (ROIs 1, 4,
  7, 9). Les ROIs 11 et 12 ne bougent pas ou baissent un peu, sauf le ROI 11 à la 4e stim
  (+0,8 mM).

**Cinétique par ROI** (`results/cinetique_rois.png`) : τ
d'activation jamais mesurable, τ de retour mesurable dans 6 somas seulement (8 à 89 s, IC très
larges), latence de 1 à 15 s : à ces photons, les constantes de temps ne se lisent que sur le
champ entier (étape 6). Réserves : ROIs grossiers sur une image reconstruite, mélange possible
entre somas voisins (4/7/8/9 se touchent), ancres par soma sur l'isotherme d'une autre tranche.

### Réserves

1. **Champ entier et ROIs reconstruits** : les étapes 2 à 7 mélangent somas, neuropile et
   fond. Les étapes 8 à 10 séparent 12 somas sur une image reconstruite sans M1/M2, aux bords
   flous, avec peu de photons par soma.
2. **Cl⁻ ou artefact ?** Moins de FRET est ce qu'une entrée de Cl⁻ par eNpHR produit, mais le
   1064 nm excite la YFP de ce champ (+9 %). Un blanchiment de l'accepteur, un effet thermique
   ou de pH donneraient la même signature. Le retour complet en ~1 min plaide pour un
   phénomène réversible (Cl⁻ entré puis extrudé), sans exclure un blanchiment compensé par la
   diffusion. Le témoin décisif : le même protocole sur une tranche SCLM seul, ou sous
   ionophores, où la [Cl⁻] est imposée.
3. **Une seule session, une tranche** ; l'expression de eNpHR dans cette tranche n'est pas
   documentée dans les données.
4. Les τ et les mM s'appuient sur l'IRF du 2026-11-09, décalée, et sur la calibration d'une
   autre tranche dans un autre mode.

### Pour les prochaines sessions

- **Par cellule** : lancer le Realtime **en images** (`[clamp] images = true` dans
  `config/spc.toml`). FLIMApp y enregistre M1/M2 et le pipeline `../sclm_flim/` découpe les
  somas (Cellpose), avec le PWM tel quel. Autre voie : le mode ROI (spirales des galvos de
  FLIMApp), polygone arrêté. L'étape 6 sépare déjà les ROIs par code de routage.
- **IRF** : en mesurer une dans la configuration du Realtime ; celle du 2026-11-09 est en retard
  de 0,13 à 0,19 ns.
- **Canal suivi par FLIMApp** : ch1 = IN1 = YFP. Pour suivre τ CFP en direct, regarder ch2, ou
  inverser l'ordre de `[verification] series`.
- **τ d'activation** : des stims de 30 à 60 s (ou plusieurs puissances) pour atteindre le
  plateau.
- **Témoins** : même protocole sans eNpHR, ou sous ionophores.

## Sorties

| Fichier | Contenu |
|---|---|
| `mesures/passes.npz` | déclins par passe (acquisition, phases opto, pauses) et par canal, profils du PWM, stims, code de chaque passe (non versionné) |
| `qc/validation_flimapp.csv`, `qc/decodage.json` | contrôle contre frames.csv, compteurs du décodage |
| `qc/commandes_creneaux.csv` | commande relue et galvos, créneau par créneau (non versionné) |
| `results/stims.csv` | les 4 stims (passes, temps sur l'horloge de la carte) |
| `results/fuite.png`, `fuite.csv`, `fuite_profil.csv` | étape 2 : photons ajoutés par le 1064 nm, profil dans la période du PWM |
| `results/series.png`, `declins.png`, `bruit.png` | étape 3 : séries, déclins de la session, bruit et FLIMApp en direct |
| `results/series_passes.csv` | lectures passe par passe (non versionné) |
| `results/declins_session.csv`, `derive.csv`, `bruit.csv` | étape 3 : ajustements de la session, dérives, bruit |
| `results/reponses.png`, `reponses_temoins.png`, `reponses_*.csv` | étape 4 : réponses, témoins, modèles, gains par stim, mM |
| `results/absolu.png`, `absolu.csv`, `absolu_temps.csv` | étape 5 : lectures absolues avec la calibration du 2026-11-09 |
| `results/cinetique.png`, `cinetique.csv` | étape 6 : latence, τ d'activation et de retour, par stim et ensemble |
| `results/hybride.png`, `hybride_ancres.csv`, `hybride_reponses.csv` | étape 7 : lecture hybride (ancrages continu et unique), diagnostic de la dérive, réponses en mM |
| `results/hybride_ancrage_unique.png` | étape 7 : le rapport seul, échelle ancrée sur τ avant la 1re stim puis fixée |
| `results/avant_apres_stim.png`, `avant_apres_stim.csv` | étape 7 : avant / après chaque stim (comme en octobre) : [Cl⁻] du rapport ancré, ΔR/R₀, Δτ CFP ; témoins dans le CSV seulement |
| `results/avant_apres_reunies.png` | étape 7 : Δ[Cl⁻] avant / après, les 4 stims réunies (+0,37 mM [0,29 ; 0,46], 4/4 ↑, p = 0,0007) |
| `mesures/images.npz` | étape 8 : images reconstruites (par canal, tous / acquisition), suivi de phase, recalage des blocs (non versionné) |
| `results/images.png`, `qc/images_suivi.png` | étape 8 : images reconstruites ; suivi de la phase du polygone et recalage |
| `masks/cellpose_masks.npy`, `cellpose_params.json`, `results/rois.png` | étape 9 : ROIs Cellpose |
| `mesures/rois_passes.npz` | étape 10 : déclins par passe, ROI et canal (non versionné) |
| `results/rois_cl.csv`, `rois_avant_apres.csv`, `rois_cinetique.csv` | étape 10 : ancres, avant / après, cinétique par ROI |
| `results/rois_series.png`, `avant_apres_rois.png`, `avant_apres_rois_reunies.png`, `cinetique_rois.png` | étape 10 : figures par ROI |
| `results/hybride_ancrage_unique_rois.png` | étape 10 : [Cl⁻] de chaque soma dans le temps, ancrage unique sur son τ (moyenne des somas en haut) |
| `results/hybride.csv` | étape 7 : les lectures passe par passe (non versionné) |
