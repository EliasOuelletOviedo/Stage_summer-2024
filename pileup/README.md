# Empilement (pile-up) — SPC-150N et SPC-QC-104

Données du 2026-10-06 (tag `Data_pileup`, `2026-10-06.zip`) : une solution de MQAE imagée
à 750 nm (objectif 16X, 1024 × 1024, 50 ns par pixel, laser à 80,07 MHz), à plusieurs
puissances, par chacune des deux cartes B&H. `low` / `high` = faible / forte concentration
de Cl⁻ : le Cl⁻ éteint la MQAE (τ ≈ 5 ns en faible Cl⁻, ≈ 1,06 ns en fort Cl⁻).

**Question** : de combien la durée de vie mesurée est-elle faussée quand le taux de
photons augmente, avec chaque carte ?

## En bref

| | SPC-QC-104 | SPC-150N |
|---|---|---|
| temps mort (mesuré) | **7,5 ns** par entrée (< période de 12,49 ns) | **≈ 115 ns** (≈ 9,2 périodes) |
| taux enregistré maximal ici | 14 MHz, sans perte visible | 5,5 MHz pour 12 MHz au détecteur (58 % perdus) |
| biais de τ, faible Cl⁻ (τ ≈ 5 ns) | **−7,9 % par 0,1 photon/impulsion** | −2,5 ± 5,5 % par 0,1 (intra-image) ; ≈ 0 entre acquisitions |
| biais de τ, fort Cl⁻ (τ ≈ 1,06 ns) | **−3,7 % par 0,1 photon/impulsion** | −5,2 ± 4,2 % par 0,1 (intra-image) ; ≈ −1 % entre acquisitions |
| biais au plus fort taux mesuré | −16 % (faible Cl⁻, 0,21 ph/imp.) ; −7 % (fort Cl⁻, 0,20) | ≈ −1 à −3 % (0,07 ph/imp.) |
| taux pour un biais < 1 % (simulation) | 1,0 MHz (τ ≈ 5 ns) ; 2,5 MHz (τ ≈ 1 ns) | ≈ 3,7 MHz enregistrés (les deux solutions) |
| taux pour un biais < 2 % (simulation) | 2,0 MHz ; 4,9 MHz | ≈ 4,6–5,3 MHz enregistrés |

(1 photon enregistré par impulsion = 80 MHz ; 0,1 = 8 MHz.)

1. **La QC-104 empile comme le prévoit la théorie.** Le premier photon d'une impulsion
   est favorisé : le déclin paraît plus rapide et τ baisse linéairement avec le taux
   local, d'autant plus que τ est long. Une simulation sans paramètre libre (photons
   poissoniens par impulsion + temps mort de 7,5 ns mesuré sur les intervalles entre
   photons) reproduit la pente des τ mesurés (à l'écart propre à chaque acquisition près,
   ~1–3 %) et la déformation des déclins canal par canal (figures 3 et 4). Le temps mort, plus court que la période, ne protège pas :
   deux photons d'une même impulsion sont presque toujours à moins de 7,5 ns l'un de
   l'autre. Pour un biais < 1 % sur la MQAE faible Cl⁻, il faut rester sous ~1 MHz par
   détecteur (0,013 photon/impulsion), la règle classique des « 1 % de la fréquence du laser ».
2. **La 150N biaise beaucoup moins, mais pour une raison fragile.** Son temps mort de
   ~115 ns plafonne le taux à ~5,5 MHz (figure 1). La simulation montre que deux effets
   s'y opposent : l'empilement (qui favorise les photons précoces) et la libération de la
   carte au milieu d'une période (qui favorise les photons tardifs). Le résultat net
   **oscille avec D modulo la période** (figure 6). Avec D ≈ 114,5 ns, la carte se libère
   près de la phase du photon qui l'a bloquée et les deux effets se compensent presque.
   Ce D redonne aussi les temps morts effectifs mesurés (115/118 ns simulés, 116/120 ns
   mesurés). Un autre taux de répétition, ou quelques ns de plus ou de moins sur D,
   donneraient un biais de −5 à −10 % : c'est une coïncidence du montage, pas une propriété
   de la carte.
3. **À taux enregistré égal (≤ 5 MHz)**, la 150N donne un τ plus juste (−1 à −2 %) que
   la QC-104 (−4 % en faible Cl⁻, −2 % en fort Cl⁻ à 5 MHz). Mais la 150N jette
   au-delà de 2 MHz une grande part des photons. **La QC-104 collecte plus de photons et
   son biais est prévisible** : on peut le corriger (même simulation, ou correction de
   Coates sur l'histogramme), ou s'en tenir à ≤ 1–2 MHz par détecteur.
4. **Détecteurs coupés** (protection contre la surcharge) : QC-104 faible Cl⁻ P25 (canal
   1 : 233 000 photons, déclin plat) et QC-104 fort Cl⁻ P70 (les deux canaux, idem pour le
   `70.spc` à la racine) ; la 1re série de la 150N n'a pas de module 1 à P25. Ces
   acquisitions sont écartées.

## Méthode

**Données brutes.** Tout part des flux photon par photon `.spc` (les `.sdt` d'SPCM ne
contiennent que la première seconde, et rien sur le temps mort). Formats : FIFO_150 pour
la 150N (un module par détecteur, `_m1` / `_m2`, macrotemps de 25 ns, ADC inversé) ;
format établi par FLIMApp (`format_fifo_qc104.jl`) pour la QC-104 (deux entrées dans un
flux, macrotemps de 2,048 ns, TDC de 4 ps). Le canal de 4 ps est tiré uniformément dans
sa largeur avant d'être rangé en 256 canaux sur la période. Sinon, on obtient un peigne
de ±4 %. Plusieurs `.spc` de la QC-104 s'arrêtent à 64 Mio : sans effet, le taux est
mesuré pendant le balayage.

**Taux local.** La solution est uniforme, mais l'image ne l'est pas (vignettage, bande
sombre en bas des images fort Cl⁻ de la QC-104). On divise la carte d'intensité, lissée
(σ = 24 px), par le nombre de passages du balayage sur chaque ligne. Cela donne le taux
local en photons enregistrés par impulsion, μ_enr. Les photons de chaque image sont
rangés en 8 tranches de μ_enr, de même nombre de photons. **Toutes les tranches d'une
même image ont été acquises en même temps, sur le même échantillon** : seul le taux
change. Contrôle : à μ_enr égal, des acquisitions de puissances différentes donnent le
même τ (QC-104 : P15 et P20 à μ = 0,107, 4,77 et 4,79 ns).

**Durée de vie.** Ajustement au **maximum de vraisemblance (Poisson)** d'une
**exponentielle incomplète** h(t) = exp(−t/τ) / (1 − exp(−P/τ)), qui somme les déclins
de toutes les impulsions précédentes. Ce modèle est **reconvolué par l'IRF mesurée** de
la carte et du canal (`irf_16x_750nm_ch1/ch2.sdt` pour la QC-104, `irf_16x_750nm_150N.sdt`
pour la 150N ; fond retiré, ramenée sur la grille des déclins), avec un fond uniforme et
un décalage de l'IRF (`ajustement.py`). Fenêtre : là où la carte a des photons (150N :
limites du TAC ; QC-104 : de 0,82 ns à P). Le décalage est ajusté sur le déclin global
de chaque acquisition, puis fixé pour ses tranches. Biais = τ / τ0 − 1, où τ0 est
l'ordonnée à l'origine de la droite τ(μ_enr) de la carte et de la solution.

**Temps morts.** QC-104 : le macrotemps dit à quelle impulsion appartient chaque photon
et le microtemps donne l'intervalle exact entre deux photons successifs d'une même
entrée. Aucun intervalle sous 7 ns, demi-montée à 7,5 ns (figure 2). 150N : le
macrotemps de 25 ns est trop grossier. Le temps mort effectif vient des taux de StopInfo,
1/r_TAC − 1/r_CFD, constant de 0,7 à 12 MHz (non paralysable) : 116,3 ± 1,1 ns en
faible Cl⁻, 120,0 ± 0,9 ns en fort Cl⁻. La différence vient de la forme du déclin, que la
simulation reproduit.

**Simulation** (`simulation.py`). Chaque impulsion envoie au détecteur un nombre de
photons de Poisson. Leur temps dans la période suit le déclin mesuré au plus bas taux.
Un photon n'est enregistré que si la carte n'est pas morte. Deux variantes pour la 150N :
temps mort compté depuis le photon (« photon »), ou conversion lancée par la synchro qui
suit le photon (« synchro », TAC en start-stop inversé). Le D de la 150N est choisi pour
redonner à la fois les temps morts effectifs et le biais intra-image mesurés.
La variante « photon » avec D = 114,5 ns est la meilleure (χ² = 7,3 pour 4 mesures ;
« synchro » : D = 106,5 ns, χ² = 8,3). Les deux variantes encadrent la prédiction de la
figure 4.

## Limites

- **150N : biais mal déterminé.** Le biais intra-image (−2,5 % et −5,2 % par 0,1) et
  celui entre acquisitions (+0,6 % et −0,9 % par 0,1) ne s'accordent qu'à 1–2 σ. Le
  domaine de μ_enr est court (≤ 0,076). Le front de montée de l'IRF ne change pas avec la
  puissance (10–90 % : 137–141 ps en fort Cl⁻, 171–179 ps en faible Cl⁻), donc ce n'est
  pas la réponse instrumentale. Pour la QC-104, les
  deux pentes concordent (−7,9/−8,3 % et −3,7/−3,9 %).
- **Le D de la 150N est déduit, pas mesuré.** La compensation dépend de D au ns près
  (figure 6). Pour l'établir, il faudrait des acquisitions à plus haut taux, ou mesurer
  directement la phase de libération de la carte.
- **τ absolus.** La mono-exponentielle n'est pas parfaite à des millions de photons
  (déviance réduite 10 à 500 à fort taux), et les deux cartes ne donnent pas le même τ0 en
  faible Cl⁻ (QC-104 : 5,23 ns, 150N : 4,81 ns ; IRF et électronique différentes). Les
  biais se comparent donc à l'intérieur d'une carte. D'une acquisition à l'autre, τ varie
  aussi de 1–3 % hors empilement. La 1re série de la 150N (`SPC-150N/*.sdt`, sans `.spc`,
  6 min plus tôt) donne 5,07–5,14 ns, sans tendance avec la puissance.
- **Canal 2** (~10 fois moins de photons, ≤ 1,5 MHz) : τ stable à ±1–2 % pour les deux
  cartes. C'est le témoin sans empilement.

## Lancer

Données : dézipper `2026-10-06.zip` (release `Data_pileup`) dans `pileup/data/`, ou donner
le dossier qui contient `SPC-150N/` et `SPC-QC-104/` en argument (ou `PILEUP_DATA`).

```bash
python3 step1_photons.py      # ~3 min : flux .spc -> déclins par tranche de taux, intervalles
python3 step2_durees.py       # ~1 min : τ (exponentielle incomplète ⊛ IRF, MLE), temps morts
python3 step3_simulation.py   # ~15 min : pentes mesurées, calage du D de la 150N, courbes simulées
python3 step4_figures.py      # ~30 s : figures
```

Dépendances : numpy, scipy, matplotlib, sdtfile.

## Fichiers

| | |
|---|---|
| `config.py` | séries, puissances, fenêtres, constantes du laser et du balayage |
| `lecture.py` | décodage des `.spc` (150N, QC-104), position des photons dans l'image, IRF |
| `ajustement.py` | exponentielle incomplète reconvoluée par l'IRF, maximum de vraisemblance de Poisson |
| `simulation.py` | Monte-Carlo : photons poissoniens par impulsion, temps mort |
| `step1` … `step4` | les étapes |
| `resultats/acquisitions.csv` | par acquisition et canal : photons, durée, taux (flux, StopInfo, CFD), μ_enr |
| `resultats/durees_de_vie.csv` | τ, erreur, fond, déviance : déclin global et 8 tranches de taux |
| `resultats/temps_morts.csv` | temps mort de la QC-104 (intervalles), temps morts effectifs de la 150N |
| `resultats/pentes.csv` | τ0 et biais par photon/impulsion, toutes tranches et intra-acquisition |
| `resultats/balayage_150N.csv` | 150N : biais et temps mort effectif simulés selon D |
| `resultats/simulation.csv` | biais simulé selon μ_det (photons au détecteur par impulsion) |
| `resultats/figures/` | 1 taux · 2 temps morts · 3 déclins · 4 biais par tranche · 5 biais au détecteur · 6 balayage de D · 7 τ par puissance |
