**Figure A.2 — Erreur sur la concentration de chlorure causée par l'autofluorescence, et prix de sa
correction.**
Le capteur simulé est SuperClomeleon en solution, tel que mesuré à 850 nm : durées de vie, réponses
instrumentales, lumière diffusée et fond du détecteur proviennent d'une titration de 0 à 131 mM de Cl⁻.
Le donneur (canal bleu) a deux populations, à 0,655 ns (FRET) et 2,57 ns (sans FRET) ; le canal jaune a
deux composantes, à 1,02 et 3,22 ns. On y ajoute une autofluorescence d'intensité indépendante du Cl⁻,
égale à 5 % des photons du canal à 20 mM, avec des valeurs typiques sous excitation à deux photons
(Kolenc et Quinn, 2019) :
- dans le canal bleu, le NAD(P)H : 75 % de l'amplitude à 0,4 ns (forme libre), 25 % à 2,0 ns (forme
  liée) ;
- dans le canal jaune, les flavines (FAD) : 75 % à 0,3 ns (forme liée, éteinte), 25 % à 2,6 ns (forme
  libre).

Trois lectures de la concentration sont comparées :
- le rapport des photons jaunes sur bleus ;
- la durée de vie du donneur, par un ajustement mono-exponentiel du seul canal bleu ;
- un ajustement simultané des deux canaux, où les durées de vie du capteur et de l'autofluorescence
  sont connues. Seules l'état du capteur, l'échelle et l'amplitude de l'autofluorescence de chaque
  canal sont libres.

Chaque lecture est calibrée sur le capteur sans autofluorescence, puis appliquée au capteur contaminé,
à 10 mM. La réponse instrumentale, l'origine des temps, la lumière diffusée et le fond du détecteur
sont connus pour toutes.

**(a) Erreur sur la [Cl⁻] lue** (mM ; positive : surestimation).
- **Rapport :** il confond l'autofluorescence avec le signal. Le NAD(P)H gonfle le canal bleu et fait
  lire +2,9 mM ; les flavines gonflent le jaune et font lire −2,3 mM.
- **Durée de vie du donneur :** insensible aux flavines, qui n'émettent pas dans le bleu ici, mais le
  NAD(P)H, plus bref que le donneur, la raccourcit (−2,2 mM).
- **Ajustement à deux canaux avec un terme d'autofluorescence :** erreur nulle quand les durées de vie
  de l'autofluorescence sont celles du modèle (troisième ligne).
- **Formes mal connues** (dernière ligne, barres en contour) : l'autofluorescence réelle diffère
  alors de celle supposée (NAD(P)H : 60 % à 0,5 ns + 40 % à 2,5 ns ; flavines : 60 % à 0,4 ns + 40 % à
  2,9 ns). L'ajustement retire encore 70 % de l'erreur due aux flavines (−0,7 mM), mais seulement 30 %
  de celle due au NAD(P)H (+2,0 mM), dont les durées de vie recouvrent celles du donneur.
- **Contrôle, non montré :** sans le terme d'autofluorescence, le même ajustement conserve 85 à 95 % de
  l'erreur du rapport. La correction vient du terme explicite, pas de l'usage des durées de vie en
  soi.

**(b) Écart-type de la [Cl⁻] lue** pour 10 000 photons dans le canal bleu, à 10 mM. C'est l'écart-type
robuste (intervalle interquartile / 1,349) de 600 simulations avec bruit de Poisson, NAD(P)H et
flavines présents à 5 %. Les deux amplitudes supplémentaires multiplient le bruit par 2,7 par rapport
au rapport (1,5 contre 0,6 mM). La durée de vie du donneur seul est 6 fois plus bruitée (3,6 mM).

**Portée.**
- Les erreurs sont proportionnelles au niveau d'autofluorescence et croissent avec la [Cl⁻].
- 5 % des photons bleus correspondent à 6 % des photons dans ±0,25 ns autour du pic ; dans nos
  tranches de cerveau, une composante additive brève y atteint 3 à 24 %.
- À 850 nm, le NAD(P)H est peu excité (excitation efficace entre 710 et 780 nm ; Huang et al., 2002),
  si bien que les flavines dominent.
- Le modèle n'inclut pas de diaphonie spectrale entre canaux.

Références : Kolenc O. I. et Quinn K. P., *Antioxid. Redox Signal.* 30, 875–889 (2019) ; Huang S.,
Heikal A. A. et Webb W. W., *Biophys. J.* 82, 2811–2825 (2002).
