# Pipeline SCLM + eNpHR — sessions FLIMApp (FLIM en temps réel, images)

Pour les jeux de données au format FLIMApp : une session « Realtime en images » (le
scanner du microscope balaie le champ, la SPC-QC-104 date chaque photon) et les images
fixes Bliq prises autour. Même démarche que `../SCLM+eNpHR_2026-09-21/` — ROIs Cellpose
triés à la main, ratio FRET en valeurs brutes, avant / après chaque stim — avec en plus
la **durée de vie du donneur**, que les données Bliq n'avaient pas.

Le code est ici, une seule fois ; chaque jeu de données garde ses données et ses
résultats dans son propre dossier.

## Format attendu

```
<jeu>/data/FLIM/        session FLIMApp, telle qu'écrite par l'appli
    run.toml            protocole, ROIs dessinées, réglages SPC
    frames.csv          ce que FLIMApp a mesuré en direct (une ligne par ROI et par bloc)
    visits.csv          commandes envoyées par la NI, slot par slot (300 ms)
    irf.csv / irf.toml  IRF par canal, 256 canaux sur 12,5 ns
    spc/<série>.spc     flux brut de chaque entrée (3T0089_IN1, _IN2) + _acquisition.ini
<jeu>/data/FOVs/        images fixes Bliq : FOV1_800nm-001, FOV1_950nm_post-001, …
```

`config.py` sait lire tout ça (`session()`, `irf()`, `fovs()`) ; le nom des images fixes
doit commencer par `<FOV>_<λ>nm`, `_post` pour celles d'après la session.

## Lancer

Depuis ce dossier, le jeu de données en argument :

```bash
python3 step1_decode.py ../<jeu>                          # ~2 min : photons, stim, dérive, recalage Bliq
~/.venvs/cellpose/bin/python step2_cellpose.py ../<jeu>   # ~30 s : somas (Cellpose requis)
python3 step3_select.py ../<jeu>                          # tri manuel (fenêtre) ; --tout = tout garder
python3 step4_declins.py ../<jeu>                         # ~2 min : déclins par ROI et par bloc
python3 step5_analyse.py ../<jeu>                         # ~30 s : durées de vie, réponses, figures
python3 step6_ratio_corrige.py ../<jeu>                   # quelques s : ratio image par image corrigé de sa dérive
python3 step7_retour.py ../<jeu>                          # ~10 s : retour à l'équilibre, temps de remontée
```

Les figures de présentation (fond noir) sont produites par un script rangé avec elles,
`<jeu>/présentation/presentation.py` : il prend son jeu de données dans le dossier parent
et le code ici. Pour un nouveau jeu, copier ce fichier dans son dossier `présentation/`.

Seule l'étape 2 a besoin de Cellpose. Après un tri manuel, il suffit de relancer
l'étape 5 : l'étape 4 mesure tous les ROIs une fois pour toutes.

## Ce que fait chaque étape

**1 — décodage.** `spc.py` relit les `.spc` (format FIFO_150, traduit de la QC-104 par
FLIMApp) et range chaque photon dans son pixel exactement comme FLIMApp
(`src/spc/rangement.jl`) : ligne au marqueur M1, trame au M2, 32 lignes ignorées en haut,
21 pixels à gauche. Le pixel vaut **25 tics = 51,2 ns** et non 50 : FLIMApp calcule en
dixièmes de ns entiers et arrondit le tic de 2,048 ns à 2,0. Un bloc = 50 images
(3,0 s), le bloc de FLIMApp.

Les fenêtres de stim viennent de la commande 1 de `visits.csv` (horloge de la NI),
recalées sur l'horloge des cartes par la fuite de la stim dans ch2, image par image
(60 ms). La dérive est mesurée bloc par bloc ; sous le demi-pixel, les masques restent
fixes. Les images Bliq sont ramenées au repère FLIM : échelle en x (pixel Bliq de 51,8
ns) puis translation, mesurées sur l'image qui ressemble le plus au FLIM et appliquées
à toutes (c'est une propriété du montage). Une image dont la ressemblance à cette
référence est sous `MIN_NCC_MEME_PLAN` (0,5) est marquée « autre plan focal ».

**2 — Cellpose** (`cpsam`) sur la somme des photons ch1 + ch2 de toute la session, dans
le repère exact des déclins. Diamètre 100 px : à 60 px (réglage du 2026-09-21), des
somas brillants étaient manqués.

**3 — tri manuel**, comme au 2026-09-21 (clic = garder/rejeter, `b` = fond, `q` =
enregistrer). Les ROIs FLIMApp ne se trient pas.

**4 — déclins.** Second passage sur les photons : un déclin de 256 canaux par ROI, par
bloc et par canal, pour les ROIs Cellpose, les ROIs FLIMApp (pixel pour pixel comme
`roi_pixel_mask`) et le neuropile (hors somas dilatés de 10 px, au-dessus du 25e
centile). Les photons de chaque ROI sont aussi comptés image par image (60 ms), pour
l'étape 6. **Contrôle intégré** : les photons des ROIs FLIMApp, bloc par bloc, doivent
être identiques à `frames.csv` (`qc/validation_flimapp.csv`).

**5 — analyse.** Ratio = photons ch2 / ch1 (YFP / CFP). Durées de vie : modèle de
FLIMApp (mono-exponentielle reconvoluée par l'IRF, vraisemblance de Poisson,
`flim_fit.py`), le décalage de l'IRF ajusté une fois par ROI sur toute la session puis
fixé ; τm = temps moyen − celui de l'IRF, sans ajustement. Réponse = 10 s après la fin
de la stim − 10 s avant son début. **Témoins** : la même mesure sur une fausse stim de
même durée au milieu de chaque intervalle sans stim — ce que la dérive et le bruit
font seuls. Comparaison stims / témoins appariée par ROI (Wilcoxon).

**6 — ratio corrigé de sa descente graduelle**, image par image, pour le champ entier
et pour les somas retenus (photons réunis). La pente est estimée hors des stims : la
base, puis de `EXCLURE_APRES_STIM_S` (30 s) après chaque stim jusqu'à la suivante, où le
ratio a fini de remonter. Modèle ln(ratio) = a_k + β·t, une pente commune et **un niveau
libre par intervalle** : les marches des stims restent hors de la pente. Une
exponentielle lisse passant par toutes les bases les aurait absorbées et aurait fabriqué
un faux retour après chaque stim. Ratio corrigé = ratio · exp(−β·t).

**7 — retour à l'équilibre.** Sur le ratio corrigé des somas (tranches de 1 s depuis la
fin de chaque stim, en % des 10 s d'avant), mono-exponentielle C + A·exp(−t/τ), équilibre
libre, sur une fenêtre commune à toutes les stims. Trois niveaux : somas réunis (un τ par
stim et un sur la moyenne des retours), chaque soma (moyenne de ses retours) et chaque
soma × stim. La distribution entre cellules est comparée à un **témoin « bruit seul »** :
chaque soma resimulé avec le τ commun et ses propres résidus décalés circulairement, puis
réajusté — ce qu'on observerait si toutes les cellules revenaient à la même vitesse.

## Canaux

ch1 = IN1 = **CFP** (donneur, τ ~1,2 ns, raccourci par le FRET) ; ch2 = IN2 = **YFP**
(accepteur, ~4 fois plus de photons, τ apparent ~2,4 ns). C'est aussi ce que suivait
FLIMApp pendant l'expérience (durée de vie de ch1, rapport ch2/ch1). Une baisse du
ratio ou une **hausse** de τ CFP = moins de FRET, ce qu'on attend d'une entrée de Cl⁻.

Deux pièges propres à ce montage :

- eNpHR est fusionnée à la YFP : **ch2 contient aussi la YFP de eNpHR**, qui ne fait
  pas de FRET. Le ratio en est dilué (surtout dans le neuropile, riche en membranes) ;
  τ CFP, lui, n'en dépend pas ;
- la colonne `concentration_ch*` de `frames.csv` est une conversion fictive de FLIMApp
  (`(9,5/τ − 1)/0,025`, Stern-Volmer de substitution) : à ignorer.

## Fichiers

| | |
|---|---|
| `config.py` | paramètres, lecture d'une session, masques des ROIs FLIMApp |
| `spc.py` | décodage des `.spc` et rangement des photons (vectorisé, un passage) |
| `flim_fit.py` | ajustement des durées de vie, temps moyen |
| `figures.py` | style commun des figures |
| `step1` … `step7` | les étapes |
| `<jeu>/présentation/presentation.py` | (dans chaque jeu) figures sur fond noir : images avec et sans ROIs, photons et ratio image par image, ratio corrigé, retour à l'équilibre, temps de remontée. Calculs importés d'ici ; couleurs, polices, tailles et tous les textes réglables en tête du script |
