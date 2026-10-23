"""Configuration — empilement (pile-up) de la SPC-150N et de la SPC-QC-104.

Données du 2026-10-06 (tag `Data_pileup`, 2026-10-06.zip) : une solution de MQAE
imagée à 750 nm (objectif 16X), à plusieurs puissances, par les deux cartes B&H tour à
tour. `low` / `high` = faible / forte concentration de Cl⁻ : le Cl⁻ éteint la MQAE, la
durée de vie est plus courte et il y a moins de photons à puissance égale.

    <données>/SPC-QC-104/Images_low|Images_high/<P>.spc|.sdt|.set
    <données>/SPC-150N/Image_low|Image_high/<P>_m1.spc, <P>_m2.spc, <P>.sdt
    <données>/SPC-150N/<P>.sdt                   1re série, sans .spc (voir SERIES)
    <données>/irf_16x_750nm_ch1|ch2.sdt          IRF de la QC-104, une entrée par fichier
    <données>/irf_16x_750nm_150N.sdt             IRF de la SPC-150N, les deux modules

Le dossier des données (celui qui contient SPC-150N/ et SPC-QC-104/) se donne en
premier argument de chaque étape, ou par la variable PILEUP_DATA ; à défaut,
`pileup/data/2026-10-06`.
"""
import os
import sys
from pathlib import Path

import numpy as np

ICI = Path(__file__).resolve().parent
RESULTATS = ICI / "resultats"
MESURES = RESULTATS / "mesures"          # intermédiaires .npz (ignorés par git)
FIGURES = RESULTATS / "figures"

# --- laser et fenêtre -----------------------------------------------------------------
F_LASER = 8.007312e7        # Hz, taux de synchro lu par les deux cartes (StopInfo)
PERIODE_NS = 1e9 / F_LASER  # 12,489 ns
N_CANAUX = 256              # déclins : 256 canaux sur une période
# Fenêtre des ajustements, en canaux [début, fin) : là où chaque carte a des photons.
# 150N : limites du TAC, 5,1 % et 94,9 % de 12,5 ns (13 canaux de chaque côté) ;
# QC-104 : limite basse à 5 % de la plage du TDC (0,82 ns, canal 17), rien au-delà de P.
FENETRE = {"SPC-150N": (12, 243), "SPC-QC-104": (17, 256)}

# --- balayage ----------------------------------------------------------------------------
PIXEL_NS = 50.0             # temps par pixel (SP_PIX_T des .set)
PIXELS_LIGNE = 1024
LISSAGE_PX = 24             # écart-type (px) de la carte d'intensité qui donne le taux local
N_TRANCHES = 8              # pixels regroupés par taux local : 8 tranches de même nombre de photons

# --- cartes ------------------------------------------------------------------------------
# 150N : flux FIFO_150, un module (une carte) par détecteur, fichiers _m1 / _m2 ; macrotemps
# de 25 ns (en-tête du .spc), ADC 12 bits, temps = 4095 − ADC.
# QC-104 : une carte, deux entrées dans le même flux (bits 28-29) ; macrotemps 2,048131 ns ;
# microtemps du TDC en canaux de 4 ps (4096 sur 16,384 ns), relatif à la synchro.
# Format établi par FLIMApp (format_fifo_qc104.jl, branche FLIMApp).
CARTES = ("SPC-150N", "SPC-QC-104")
NOM_COURT = {"SPC-150N": "150N", "SPC-QC-104": "QC-104"}
TIC_QC_S = 2.048131e-9
CANAL_QC_NS = 16.384 / 4096

# Canal 1 = détecteur principal (~10 fois plus de photons), canal 2 = second détecteur.
CANAUX = (1, 2)

# --- séries ------------------------------------------------------------------------------
# (carte, solution, sous-dossier, puissances) ; les puissances sont celles des noms de
# fichiers (consigne du laser).
SOLUTIONS = {"low": "faible Cl⁻", "high": "fort Cl⁻"}
SERIES = [
    ("SPC-QC-104", "low", "SPC-QC-104/Images_low", (5, 10, 15, 20, 25)),
    ("SPC-QC-104", "high", "SPC-QC-104/Images_high", (10, 20, 30, 40, 50, 60, 70)),
    ("SPC-150N", "low", "SPC-150N/Image_low", (5, 10, 15, 20)),
    ("SPC-150N", "high", "SPC-150N/Image_high", (10, 20, 30, 40)),
]
# 1re série de la 150N (17 h 29 – 17 h 31, avant Image_high et Image_low), au dossier racine
# SPC-150N/ : .sdt seulement. Mêmes taux et même forme de déclin que Image_low :
# c'est la solution faible Cl⁻. Elle ne sert que de répétition, sur ses déclins globaux.
SERIE_150N_SDT = ("SPC-150N", "low", "SPC-150N", (5, 10, 15, 20, 25))

IRF = {"SPC-QC-104": {1: ("irf_16x_750nm_ch1.sdt", 0, 0), 2: ("irf_16x_750nm_ch2.sdt", 0, 1)},
       "SPC-150N": {1: ("irf_16x_750nm_150N.sdt", 0, 0), 2: ("irf_16x_750nm_150N.sdt", 1, 0)}}

# Un déclin n'est gardé que si ses photons montrent de la fluorescence : à pleine
# extinction du détecteur (QC-104 low/25, ch1 : 207 photons) ou sans échantillon
# (les 70.spc), le déclin est plat : bruit de fond seul.
MIN_PHOTONS = 20_000


def dossier_donnees(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    for a in argv:
        if not a.startswith("-") and (Path(a) / "SPC-150N").is_dir():
            return Path(a).resolve()
    env = os.environ.get("PILEUP_DATA")
    if env:
        return Path(env).resolve()
    d = ICI / "data" / "2026-10-06"
    if (d / "SPC-150N").is_dir():
        return d
    sys.exit("Données introuvables : donne le dossier qui contient SPC-150N/ et SPC-QC-104/ "
             "(ex. python3 step1_photons.py data/2026-10-06)")


def preparer():
    for p in (RESULTATS, MESURES, FIGURES):
        p.mkdir(parents=True, exist_ok=True)


def temps_canaux():
    """Centres des 256 canaux d'une période, en ns."""
    return (np.arange(N_CANAUX) + 0.5) * PERIODE_NS / N_CANAUX


def cle(carte, sol, p):
    return f"{NOM_COURT[carte]}_{sol}_{p}"
