"""Configuration — calibration de SCLM en tranche, [Cl⁻] imposée par ionophores.

Un jeu de données est un dossier qui contient :

    data/FLIM/   fichiers du logiciel SPCM (Becker & Hickl, SPC-QC-104 en FIFO_IMAGE) :
                 <condition>/fov<k>_c<n>.spc   un cycle de 10 s par fichier (c1…c6, ou c01…)
                 <condition>/fov<k>_c<n>.sdt/.set   images du dernier cycle, réglages
                 <transition>/<transition>_c<nn>.spc   time-lapse : un cycle toutes les ~41 s
    data/Bliq/   images Bliq VMS (Nirvana) prises en même temps, moyenne de 333 trames
                 (10 s) : <condition>/fov<k>-001/Bliq VMS/C1|C2|C3/*.tif

Les conditions fixes sont « basal » (ACSF, sans ionophore) et les [Cl⁻] imposées (« 0 mM »,
« 20 mM », « 131 mM ») ; une transition « a_to_b » est un time-lapse lancé au changement
de solution (t = 0 au premier cycle, d'après l'expérimentateur). Une transition peut être
coupée en plusieurs dossiers (« basal_to_0_1 », « _2 », « _3 ») : ils sont mis bout à bout.

Tout ce que produit le pipeline est écrit à côté de `data/`, dans le dossier du jeu.
Le dossier se donne en premier argument de chaque étape, ou par la variable
SCLM_CALIB_DATA ; à défaut, le dossier courant.
"""
import os
import re
import sys
import tomllib
from pathlib import Path

import numpy as np

# --- canaux ----------------------------------------------------------------------------
# Champ « entrée » du flux FIFO de la QC-104 (bits 28-29) : 0 = IN1, 1 = IN2 (FLIMApp,
# src/spc/SPCLite.jl). Le 2026-11-09, IN1 = YFP et IN2 = CFP, établi par les données :
# IN1 a 4 à 5 fois plus de photons et un déclin plus lent, IN1/IN2 baisse quand le Cl⁻
# monte (5,1 à 0 mM, 1,7 à 131 mM) et le déclin de IN2 s'allonge (le donneur perd son FRET).
# C'est l'inverse des sessions FLIMApp du 2026-10-08/09 (IN1 = CFP) : câblage changé.
# Les images Bliq suivent les mêmes détecteurs : C1 = IN1 (NCC 0,96), C2 = IN2 (0,89) ;
# C3 (gain 0,2) ne montre que les rayures du polygone, pas de cellules : ignoré.
ENTREE = {"yfp": 0, "cfp": 1}
NOM = {"yfp": "YFP (IN1, Bliq C1)", "cfp": "CFP (IN2, Bliq C2)"}
CANAUX = ("cfp", "yfp")                       # ordre des tableaux : 0 = CFP, 1 = YFP

# --- temps -----------------------------------------------------------------------------
# Macrotemps : 2,048131 ns (en-tête des .spc, en fs). Microtemps : canal du TDC de 4 ps,
# croissant avec le temps depuis la synchro. Laser : 80,2475 MHz (StopInfo des .sdt).
TIC_NS = 2.048131
CANAL_TDC_NS = 16.384 / 4096
F_LASER_HZ = 8.02475e7
PERIODE_NS = 1e9 / F_LASER_HZ                 # 12,461 ns
N_CANAUX = 256                                # déclins : 256 canaux sur une période
DT_NS = PERIODE_NS / N_CANAUX                 # 48,7 ps
# Le TDC n'enregistre rien sous 0,82 ns (canal 205 de 4 ps) : fenêtre des déclins
# [17, 256), de 0,83 ns à la fin de la période (comme pileup/ pour la QC-104).
FENETRE = (17, 256)

# --- géométrie ---------------------------------------------------------------------------
# Scanner à polygone : une ligne au marqueur M1 toutes les 55,56 µs (500 tours/s × 36
# facettes), une trame au M2 toutes les 30 ms (540 lignes). Pixel de 51,82 ns (« Dwell
# time » des images Bliq). L'image Bliq (512 × 1024) est le rectangle des lignes 14 à 525
# de la trame et des pixels 28 à 1051 de la ligne : décalages trouvés par corrélation des
# photons de IN1 avec Bliq C1 (NCC 0,96), vérifiés acquisition par acquisition (étape 1).
DWELL_NS = 51.82421227197347
LIGNES_IMAGE, PIXELS_LIGNE = 512, 1024
DECALAGE_LIGNES, DECALAGE_PIXELS = 14, 28
IMPULSIONS_PAR_PIXEL = DWELL_NS * 1e-9 * F_LASER_HZ   # 4,16 impulsions par passage
PIXEL_UM = 0.5087830142943353                 # métadonnée OME des images Bliq (16X)

# --- Cellpose ------------------------------------------------------------------------------
# Somas de ~15 µm à 0,509 µm/px : 30 px, comme Tranches/sclm (même objectif, même pixel).
DIAMETRE_PX = 30
MIN_TAILLE_PX = 120

# --- dérive ---------------------------------------------------------------------------------
# Décalage de chaque cycle sur la somme de l'acquisition (corrélation de phase) ; arrondi
# au pixel et appliqué aux masques à l'étape 4 (sous le demi-pixel : rien).
# Neuropile : hors des somas dilatés de 10 px, au-dessus du 25e centile (comme sclm_flim).
DILATATION_NEUROPILE = 10
CENTILE_NEUROPILE = 25

# --- IRF ------------------------------------------------------------------------------------
# IRF mesurée par SPCM pour le jeu (data/FLIM/irf*.sdt ; 2026-11-09 : irf_850nm_16X.sdt, prise le
# 2026-11-10 avec les mêmes décalages du TDC que les données, donc dans le même repère de
# temps). Utilisée telle quelle : seul le fond plat (médiane avant le pic) est retiré, la
# traîne est gardée. Le .sdt a deux courbes non vides, lignes 0 et 2 : attribuées dans l'ordre
# des entrées, ligne 0 = IN1 (YFP), ligne 2 = IN2 (CFP) — à confirmer.
IRF_LIGNE = {"yfp": 0, "cfp": 2}

# --- seuils ----------------------------------------------------------------------------------
MIN_PHOTONS_FIT = 2000        # sous ce nombre de photons CFP, pas d'ajustement de τ
MIN_DUREE_CYCLE_S = 3.0       # cycle interrompu plus court : écarté

CONDITIONS_FIXES = ("basal", "0 mM", "20 mM", "131 mM")
CL_MM = {"basal": np.nan, "0 mM": 0.0, "20 mM": 20.0, "131 mM": 131.0}
COULEUR_CONDITION = {"basal": "#52514e", "0 mM": "#2a78d6", "20 mM": "#1baf7a", "131 mM": "#eb6834"}


# ------------------------------------------------------------------ jeu de données
def racine(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    for a in argv:
        if not a.startswith("-") and (Path(a) / "data" / "FLIM").is_dir():
            return Path(a).resolve()
    env = os.environ.get("SCLM_CALIB_DATA")
    if env:
        return Path(env).resolve()
    if (Path.cwd() / "data" / "FLIM").is_dir():
        return Path.cwd()
    sys.exit("Jeu de données introuvable : donne le dossier qui contient data/FLIM "
             "(ex. python3 step1_decode.py ../SCLM_calibration_2026-11-09)")


def dirs(root):
    d = {k: root / k for k in ("projections", "masks", "selection", "mesures", "results", "qc")}
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d


def _num_cycle(p):
    return int(re.search(r"_c(\d+)\.spc$", p.name).group(1))


def _num_T(p):
    return int(re.search(r"-T(\d+)\.tif$", p.name).group(1))


def _cl_de(nom):
    """« 0 mM » -> 0.0 ; « basal » -> nan."""
    m = re.match(r"^(\d+(?:\.\d+)?) ?mM$", nom)
    return float(m.group(1)) if m else np.nan


def acquisitions(root):
    """Toutes les acquisitions du jeu : FOVs des conditions fixes, puis transitions.

    Chaque acquisition : {nom, type ('fov' | 'transition'), condition, cl_mM, avant, apres,
    cycles: [Path .spc triés], sdt: Path | None, bliq: {C1: [Path], C2: [Path]}}.
    Pour une transition coupée en plusieurs dossiers, les cycles sont mis bout à bout.
    """
    F, B = root / "data" / "FLIM", root / "data" / "Bliq"
    out = []
    for cond in sorted((p.name for p in F.iterdir() if p.is_dir() and "_to_" not in p.name),
                       key=lambda c: (c != "basal", _cl_de(c) if c != "basal" else -1)):
        fovs = sorted({re.match(r"(fov\d+)_c", p.name).group(1) for p in (F / cond).glob("fov*_c*.spc")},
                      key=lambda s: int(s[3:]))
        for fov in fovs:
            # Par heure d'écriture, pas par numéro : un FOV relancé écrase ses premiers cycles
            # (basal/fov2 du 2026-11-09 : c5-c6 du 1er passage, puis c1-c4 du second).
            cycles = sorted((F / cond).glob(f"{fov}_c*.spc"), key=lambda p: p.stat().st_mtime)
            sdt = sorted((F / cond).glob(f"{fov}_c*.sdt"))
            bliq = {c: sorted((B / cond / f"{fov}-001" / "Bliq VMS" / c).glob("*.tif"), key=_num_T)
                    for c in ("C1", "C2")}
            out.append(dict(nom=f"{cond.replace(' ', '')}_{fov}", type="fov", condition=cond,
                            cl_mM=_cl_de(cond), cycles=cycles, sdt=sdt[-1] if sdt else None, bliq=bliq))
    groupes = {}
    for p in sorted(F.iterdir()):
        if p.is_dir() and "_to_" in p.name:
            m = re.match(r"^(.+_to_[^_]+)_(\d+)$", p.name)   # basal_to_0_1 -> basal_to_0, partie 1
            base, partie = (m.group(1), int(m.group(2))) if m else (p.name, 0)
            groupes.setdefault(base, []).append((partie, p))
    trans = []
    for base, parts in groupes.items():
        parts = [p for _, p in sorted(parts)]
        cycles, bliq = [], {"C1": [], "C2": []}
        for p in parts:
            cycles += sorted(p.glob("*.spc"), key=_num_cycle)
            for c in bliq:
                bliq[c] += sorted((B / f"{p.name}-001" / "Bliq VMS" / c).glob("*.tif"), key=_num_T)
        avant, apres = base.split("_to_")
        trans.append(dict(nom=base, type="transition", condition=f"{avant}→{apres}",
                          avant=avant if avant == "basal" else f"{avant} mM", apres=f"{apres} mM",
                          cl_mM=_cl_de(f"{apres} mM"), cycles=cycles, sdt=None, bliq=bliq))
    # transitions dans l'ordre de l'expérience
    out += sorted(trans, key=lambda a: a["cycles"][0].stat().st_mtime if a["cycles"] else 0)
    return out


def acquisition(root, nom):
    for a in acquisitions(root):
        if a["nom"] == nom:
            return a
    sys.exit(f"acquisition inconnue : {nom}")


def irf_sdt(root):
    """Fichier IRF du jeu (data/FLIM/irf*.sdt)."""
    f = sorted((root / "data" / "FLIM").glob("irf*.sdt"))
    if not f:
        sys.exit("Aucune IRF (data/FLIM/irf*.sdt) dans le jeu")
    return f[0]


def _sur_grille(t, y):
    """Histogramme à pas constant (bords t[0..n]) -> 256 canaux de P/256, somme 1, par la
    fonction cumulée ; ce qui dépasse la période est replié au début."""
    dt = t[1] - t[0]
    bords = np.r_[t, t[-1] + dt]
    cum = np.r_[0, np.cumsum(y)]
    g = np.arange(N_CANAUX + 1) * DT_NS
    z = np.diff(np.interp(g, bords, cum)) + np.diff(np.interp(g + PERIODE_NS, bords, cum, right=cum[-1]))
    return z / z.sum()


def charger_irf(canal, root):
    """IRF d'un canal (« cfp » ou « yfp ») sur la grille des déclins : fond plat retiré (médiane
    des canaux avant le pic, au-dessus de la limite basse du TDC), somme 1."""
    import sdtfile
    s = sdtfile.SdtFile(irf_sdt(root))
    y = s.data[0][IRF_LIGNE[canal]].astype(float)
    t = s.times[0] * 1e9
    pic = int(np.argmax(y))
    k = np.arange(len(y))
    fond = np.median(y[(t >= 0.9) & (k < pic - 3)])
    y = np.clip(y - fond, 0, None)
    y[(k < pic - 6) | (t >= PERIODE_NS)] = 0
    return _sur_grille(t, y)


def temps_canaux():
    """Centres des 256 canaux d'une période, en ns."""
    return (np.arange(N_CANAUX) + 0.5) * DT_NS
