"""Configuration partagée — SCLM + eNpHR, sessions FLIMApp (Realtime en images).

Un jeu de données est un dossier qui contient `data/FLIM/` (une session FLIMApp :
run.toml, frames.csv, visits.csv, irf.csv, spc/*.spc) et `data/FOVs/` (images fixes
Bliq, ex. FOV1_800nm-001, FOV1_950nm_post-001). Tout ce que produit le pipeline est
écrit à côté de `data/`, dans ce même dossier.

Le dossier se donne en premier argument de chaque étape, ou par la variable
SCLM_FLIM_DATA ; à défaut, le dossier courant.
"""
import csv
import math
import os
import re
import sys
import tomllib
from pathlib import Path

import numpy as np

# --- canaux ----------------------------------------------------------------------
# Les deux entrées de la QC-104 (IN1, IN2) deviennent les canaux 1 et 2 de FLIMApp.
# ch1 = CFP (donneur) : durée de vie ~1,3 ns, raccourcie par le FRET ; ch2 = YFP
# (accepteur) : ~3,7 fois plus de photons, durée de vie apparente ~2,4 ns. C'est aussi
# le choix de FLIMApp pendant l'expérience (durée de vie suivie sur ch1, rapport
# ch2/ch1 affiché). Ratio FRET = ch2 / ch1, comme C2/C3 (YFP/CFP) des données Bliq :
# une baisse est une entrée de Cl⁻.
CANAUX = (1, 2)
NOMS_CANAUX = {1: "CFP", 2: "YFP"}

CANAUX_TEMPS = 256          # canaux des déclins (ceux du Realtime et de l'IRF)
FENETRE_NS = 12.5           # période du laser, fenêtre des déclins

# --- géométrie ---------------------------------------------------------------------
PIXEL_SIZE_UM = 0.5087830142943353   # métadonnée OME des images Bliq (objectif 16X)
CELL_DIAMETER_PX = 100               # 2026-10-08 : à 60 px (réglage du 09-21) Cellpose
                                     # manquait des somas brillants ; aire médiane des ROIs
                                     # 8 800 px, soit ~106 px de diamètre équivalent
BIN_DERIVE = 4                       # images de bloc réduites 4 x 4 pour la dérive

# --- recalage ------------------------------------------------------------------------
ECHELLES_X = np.arange(0.97, 1.0501, 0.0025)   # échelle x Bliq -> FLIM essayée
# Le pixel FLIMApp vaut 25 tics = 51,2 ns ; celui de Bliq 51,8 ns (OME : dwell time).
# Même scanner, même ligne : l'image Bliq est un peu plus large, d'où l'échelle en x.

# Une image fixe n'est au plan de la session que si elle ressemble à l'image de référence
# (celle qui ressemble le plus au FLIM). Le 2026-10-08 : même plan à deux longueurs
# d'onde, NCC 0,58–0,60 ; plans différents à la même longueur d'onde, 0,33–0,38.
MIN_NCC_MEME_PLAN = 0.5

MIN_VALID_FRAC = 0.90       # part du ROI qui doit rester dans le champ après recalage

# --- stimulation ------------------------------------------------------------------
# La commande 1 du protocole (visits.csv, horloge de la NI) donne la fenêtre de stim ;
# elle est recalée sur l'horloge des cartes SPC par la fuite de la stim dans ch2
# (moyenne par image), qui la montre à 60 ms près.
STIM_GARDE_S = 0.0          # garde ajoutée de chaque côté de la fenêtre
RESPONSE_WINDOW_S = 10.0    # réponse : moyenne des 10 s après la fin de la stim
PRE_WINDOW_S = 10.0         # base : moyenne des 10 s avant le début de la stim
BASE_SKIP = 1               # premiers blocs ignorés dans R0 (transitoire de départ)
# Tranche d'analyse de l'étape 5 : les blocs FLIMApp (images_par_bloc, réglage d'affichage
# du direct) sont regroupés pour faire ~3 s, assez de photons pour ajuster τ bloc par bloc.
# 2026-10-08 : blocs de 50 images (3,0 s), gardés tels quels ; 2026-10-09 : 17 images
# (1,0 s), regroupés par 3.
ANALYSE_BLOC_S = 3.0

# Descente graduelle du ratio (étape 6) : estimée entre les stims seulement, en écartant
# les EXCLURE_APRES_STIM_S secondes qui suivent chaque stim, où le ratio remonte en
# partie. Le 2026-10-08 (champ entier), la pente vaut −0,2 %/100 s sans exclusion, puis
# −0,65, −0,93, −1,06 et −1,14 %/100 s en excluant 10, 20, 30 et 45 s ; elle se stabilise
# vers 30 s, et rejoint celle de la base avant la 1re stim (−0,90 %/100 s, sans stim).
EXCLURE_APRES_STIM_S = 30.0

# --- ajustement des durées de vie ---------------------------------------------------
# Mono-exponentielle reconvoluée par l'IRF, maximum de vraisemblance de Poisson : le
# modèle de FLIMApp (src/analysis/lifetime_analysis.jl). Le décalage de l'IRF est
# ajusté une fois par ROI et par canal sur le déclin de toute l'acquisition, puis fixé
# bloc par bloc (il est instrumental) : seuls τ et le fond restent libres.
COUPE_BAS = 13              # canaux ignorés au début (TAC limit low, 5,1 %)
COUPE_HAUT = 13             # ... et à la fin (TAC limit high, 94,9 %)
MIN_PHOTONS_FIT = 500       # en dessous, pas d'ajustement (comme FLIMApp)


# ------------------------------------------------------------------ jeu de données
def racine(argv=None):
    """Dossier du jeu de données : 1er argument qui est un dossier, SCLM_FLIM_DATA, ou cwd."""
    argv = sys.argv[1:] if argv is None else argv
    for a in argv:
        if not a.startswith("-") and (Path(a) / "data").is_dir():
            return Path(a).resolve()
    env = os.environ.get("SCLM_FLIM_DATA")
    if env:
        return Path(env).resolve()
    if (Path.cwd() / "data" / "FLIM").is_dir():
        return Path.cwd()
    sys.exit("Jeu de données introuvable : donne le dossier qui contient data/FLIM "
             "(ex. python3 step1_decode.py ../SCLM+eNpHR_2026-10-08)")


def dirs(root):
    d = {k: root / k for k in ("projections", "masks", "selection", "mesures", "results", "qc")}
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d


def _csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def session(root):
    """Ce que dit la session FLIMApp : protocole, ROIs dessinées, fichiers .spc,
    fenêtres de commande (horloge NI)."""
    flim = root / "data" / "FLIM"
    run = tomllib.load(open(flim / "run.toml", "rb"))
    series = run["spc"]["series"]                       # ["3T0089/IN1", "3T0089/IN2"]
    spc = {}
    for k, s in enumerate(series, start=1):
        nom = s.replace("/", "_")
        spc[k] = {"serie": s, "spc": flim / "spc" / f"{nom}.spc",
                  "ini": flim / "spc" / f"{nom}_acquisition.ini"}

    visits = _csv(flim / "visits.csv")
    fs = float(run["daq"]["sample_rate_hz"])
    slot = (float(run["protocol"]["scan_time"]) + float(run["protocol"]["shift_time"])) / 1000
    on = np.array([_float(v["command1_v"]) > 0 for v in visits])
    t_slot = np.array([int(v["first_sample"]) / fs for v in visits])
    commandes = []
    if on.any():
        d = np.diff(np.r_[0, on.astype(int), 0])
        for a, b in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1):
            commandes.append((float(t_slot[a]), float(t_slot[b] + slot),
                              float(_float(visits[a]["command1_v"]))))

    rois = [{"index": r["index"], "nom": r["name"], "xs": r["xs"], "ys": r["ys"]}
            for r in run["roi"]["list"]] if run.get("roi", {}).get("active") else []
    return {"run": run, "spc": spc, "commandes": commandes, "rois_flimapp": rois,
            "images_par_bloc": int(run["spc"]["images_par_bloc"]),
            "image_size": tuple(run["roi"]["image_size"]),
            "debut": run.get("started", ""), "protocole": run["protocol"]}


def irf(root):
    """IRF par canal (256 canaux, temps croissant) et pas en ns, depuis irf.csv."""
    rows = _csv(root / "data" / "FLIM" / "irf.csv")
    t = np.array([float(r["time_ns"]) for r in rows])
    return {k: np.array([float(r[f"ch{k}"]) for r in rows]) for k in CANAUX}, float(t[1] - t[0])


def frames_flimapp(root):
    """frames.csv : ce que FLIMApp a mesuré en direct, une ligne par ROI et par bloc."""
    return _csv(root / "data" / "FLIM" / "frames.csv")


# ------------------------------------------------------------------ images fixes Bliq
# Deux façons de nommer les images fixes ont servi : FOV1_800nm-001 / FOV1_800nm_post-001
# (2026-10-08) et FOV2_800_pre-001 / FOV3_950_post-001 (2026-10-09). Sans « _post »,
# l'image est d'avant la session.
_FOV_RE = re.compile(r"^(?P<fov>[^_]+)_(?P<nm>\d+)\s*(?:nm)?(?:_(?P<moment>pre|post))?(?:-|$)", re.I)


def fovs(root):
    """Images fixes Bliq : [{nom, fov, nm, moment ('avant'/'après'), canaux: {C1: path}}]."""
    out = []
    base = root / "data" / "FOVs"
    if not base.is_dir():
        return out
    for d in sorted(p for p in base.iterdir() if p.is_dir()):
        m = _FOV_RE.match(d.name)
        if not m:
            continue
        ch = {}
        for c in sorted((d / "Bliq VMS").glob("C*")):
            tifs = sorted(c.glob("*.tif"))
            if tifs:
                ch[c.name] = tifs[0]
        out.append({"nom": d.name, "fov": m["fov"], "nm": int(m["nm"]),
                    "moment": "après" if (m["moment"] or "").lower() == "post" else "avant",
                    "canaux": ch})
    return out


def ome_date(path):
    import tifffile
    with tifffile.TiffFile(path) as tf:
        xml = tf.ome_metadata or ""
    m = re.search(r"<AcquisitionDate>(.*?)</", xml)
    return m.group(1) if m else ""


# ------------------------------------------------------------------ ROIs FLIMApp
def masque_polygone(xs, ys, nx=1024, ny=1024):
    """Pixels d'une ROI FLIMApp : ceux dont le centre (x, y), en indices de 0, est dans
    le polygone, règle pair-impair — exactement `roi_pixel_mask` (gui/roi_popup.jl)."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    m = np.zeros((ny, nx), bool)
    x0, x1 = max(int(np.floor(xs.min())), 0), min(int(np.ceil(xs.max())), nx - 1)
    y0, y1 = max(int(np.floor(ys.min())), 0), min(int(np.ceil(ys.max())), ny - 1)
    px, py = np.meshgrid(np.arange(x0, x1 + 1, dtype=float), np.arange(y0, y1 + 1, dtype=float))
    dedans = np.zeros(px.shape, bool)
    j = len(xs) - 1
    for i in range(len(xs)):
        xi, yi, xj, yj = xs[i], ys[i], xs[j], ys[j]
        with np.errstate(divide="ignore", invalid="ignore"):
            dedans ^= ((yi > py) != (yj > py)) & (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
        j = i
    m[y0:y1 + 1, x0:x1 + 1] = dedans
    return m
