"""Configuration partagée — sessions FLIMApp en mode Realtime (passes), sans images.

Un jeu de données est un dossier qui contient `data/FLIM/` : une session FLIMApp telle
qu'écrite par l'appli (run.toml, frames.csv, visits.csv, readback.bin, irf.csv,
spc/<série>.spc + _acquisition.ini). Tout ce que produit le pipeline est écrit à côté de
`data/`, dans ce même dossier (mesures/, qc/, results/).

Le dossier se donne en premier argument de chaque étape ; à défaut, le dossier courant.
"""
import sys
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

# --- canaux ------------------------------------------------------------------------
# Canal FLIMApp k = entrée IN<k> de la QC-104 ([spc] series de run.toml). Le 2026-11-11,
# IN1 = YFP et IN2 = CFP, comme le 2026-11-09 et à l'inverse des sessions d'octobre
# (ch1 = CFP) : IN1 a ~6 fois plus de photons et une durée de vie de ~2,9 ns, IN2 ~1,5 ns
# (donneur raccourci par le FRET). L'étape 1 vérifie cette attribution sur les déclins
# (le CFP a le temps moyen le plus court) et s'arrête si elle ne tient pas.
CANAL = {"yfp": 1, "cfp": 2}
NOMS = {"cfp": "CFP", "yfp": "YFP"}

# --- temps ---------------------------------------------------------------------------
# FLIMApp rééchantillonne les 4096 canaux du TDC (4 ps, 16,384 ns) sur [qc] fenetre_ns =
# 12,5 ns, puis les déclins du Realtime en ont 256 : 48,83 ps par canal. La période du
# laser (80,25 MHz) est plus courte, 12,461 ns : les canaux au-delà de 255,2 sont vides, et
# la limite basse du TDC (5 % de 16,384 ns = 0,82 ns) vide les canaux 0 à 16.
N_CANAUX = 256
FENETRE_NS = 12.5
DT_NS = FENETRE_NS / N_CANAUX
F_LASER_HZ = 80.2475e6                    # SYNC mesuré par la carte (debug_report : 8,02e7 /s)
PERIODE_NS = 1e9 / F_LASER_HZ             # 12,461 ns
LIMITE_BASSE_NS = 0.05 * 16.384           # limite basse du TDC
FENETRE = (17, 255)                       # canaux pleins : [17, 255) — 16 et 255 sont entamés

# --- découpage du temps --------------------------------------------------------------
GROUPE_S = 2.0          # tranches d'analyse (passes regroupées) pour les durées de vie
PRE_S = 10.0            # base : les 10 s avant le début de la stim
POST_S = 10.0           # réponse : les 10 s après la fin de la stim (convention des tranches)
EXCLURE_APRES_STIM_S = 60.0   # dérive estimée hors des stims et de ce qui les suit


# ------------------------------------------------------------------ jeu de données
def racine(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    for a in argv:
        if not a.startswith("-") and (Path(a) / "data" / "FLIM").is_dir():
            return Path(a).resolve()
    if (Path.cwd() / "data" / "FLIM").is_dir():
        return Path.cwd()
    sys.exit("Jeu de données introuvable : donne le dossier qui contient data/FLIM "
             "(ex. python3 step1_decode.py ../SCLM_realtime_2026-11-11)")


def dirs(root):
    d = {k: root / k for k in ("mesures", "results", "qc")}
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d


def lire_ini(path):
    out, sec = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            sec = line[1:-1]
            out[sec] = {}
        elif "=" in line and sec is not None:
            k, v = line.split("=", 1)
            out[sec][k.strip()] = v.strip()
    return out


def session(root):
    """Ce que dit la session : protocole, fichiers .spc par canal FLIMApp, rythme des créneaux."""
    flim = root / "data" / "FLIM"
    run = tomllib.load(open(flim / "run.toml", "rb"))
    spc = {}
    for k, s in enumerate(run["spc"]["series"], start=1):
        nom = s.replace("/", "_")
        spc[k] = {"serie": s, "spc": flim / "spc" / f"{nom}.spc",
                  "ini": lire_ini(flim / "spc" / f"{nom}_acquisition.ini")}
    p, d = run["protocol"], run["daq"]
    fs = float(d["sample_rate_hz"])
    return dict(run=run, flim=flim, spc=spc, fs=fs,
                scan=int(d["scan_samples"]), pause=int(d["shift_samples"]),
                pwm=tuple(int(x) for x in d.get("pwm_samples", (0, 0))),
                tic_s=float(spc[1]["ini"]["acquisition"]["tic_s"]),
                rois_actives=bool(run["roi"]["active"]) and len(run["roi"]["list"]) > 0)


def readback(root):
    """Échantillons relus par la NI (10 kHz) : tableau (n, 12) et noms des signaux."""
    flim = root / "data" / "FLIM"
    noms = None
    for line in open(flim / "readback.txt", encoding="utf-8"):
        if line.startswith("signals"):
            noms = [s.strip() for s in line.split("=", 1)[1].split(",")]
    r = np.fromfile(flim / "readback.bin", dtype="<f4").reshape(-1, len(noms))
    return r, noms


def frames(root):
    return pd.read_csv(root / "data" / "FLIM" / "frames.csv")


def visits(root):
    return pd.read_csv(root / "data" / "FLIM" / "visits.csv")


def irf(root, canal):
    """IRF de FLIMApp pour ce canal (irf.csv, 256 canaux de 12,5/256 ns), somme 1."""
    a = np.loadtxt(root / "data" / "FLIM" / "irf.csv", delimiter=",", skiprows=1)
    y = np.clip(a[:, CANAL[canal]], 0, None)
    return y / y.sum()


def temps_canaux():
    """Centres des 256 canaux, en ns."""
    return (np.arange(N_CANAUX) + 0.5) * DT_NS


def charger_passes(D):
    """Ce que l'étape 1 a rangé (mesures/passes.npz), en dictionnaire."""
    z = np.load(D["mesures"] / "passes.npz")
    return {k: z[k] for k in z.files}
