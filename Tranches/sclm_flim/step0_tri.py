"""Étape 0 — tri d'une journée d'acquisitions : ne garder que les trains de stimulation.

Usage :  python3 step0_tri.py <dossier de la journée> <dossier d'analyse>
ex.      python3 step0_tri.py ../SCLM+eNpHR_2026-10-08/data/2026-10-09 ../SCLM+eNpHR_2026-10-09

Une journée = des sessions FLIMApp (dossiers <date>_<heure> : run.toml, visits.csv,
spc/…) et des images fixes Bliq (dossiers FOVk_<λ>[nm][_pre|_post]-NNN).

Une session est **gardée** si la commande de stimulation (visits.csv, commande 1, horloge
de la NI) a produit au moins MIN_STIMS stims complètes : chacune dure au moins la moitié
de la durée prévue par le protocole et se termine avant la fin de l'enregistrement. Les
autres — essais interrompus, sessions sans stim — sont **écartées**. Rien n'est effacé ni
déplacé : tri.csv dit ce qui a été écarté et pourquoi.

Chaque session gardée devient un jeu de données du pipeline, <dossier d'analyse>/<nom>/,
avec deux liens symboliques : data/FLIM vers la session, et data/FOVs/ vers les images
fixes de son champ. Un champ (FOVk) reçoit les sessions qui commencent entre sa première
image et la première image du champ suivant ; une session prise avant toute image fixe
n'en a pas.

Sortie : <dossier d'analyse>/tri.csv — toutes les sessions, gardées ou non.
"""
import csv
import math
import os
import sys
import tomllib
from datetime import datetime
from pathlib import Path

import config
import spc

MIN_STIMS = 2          # un train : au moins deux stims complètes
FRACTION_MIN = 0.5     # une stim complète dure au moins la moitié de la durée prévue
MARGE_FIN_S = 0.5      # ... et finit au moins 0,5 s avant la fin de l'enregistrement


def stims(session):
    """Fenêtres de la commande 1 [(début, fin, valeur)], en s depuis le début de la session."""
    run = tomllib.load(open(session / "run.toml", "rb"))
    fs = float(run["daq"]["sample_rate_hz"])
    slot = (float(run["protocol"]["scan_time"]) + float(run["protocol"]["shift_time"])) / 1000
    rows = list(csv.DictReader(open(session / "visits.csv")))
    out, debut = [], None
    for i, r in enumerate(rows + [None]):
        on = r is not None and float(r["command1_v"]) > 0
        if on and debut is None:
            debut = i
        elif not on and debut is not None:
            a = int(rows[debut]["first_sample"]) / fs
            b = int(rows[i - 1]["first_sample"]) / fs + slot
            out.append((a, b, float(rows[debut]["command1_v"])))
            debut = None
    return run, out


def classer(session):
    run, w = stims(session)
    ini = spc.lire_ini(session / "spc" / "3T0089_IN1_acquisition.ini")
    duree = float(ini["acquisition"]["duree_s"])
    prevue = next((t for t in run["protocol"].get("times", []) if not math.isnan(t)), float("nan"))
    completes = [(a, b, v) for a, b, v in w
                 if b - a >= FRACTION_MIN * prevue and b <= duree - MARGE_FIN_S]
    debut = datetime.fromisoformat(run["started"])
    info = dict(session=session.name, debut=debut.strftime("%H:%M:%S"), duree_s=round(duree, 1),
                n_rois_flimapp=len(run.get("roi", {}).get("list", [])),
                images_par_bloc=run["spc"]["images_par_bloc"],
                n_commandes=len(w), n_stims_completes=len(completes))
    if completes:
        d = [b - a for a, b, _ in completes]
        info.update(duree_stim_s=round(sum(d) / len(d), 1),
                    puissance_pct=round(100 * completes[0][2]),
                    intervalle_s=round(completes[1][0] - completes[0][0], 1) if len(completes) > 1 else "",
                    base_avant_1re_s=round(completes[0][0], 1),
                    apres_derniere_s=round(duree - completes[-1][1], 1))
    if len(completes) >= MIN_STIMS:
        info.update(gardee=True, raison=f"train de {len(completes)} stims")
        if completes[0][0] < config.PRE_WINDOW_S:
            info["raison"] += f" ; pas de ligne de base ({completes[0][0]:.1f} s avant la 1re)"
        if len(w) > len(completes):
            info["raison"] += f" ; {len(w) - len(completes)} stim incomplète ignorée"
    elif not w:
        info.update(gardee=False, raison=f"aucune stim envoyée (enregistrement de {duree:.0f} s)")
    else:
        info.update(gardee=False, raison=f"{len(completes)} stim complète sur {len(w)} "
                                         f"(enregistrement de {duree:.0f} s)")
    return info, debut


def champs(jour):
    """{FOVk: (heure de la 1re image, [dossiers d'images])}, triés par heure."""
    out = {}
    for d in sorted(p for p in jour.iterdir() if p.is_dir()):
        m = config._FOV_RE.match(d.name)
        tifs = sorted(d.glob("Bliq VMS/C*/*.tif"))
        if not m or not tifs:
            continue
        t = datetime.fromisoformat(config.ome_date(tifs[0]))
        fov = out.setdefault(m["fov"], [t, []])
        fov[0] = min(fov[0], t)
        fov[1].append(d)
    return dict(sorted(out.items(), key=lambda kv: kv[1][0]))


def lien(cible, lien_):
    lien_.parent.mkdir(parents=True, exist_ok=True)
    if lien_.is_symlink() or lien_.exists():
        lien_.unlink()
    lien_.symlink_to(os.path.relpath(cible, lien_.parent))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if len(args) != 2:
        sys.exit(__doc__)
    jour, sortie = Path(args[0]).resolve(), Path(args[1]).resolve()
    sessions = sorted(p for p in jour.iterdir() if p.is_dir() and (p / "run.toml").exists())
    fovs = champs(jour)
    rows = []
    for s in sessions:
        info, debut = classer(s)
        # le champ de la session : le dernier dont la 1re image précède son début
        fov = [k for k, (t, _) in fovs.items() if t <= debut.replace(tzinfo=t.tzinfo)]
        info["fov"] = fov[-1] if fov else ""
        if info["gardee"]:
            nom = (f"{info['fov'] + '_' if info['fov'] else ''}{debut:%H%M}_"
                   f"{info['n_stims_completes']}x{info['duree_stim_s']:.0f}s_{info['puissance_pct']}pct")
            info["jeu"] = nom
            lien(s, sortie / nom / "data" / "FLIM")
            for d in (fovs[info["fov"]][1] if info["fov"] else []):
                lien(d, sortie / nom / "data" / "FOVs" / d.name)
        rows.append(info)
        print(f"  {'GARDÉE ' if info['gardee'] else 'écartée'}  {s.name}  {info['debut']}  "
              f"{info['duree_s']:6.1f} s  {info['fov'] or '-':5s} {info['raison']}"
              + (f"  -> {info['jeu']}" if info.get("jeu") else ""))
    cles = ["session", "debut", "duree_s", "gardee", "raison", "jeu", "fov", "n_commandes",
            "n_stims_completes", "duree_stim_s", "puissance_pct", "intervalle_s",
            "base_avant_1re_s", "apres_derniere_s", "n_rois_flimapp", "images_par_bloc"]
    sortie.mkdir(parents=True, exist_ok=True)
    with open(sortie / "tri.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cles, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\n{sum(r['gardee'] for r in rows)} sessions gardées sur {len(rows)} -> {sortie}")


if __name__ == "__main__":
    main()
