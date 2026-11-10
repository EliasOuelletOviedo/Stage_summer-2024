"""Étape 2 — segmentation des somas par Cellpose-SAM, une fois par acquisition.

Usage :  ~/.venvs/cellpose/bin/python step2_cellpose.py <dossier> [acquisition …] [--diameter 30]
Sorties : masks/<acq>_masks.npy, masks/<acq>_overlay.png, masks/cellpose_params.json

Image segmentée : racine de la somme des photons CFP + YFP de toute l'acquisition (cycles
recalés, étape 1), dans le repère exact des déclins. La racine resserre la dynamique (les
somas les plus brillants ne saturent pas la normalisation) et stabilise la variance d'un
comptage de photons (sclm_flim, 2026-10-09). Les acquisitions d'une condition fixe sont
des champs différents (vérifié le 2026-11-09, NCC < 0,3 entre conditions) : pas de reprise
des ROIs d'une acquisition à l'autre.
"""
import argparse
import json

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cellpose import models

import config as C

# Réglages des tranches 2026-08-17 (même objectif, même pixel) : somas compacts et
# contrastés ; min_size écarte débris et taches de neuropile ; max_size_fraction évite
# qu'un ROI avale le champ.
EVAL_KW = dict(flow_threshold=0.4, cellprob_threshold=0.0, min_size=C.MIN_TAILLE_PX,
               max_size_fraction=0.25, normalize=True)


def overlay(img, masks, title, out):
    lo, hi = np.percentile(img, (1, 99.7))
    fig, ax = plt.subplots(figsize=(14, 7.4))
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
    for lab in range(1, int(masks.max()) + 1):
        m = masks == lab
        if not m.any():
            continue
        ax.contour(m, levels=[0.5], colors="#ff3b30", linewidths=0.9)
        ys, xs = np.nonzero(m)
        ax.text(xs.mean(), ys.mean(), str(lab), color="#ffd60a", fontsize=6.5, ha="center",
                va="center", fontweight="bold")
    ax.set_title(title, fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("acquisitions", nargs="*")
    ap.add_argument("--diameter", type=float, default=C.DIAMETRE_PX)
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    A = [a for a in C.acquisitions(root) if not args.acquisitions or a["nom"] in args.acquisitions]
    gpu = torch.backends.mps.is_available()
    print(f"Cellpose-SAM (MPS={gpu}), diamètre={args.diameter:g} px, {len(A)} acquisitions")
    model = models.CellposeModel(gpu=gpu)
    p_json = D["masks"] / "cellpose_params.json"
    params = json.loads(p_json.read_text()) if p_json.exists() else {}
    for a in A:
        z = np.load(D["projections"] / f"{a['nom']}.npz")
        img = np.sqrt(z["somme"].sum(0)).astype(np.float32)
        # eval() normalise son entrée en place : on lui donne une copie.
        masks, _, _ = model.eval(img.copy(), diameter=args.diameter, **EVAL_KW)
        masks = masks.astype(np.int32)
        n = int(masks.max())
        ancien = D["masks"] / f"{a['nom']}_masks.npy"
        sel = D["selection"] / f"{a['nom']}_selection.json"
        if sel.exists() and (not ancien.exists() or not np.array_equal(np.load(ancien), masks)):
            sel.rename(sel.with_name(f"{a['nom']}_selection_ancienne.json"))
            print(f"  {a['nom']} : ROIs changés, ancienne sélection mise de côté")
        np.save(ancien, masks)
        params[a["nom"]] = {"modele": "cpsam", "diametre_px": args.diameter, "image": "racine",
                            **EVAL_KW, "n_rois": n}
        overlay(img, masks, f"{a['nom']} — Cellpose-SAM sur √photons (CFP + YFP), diamètre "
                f"{args.diameter:g} px (~{args.diameter * C.PIXEL_UM:.0f} µm) — {n} ROIs",
                D["masks"] / f"{a['nom']}_overlay.png")
        print(f"  {a['nom']:>12} : {n} ROIs")
    p_json.write_text(json.dumps(params, indent=1))


if __name__ == "__main__":
    main()
