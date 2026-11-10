"""Étape 2 — segmentation des somas par Cellpose-SAM sur l'image FLIM (ch1 + ch2).

Usage :  ~/.venvs/cellpose/bin/python step2_cellpose.py <dossier> [--diameter 100] [--image racine|brut]
         ~/.venvs/cellpose/bin/python step2_cellpose.py <dossier> --masques-de <autre jeu du même champ>
Sortie : masks/cellpose_masks.npy, masks/cellpose_overlay.png, masks/cellpose_params.json

L'image segmentée est la somme des photons des deux canaux sur toute la session
(~450 photons par pixel) : c'est le repère exact des déclins, sans recalage. Par défaut,
Cellpose voit sa racine carrée (`--image racine`) : sur les photons bruts, les somas les
plus brillants saturent après la normalisation et sont manqués (2026-10-09 : 25 ROIs au
lieu de 31 sur un champ). La racine resserre la dynamique et stabilise la variance d'un
comptage de photons. Le 2026-10-08 a été segmenté sur les photons bruts (`--image brut`,
34 ROIs ; 32 avec la racine). Les ROIs
dessinées dans FLIMApp pendant l'expérience sont tracées en cyan pour comparaison ;
elles sont mesurées à part (étape 4), telles quelles.

Plusieurs sessions sur le même champ : `--masques-de` reprend les ROIs d'une autre session
au lieu de relancer Cellpose, déplacés du décalage entre les deux images FLIM. Le même
numéro désigne alors la même cellule d'une session à l'autre, comme au 2026-09-21.
"""
import argparse
import json
import sys

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cellpose import models

import config
from config import CELL_DIAMETER_PX, PIXEL_SIZE_UM

# Mêmes réglages que les tranches 2026-09-21 : somas compacts et contrastés. min_size
# écarte débris et taches de neuropile ; max_size_fraction évite qu'un ROI avale le champ.
EVAL_KW = dict(flow_threshold=0.4, cellprob_threshold=0.0,
               min_size=300, max_size_fraction=0.25, normalize=True)


def overlay(img, masks, polys, title, out):
    lo, hi = np.percentile(img, (1, 99.7))
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
    for lab in range(1, int(masks.max()) + 1):
        m = masks == lab
        ax.contour(m, levels=[0.5], colors="#ff3b30", linewidths=1.0)
        ys, xs = np.nonzero(m)
        ax.text(xs.mean(), ys.mean(), str(lab), color="#ffd60a",
                fontsize=7, ha="center", va="center", fontweight="bold")
    for r in polys:
        ax.plot(r["xs"], r["ys"], color="#00e5ff", lw=1.2, ls="--")
        ax.text(np.mean(r["xs"][:-1]), np.mean(r["ys"][:-1]), f"F{r['index']}", color="#00e5ff",
                fontsize=8, ha="center", va="center", fontweight="bold")
    ax.set_title(title, fontsize=12)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=115)
    plt.close(fig)


def reprendre(root, source, img_ici):
    """ROIs d'une autre session du même champ, déplacés du décalage entre les deux images
    FLIM (corrélation de phase, au pixel près) : même numéro = même cellule."""
    from scipy import ndimage as ndi
    from recalage import _prep, correlation_phase, ncc
    src = config.racine([source])
    ms = np.load(src / "masks" / "cellpose_masks.npy")
    ref = np.load(src / "projections" / "flim.npz")["somme"].astype(float)
    A, B = _prep(img_ici), _prep(ref)
    dy, dx = correlation_phase(A, B, sous_pixel=True)
    s = (int(round(dy)), int(round(dx)))
    m = np.zeros(A.shape, bool)
    m[60:-60, 60:-60] = True
    score = ncc(A, ndi.shift(B, s, order=0), m)
    if score < 0.5:
        sys.exit(f"{src.name} ne ressemble pas assez à {root.name} (NCC {score:.2f}) : "
                 "pas le même champ ?")
    masks = ndi.shift(ms, s, order=0, cval=0).astype(np.int32)
    return masks, dict(origine=src.name, decalage_px=[dy, dx], decalage_applique=list(s), ncc=score)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--diameter", type=float, default=CELL_DIAMETER_PX)
    ap.add_argument("--image", choices=("racine", "brut"), default="racine")
    ap.add_argument("--masques-de", dest="source",
                    help="autre jeu du même champ dont reprendre les ROIs (recalés), au lieu de Cellpose")
    args = ap.parse_args()
    root = config.racine([args.dossier] if args.dossier else [])
    D = config.dirs(root)
    S = config.session(root)
    p = D["projections"] / "flim.npz"
    if not p.exists():
        sys.exit("projections/flim.npz absent : lance d'abord step1_decode.py")
    brut = np.load(p)["somme"].astype(np.float32)
    img = np.sqrt(brut) if args.image == "racine" else brut

    if args.source:
        masks, params = reprendre(root, args.source, brut)
        titre = (f"{root.name} — ROIs repris de {params['origine']} (décalés de "
                 f"{params['decalage_applique'][0]:+d}, {params['decalage_applique'][1]:+d} px ; "
                 f"NCC {params['ncc']:.2f})")
    else:
        gpu = torch.backends.mps.is_available()
        print(f"Cellpose-SAM (MPS={gpu}), diamètre={args.diameter:g} px")
        model = models.CellposeModel(gpu=gpu)
        # eval() normalise son entrée en place : on lui donne une copie.
        masks, _, _ = model.eval(img.copy(), diameter=args.diameter, **EVAL_KW)
        masks = masks.astype(np.int32)
        params = {"modele": "cpsam", "diametre_px": args.diameter, "image": args.image, **EVAL_KW}
        titre = (f"{root.name} — Cellpose-SAM sur {'√photons' if args.image == 'racine' else 'photons'}, "
                 f"diamètre {args.diameter:g} px (~{args.diameter * PIXEL_SIZE_UM:.0f} µm)")
    n = int(masks.max())

    # De nouveaux ROIs rendent caduque une sélection faite sur les anciens : elle est mise
    # de côté, pas effacée.
    ancien = D["masks"] / "cellpose_masks.npy"
    sel = D["selection"] / "cellpose_selection.json"
    if sel.exists() and (not ancien.exists() or not np.array_equal(np.load(ancien), masks)):
        sel.rename(sel.with_name("cellpose_selection_ancienne.json"))
        print(f"  ROIs changés : l'ancienne sélection est mise de côté ({sel.parent.name}/"
              "cellpose_selection_ancienne.json) ; refais le tri (étape 3)")
    np.save(ancien, masks)
    (D["masks"] / "cellpose_params.json").write_text(json.dumps({**params, "n_rois": n}, indent=1))
    overlay(img, masks, S["rois_flimapp"], f"{titre} — {n} ROIs ; cyan : ROIs FLIMApp",
            D["masks"] / "cellpose_overlay.png")
    print(f"{n} ROIs -> {D['masks']}")


if __name__ == "__main__":
    main()
