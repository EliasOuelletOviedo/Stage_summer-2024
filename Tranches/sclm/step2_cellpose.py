"""Étape 2 — segmentation des somas par Cellpose-SAM sur les images somme.

Usage :  ~/.venvs/cellpose/bin/python step2_cellpose.py [--diameter 30] [Slice1_FOV3 ...]
Sortie : masks/<FOV>_masks.npy  +  masks/<FOV>_overlay.png
"""
import argparse
import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cellpose import models

from config import CELL_DIAMETER_PX, PIXEL_SIZE_UM, SLICES, N_FOV, fov_key, dirs

# Somas de neurones en tranche : compacts, contrastés, ~15 µm.
# min_size écarte les débris ; max_size_fraction évite qu'un ROI avale le champ.
EVAL_KW = dict(flow_threshold=0.4, cellprob_threshold=0.0,
               min_size=120, max_size_fraction=0.25, normalize=True)


def overlay(img, masks, key, out, diameter):
    lo, hi = np.percentile(img, (1, 99.7))
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
    n = int(masks.max())
    for lab in range(1, n + 1):
        m = masks == lab
        ax.contour(m, levels=[0.5], colors="#ff3b30", linewidths=1.0)
        ys, xs = np.nonzero(m)
        ax.text(xs.mean(), ys.mean(), str(lab), color="#ffd60a",
                fontsize=7, ha="center", va="center", fontweight="bold")
    ax.set_title(f"{key} — Cellpose-SAM, diamètre {diameter} px "
                 f"(~{diameter*PIXEL_SIZE_UM:.0f} µm) — {n} ROIs", fontsize=12)
    ax.axis("off")
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diameter", type=float, default=CELL_DIAMETER_PX)
    ap.add_argument("fovs", nargs="*")
    args = ap.parse_args()

    D = dirs()
    keys = args.fovs or [fov_key(s, i) for s in SLICES for i in range(N_FOV)]
    gpu = torch.backends.mps.is_available()
    print(f"Cellpose-SAM (MPS={gpu}), diamètre={args.diameter:g} px")
    model = models.CellposeModel(gpu=gpu)

    total = 0
    for key in keys:
        img = np.load(D["registered"] / f"{key}.npz")["summed"].astype(np.float32)
        # eval() normalise son entrée en place : on lui donne une copie.
        masks, _, _ = model.eval(img.copy(), diameter=args.diameter, **EVAL_KW)
        masks = masks.astype(np.int32)
        np.save(D["masks"] / f"{key}_masks.npy", masks)
        overlay(img, masks, key, D["masks"] / f"{key}_overlay.png", args.diameter)
        n = int(masks.max()); total += n
        print(f"  {key:14s} {n:4d} ROIs")
    print(f"\n{total} ROIs sur {len(keys)} FOVs -> {D['masks']}")


if __name__ == "__main__":
    main()
