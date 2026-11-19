"""Étape 9 — ROIs grossiers par Cellpose-SAM sur l'image reconstruite (étape 8).

Usage :  ~/.venvs/cellpose/bin/python step9_cellpose.py <dossier> [--diameter 50]

L'image est celle de tous les photons de la session, blocs recalés (mesures/images.npz,
« recale »), en racine carrée des photons comme pour les tranches d'octobre. Les lignes du
retour du galvanomètre (en bas) et les premières lignes après lui (en haut, où le balayage
n'est pas encore stable) sont écartées : un ROI dont le centre y tombe, ou dont plus de 20 %
y tombe, est retiré. Pixels
de ~51,8 ns en x (le pixel du Bliq, ~0,51 µm) et d'une ligne en y : les somas font ~50 px.

Sorties : masks/cellpose_masks.npy, masks/cellpose_params.json, results/rois.png
"""
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402
from cellpose import models  # noqa: E402

import config as C  # noqa: E402

HAUT, BAS = 25, 22          # lignes écartées en haut (balayage pas encore stable) et en bas (retour)
EVAL_KW = dict(flow_threshold=0.4, cellprob_threshold=0.0, min_size=300, max_size_fraction=0.25, normalize=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--diameter", type=float, default=50.0)
    a = ap.parse_args()
    root = C.racine([a.dossier] if a.dossier else [])
    D = C.dirs(root)
    (root / "masks").mkdir(exist_ok=True)
    z = np.load(D["mesures"] / "images.npz")
    brut = z["recale"].astype(np.float32)
    img = np.sqrt(np.maximum(brut, 0))
    gpu = torch.backends.mps.is_available()
    print(f"Cellpose-SAM (MPS={gpu}), diamètre {a.diameter:g} px")
    model = models.CellposeModel(gpu=gpu)
    masks, _, _ = model.eval(img.copy(), diameter=a.diameter, **EVAL_KW)
    masks = masks.astype(np.int32)
    # ROIs dans les bandes écartées : retirés ; renumérotés de 1 à n
    ny = masks.shape[0]
    bande = np.zeros_like(masks, bool)
    bande[:HAUT] = True
    bande[ny - BAS:] = True
    yy = np.arange(ny)[:, None]
    garde = [r for r in range(1, masks.max() + 1)
             if (masks == r).any() and (bande & (masks == r)).sum() < 0.2 * (masks == r).sum()
             and HAUT <= (yy * (masks == r)).sum() / (masks == r).sum() < ny - BAS]
    out = np.zeros_like(masks)
    for k, r in enumerate(garde, start=1):
        out[masks == r] = k
    n = int(out.max())
    np.save(root / "masks" / "cellpose_masks.npy", out)
    (root / "masks" / "cellpose_params.json").write_text(json.dumps(
        {"modele": "cpsam", "diametre_px": a.diameter, "image": "racine des photons, blocs recalés",
         "lignes_ecartees": [HAUT, BAS], "n_rois": n, **EVAL_KW}, indent=1))

    fig, ax = plt.subplots(figsize=(16, 8.6))
    ax.imshow(img, cmap="gray", vmin=0, vmax=np.percentile(img, 99.7))
    from scipy import ndimage as ndi
    for r in range(1, n + 1):
        m = out == r
        ax.contour(m, levels=[0.5], colors=["#eb6834"], linewidths=1.0)
        cy, cx = ndi.center_of_mass(m)
        ax.text(cx, cy, str(r), color="#ffd166", fontsize=10, ha="center", va="center", fontweight="bold")
    ax.axhline(HAUT, color="#2a78d6", lw=0.8, ls="--")
    ax.axhline(ny - BAS, color="#2a78d6", lw=0.8, ls="--")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{root.name} — image reconstruite sans M1/M2 (√photons, toute la session) ; Cellpose-SAM, "
                 f"diamètre {a.diameter:g} px : {n} ROIs ; tirets : bandes écartées", loc="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(D["results"] / "rois.png", dpi=110, facecolor="#fcfcfb")
    print(f"{n} ROIs -> masks/cellpose_masks.npy, results/rois.png")


if __name__ == "__main__":
    main()
