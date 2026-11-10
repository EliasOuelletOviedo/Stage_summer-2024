"""Étape 4 — déclins par ROI, par cycle et par canal (second passage sur les photons).

Usage :  python3 step4_declins.py <dossier> [--processus 4]

Pour chaque acquisition : un déclin de 256 canaux par ROI Cellpose (tous, triés ou non :
le tri s'applique à l'étape 6), par cycle et par canal, plus le neuropile ; et les photons
de chaque ROI par seconde de cycle (cinétiques rapides des transitions, étape 8).

Les masques suivent la dérive mesurée à l'étape 1 : déplacés, à chaque cycle, du décalage
arrondi au pixel (sous le demi-pixel : rien). Neuropile : pixels à plus de 10 px de tout
soma et au-dessus du 25e centile de l'image somme (comme sclm_flim).

Sortie : mesures/<acq>.npz — declins (cycles, rois + 2, 2, 256), secondes (cycles, rois + 2,
2, 11), etiquettes (0 = hors ROIs, 1..n = somas, n + 1 = neuropile), decalages appliqués.
"""
import argparse
import sys
from multiprocessing import Pool

import numpy as np
from scipy import ndimage as ndi

import config as C
import recalage
import spcm

N_SECONDES = 11


def etiquettes(masks, somme):
    """Somas 1..n, neuropile n + 1."""
    n = int(masks.max())
    lab = masks.astype(np.int32).copy()
    loin = ~ndi.binary_dilation(masks > 0, iterations=C.DILATATION_NEUROPILE)
    s = somme.sum(0)
    np_ = loin & (s > np.percentile(s, C.CENTILE_NEUROPILE))
    lab[np_] = n + 1
    return lab, n


def traiter(args):
    path, lab, n_lab = args
    ph = spcm.lire(path)
    return spcm.declins_rois(ph, lab, n_lab, bins_s=N_SECONDES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--processus", type=int, default=4)
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    A = [a for a in C.acquisitions(root) if (D["masks"] / f"{a['nom']}_masks.npy").exists()]
    if not A:
        sys.exit("Aucun masque — lance d'abord step2_cellpose.py")
    with Pool(args.processus) as pool:
        for a in A:
            z = np.load(D["projections"] / f"{a['nom']}.npz")
            masks = np.load(D["masks"] / f"{a['nom']}_masks.npy")
            lab, n = etiquettes(masks, z["somme"])
            d = np.nan_to_num(z["derive"])
            taches, decal = [], []
            for k, p in enumerate(a["cycles"]):
                dk = (int(round(d[k, 0])), int(round(d[k, 1])))
                decal.append(dk)
                lk = recalage.deplacer(lab, dk) if dk != (0, 0) else lab
                taches.append((p, lk.astype(np.uint16), n + 1))
            res = pool.map(traiter, taches, chunksize=1)
            dec = np.stack([r[0] for r in res]).astype(np.int32)
            sec = np.stack([r[1] for r in res]).astype(np.int32)
            np.savez_compressed(D["mesures"] / f"{a['nom']}.npz", declins=dec, secondes=sec,
                                etiquettes=lab, n_rois=n, decalages=np.array(decal),
                                garde=z["garde"], debut=z["debut"], duree_s=z["duree_s"])
            tot = dec[z["garde"]].sum((0, 3))
            print(f"  {a['nom']:>12} : {n} somas ; photons somas CFP {tot[1:n + 1, 0].sum() / 1e6:.2f} M, "
                  f"YFP {tot[1:n + 1, 1].sum() / 1e6:.2f} M ; neuropile CFP {tot[n + 1, 0] / 1e6:.2f} M")


if __name__ == "__main__":
    main()
