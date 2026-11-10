"""Étape 4 — déclins par ROI et par bloc, pour les deux canaux.

Usage :  python3 step4_declins.py <dossier>

Second passage sur les .spc : chaque photon va dans le déclin (256 canaux) de la ROI
qui contient son pixel, pour son bloc de `images_par_bloc` images. Trois familles :

- `cellpose` : tous les ROIs de l'étape 2 — le tri de l'étape 3 s'applique à l'étape 5,
  sans repasser sur les photons ;
- `flimapp` : les ROIs dessinées dans FLIMApp pendant l'expérience, pixel pour pixel
  comme FLIMApp (`roi_pixel_mask`), qui peuvent se chevaucher. Leurs photons par bloc
  sont comparés à frames.csv : ils doivent être identiques ;
- `neuropile` : les pixels hors somas (masques Cellpose dilatés de NEUROPILE_MARGE px),
  au-dessus du NEUROPILE_CENTILE e centile de l'image — on écarte les vaisseaux et les
  bords sombres. C'est le témoin « hors cellules » de la réponse.

Si la dérive (étape 1) dépasse un demi-pixel, les masques Cellpose et du neuropile sont
déplacés bloc par bloc de l'inverse de la dérive, comme dans les tranches 2026-09-21 ; un
ROI dont moins de MIN_VALID_FRAC des pixels restent dans le champ est marqué pour ce bloc.
Les ROIs FLIMApp, elles, restent fixes, comme en direct.

Les photons de chaque ROI sont aussi comptés image par image (60 ms), pour les figures
à la résolution de l'image (étape 6).

Sortie : mesures/declins.npz, qc/validation_flimapp.csv
"""
import csv
import json
import time

import numpy as np
from scipy import ndimage as ndi

import config
import spc
from config import CANAUX, CANAUX_TEMPS, MIN_VALID_FRAC

NEUROPILE_MARGE = 10        # px autour des somas exclus du neuropile
NEUROPILE_CENTILE = 25      # en dessous : vaisseaux, bords sombres


def neuropile(masks, somme):
    hors = ~ndi.binary_dilation(masks > 0, iterations=NEUROPILE_MARGE)
    return hors & (somme > np.percentile(somme, NEUROPILE_CENTILE))


def decalages(info, nb):
    """Décalage entier (dy, dx) du masque pour chaque bloc : l'inverse de la dérive,
    ou zéro partout si la dérive reste sous le demi-pixel."""
    d = np.array([(x["dy"], x["dx"]) for x in info["derive"]][:nb])
    if np.abs(d).max() < 0.5:
        return np.zeros((nb, 2), int)
    return -np.round(d).astype(int)


def mesurer(path, geo, fpb, nb, n_tr, couches, shifts):
    """`couches` : {nom: (étiquettes int32 [ny, nx], n)} — un photon compte dans
    l'étiquette de son pixel. Rend les déclins par bloc {nom: [n+1, nb, canaux]} et les
    photons image par image {nom: [n+1, n_tr]} — toutes les images, y compris celles
    qui suivent le dernier bloc complet (elles prennent le décalage de ce bloc)."""
    ny, nx = geo.lignes_image, geo.pixels_par_ligne
    groupes = {}
    for b, s in enumerate(map(tuple, shifts)):
        groupes.setdefault(s, []).append(b)
    plats = {}
    for s, blocs in groupes.items():
        sel = np.zeros(nb, bool)
        sel[blocs] = True
        # Les ROIs FLIMApp restent où elles ont été dessinées, comme en direct : elles ne
        # suivent pas la dérive (sinon elles ne reproduiraient plus frames.csv).
        plats[s] = (sel, {k: (ndi.shift(e, s, order=0, cval=0) if any(s) and not k.startswith("flimapp")
                              else e).ravel()
                          for k, (e, n) in couches.items()})
    H = {k: np.zeros((n + 1) * nb * CANAUX_TEMPS, np.int64) for k, (e, n) in couches.items()}
    P = {k: np.zeros((n + 1) * n_tr, np.int64) for k, (e, n) in couches.items()}
    for ph in spc.Rangement(path, geo, canaux=CANAUX_TEMPS):
        tr = ph["trame"]
        complet = tr // fpb < nb
        b = np.minimum(tr // fpb, nb - 1)
        pix, c = ph["y"] * nx + ph["x"], ph["canal"]
        for s, (sel, et) in plats.items():
            m = sel[b] if len(plats) > 1 else np.ones(len(b), bool)
            trm, bm, pm, cm, okm = tr[m], b[m], pix[m], c[m], complet[m]
            for k, (e, n) in couches.items():
                lab = et[k][pm]
                g = lab > 0
                P[k] += np.bincount(lab[g] * n_tr + trm[g], minlength=(n + 1) * n_tr)
                h = g & okm
                H[k] += np.bincount((lab[h] * nb + bm[h]) * CANAUX_TEMPS + cm[h],
                                    minlength=(n + 1) * nb * CANAUX_TEMPS)
    return ({k: H[k].reshape(couches[k][1] + 1, nb, CANAUX_TEMPS) for k in H},
            {k: P[k].reshape(couches[k][1] + 1, n_tr) for k in P})


def frac_valide(etiquettes, n, shifts):
    """Part des pixels de chaque ROI encore dans le champ, bloc par bloc [nb, n+1]."""
    npx = np.bincount(etiquettes.ravel(), minlength=n + 1).astype(float)
    out = np.ones((len(shifts), n + 1))
    for b, s in enumerate(map(tuple, shifts)):
        if any(s):
            c = np.bincount(ndi.shift(etiquettes, s, order=0, cval=0).ravel(), minlength=n + 1)
            with np.errstate(invalid="ignore", divide="ignore"):
                out[b] = np.where(npx > 0, c / npx, 0)
    return out


def main():
    root = config.racine()
    D = config.dirs(root)
    S = config.session(root)
    info = json.loads((D["projections"] / "session.json").read_text())
    fpb, nb = info["images_par_bloc"], info["n_blocs"]
    f = np.load(D["projections"] / "flim.npz")
    masks = np.load(D["masks"] / "cellpose_masks.npy").astype(np.int32)
    n_cp = int(masks.max())
    ny, nx = masks.shape
    np_mask = neuropile(masks, f["somme"])
    shifts = decalages(info, nb)
    print(f"{root.name} — {n_cp} ROIs Cellpose, {len(S['rois_flimapp'])} ROIs FLIMApp, "
          f"neuropile {np_mask.sum():,} px ; masques "
          f"{'déplacés selon la dérive' if np.abs(shifts).max() else 'fixes (dérive < 0,5 px)'}")

    # Les ROIs FLIMApp peuvent se chevaucher : une couche par ROI.
    polys = [config.masque_polygone(r["xs"], r["ys"], nx, ny) for r in S["rois_flimapp"]]
    couches = {"cellpose": (masks, n_cp), "neuropile": (np_mask.astype(np.int32), 1)}
    for i, m in enumerate(polys):
        couches[f"flimapp_{i + 1}"] = (m.astype(np.int32), 1)

    n_tr = len(f["trames_t"])
    H, P = {}, {}
    for k in CANAUX:
        t0 = time.time()
        geo = spc.geometrie(S["spc"][k]["ini"])
        H[k], P[k] = mesurer(S["spc"][k]["spc"], geo, fpb, nb, n_tr, couches, shifts)
        print(f"  ch{k} : {time.time() - t0:.0f} s")

    stack = lambda key: np.stack([H[k][key] for k in CANAUX], axis=2)     # [n+1, nb, 2, 256]
    cp = stack("cellpose")
    fa = np.concatenate([stack(f"flimapp_{i + 1}")[1:] for i in range(len(polys))]) \
        if polys else np.zeros((0, nb, 2, CANAUX_TEMPS), np.int64)
    npl = stack("neuropile")[1]
    par_image = lambda key: np.stack([P[k][key] for k in CANAUX], axis=2)  # [n+1, n_tr, 2]
    np.savez_compressed(
        D["mesures"] / "declins.npz",
        cellpose=cp.astype(np.int32), flimapp=fa.astype(np.int32), neuropile=npl.astype(np.int32),
        npx_cellpose=np.bincount(masks.ravel(), minlength=n_cp + 1),
        npx_flimapp=np.array([m.sum() for m in polys]), npx_neuropile=np_mask.sum(),
        frac_cellpose=frac_valide(masks, n_cp, shifts),
        masque_neuropile=np_mask, shifts=shifts,
        par_image_cellpose=par_image("cellpose").astype(np.int32),
        par_image_flimapp=(np.concatenate([par_image(f"flimapp_{i + 1}")[1:] for i in range(len(polys))])
                           if polys else np.zeros((0, n_tr, 2), np.int64)).astype(np.int32),
        par_image_neuropile=par_image("neuropile")[1].astype(np.int32))

    # --- validation contre ce que FLIMApp a mesuré en direct
    fr = config.frames_flimapp(root)
    # Une valeur non entière dans frames.csv est une moyenne glissante de FLIMApp : le
    # réglage « binning » de l'affichage, s'il est changé pendant la session, lisse aussi les
    # photons écrits (2026-10-09, FOV3 : blocs 47 à 63). Le lissage vaut pour tout le bloc
    # (une moyenne peut tomber juste sur un entier) : un bloc dont une valeur n'est pas
    # entière n'est pas comparé.
    lisses = {r["pass"] for r in fr for k in CANAUX
              if abs(float(r[f"photons_ch{k}"]) - round(float(r[f"photons_ch{k}"]))) > 1e-6}
    rows, n_diff, n_lisse = [], 0, 0
    for r in fr:
        i, b = int(r["roi_index"]) - 1, int(r["pass"]) - 1
        if i >= len(polys) or b >= nb:
            continue
        row = {"roi": i + 1, "bloc": b + 1}
        for k in CANAUX:
            moi, eux = int(fa[i, b, k - 1].sum()), float(r[f"photons_ch{k}"])
            row[f"photons_ch{k}"] = moi
            row[f"photons_ch{k}_flimapp"] = eux
            if r["pass"] in lisses:
                n_lisse += 1
                row["lisse_par_flimapp"] = True
            else:
                n_diff += moi != round(eux)
        rows.append(row)
    if rows:
        with open(D["qc"] / "validation_flimapp.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) + ["lisse_par_flimapp"], restval=False)
            w.writeheader()
            w.writerows(rows)
        print(f"  ROIs FLIMApp contre frames.csv : {len(rows)} blocs x ROIs, "
              f"{n_diff} comptes différents sur {len(rows) * len(CANAUX) - n_lisse}"
              + (f" ; {len(lisses)} blocs lissés par FLIMApp (binning de l'affichage) non comparés"
                 if n_lisse else ""))
    print(f"\nSorties -> {D['mesures'] / 'declins.npz'}")


if __name__ == "__main__":
    main()
