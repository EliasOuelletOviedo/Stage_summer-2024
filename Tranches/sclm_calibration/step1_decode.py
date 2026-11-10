"""Étape 1 — décodage des .spc : images de photons cycle par cycle, dérive, contrôles.

Usage :  python3 step1_decode.py <dossier> [--processus 4]
         python3 step1_decode.py <dossier> --sans-decodage   (refait dérive, plan focal, figures)

Chaque cycle (10 s) de chaque acquisition est décodé une fois (`spcm.lire`) : photons de
chaque canal par pixel (512 × 1024, repère des images Bliq), nombre de balayages de chaque
ligne (pour le taux de photons par impulsion, étape 6), déclin du champ entier, histogramme
brut du TDC (contrôle contre le .sdt d'SPCM).

**Plan focal.** Une retouche de la mise au point pendant une acquisition change les cellules
vues sans translation (2026-11-09 : basal/fov2, relancé, et le time-lapse basal→0 vers 15:35
et 15:52). Chaque cycle est comparé, images lissées (σ = 2 px) et recalées, à la somme des
autres cycles gardés (FOV fixes, retraits successifs) ou à ses 3 voisins de chaque côté
(transitions, où le champ évolue lentement) ; sous une ressemblance (NCC) de 0,95 il est
écarté. Deux cycles du même plan se ressemblent à ≥ 0,99 ; les retouches tombent à 0,16–0,95.

Sorties
  projections/<acq>.npz   images (cycles, 2, 512, 1024) [CFP, YFP], passages (cycles, 512),
                          debut (s, horloge du PC), duree_s, n_trames, derive (cycles, 2) en
                          pixels, ncc_plan, plan, garde, somme (2, 512, 1024) recalée, declin (cycles, 2, 256)
  qc/inventaire.csv       un cycle par ligne : heure, durée, trames, photons, plan, gardé ou non
  qc/validation_sdt.csv   dernier cycle de chaque FOV contre le .sdt écrit par SPCM
  qc/recalage_bliq.csv    ressemblance FLIM / Bliq au décalage nominal, et meilleur décalage
  qc/derive.png, qc/apercu_acquisitions.png
"""
import argparse
import csv
import datetime as dt
import re
import sys
from multiprocessing import Pool

import numpy as np
import tifffile

import config as C
import figures as G
import recalage
import spcm
from figures import plt


def traiter(path):
    ph = spcm.lire(path)
    im = spcm.image(ph)
    dec = np.stack([spcm.declin_total(ph, n) for n in C.CANAUX])
    # histogramme brut du TDC (4096 canaux de 4 ps) par entrée, pour le contrôle .sdt
    raw = np.fromfile(path, dtype="<u4")[1:]
    p = (raw & 0xC0000000) == 0
    e = (raw[p] >> 28) & 3
    tdc = (raw[p] >> 16) & 0xFFF
    h = np.stack([np.bincount(tdc[e == k], minlength=4096) for k in (0, 1)])
    return dict(im=im, passages=ph["passages"], n_trames=ph["n_trames"], duree_s=ph["duree_s"],
                n_entree=ph["n_entree"], declin=dec, tdc=h, mtime=path.stat().st_mtime)


def date_bliq(f):
    with tifffile.TiffFile(f) as t:
        d = t.pages[0].tags["ImageDescription"].value
    s = re.search(r"<AcquisitionDate>([^<]+)<", d).group(1)
    return dt.datetime.fromisoformat(s).timestamp()


def lire_bliq(f):
    return np.squeeze(tifffile.imread(f)).astype(float)


def controle_bliq(a, res, debut):
    """NCC entre l'image Bliq la plus proche dans le temps du premier cycle et l'image FLIM
    de ce cycle, au décalage nominal, puis meilleur décalage entier dans ±6 px."""
    out = {}
    if not a["bliq"]["C1"]:
        return out
    tb = np.array([date_bliq(f) for f in a["bliq"]["C1"]])
    k_bliq = int(np.argmin(np.abs(tb - debut[0])))
    k_flim = int(np.argmin(np.abs(debut - tb[k_bliq])))
    for c, i, nom in (("C1", 1, "yfp"), ("C2", 0, "cfp")):
        b = lire_bliq(a["bliq"][c][k_bliq])
        f = res[k_flim]["im"][i].astype(float)
        m = np.zeros(b.shape, bool)
        m[10:-10, 10:-10] = True
        best = max((recalage.ncc(np.roll(f, (dy, dx), (0, 1)), b, m), dy, dx)
                   for dy in range(-6, 7) for dx in range(-6, 7))
        out[nom] = dict(bliq=c, ncc_nominal=recalage.ncc(f, b, m), ncc_meilleur=best[0],
                        dy=best[1], dx=best[2], ecart_s=float(tb[k_bliq] - debut[k_flim]))
    return out


SEUIL_PLAN = 0.95          # NCC (images lissées σ = 2 px) sous lequel un cycle est d'un autre plan
VOISINS_PLAN = 3


def deriver(tot, iterations=2):
    """Décalage de chaque cycle (images CFP + YFP) sur la somme recalée (position médiane = 0)."""
    d = np.zeros((len(tot), 2))
    ref = tot.sum(0)
    for _ in range(iterations):
        d = np.array([recalage.decalage(t, ref) for t in tot])
        d -= np.median(d, axis=0)
        ref = sum(recalage.deplacer(t, -d[k]) for k, t in enumerate(tot))
    return d


def plan_focal(reg, transition):
    """NCC de chaque cycle (images recalées, lissées) : aux 3 voisins de chaque côté pour une
    transition ; à la somme des autres cycles gardés pour un FOV, en retirant un à un les
    cycles sous le seuil. Rend (ncc, même_plan)."""
    from scipy import ndimage as ndi
    lis = np.stack([ndi.gaussian_filter(r, 2) for r in reg])
    m = np.zeros(lis.shape[1:], bool)
    m[20:-20, 20:-20] = True
    n = len(lis)
    if transition:
        nc = np.array([recalage.ncc(lis[k], lis[[j for j in range(max(0, k - VOISINS_PLAN), min(n, k + VOISINS_PLAN + 1))
                                                   if j != k]].sum(0), m) for k in range(n)])
        return nc, nc >= SEUIL_PLAN
    ok = np.ones(n, bool)
    while True:
        nc = np.array([recalage.ncc(lis[k], lis[ok & (np.arange(n) != k)].sum(0), m) for k in range(n)])
        cand = np.nonzero(ok & (nc < SEUIL_PLAN))[0]
        if not len(cand) or ok.sum() <= 2:
            return nc, ok
        ok[cand[np.argmin(nc[cand])]] = False


def decoder_tout(root, D, A, processus):
    taches = [(i, k, p) for i, a in enumerate(A) for k, p in enumerate(a["cycles"])]
    print(f"{len(A)} acquisitions, {len(taches)} cycles ; {processus} processus")
    sdt_rows = []
    par_acq = {i: [None] * len(a["cycles"]) for i, a in enumerate(A)}
    restants = {i: len(a["cycles"]) for i, a in enumerate(A)}
    with Pool(processus) as pool:
        for (i, k, p), r in zip(taches, pool.imap(traiter, [t[2] for t in taches], chunksize=1)):
            par_acq[i][k] = r
            restants[i] -= 1
            if restants[i]:
                continue
            a, res = A[i], par_acq.pop(i)
            # contrôle .sdt : dernier cycle d'un FOV ; bloc b = entrée b (IN1, IN2)
            if a["sdt"] is not None:
                import sdtfile
                s = sdtfile.SdtFile(a["sdt"])
                cyc = [k2 for k2, p2 in enumerate(a["cycles"]) if p2.stem == a["sdt"].stem]
                if cyc:
                    r = res[cyc[0]]
                    for b, blk in enumerate(s.data[:2]):
                        dsdt = blk.reshape(-1, blk.shape[-1]).sum(0).astype(float)
                        tot = r["tdc"][b].sum()
                        mien = r["tdc"][b].reshape(256, 16).sum(1).astype(float)
                        rel = np.max(np.abs(dsdt / dsdt.sum() - mien / mien.sum())) / (mien / mien.sum()).max()
                        sdt_rows.append(dict(acquisition=a["nom"], fichier=a["sdt"].name, bloc=b,
                                             entree=f"IN{b + 1}", photons_sdt=int(dsdt.sum()),
                                             photons_spc=int(tot), rapport=round(dsdt.sum() / tot, 5),
                                             ecart_max_forme=round(float(rel), 5)))
            np.savez_compressed(D["projections"] / f"{a['nom']}.npz",
                                images=np.stack([r["im"] for r in res]),
                                passages=np.stack([r["passages"] for r in res]),
                                debut=np.array([r["mtime"] - r["duree_s"] for r in res]),
                                duree_s=np.array([r["duree_s"] for r in res]),
                                n_trames=np.array([r["n_trames"] for r in res]),
                                n_entree=np.stack([r["n_entree"] for r in res]),
                                declin=np.stack([r["declin"] for r in res]))
            print(f"  {a['nom']:>12} : {len(res)} cycles décodés")
    ecrire(D, "validation_sdt.csv", sdt_rows)


def finaliser(root, D, a):
    """Dérive, plan focal, somme recalée, contrôle Bliq : complète projections/<acq>.npz."""
    p = D["projections"] / f"{a['nom']}.npz"
    z = dict(np.load(p))
    ims = z["images"]
    duree_ok = z["duree_s"] >= C.MIN_DUREE_CYCLE_S
    tot = ims.sum(1).astype(np.float32)
    d = np.full((len(ims), 2), np.nan)
    d[duree_ok] = deriver(tot[duree_ok])
    reg = np.stack([recalage.deplacer(t, -np.nan_to_num(d[k])) for k, t in enumerate(tot)])
    nc = np.full(len(ims), np.nan)
    plan = np.zeros(len(ims), bool)
    nc[duree_ok], plan[duree_ok] = plan_focal(reg[duree_ok], a["type"] == "transition")
    garde = duree_ok & plan
    if garde.sum() and garde.sum() < duree_ok.sum():          # dérive refaite sur les seuls cycles gardés
        d[garde] = deriver(tot[garde])
        d[duree_ok & ~garde] = [recalage.decalage(tot[k], reg[garde].sum(0)) for k in np.nonzero(duree_ok & ~garde)[0]]
    somme = np.zeros((2,) + ims.shape[2:], np.float32)
    for k in np.nonzero(garde)[0]:
        for c in range(2):
            somme[c] += recalage.deplacer(ims[k, c].astype(np.float32), -d[k])
    z.update(derive=d, ncc_plan=nc, plan=plan, garde=garde, somme=somme)
    np.savez_compressed(p, **z)
    bl = controle_bliq(a, [dict(im=im) for im in ims], z["debut"])
    inv = []
    for k, p2 in enumerate(a["cycles"]):
        n_img = ims[k].reshape(2, -1).sum(1)
        inv.append(dict(acquisition=a["nom"], type=a["type"], condition=a["condition"], cycle=k,
                        fichier=str(p2.relative_to(root / "data" / "FLIM")),
                        debut=dt.datetime.fromtimestamp(z["debut"][k]).strftime("%H:%M:%S"),
                        duree_s=round(float(z["duree_s"][k]), 3), trames=int(z["n_trames"][k]),
                        photons_cfp=int(n_img[0]), photons_yfp=int(n_img[1]),
                        yfp_sur_cfp=round(n_img[1] / max(n_img[0], 1), 4),
                        MHz_cfp=round(z["n_entree"][k][C.ENTREE["cfp"]] / z["duree_s"][k] / 1e6, 4),
                        MHz_yfp=round(z["n_entree"][k][C.ENTREE["yfp"]] / z["duree_s"][k] / 1e6, 4),
                        derive_dy=round(float(d[k, 0]), 2), derive_dx=round(float(d[k, 1]), 2),
                        ncc_plan=round(float(nc[k]), 4), meme_plan=bool(plan[k]), garde=bool(garde[k])))
    dmax = np.nanmax(np.hypot(d[garde, 0], d[garde, 1])) if garde.any() else np.nan
    ecartes = [k for k in range(len(ims)) if not garde[k]]
    print(f"  {a['nom']:>12} : {garde.sum()}/{len(ims)} cycles gardés"
          + (f" (écartés : {ecartes})" if ecartes else "") + f", dérive max {dmax:.2f} px"
          + (f", Bliq C1 {bl['yfp']['ncc_nominal']:.2f} (meilleur {bl['yfp']['dy']:+d},{bl['yfp']['dx']:+d})" if bl else ""))
    return inv, [dict(acquisition=a["nom"], canal=nom, **{k: (round(x, 4) if isinstance(x, float) else x)
                                                         for k, x in v.items()}) for nom, v in bl.items()]


def ecrire(D, nom, rows):
    if rows:
        with open(D["qc"] / nom, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"  -> qc/{nom}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--processus", type=int, default=4)
    ap.add_argument("--sans-decodage", action="store_true")
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    A = C.acquisitions(root)
    if not args.sans_decodage:
        decoder_tout(root, D, A, args.processus)
    inventaire, bliq_rows = [], []
    for a in A:
        inv, bl = finaliser(root, D, a)
        inventaire += inv
        bliq_rows += bl
    ecrire(D, "inventaire.csv", inventaire)
    ecrire(D, "recalage_bliq.csv", bliq_rows)
    figures(root, D, A)


def figures(root, D, A):
    # dérive des transitions
    trans = [a for a in A if a["type"] == "transition"]
    fig, axes = plt.subplots(len(trans), 1, figsize=(10, 2.6 * len(trans)), squeeze=False)
    for ax, a in zip(axes[:, 0], trans):
        z = np.load(D["projections"] / f"{a['nom']}.npz")
        t = (z["debut"] - z["debut"][0]) / 60
        G.style(ax)
        ax.plot(t, z["derive"][:, 0], "o-", ms=3, color=G.PALETTE[0], label="dy (lignes)")
        ax.plot(t, z["derive"][:, 1], "o-", ms=3, color=G.PALETTE[1], label="dx (pixels)")
        ax2 = ax.twinx()
        ax2.plot(t, z["ncc_plan"], color=G.TRACE, lw=1)
        ko = ~z["garde"]
        ax2.plot(t[ko], z["ncc_plan"][ko], "x", color=G.PALETTE[4], ms=6, mew=1.5)
        ax2.axhline(SEUIL_PLAN, color=G.PALETTE[4], lw=0.8, ls=":")
        ax2.set_ylabel("NCC plan focal", color=G.INK_2, fontsize=9)
        ax2.tick_params(colors=G.INK_2, length=0)
        for s in ax2.spines.values():
            s.set_visible(False)
        ax.set_ylabel("décalage (px)", color=G.INK, fontsize=10)
        ax.set_title(f"{a['nom']} : dérive cycle par cycle (gris : ressemblance aux voisins ; × : cycle écarté)",
                     color=G.INK, fontsize=11, loc="left")
        ax.legend(frameon=False, fontsize=9, loc="upper left")
    axes[-1, 0].set_xlabel("minutes depuis le début du time-lapse", color=G.INK)
    fig.tight_layout()
    G.sauver(fig, D["qc"] / "derive.png")

    # aperçu de toutes les acquisitions
    n = len(A)
    nc = 5
    nr = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 2.4 * nr), facecolor=G.SURFACE)
    for ax in axes.flat:
        ax.set_axis_off()
    for ax, a in zip(axes.flat, A):
        z = np.load(D["projections"] / f"{a['nom']}.npz")
        s = np.sqrt(z["somme"].sum(0))
        ax.imshow(s, cmap="gray", **G.contraste(s))
        tot = z["images"][z["garde"]].reshape(z["garde"].sum(), 2, -1).sum(2).sum(0)
        ax.set_title(f"{a['nom']}  —  YFP/CFP {tot[1] / tot[0]:.2f}", fontsize=9, color=G.INK)
    fig.suptitle("Somme des photons (CFP + YFP, racine) de chaque acquisition, cycles recalés",
                 color=G.INK, fontsize=12)
    fig.tight_layout()
    G.sauver(fig, D["qc"] / "apercu_acquisitions.png", dpi=110)


if __name__ == "__main__":
    main()
