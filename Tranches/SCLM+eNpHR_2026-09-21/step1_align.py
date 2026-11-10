"""Étape 1 — fenêtre de stim, recalage, projections moyennes et image somme par FOV.

Pour chaque acquisition de `Data/<FOV>/` :
  1. la fenêtre de stimulation est détectée sur la moyenne pleine image de C2
     (la lumière de stim fuit dans ce canal et sature le signal FRET) ;
  2. la dérive de platine est estimée par blocs de `BLOCK` images hors stim
     (une image isolée est trop bruitée pour un pic de corrélation fiable),
     puis interpolée sur toute l'acquisition ;
  3. la projection moyenne hors stim (C2+C3) est recalée sur l'acquisition de
     référence du FOV — les ROIs sont ainsi communs à toutes les acquisitions.

Les acquisitions à une seule image n'ont pas de décours. Elles servent au recalage et
entrent dans l'image somme ; seul leur C2 à 800 et 950 nm est affiché, dans
`<FOV>_images_fixes.png`, et mesuré par l'étape 4.

Sorties : `projections/<FOV>.npz` (projections recalées + somme), `<FOV>_sum.tif/.png`,
`<FOV>_projections.png`, `<FOV>_images_fixes.png`, et `qc/alignement.tsv` +
`projections/<FOV>_align.json` (temps, fenêtre de stim et décalage de chaque image,
relus par l'étape 4).
"""
import json

import numpy as np
import tifffile
from scipy.ndimage import binary_dilation, gaussian_filter, label
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from config import (BLOCK, CHANNELS, MAX_SHIFT, MIN_FRAMES_TIMELAPSE, MIN_NCC,
                    STIM_EDGE_W, STIM_MIN_FRAMES, STIM_MIN_RISE, STIM_PAD,
                    STIM_STEP_K, dirs, frame_files, slug, survey)


# ------------------------------------------------------------------ entrées/sorties
def load_stack(fov, rec, ch):
    """Toutes les images d'un canal, en uint8 (300 Mo max par acquisition)."""
    return np.stack([tifffile.imread(p).squeeze() for p in frame_files(fov, rec, ch)])


def read_times(fov, rec):
    """Axe des temps (s) : DeltaT de l'en-tête OME, identique pour C2 et C3."""
    return np.array([config.delta_t(p) for p in frame_files(fov, rec, "C2")])


# ------------------------------------------------------------------ stimulation
def edge_response(m, w=STIM_EDGE_W):
    """Contraste de bord : moyenne des w images suivantes moins celle des w précédentes.

    Une marche franche y répond de toute sa hauteur, une dérive lente de sa pente fois
    w seulement, et le bruit est divisé par sqrt(w/2) par rapport à une différence
    entre images voisines. Indéfini sur les w bords du décours.
    """
    n = len(m)
    g = np.full(n, np.nan)
    if n < 2 * w:
        return g
    c = np.concatenate([[0.0], np.cumsum(m)])
    i = np.arange(w, n - w + 1)
    g[i] = (c[i + w] - c[i]) / w - (c[i] - c[i - w]) / w
    return g


def stim_window(c2):
    """Images contaminées par la lumière de stim, repérées par les bords du plateau.

    La stim illumine tout le champ : elle se voit comme un plateau rectangulaire sur
    la moyenne pleine image de C2. On la cherche par ses **bords** — un allumage et
    une extinction francs — et non par son niveau : le bruit des contrastes de bord
    n'est pas gonflé par le plateau, alors qu'un MAD sur les niveaux l'est d'autant
    plus que la stim est longue et forte. Le photoblanchiment et la dérive, eux, sont
    progressifs : aucun bord ne ressort.

    Les deux bords trouvés, la fenêtre est délimitée par un seuil à mi-chemin entre le
    niveau du plateau et celui du reste, ramenée au plus long bloc contigu, puis
    élargie de `STIM_PAD` image de chaque côté — l'obturateur s'ouvre et se ferme en
    milieu d'image, qui est alors à moitié éclairée et tombe pile sur le seuil.

    Retourne (masque, moyenne pleine image, diagnostic).
    """
    m = c2.reshape(len(c2), -1).mean(axis=1)
    none = np.zeros(len(m), bool)
    diag = {"bord_sur_bruit": 0.0, "niveau": None, "hausse": None}
    if len(m) < max(4 * STIM_MIN_FRAMES, 4 * STIM_EDGE_W):
        return none, m, diag

    g = edge_response(m)
    ok = np.isfinite(g)
    scale = 1.4826 * np.median(np.abs(g[ok] - np.median(g[ok])))
    # g[i] compare m[i-w:i] et m[i:i+w] : un bord montant culmine sur la première
    # image éclairée, un bord descendant sur la première image revenue à la base.
    i_up = int(np.nanargmax(np.where(ok, g, -np.inf)))
    i_dn = int(np.nanargmin(np.where(ok, g, np.inf)))
    diag["bord_sur_bruit"] = float(min(g[i_up], -g[i_dn]) / max(scale, 1e-12))
    if not (g[i_up] > STIM_STEP_K * scale and -g[i_dn] > STIM_STEP_K * scale
            and i_dn - i_up >= STIM_MIN_FRAMES):
        return none, m, diag

    # Niveaux comparés hors images de transition (à moitié éclairées).
    inside = np.zeros(len(m), bool); inside[i_up + 1:i_dn - 1] = True
    outside = np.ones(len(m), bool); outside[max(i_up - 1, 0):i_dn + 1] = False
    if not inside.any() or not outside.any():
        return none, m, diag
    hi, lo = float(np.median(m[inside])), float(np.median(m[outside]))
    diag["hausse"] = hi / lo - 1 if lo > 0 else None
    if hi <= lo or diag["hausse"] < STIM_MIN_RISE:
        return none, m, diag

    level = (hi + lo) / 2
    diag["niveau"] = level
    hit = m > level
    lab, n = label(hit)
    if n > 1:
        sizes = [int((lab == j).sum()) for j in range(1, n + 1)]
        hit = lab == 1 + int(np.argmax(sizes))
    if hit.sum() < STIM_MIN_FRAMES:
        return none, m, diag
    return binary_dilation(hit, np.ones(2 * STIM_PAD + 1, bool)), m, diag


# ------------------------------------------------------------------ recalage
def _prep(im):
    """Passe-bande + z-score : rend comparables des images d'intensités très différentes."""
    im = gaussian_filter(im, 1.5) - gaussian_filter(im, 20.0)
    return (im - im.mean()) / (im.std() + 1e-9)


def estimate_shift(moving, reference):
    """Translation entière (dy, dx) qui amène `moving` sur `reference`, + score NCC."""
    a, b = _prep(reference), _prep(moving)
    cc = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real / a.size)
    c = np.array(cc.shape) // 2
    w = cc[c[0] - MAX_SHIFT:c[0] + MAX_SHIFT + 1, c[1] - MAX_SHIFT:c[1] + MAX_SHIFT + 1]
    i = np.unravel_index(np.argmax(w), w.shape)
    return (int(i[0] - MAX_SHIFT), int(i[1] - MAX_SHIFT)), float(w.max())


def apply_shift(im, dy, dx, fill=np.nan):
    """Translation entière sans interpolation ; les bords entrants valent `fill`."""
    out = np.full_like(im, fill, dtype=np.float64)
    ys, xs = im.shape
    y0d, y1d = max(dy, 0), min(ys + dy, ys)
    y0s, y1s = max(-dy, 0), min(ys - dy, ys)
    x0d, x1d = max(dx, 0), min(xs + dx, xs)
    x0s, x1s = max(-dx, 0), min(xs - dx, xs)
    out[y0d:y1d, x0d:x1d] = im[y0s:y1s, x0s:x1s]
    return out


def drift_shifts(total, stim, idx_ok, log_rows, tag):
    """Dérive intra-acquisition : décalage de chaque image vers le premier bloc hors stim.

    Estimée par blocs de `BLOCK` images consécutives hors stim, puis interpolée
    linéairement sur l'axe du temps. Les images de stim n'entrent pas dans
    l'estimation : leur fond uniforme fausserait la corrélation.
    """
    n = len(total)
    if len(idx_ok) < BLOCK:
        return np.zeros((n, 2), int)

    blocks = [idx_ok[i:i + BLOCK] for i in range(0, len(idx_ok), BLOCK)]
    if len(blocks) > 1 and len(blocks[-1]) < BLOCK // 2:
        # Dernier bloc trop court : fusionné avec le précédent. En deux temps — écrit
        # `blocks[-2] = concat(blocks[-2], blocks.pop())`, l'index -2 serait évalué
        # après le pop et viserait le mauvais bloc.
        last = blocks.pop()
        blocks[-1] = np.concatenate([blocks[-1], last])

    ref = total[blocks[0]].mean(axis=0)
    centers, shifts = [], []
    for b in blocks:
        (dy, dx), ncc = estimate_shift(total[b].mean(axis=0), ref)
        if ncc < MIN_NCC:
            dy = dx = 0
            note = "NON RECALE (NCC faible)"
        else:
            note = "ok"
        centers.append(float(np.mean(b)))
        shifts.append((dy, dx))
        log_rows.append(f"{tag}\tbloc {b[0]+1}-{b[-1]+1}\tdy={dy:+d}\tdx={dx:+d}\t"
                        f"NCC={ncc:.3f}\t{note}")

    x = np.arange(n, dtype=float)
    s = np.array(shifts, float)
    return np.round(np.column_stack([np.interp(x, centers, s[:, 0]),
                                     np.interp(x, centers, s[:, 1])])).astype(int)


# ------------------------------------------------------------------ figures
def _show(ax, img, title):
    img = np.nan_to_num(img, nan=0.0)
    lo, hi = np.percentile(img, (1, 99.7))
    ax.imshow(img, cmap="gray", vmin=0, vmax=60)
    ax.set_title(title, fontsize=10)
    ax.axis("off")


def fixed_images_figure(fov, recs, D):
    """Images fixes retenues (800 / 950 nm) : C2 seul, tel qu'acquis."""
    fixed = config.fixed_recordings(recs)
    if not fixed:
        return
    chs = ["C2", "C3"]
    fig, axes = plt.subplots(1, len(fixed), figsize=(4.4 * len(fixed), 4.9), squeeze=False)
    for ax, rec in zip(axes[0], chs):
        print(rec)
        im = tifffile.imread(frame_files(fov, "800-001", rec)[0]).squeeze().astype(np.float64)
        _show(ax, im, f"{"950-001"} — {rec}")
    fig.suptitle(f"{fov} — images fixes, {"800-001"} seul (telles qu'acquises)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(D["projections"] / f"{fov}_images_fixes.png", dpi=120)
    plt.close(fig)


def figures(fov, proj, total, D):
    lo, hi = np.percentile(total, (1, 99.7))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(total, cmap="gray", vmin=lo, vmax=hi)
    ax.set_title(f"{fov} — somme des projections C2+C3 ({len(proj)} acquisitions)", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(D["projections"] / f"{fov}_sum.png", dpi=130)
    plt.close(fig)

    ncol = min(4, len(proj))
    nrow = int(np.ceil(len(proj) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 4.3 * nrow), squeeze=False)
    for ax, (rec, im) in zip(axes.ravel(), proj.items()):
        _show(ax, im, rec)
    for ax in axes.ravel()[len(proj):]:
        ax.axis("off")
    fig.suptitle(f"{fov} — projections moyennes hors stim, recalées", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(D["projections"] / f"{fov}_projections.png", dpi=120)
    plt.close(fig)


# ------------------------------------------------------------------ traitement
def process_recording(fov, rec, n_frames, D, log):
    """Retourne (projection non recalée sur le FOV, infos par image)."""
    tag = f"{fov}/{rec}"
    stacks = {ch: load_stack(fov, rec, ch) for ch in CHANNELS}
    times = read_times(fov, rec)
    total = stacks["C2"].astype(np.float32) + stacks["C3"]

    if n_frames >= MIN_FRAMES_TIMELAPSE:
        stim, c2_mean, diag = stim_window(stacks["C2"])
        idx_ok = np.flatnonzero(~stim)
        shifts = drift_shifts(total, stim, idx_ok, log, tag)
    else:
        stim = np.zeros(n_frames, bool)
        c2_mean = stacks["C2"].reshape(n_frames, -1).mean(axis=1)
        diag = {}
        idx_ok = np.arange(n_frames)
        shifts = np.zeros((n_frames, 2), int)

    # Projection de référence de l'acquisition : moyenne des images hors stim,
    # chacune ramenée dans le repère du premier bloc.
    acc = np.zeros(total.shape[1:], np.float64)
    for i in idx_ok:
        dy, dx = shifts[i]
        acc += total[i] if (dy == 0 and dx == 0) else apply_shift(total[i], dy, dx, fill=0.0)
    proj = acc / len(idx_ok)

    info = dict(n_frames=int(n_frames),
                times=[round(float(t), 4) for t in times],
                t_index=[config._t_index(p) for p in frame_files(fov, rec, "C2")],
                stim=[bool(s) for s in stim],
                c2_full_frame_mean=[round(float(v), 4) for v in c2_mean],
                stim_detection=diag,
                drift=[[int(a), int(b)] for a, b in shifts])
    if stim.any():
        j = np.flatnonzero(stim)
        info["stim_frames"] = [int(j[0]) + 1, int(j[-1]) + 1]
        info["stim_times"] = [round(float(times[j[0]]), 2), round(float(times[j[-1]]), 2)]
        log.append(f"{tag}\tstim\timages {j[0]+1}-{j[-1]+1}\t"
                   f"t={times[j[0]]:.1f}-{times[j[-1]]:.1f} s\t"
                   f"C2 x{c2_mean[j].mean()/np.median(c2_mean[~stim]):.1f}\t"
                   f"bord/bruit={diag['bord_sur_bruit']:.0f}")
    else:
        info["stim_frames"] = None
        why = (f"bord/bruit={diag['bord_sur_bruit']:.1f} < {STIM_STEP_K:g}"
               if diag.get("bord_sur_bruit") is not None else "")
        log.append(f"{tag}\tstim\taucune\t\t\t{why}")
    return proj, info


def process_fov(fov, recs, D, log):
    ref_rec = config.reference_recording(recs)
    raw, infos = {}, {}
    for rec, n, date in recs:
        raw[rec], infos[rec] = process_recording(fov, rec, n, D, log)
        infos[rec]["acquisition_date"] = date
        print(f"    {rec:16s} {n:4d} images   "
              f"stim {infos[rec]['stim_frames'] or '--'}")

    # Recalage inter-acquisitions sur la référence du FOV.
    proj, valid = {}, {}
    for rec in raw:
        if rec == ref_rec:
            (dy, dx), ncc, note = (0, 0), 1.0, "ref"
        else:
            (dy, dx), ncc = estimate_shift(raw[rec], raw[ref_rec])
            if ncc < MIN_NCC:
                dy = dx = 0
                note = "NON RECALE (NCC faible)"
            else:
                note = "ok"
        proj[rec] = apply_shift(raw[rec], dy, dx)
        valid[rec] = apply_shift(np.ones_like(raw[rec]), dy, dx, fill=0.0) > 0.5
        infos[rec]["fov_shift"] = [dy, dx]
        infos[rec]["fov_ncc"] = round(ncc, 4)
        log.append(f"{fov}/{rec}\tFOV -> {ref_rec}\tdy={dy:+d}\tdx={dx:+d}\tNCC={ncc:.3f}\t{note}")

    total = sum(np.nan_to_num(p, nan=0.0) for p in proj.values())
    np.savez_compressed(
        D["projections"] / f"{fov}.npz",
        **{f"proj_{slug(r)}": proj[r].astype(np.float32) for r in proj},
        **{f"valid_{slug(r)}": valid[r] for r in valid},
        summed=total.astype(np.float32))
    tifffile.imwrite(D["projections"] / f"{fov}_sum.tif", total.astype(np.float32))
    (D["projections"] / f"{fov}_align.json").write_text(json.dumps(
        {"fov": fov, "reference": ref_rec,
         "recordings": [r for r, _, _ in recs], "info": infos}, indent=1))
    figures(fov, proj, total, D)
    fixed_images_figure(fov, recs, D)


def main():
    D = dirs()
    data = survey()
    if not data:
        raise SystemExit(f"Aucune acquisition exploitable dans {config.DATA_ROOT}")
    log = ["acquisition\tetape\tdy / detail\tdx\tNCC\tstatut"]
    for fov, recs in data.items():
        print(f"  {fov} ({len(recs)} acquisitions, référence {config.reference_recording(recs)})")
        process_fov(fov, recs, D, log)
    (D["qc"] / "alignement.tsv").write_text("\n".join(log) + "\n")
    print(f"\n{len(data)} FOVs -> {D['projections']}")


if __name__ == "__main__":
    main()
