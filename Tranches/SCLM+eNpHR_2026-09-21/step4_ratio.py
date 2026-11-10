"""Étape 4 — ratio FRET C2/C3 dans le temps, par ROI et par acquisition.

ratio = moyenne(C2) / moyenne(C3) sur les pixels du ROI, **valeurs brutes**
(aucune soustraction de fond), exactement comme dans `../sclm/`, mais calculé
image par image le long du décours.

La stimulation illumine tout le champ et fuit dans C2 (jusqu'à x25 le niveau de base) :
le ratio brut y est ininterprétable. Ces images restent dans le CSV avec `stim=True`,
mais le trait s'interrompt dans les figures et elles ne comptent pas dans R0.

Les colonnes `fond_C2` / `fond_C3` (5e centile de l'image) et `ratio_fond_soustrait`
sont des diagnostics : hors stim le fond vaut 0 et le ratio corrigé est identique au
ratio brut ; pendant la stim la correction ne **rétablit pas** la valeur de base
(54 -> 18 alors que la base est 5 pour stim_max), car la fuite n'est pas un simple
décalage uniforme. Ne pas s'en servir comme d'une mesure du FRET pendant la stim.

Le ratio n'est calculé que sur les acquisitions d'au moins `MIN_FRAMES_TIMELAPSE` images.
Les images fixes, elles, ne sont mesurées que sur C2 et qu'à 800 et 950 nm : C2 moyen
par ROI, sans ratio (`images_fixes_C2.csv`).

Usage :  python3 step4_ratio.py   (après le tri manuel de step3_select.py)
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config
from config import (CHANNELS, FIXED_CHANNEL, MIN_FRAMES_TIMELAPSE, MIN_VALID_FRAC,
                    PRE_WINDOW_S, RESPONSE_WINDOW_S, dirs, frame_files, slug, survey)
from step1_align import apply_shift

BG_PERCENTILE = 5      # estimation du fond : 5e centile de l'image
BASE_SKIP = 2          # premières images ignorées dans la ligne de base (transitoire PMT)
SMOOTH = 5             # fenêtre du lissage (images) — figures seulement, jamais le CSV

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
TRACE = "#b4b4ac"
STIM_BAND = "#f4d35e"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#9b51e0", "#d64550", "#00867d"]


# ------------------------------------------------------------------ mesures
def _groups(shifts):
    """Indices d'images regroupés par décalage identique (le masque n'est déplacé qu'une fois)."""
    out = {}
    for i, (dy, dx) in enumerate(shifts):
        out.setdefault((int(dy), int(dx)), []).append(i)
    return out


def measure_recording(fov, rec, info, labels, n_lab, channels=CHANNELS):
    """Moyennes par ROI et par image des `channels`. Tableaux [n_images, n_lab+1]."""
    n = info["n_frames"]
    shifts = np.asarray(info["drift"]) + np.asarray(info["fov_shift"])
    files = {ch: frame_files(fov, rec, ch) for ch in channels}

    npx = np.bincount(labels.ravel(), minlength=n_lab + 1).astype(float)
    mean = {ch: np.full((n, n_lab + 1), np.nan) for ch in channels}
    bg = {ch: np.full(n, np.nan) for ch in channels}
    frac = np.zeros((n, n_lab + 1))

    for (dy, dx), idx in _groups(shifts).items():
        # Le ROI est défini dans le repère du FOV ; l'image brute, elle, n'est pas
        # recalée : on déplace le masque de l'inverse plutôt que d'interpoler l'image.
        lab_raw = apply_shift(labels.astype(np.float64), -dy, -dx, fill=0.0).astype(np.int32)
        flat = lab_raw.ravel()
        counts = np.bincount(flat, minlength=n_lab + 1).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(npx > 0, counts / npx, 0.0)
        for i in idx:
            frac[i] = f
            for ch in channels:
                im = tifffile.imread(files[ch][i]).squeeze().astype(np.float64)
                s = np.bincount(flat, weights=im.ravel(), minlength=n_lab + 1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    mean[ch][i] = np.where(counts > 0, s / counts, np.nan)
                bg[ch][i] = np.percentile(im[::2, ::2], BG_PERCENTILE)
    return mean, bg, frac, npx


def measure():
    D = dirs()
    data = survey()
    rows, missing = [], []
    for fov, recs in data.items():
        sel_path = D["selection"] / f"{fov}_selection.json"
        align_path = D["projections"] / f"{fov}_align.json"
        if not sel_path.exists() or not align_path.exists():
            missing.append(fov)
            continue
        kept = json.loads(sel_path.read_text())["kept"]
        align = json.loads(align_path.read_text())
        masks = np.load(D["masks"] / f"{fov}_masks.npy")
        n_lab = int(masks.max())
        labels = np.where(np.isin(masks, kept), masks, 0).astype(np.int32)

        for rec, n_frames, date in recs:
            if n_frames < MIN_FRAMES_TIMELAPSE:
                continue          # image fixe : affichée par l'étape 1, pas mesurée
            info = align["info"][rec]
            mean, bg, frac, npx = measure_recording(fov, rec, info, labels, n_lab)
            times = np.asarray(info["times"])
            stim = np.asarray(info["stim"], bool)
            for lab in kept:
                c2, c3 = mean["C2"][:, lab], mean["C3"][:, lab]
                ok = frac[:, lab] >= MIN_VALID_FRAC
                with np.errstate(invalid="ignore", divide="ignore"):
                    ratio = np.where(ok & (c3 > 0), c2 / c3, np.nan)
                    num, den = c2 - bg["C2"], c3 - bg["C3"]
                    ratio_bg = np.where(ok & (den > 0), num / den, np.nan)
                for i in range(n_frames):
                    note = ("hors champ après recalage" if not ok[i]
                            else "stim : C2 contaminé" if stim[i] else "")
                    rows.append(dict(
                        fov=fov, acquisition=rec, date=date, roi=int(lab),
                        n_pixels=int(npx[lab]), image=i + 1, t_s=round(float(times[i]), 4),
                        mean_C2=c2[i], mean_C3=c3[i], ratio=ratio[i],
                        fond_C2=bg["C2"][i], fond_C3=bg["C3"][i],
                        ratio_fond_soustrait=ratio_bg[i],
                        stim=bool(stim[i]), frac_valide=round(float(frac[i, lab]), 4),
                        note=note))
            print(f"    {fov}/{rec:16s} {n_frames:4d} images x {len(kept)} ROIs")

    if not rows:
        sys.exit(
            f"Aucune sélection de ROIs trouvée dans {D['selection']}\n\n"
            f"Le tri manuel (étape 3) doit être fait avant de calculer les ratios :\n"
            f"    cd {Path(__file__).resolve().parent}\n"
            f"    python3 step3_select.py\n\n"
            f"Clic = garder/rejeter, flèches = changer de FOV, q = enregistrer et quitter.")
    if missing:
        print(f"  ! FOV(s) sans sélection, ignorés : {', '.join(missing)}")

    df = pd.DataFrame(rows)
    df.to_csv(D["results"] / "ratio_temps.csv", index=False)
    return df, D


def measure_fixed(D):
    """Images fixes : C2 moyen par ROI conservé, à 800 et 950 nm seulement."""
    col = f"mean_{FIXED_CHANNEL}"
    rows = []
    for fov, recs in survey().items():
        sel_path = D["selection"] / f"{fov}_selection.json"
        align_path = D["projections"] / f"{fov}_align.json"
        if not sel_path.exists() or not align_path.exists():
            continue
        kept = json.loads(sel_path.read_text())["kept"]
        info = json.loads(align_path.read_text())["info"]
        masks = np.load(D["masks"] / f"{fov}_masks.npy")
        n_lab = int(masks.max())
        labels = np.where(np.isin(masks, kept), masks, 0).astype(np.int32)
        for rec in config.fixed_recordings(recs):
            mean, _, frac, npx = measure_recording(fov, rec, info[rec], labels, n_lab,
                                                   channels=(FIXED_CHANNEL,))
            for lab in kept:
                ok = frac[0, lab] >= MIN_VALID_FRAC
                rows.append({"fov": fov, "acquisition": rec,
                             "longueur_onde_nm": config.wavelength(rec), "roi": int(lab),
                             "n_pixels": int(npx[lab]),
                             col: mean[FIXED_CHANNEL][0, lab] if ok else np.nan,
                             "frac_valide": round(float(frac[0, lab]), 4),
                             "note": "" if ok else "hors champ après recalage"})
    fx = pd.DataFrame(rows)
    if len(fx):
        fx.to_csv(D["results"] / f"images_fixes_{FIXED_CHANNEL}.csv", index=False)
    return fx


# ------------------------------------------------------------------ dérivés
def add_normalised(df):
    """R0 = médiane du ratio brut avant la stim ; ajoute dR/R0 en % et la phase."""
    by_rec = ["fov", "acquisition"]
    lim = df[df["stim"]].groupby(by_rec)["image"].agg(["min", "max"])
    df = df.merge(lim.rename(columns={"min": "stim_i0", "max": "stim_i1"}),
                  on=by_rec, how="left")

    avant = df["stim_i0"].isna() | (df["image"] < df["stim_i0"])
    base = df["ratio"].where(avant & (df["image"] > BASE_SKIP))
    df["R0"] = base.groupby([df[k] for k in by_rec + ["roi"]]).transform("median")
    with np.errstate(invalid="ignore", divide="ignore"):
        df["dR_R0_pct"] = np.where(df["R0"] > 0, (df["ratio"] / df["R0"] - 1) * 100, np.nan)

    df["phase"] = np.select(
        [df["stim"], df["stim_i1"].notna() & (df["image"] > df["stim_i1"])],
        ["pendant", "apres"], default="avant")
    return df.drop(columns=["stim_i0", "stim_i1"])


def immediate_response(df):
    """dR/R0 moyen sur les premières `RESPONSE_WINDOW_S` secondes après la stim.

    La moyenne sur toute la période d'après suit la dérive lente de la ligne de base ;
    les premières secondes la suivent bien moins, et c'est là que la réponse est la
    plus franche.
    """
    ap = df[df["phase"] == "apres"]
    t0 = ap.groupby(["fov", "acquisition"])["t_s"].transform("min")
    tot = ap[ap["t_s"] <= t0 + RESPONSE_WINDOW_S]
    return tot.groupby(["fov", "acquisition", "roi"])["dR_R0_pct"].mean()


def pre_stim_slope(df):
    """Pente de la ligne de base avant la stim, en % de R0 par 10 s.

    Diagnostic : une pente non nulle décale mécaniquement la moyenne d'après, puisque
    R0 est pris au milieu d'une base qui dérive déjà.
    """
    def one(g):
        a = g[g["phase"] == "avant"]
        if a["t_s"].nunique() < 5:
            return np.nan
        m = a.groupby("t_s")["dR_R0_pct"].median()
        return float(np.polyfit(m.index.to_numpy(), m.to_numpy(), 1)[0] * 10)
    s = df.groupby(["fov", "acquisition"]).apply(one)
    return s.reindex(df.groupby(["fov", "acquisition", "phase"]).size()
                       .index.droplevel("phase")).to_numpy()


# ------------------------------------------------------------------ figures
def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID); ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, length=0)


def _matrix(sub, col):
    """[n_ROIs, n_images] du champ `col`, plus l'axe des temps."""
    w = sub.pivot(index="roi", columns="image", values=col).sort_index()
    t = sub.drop_duplicates("image").set_index("image")["t_s"].sort_index()
    return w.to_numpy(float), t.to_numpy(float), w.index.to_numpy()


def _stim_blocks(sub):
    """Blocs contigus d'images de stim, en secondes : [(début, fin), ...].
    Un seul bloc pour une stim unique, plusieurs pour un train."""
    im = (sub.loc[sub["stim"], ["image", "t_s"]].drop_duplicates("image")
             .sort_values("image"))
    if im.empty:
        return []
    cut = np.flatnonzero(np.diff(im["image"].to_numpy()) > 1) + 1
    return [(p.min(), p.max()) for p in np.split(im["t_s"].to_numpy(), cut)]


def _stim_band(ax, sub):
    for a, b in _stim_blocks(sub):
        ax.axvspan(a, b, color=STIM_BAND, alpha=0.30, lw=0, zorder=0)


def _smooth(M, w=SMOOTH):
    """Médiane glissante le long du temps (axe 1) — pour l'affichage seulement.

    `rolling` travaille sur l'index, donc sur M transposée : les lignes de M sont
    des ROIs, il ne faut surtout pas mélanger des ROIs voisins.
    """
    if w < 2 or M.shape[1] < w:
        return M
    return (pd.DataFrame(M.T).rolling(w, center=True, min_periods=max(2, w // 2))
              .median().to_numpy().T)


def _robust_ylim(ax, M, pad=0.12):
    v = M[np.isfinite(M)]
    if v.size < 10:
        return
    lo, hi = np.percentile(v, (0.5, 99.5))
    m = (hi - lo) * pad or 1.0
    ax.set_ylim(lo - m, hi + m)


def _traces(ax, sub, col, color, smooth=True):
    """Un trait mince par ROI + médiane et intervalle interquartile. Stim masquée."""
    M, t, _ = _matrix(sub, col)
    M = np.where(_matrix(sub, "stim")[0].astype(bool), np.nan, M)
    for row in M:
        ax.plot(t, row, color=TRACE, linewidth=0.5, alpha=0.5, zorder=1)
    if not M.size or not np.isfinite(M).any():
        return M, t
    S = _smooth(M) if smooth else M
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # colonnes entièrement masquées
        q1, med, q3 = (np.nanpercentile(S, q, axis=0) for q in (25, 50, 75))
    ax.fill_between(t, q1, q3, color=color, alpha=0.22, lw=0, zorder=2)
    ax.plot(t, med, color=color, linewidth=2.2, zorder=3)
    _robust_ylim(ax, M)
    return M, t


def trace_figure(sub, title, label, color, out, extra=None):
    """Ratio brut en haut, dR/R0 en bas : un trait par ROI, médiane et écart interquartile
    lissés, bande de stim masquée. `sub` doit avoir une colonne `roi` unique par trait ;
    `extra(fig, axes)` permet d'ajouter des éléments avant l'enregistrement."""
    fig, axes = plt.subplots(2, 1, figsize=(9.6, 7.0), sharex=True, facecolor=SURFACE)
    for ax in axes:
        _style(ax)
        _stim_band(ax, sub)
    _traces(axes[0], sub, "ratio", color)
    _traces(axes[1], sub, "dR_R0_pct", color)
    axes[1].axhline(0, color=INK_2, linewidth=0.9, linestyle="--", zorder=2)

    axes[0].set_ylabel("Ratio C2 / C3", color=INK, fontsize=11)
    axes[1].set_ylabel("$\\Delta$R / R$_0$  (%)", color=INK, fontsize=11)
    axes[1].set_xlabel("temps (s)", color=INK, fontsize=11)
    axes[0].set_title(title, color=INK, fontsize=13, pad=40, loc="left")
    blocks = _stim_blocks(sub)
    if len(blocks) == 1:
        bande = f"bande jaune = stimulation {blocks[0][0]:.0f}–{blocks[0][1]:.0f} s"
    else:
        bande = f"bandes jaunes = {len(blocks)} stimulations"
    note = (f"{label}, valeurs brutes   ·   médiane et écart interquartile "
            f"lissés sur {SMOOTH} images\n{bande}, images masquées : C2 contaminé"
            if blocks else f"{label}   ·   aucune stimulation détectée")
    axes[0].text(0, 1.015, note, transform=axes[0].transAxes,
                 color=INK_2, fontsize=9.5, va="bottom", linespacing=1.5)
    if extra:
        extra(fig, axes)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


def recording_figure(sub, fov, rec, color, D):
    trace_figure(sub, f"{fov} — {rec} : ratio FRET dans le temps",
                 f"{sub['roi'].nunique()} ROIs", color,
                 D["results"] / f"trace_{fov}_{slug(rec)}.png")


def heatmap_figure(sub, fov, rec, D):
    M, t, rois = _matrix(sub, "dR_R0_pct")
    M = _smooth(np.where(_matrix(sub, "stim")[0].astype(bool), np.nan, M))
    s = sub.loc[sub["stim"], "t_s"]
    cols = t > s.max() if len(s) else np.ones_like(t, bool)   # tri sur l'après-stim
    if not cols.any():
        cols = np.ones_like(t, bool)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        score = np.nanmean(np.where(cols, M, np.nan), axis=1)
    order = np.argsort(-np.nan_to_num(score, nan=-1e9))
    M, rois = M[order], rois[order]

    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        moy = np.nanmean(M, axis=0)          # moyenne des ROIs, mêmes valeurs lissées

    lim = np.nanpercentile(np.abs(M), 98) if np.isfinite(M).any() else 1.0
    h = 0.115 * len(rois) + 0.6
    fig, (top, ax) = plt.subplots(
        2, 1, sharex=True, figsize=(9.6, h + 1.9), facecolor=SURFACE,
        gridspec_kw=dict(height_ratios=[0.34, h], hspace=0.02))

    kw = dict(aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")

    # La bande du haut a sa propre échelle, étirée sur son min et son max, pour
    # occuper toute la gamme de couleurs : ses amplitudes sont bien plus petites que
    # celles des ROIs pris un à un. Elle ne se lit donc pas sur la barre de couleur,
    # et le blanc n'y est plus le zéro — ses bornes sont rappelées à gauche.
    if np.isfinite(moy).any():
        mlo, mhi = float(np.nanmin(moy)), float(np.nanmax(moy))
    else:
        mlo, mhi = -lim, lim
    if mhi - mlo < 1e-9:
        mlo, mhi = mlo - 1.0, mhi + 1.0
    top.imshow(moy[None, :], extent=(t[0], t[-1], 0.5, -0.5),
               **{**kw, "vmin": mlo, "vmax": mhi})
    top.set_yticks([0])
    top.set_yticklabels([f"moyenne\n{mlo:+.1f} … {mhi:+.1f} %"], fontsize=7.5, color=INK)
    top.tick_params(colors=INK_2, length=0)

    im = ax.imshow(M, extent=(t[0], t[-1], len(rois) - 0.5, -0.5), **kw)
    if len(s):
        for a in (top, ax):
            for x in (s.min(), s.max()):
                a.axvline(x, color=INK, linewidth=1.0, linestyle="--", zorder=4)
    ax.set_yticks(range(len(rois)))
    ax.set_yticklabels(rois, fontsize=6, color=INK_2)
    ax.set_xlabel("temps (s)", color=INK, fontsize=11)
    ax.set_ylabel("ROI (trié par réponse)", color=INK, fontsize=10)
    top.set_title(f"{fov} — {rec} : $\\Delta$R/R$_0$ par ROI, lissé sur {SMOOTH} images  "
                  f"(blanc = stim, C2 contaminé)", color=INK, fontsize=12, pad=10, loc="left")
    ax.tick_params(colors=INK_2, length=0)
    cb = fig.colorbar(im, ax=(top, ax), pad=0.015, fraction=0.03)
    cb.set_label("$\\Delta$R / R$_0$  (%)", color=INK_2, fontsize=9.5)
    cb.ax.tick_params(colors=INK_2, length=0)
    out = D["results"] / f"heatmap_{fov}_{slug(rec)}.png"
    fig.savefig(out, dpi=190, facecolor=SURFACE, bbox_inches="tight"); plt.close(fig)
    print(f"  -> {out.name}")


def fov_figure(df, fov, recs, D):
    """Petits multiples : une acquisition par panneau, axe y commun."""
    sub_all = df[(df["fov"] == fov) & df["acquisition"].isin(recs)]
    fig, axes = plt.subplots(1, len(recs), figsize=(5.0 * len(recs), 4.9),
                             sharey=True, facecolor=SURFACE, squeeze=False)
    for k, (ax, rec) in enumerate(zip(axes[0], recs)):
        _style(ax)
        sub = sub_all[sub_all["acquisition"] == rec]
        _stim_band(ax, sub)
        _traces(ax, sub, "dR_R0_pct", PALETTE[k % len(PALETTE)])
        ax.axhline(0, color=INK_2, linewidth=0.9, linestyle="--", zorder=2)
        ax.set_title(rec, color=INK, fontsize=12, pad=8)
        ax.set_xlabel("temps (s)", color=INK, fontsize=10.5)
    lims = [ax.get_ylim() for ax in axes[0]]
    axes[0][0].set_ylim(min(a for a, _ in lims), max(b for _, b in lims))
    axes[0][0].set_ylabel("$\\Delta$R / R$_0$  (%)", color=INK, fontsize=11)
    fig.suptitle(f"{fov} — réponse FRET de chaque ROI, toutes les acquisitions",
                 color=INK, fontsize=13, x=0.006, ha="left", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = D["results"] / f"trace_{fov}_toutes.png"
    fig.savefig(out, dpi=190, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


def _paired(ax, w, labels):
    """Un trait par ROI d'une colonne de `w` à l'autre, points colorés, médiane en gras."""
    _style(ax)
    x = np.arange(w.shape[1])
    # Au-delà de quelques dizaines de traits, une opacité fixe ferait un aplat gris.
    alpha = min(0.6, max(0.06, 25 / max(len(w), 1)))
    for row in w.to_numpy():
        ax.plot(x, row, color=TRACE, linewidth=0.6, alpha=alpha, zorder=1)
    for j, c in enumerate(w.columns):
        v = w[c].dropna().to_numpy()
        ax.scatter(np.full(len(v), j), v, s=12, color=PALETTE[j % len(PALETTE)],
                   alpha=0.85, edgecolors=SURFACE, linewidths=0.6, zorder=3)
    ax.plot(x, w.median().to_numpy(), color=INK, linewidth=2.2, marker="o",
            markersize=6, markerfacecolor=INK, markeredgecolor=SURFACE, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=INK, fontsize=10)
    ax.set_xlim(-0.4, len(x) - 0.6)


def fixed_figure(fx, D):
    """C2 moyen par ROI aux deux longueurs d'onde : un trait par ROI, médiane en gras."""
    col = f"mean_{FIXED_CHANNEL}"
    fovs = list(dict.fromkeys(fx["fov"]))
    wl = sorted(fx["longueur_onde_nm"].unique())
    fig, axes = plt.subplots(1, len(fovs), figsize=(2.2 * len(fovs) + 1.4, 4.9),
                             facecolor=SURFACE, squeeze=False)
    for ax, fov in zip(axes[0], fovs):
        w = (fx[fx["fov"] == fov]
             .pivot(index="roi", columns="longueur_onde_nm", values=col).reindex(columns=wl))
        _paired(ax, w, [f"{l} nm" for l in wl])
        ax.set_title(f"{fov}\nn = {len(w)}", color=INK, fontsize=11, pad=6)
    axes[0][0].set_ylabel(f"{FIXED_CHANNEL} moyen par ROI (counts)", color=INK, fontsize=11)
    fig.suptitle(f"Images fixes — {FIXED_CHANNEL} moyen de chaque ROI à "
                 f"{' et '.join(f'{l} nm' for l in wl)}   (trait noir = médiane ; "
                 f"échelle propre à chaque FOV)", color=INK, fontsize=12.5,
                 x=0.006, ha="left", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = D["results"] / f"images_fixes_{FIXED_CHANNEL}.png"
    fig.savefig(out, dpi=190, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


def pre_post(df):
    """dR/R0 moyen par ROI et par stim : `PRE_WINDOW_S` avant le début de la stim, et
    `RESPONSE_WINDOW_S` après sa fin, les images de stim exclues des deux côtés.

    Le pourcentage de changement rend les ROIs comparables d'un FOV à l'autre, ce que
    le ratio brut ne permet pas (sa valeur de base va de 1,5 à 5 selon la tranche) : on
    peut donc regrouper toutes les stims VIS d'un côté, toutes les IR de l'autre. Les
    stims sans étiquette (FOV1/stim-00x, FOV2/stim_max) n'entrent dans aucun groupe.
    """
    rows = []
    for (fov, rec), g in df.groupby(["fov", "acquisition"], sort=False):
        typ = config.stim_type(rec)
        if typ is None or not g["stim"].any():
            continue
        t0, t1 = g.loc[g["stim"], "t_s"].min(), g.loc[g["stim"], "t_s"].max()
        av = g[(g["phase"] == "avant") & (g["t_s"] >= t0 - PRE_WINDOW_S)]
        ap = g[(g["phase"] == "apres") & (g["t_s"] <= t1 + RESPONSE_WINDOW_S)]
        w = pd.DataFrame({"dR_R0_avant_pct": av.groupby("roi")["dR_R0_pct"].mean(),
                          "dR_R0_apres_pct": ap.groupby("roi")["dR_R0_pct"].mean()})
        rows.append(w.reset_index().assign(
            fov=fov, acquisition=rec, type=typ,
            n_images_avant=av["image"].nunique(), n_images_apres=ap["image"].nunique()))
    pp = pd.concat(rows, ignore_index=True)
    return pp[["type", "fov", "acquisition", "roi", "dR_R0_avant_pct", "dR_R0_apres_pct",
               "n_images_avant", "n_images_apres"]]


def pre_post_figure(pp, D):
    """Toutes les stims VIS regroupées, puis toutes les IR : un trait par ROI et par stim."""
    types = [t for t in ("VIS", "IR") if (pp["type"] == t).any()]
    labels = [f"{PRE_WINDOW_S:.0f} s\navant", f"{RESPONSE_WINDOW_S:.0f} s\naprès"]
    cols = ["dR_R0_avant_pct", "dR_R0_apres_pct"]
    fig, axes = plt.subplots(1, len(types), figsize=(3.4 * len(types) + 1.6, 5.4),
                             sharey=True, facecolor=SURFACE, squeeze=False)
    for ax, t in zip(axes[0], types):
        sub = pp[pp["type"] == t]
        w = sub[cols].dropna().reset_index(drop=True)
        _paired(ax, w, labels)
        ax.axhline(0, color=INK_2, linewidth=0.9, linestyle="--", zorder=2)
        med = w.median()
        ax.annotate(f"{med.iloc[1]:+.1f} %", (1, med.iloc[1]), textcoords="offset points",
                    xytext=(10, 0), va="center", color=INK, fontsize=10.5, zorder=5)
        n_stim = sub[["fov", "acquisition"]].drop_duplicates().shape[0]
        ax.set_title(f"{t}\n{n_stim} stims, {sub['fov'].nunique()} FOVs, n = {len(w)}",
                     color=INK, fontsize=11, pad=6)
    _robust_ylim(axes[0][0], pp[cols].to_numpy())
    axes[0][0].set_ylabel("$\\Delta$R / R$_0$ moyen  (%)", color=INK, fontsize=11)
    fig.suptitle(f"Changement du ratio C2/C3, {PRE_WINDOW_S:.0f} s avant la stim et "
                 f"{RESPONSE_WINDOW_S:.0f} s après\nun trait par ROI et par stim ; "
                 f"trait noir = médiane", color=INK, fontsize=12, x=0.01, ha="left",
                 y=0.985, linespacing=1.4)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = D["results"] / "avant_apres_stim.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


# ------------------------------------------------------------------ main
def main():
    df, D = measure()
    df = add_normalised(df)
    df.to_csv(D["results"] / "ratio_temps.csv", index=False)

    ok = df["ratio"].notna()
    print(f"\n{df[['fov','roi']].drop_duplicates().shape[0]} ROIs conservés, "
          f"{ok.sum()}/{len(df)} mesures valides")
    for k, v in df.loc[~ok & (df['note'] != ''), 'note'].value_counts().items():
        print(f"  écarté : {k} ({v})")

    # --- tableaux par ROI et par acquisition ---
    piv = df.pivot_table(index=["fov", "acquisition", "roi"], columns="phase",
                         values=["ratio", "dR_R0_pct"], aggfunc="mean")
    per_roi = pd.DataFrame({
        "R0": df.groupby(["fov", "acquisition", "roi"])["R0"].first(),
        "ratio_avant": piv.get(("ratio", "avant")),
        "ratio_apres": piv.get(("ratio", "apres")),
        "dR_R0_apres": piv.get(("dR_R0_pct", "apres")),
        "dR_R0_apres_immediat": immediate_response(df),
    }).reset_index()
    per_roi.to_csv(D["results"] / "resume_par_roi.csv", index=False)

    summary = (df[df["ratio"].notna()]
               .groupby(["fov", "acquisition", "phase"])
               .agg(n_mesures=("ratio", "size"), n_rois=("roi", "nunique"),
                    ratio_moyen=("ratio", "mean"), ratio_median=("ratio", "median"),
                    ratio_ecart_type=("ratio", "std"),
                    dR_R0_median=("dR_R0_pct", "median")).round(4))
    summary["pente_avant_pct_par_10s"] = pre_stim_slope(df).round(3)
    summary.to_csv(D["results"] / "resume_acquisitions.csv")
    print("\n" + summary.to_string())
    print()

    for fov, recs_info in survey().items():
        recs = [r for r, n, _ in recs_info if n >= MIN_FRAMES_TIMELAPSE
                and ((df["fov"] == fov) & (df["acquisition"] == r)).any()]
        for k, rec in enumerate(recs):
            sub = df[(df["fov"] == fov) & (df["acquisition"] == rec)]
            color = PALETTE[k % len(PALETTE)]
            recording_figure(sub, fov, rec, color, D)
            heatmap_figure(sub, fov, rec, D)
        if recs:
            fov_figure(df, fov, recs, D)

    pp = pre_post(df)
    pp.to_csv(D["results"] / "avant_apres_stim.csv", index=False)
    print(f"\ndR/R0 moyen par ROI, {PRE_WINDOW_S:.0f} s avant / {RESPONSE_WINDOW_S:.0f} s "
          f"après la stim (médianes ; un ROI compte une fois par stim) :")
    ok = pp.dropna(subset=["dR_R0_avant_pct", "dR_R0_apres_pct"])
    for t, g in ok.groupby("type", sort=False):
        d = g["dR_R0_apres_pct"] - g["dR_R0_avant_pct"]
        print(f"  {t:3s}  {g['dR_R0_avant_pct'].median():+.2f} % -> "
              f"{g['dR_R0_apres_pct'].median():+.2f} %   "
              f"({int((d < 0).sum())}/{len(d)} en baisse)")
        if t == "VIS":
            pw = g["acquisition"].str.extract(r"VIS_(\d+)", expand=False).fillna("?")
            for k, gg in g.groupby(pw, sort=False):
                print(f"         VIS {k:>3s} % : {gg['dR_R0_apres_pct'].median():+.2f} %  "
                      f"(n = {len(gg)})")
    pre_post_figure(pp, D)

    fx = measure_fixed(D)
    if len(fx):
        col = f"mean_{FIXED_CHANNEL}"
        med = fx.pivot_table(index="fov", columns="longueur_onde_nm", values=col,
                             aggfunc="median")
        print(f"\nImages fixes — {FIXED_CHANNEL} médian par ROI (counts) :\n"
              + med.round(2).to_string())
        fixed_figure(fx, D)

    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
