"""Étape 4 — ratio C2/C3 par ROI conservé, export CSV, boxplots et tracés appariés.

Ratio = moyenne(C2) / moyenne(C3) sur les pixels du ROI, valeurs brutes
(aucune soustraction de fond), calculé séparément dans chaque condition.

Usage :  python3 step4_ratio.py   (à lancer après le tri manuel de step3_select.py)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import CONDITIONS, SLICES, N_FOV, fov_key, dirs

MIN_VALID_FRAC = 0.90     # part des pixels du ROI qui doit rester dans le champ après recalage

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
TRACE = "#9a9a93"        # lignes reliant les points d'un même ROI
COLORS = {"Basal": "#2a78d6", "High_K": "#eb6834", "Return": "#1baf7a"}
LABELS = {"Basal": "Basal", "High_K": "High K$^+$", "Return": "Retour"}


# ------------------------------------------------------------------ mesures
def measure():
    D = dirs()
    rows, missing = [], []
    for slice_ in SLICES:
        for i in range(N_FOV):
            key = fov_key(slice_, i)
            sel_path = D["selection"] / f"{key}_selection.json"
            if not sel_path.exists():
                missing.append(key)
                continue
            kept = json.loads(sel_path.read_text())["kept"]
            npz = np.load(D["registered"] / f"{key}.npz")
            masks = np.load(D["masks"] / f"{key}_masks.npy")

            for lab in kept:
                roi = masks == lab
                npx = int(roi.sum())
                for cond in CONDITIONS:
                    rec = dict(slice=slice_, fov=key, roi=int(lab), condition=cond,
                               n_pixels=npx, mean_C2=np.nan, mean_C3=np.nan, ratio=np.nan)
                    if f"{cond}_C2" not in npz.files:
                        rec["note"] = "condition absente"
                        rows.append(rec); continue
                    valid = npz[f"valid_{cond}"]
                    frac = float((roi & valid).sum()) / npx
                    if frac < MIN_VALID_FRAC:
                        rec["note"] = f"hors champ après recalage ({frac:.0%} valide)"
                        rows.append(rec); continue
                    sel = roi & valid
                    c2 = float(np.nanmean(npz[f"{cond}_C2"][sel]))
                    c3 = float(np.nanmean(npz[f"{cond}_C3"][sel]))
                    rec.update(mean_C2=c2, mean_C3=c3,
                               ratio=c2 / c3 if c3 > 0 else np.nan,
                               note="" if c3 > 0 else "C3 nul")
                    rows.append(rec)
    if not rows:
        sys.exit(
            f"Aucune sélection de ROIs trouvée dans {D['selection']}\n\n"
            f"Le tri manuel (étape 3) doit être fait avant de calculer les ratios :\n"
            f"    cd {Path(__file__).resolve().parent}\n"
            f"    python3 step3_select.py\n\n"
            f"Clic = garder/rejeter, flèches = changer de FOV, q = enregistrer et quitter.")
    if missing:
        print(f"  ! {len(missing)} FOV(s) sans sélection, ignorés : {', '.join(missing)}")

    df = pd.DataFrame(rows)
    df["note"] = df["note"].fillna("")
    df.to_csv(D["results"] / "ratios_C2_C3.csv", index=False)
    return df, D


# ------------------------------------------------------------------ figures
def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID); ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, length=0)


def _box(ax, data, pos, color, width=0.5):
    bp = ax.boxplot([data], positions=[pos], widths=width, patch_artist=True,
                    showfliers=False, whis=1.5, zorder=2)
    bp["boxes"][0].set(facecolor=color + "26", edgecolor=color, linewidth=1.4)
    bp["medians"][0].set(color=color, linewidth=2.2)
    for part in ("whiskers", "caps"):
        for a in bp[part]:
            a.set(color=color, linewidth=1.2)
    rng = np.random.default_rng(0)
    x = pos + rng.uniform(-width * 0.28, width * 0.28, len(data))
    ax.scatter(x, data, s=16, color=color, alpha=0.75,
               edgecolors=SURFACE, linewidths=0.8, zorder=3)


def per_slice_plots(df, D):
    for slice_ in SLICES:
        sub = df[df["slice"] == slice_]
        fig, ax = plt.subplots(figsize=(6.6, 5.4), facecolor=SURFACE)
        _style(ax)
        ticks = []
        for p, cond in enumerate(CONDITIONS):
            vals = sub.loc[sub["condition"] == cond, "ratio"].dropna().values
            if len(vals):
                _box(ax, vals, p, COLORS[cond])
            ticks.append(f"{LABELS[cond]}\nn = {len(vals)}")
        if not any(sub["ratio"].notna()):
            ax.text(0.5, 0.5, "aucun ROI trié pour cette tranche\n"
                              "(lance step3_select.py)", transform=ax.transAxes,
                    ha="center", va="center", color=INK_2, fontsize=12)
        ax.set_xticks(range(len(CONDITIONS)))
        ax.set_xticklabels(ticks, color=INK, fontsize=11)
        ax.set_xlim(-0.6, len(CONDITIONS) - 0.4)
        ax.set_ylabel("Ratio C2 / C3", color=INK, fontsize=11)
        n_fov = sub["fov"].nunique()
        n_roi = sub[["fov", "roi"]].drop_duplicates().shape[0]
        ax.set_title(f"{slice_} — ratio C2/C3 par ROI",
                     color=INK, fontsize=13, pad=34, loc="left")
        ax.text(0, 1.012, f"{n_roi} ROIs sur {n_fov} FOVs", transform=ax.transAxes,
                color=INK_2, fontsize=9.5, va="bottom")
        fig.tight_layout()
        out = D["results"] / f"boxplot_{slice_.replace(' ', '')}.png"
        fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
        print(f"  -> {out.name}")


def combined_plot(df, D):
    fig, ax = plt.subplots(figsize=(9.4, 5.6), facecolor=SURFACE)
    _style(ax)
    w, centers = 0.26, []
    for gi, slice_ in enumerate(SLICES):
        centers.append(gi)
        for ci, cond in enumerate(CONDITIONS):
            vals = df.loc[(df["slice"] == slice_) & (df["condition"] == cond),
                          "ratio"].dropna().values
            if len(vals):
                _box(ax, vals, gi + (ci - 1) * w, COLORS[cond], width=w * 0.78)
    ax.set_xticks(centers)
    ax.set_xticklabels(SLICES, color=INK, fontsize=11.5)
    ax.set_xlim(-0.6, len(SLICES) - 0.4)
    ax.set_ylabel("Ratio C2 / C3", color=INK, fontsize=11)
    ax.set_title("Ratio C2/C3 par ROI — les 3 tranches", color=INK, fontsize=13,
                 pad=26, loc="left")
    handles = [plt.Line2D([], [], marker="s", linestyle="", markersize=9,
                          markerfacecolor=COLORS[c] + "40", markeredgecolor=COLORS[c],
                          label=LABELS[c]) for c in CONDITIONS]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.005), ncol=3,
              frameon=False, fontsize=10.5, labelcolor=INK_2, handletextpad=0.5,
              columnspacing=1.6)
    fig.tight_layout()
    out = D["results"] / "boxplot_combine.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


def _wide(sub):
    """Une ligne par ROI (identifié par FOV + numéro), une colonne par condition."""
    if sub.empty:
        return pd.DataFrame(columns=CONDITIONS)
    return (sub.pivot(index=["fov", "roi"], columns="condition", values="ratio")
               .reindex(columns=CONDITIONS))


def _traces(ax, sub, label_median=True):
    """Trace apparié : un trait mince par ROI, points colorés par condition."""
    x = np.arange(len(CONDITIONS))
    wide = _wide(sub)
    for vals in wide.to_numpy():
        # Un NaN interrompt le trait : les ROIs sans Return s'arrêtent d'eux-mêmes.
        ax.plot(x, vals, color=TRACE, linewidth=0.6, alpha=0.5, zorder=1)
    for j, cond in enumerate(CONDITIONS):
        v = wide[cond].dropna().to_numpy() if cond in wide else np.array([])
        ax.scatter(np.full(len(v), j), v, s=15, color=COLORS[cond], alpha=0.85,
                   edgecolors=SURFACE, linewidths=0.7, zorder=3)
    if len(wide):
        med = wide.median().to_numpy()
        ax.plot(x, med, color=INK, linewidth=2.4, zorder=4, marker="o", markersize=7,
                markerfacecolor=INK, markeredgecolor=SURFACE, markeredgewidth=1.6)
        last = np.flatnonzero(~np.isnan(med))
        if label_median and len(last):
            k = last[-1]
            ax.annotate("médiane", (k, med[k]), textcoords="offset points",
                        xytext=(11, 0), va="center", color=INK, fontsize=10)
    return wide


def _cond_ticks(ax, wide):
    ticks = [f"{LABELS[c]}\nn = {int(wide[c].notna().sum()) if c in wide else 0}"
             for c in CONDITIONS]
    ax.set_xticks(range(len(CONDITIONS)))
    ax.set_xticklabels(ticks, color=INK, fontsize=11)
    ax.set_xlim(-0.45, len(CONDITIONS) - 0.35)


def per_slice_traces(df, D):
    for slice_ in SLICES:
        sub = df[df["slice"] == slice_]
        fig, ax = plt.subplots(figsize=(6.8, 5.4), facecolor=SURFACE)
        _style(ax)
        wide = _traces(ax, sub)
        if not len(wide) or not wide.notna().to_numpy().any():
            ax.text(0.5, 0.5, "aucun ROI trié pour cette tranche\n"
                              "(lance step3_select.py)", transform=ax.transAxes,
                    ha="center", va="center", color=INK_2, fontsize=12)
        _cond_ticks(ax, wide)
        ax.set_ylabel("Ratio C2 / C3", color=INK, fontsize=11)
        ax.set_title(f"{slice_} — trajectoire de chaque ROI",
                     color=INK, fontsize=13, pad=34, loc="left")
        n_complet = int(wide.notna().all(axis=1).sum()) if len(wide) else 0
        ax.text(0, 1.012, f"{len(wide)} ROIs, dont {n_complet} suivis dans les 3 conditions",
                transform=ax.transAxes, color=INK_2, fontsize=9.5, va="bottom")
        fig.tight_layout()
        out = D["results"] / f"trace_{slice_.replace(' ', '')}.png"
        fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
        print(f"  -> {out.name}")


def combined_traces(df, D):
    """Petits multiples : une tranche par panneau, axe y commun pour comparer."""
    fig, axes = plt.subplots(1, len(SLICES), figsize=(14.5, 5.2), sharey=True,
                             facecolor=SURFACE)
    for k, (ax, slice_) in enumerate(zip(axes, SLICES)):
        _style(ax)
        # Les panneaux sont construits à l'identique : une seule étiquette suffit.
        wide = _traces(ax, df[df["slice"] == slice_], label_median=(k == 0))
        if not len(wide) or not wide.notna().to_numpy().any():
            ax.text(0.5, 0.5, "aucun ROI trié", transform=ax.transAxes,
                    ha="center", va="center", color=INK_2, fontsize=11)
        _cond_ticks(ax, wide)
        ax.set_title(slice_, color=INK, fontsize=12, pad=8)
    axes[0].set_ylabel("Ratio C2 / C3", color=INK, fontsize=11)
    fig.suptitle("Trajectoire de chaque ROI — les 3 tranches", color=INK, fontsize=13,
                 x=0.006, ha="left", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = D["results"] / "trace_combine.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


def main():
    df, D = measure()
    ok = df["ratio"].notna()
    print(f"\n{df[['fov','roi']].drop_duplicates().shape[0]} ROIs conservés, "
          f"{ok.sum()}/{len(df)} mesures valides")
    if excl := df.loc[~ok & (df["note"] != ""), "note"].value_counts().to_dict():
        for k, v in excl.items():
            print(f"  écarté : {k} ({v})")
    summary = (df.dropna(subset=["ratio"]).groupby(["slice", "condition"])["ratio"]
                 .agg(n="size", moyenne="mean", mediane="median", ecart_type="std",
                      q1=lambda s: s.quantile(.25), q3=lambda s: s.quantile(.75))
                 .reindex(pd.MultiIndex.from_product([SLICES, CONDITIONS])).round(4))
    summary.to_csv(D["results"] / "resume_stats.csv")
    print("\n" + summary.to_string())
    print()
    per_slice_plots(df, D)
    combined_plot(df, D)
    per_slice_traces(df, D)
    combined_traces(df, D)
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
