"""Slide des résultats préliminaires — SCLm + eNpHR en tranches organotypiques.

Produit présentation/slide_resultats_FR.png et présentation/slide_resultats_EN.png, leurs
versions sur fond noir (_FR_dark.png, _EN_dark.png),
et FOV11_train/results/tau_par_cellule.csv.

Tout porte sur FOV11 :
  A  images à 800 nm et à 950 nm côte à côte, canal C2 (YFP) : SCLm remplit
     les corps cellulaires, eNpHR marque les membranes ;
  B  premier train : ligne de base, première stimulation et retour ;
  C  premier train complet : quatre stimulations ;
  D  temps de retour tau de chaque cellule, pour les 8 stimulations des deux trains.

Choix convenus :
  - ratio YFP/CFP = C2/C3 ; une baisse du ratio = entrée de Cl- ;
  - cellules retenues : au moins 10 % plus brillantes que le neuropile qui les entoure
    (MIN_CONTRASTE) — les cellules écartées répondent moins, ce qui confirme qu'elles
    expriment peu le capteur ;
  - ratio brut (les sauts instrumentaux de CFP du second train tombent hors des
    fenêtres de mesure ; voir le README) ;
  - tau : pour les N_TAU cellules à la réponse médiane la plus forte sur les 8
    stimulations des deux trains, ajustement base + A·exp(-t/tau) du retour après
    chaque stimulation (base = niveau des 10 s avant) ; les tau > TAU_MAX sont rejetés
    (retour trop lent pour être mesuré dans la fenêtre).

Usage :  python3 slide_resultats.py   (après train_FOV11.py)
"""
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("SCLM_FOVS", "FOV11")      # pour lister les images de FOV11

import json                       # noqa: E402
import warnings                   # noqa: E402
import textwrap                   # noqa: E402
import numpy as np                # noqa: E402
import pandas as pd               # noqa: E402
import tifffile                   # noqa: E402
import matplotlib                 # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
from matplotlib.colors import LinearSegmentedColormap   # noqa: E402
from matplotlib.ticker import FuncFormatter      # noqa: E402
from scipy.ndimage import binary_dilation, label as ndi_label   # noqa: E402
from scipy.optimize import curve_fit             # noqa: E402

import config                                    # noqa: E402
from step1_align import apply_shift              # noqa: E402

OUT = HERE / "présentation"
TR_ROOT, FOV = HERE / "FOV11_train", "FOV11"
TRAINS = ("train_stim-002", "train_stim-003")
FIXED = (("800-001", "C2"), ("950-001", "C2"))
MIN_CONTRASTE = 1.1
N_TAU = 10            # meilleurs ROIs dont on montre le tau à chaque stimulation
ZOOM_ROI = 6          # ROI agrandi dans un encart en bas à droite des deux images de A
ZOOM_COTE = 150       # côté du carré agrandi, en px, centré sur le ROI
# Fausses couleurs des images de A : noir -> cyan (800 nm), noir -> jaune (950 nm).
CMAPS_A = (LinearSegmentedColormap.from_list("noir_cyan", ["#000000", "#00ffff"]),
           LinearSegmentedColormap.from_list("noir_jaune", ["#000000", "#ffff00"]))
B_TRAIN, B_PULSE = TRAINS[1], 2   # panneau B : ligne de base puis cette stimulation
B_AVANT = 0.0                      # s montrées avant le début de la stimulation
# Remontée de B ajustée par base + A·exp(−t/τ), t depuis la fin de la stimulation, sur la
# médiane non lissée des cellules. La base est fixée à la ligne de base affichée (0), comme
# pour les τ de D (base fixée) : sur cette fenêtre de 49 s, une base libre absorbe le
# retour (τ = 46 ± 33 s, contre 23,8 ± 2,9 s base fixée). None : base libre.
B_FIT_BASE = 0.0
TAU_MAX = 60          # tau au-delà rejetés : retour trop lent pour la fenêtre (42 à 53 s)
WIN, SMOOTH = 10.0, 5

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
ACCENT, ACCENT_2, STIM = "#2a78d6", "#eb6834", "#f4d35e"

# Deux thèmes : chaque slide est produite en clair (slide_resultats_FR.png) et en sombre,
# sur fond noir (slide_resultats_FR_dark.png). Les teintes sombres sont celles de la
# palette de référence pour un fond sombre, validées sur #000000 (contraste ≥ 3:1).
THEMES = {
    "": dict(SURFACE=SURFACE, INK=INK, INK_2=INK_2, GRID=GRID, ACCENT=ACCENT,
             ACCENT_2=ACCENT_2, STIM=STIM, STIM_ALPHA=0.45),
    "_dark": dict(SURFACE="#000000", INK="#ffffff", INK_2="#c3c2b7", GRID="#2f2f2c",
                  ACCENT="#3987e5", ACCENT_2="#d95926", STIM="#c9a227", STIM_ALPHA=0.35),
}
STIM_ALPHA = 0.45


def appliquer_theme(suffixe):
    globals().update(THEMES[suffixe])

TEXT = {
    "FR": dict(
        title="SCLm + eNpHR en tranches organotypiques : résultats préliminaires",
        takeaway="L'activation d'eNpHR fait baisser le ratio YFP/CFP — une entrée de Cl⁻ — à chaque "
                 "stimulation, avec un retour vers la base en ~{tau:.0f} s",
        A="Expression des deux constructions", A_sub="même champ, canal YFP",
        A1="800 nm : SCLm", A2="950 nm : eNpHR",
        B="Une stimulation, puis le retour", B_sub="essai 2 : ligne de base, puis stimulation 2 — "
                                                  "médiane ± écart interquartile, n = {n} cellules",
        C="Quatre stimulations, quatre réponses",
        C_sub="médiane ± écart interquartile, n = {n} cellules",
        D="Temps de retour τ", D_sub="2 essais, un point par cellule et par essai\n"
                                    "barre = médiane ± écart interquartile",
        stim_k="stim {k}",
        time="temps (s)", y="ΔR/R₀ YFP/CFP  (%)", stim="stim", tau="τ  (s)",
        trial="essai {i}", tick="stim {k}\n{d} s",
        foot=("Ratio YFP/CFP (canaux C2/C3) ; ΔR/R₀ par rapport à la ligne de base avant la première "
              "stimulation ; courbes lissées sur 5 images. Cellules retenues : au moins 10 % plus brillantes "
              "que le neuropile qui les entoure ({n}/{N}). Images acquises pendant la stimulation masquées "
              "(fuite de la lumière de stimulation dans le canal YFP).\n"
              "A : une image par longueur d'onde, contraste propre à chaque image. B–C : essai 1, quatre "
              "stimulations de 8 à 10 s. D : essais 1 et 2 (essai 2 : 6, 10, 15 et 20 s) ; τ ajusté par "
              "A·exp(−t/τ) + C sur le retour de chaque cellule jusqu'à la stimulation suivante, retenu si "
              "la réponse (10 s après − 10 s avant) dépasse 3 fois le bruit de la cellule, R² ≥ 0,3 et τ "
              "plus court que la fenêtre d'observation (42 à 53 s)."),
        scale="50 µm"),
    "EN": dict(
        title="SCLm + eNpHR in organotypic slices: preliminary results",
        takeaway="eNpHR activation lowers the YFP/CFP ratio — a Cl⁻ influx — at every stimulation, "
                 "recovering towards baseline in ~{tau:.0f} s",
        A="Expression of both constructs", A_sub="same field, YFP channel",
        A1="800 nm: SCLm", A2="950 nm: eNpHR",
        B="One stimulation, then recovery", B_sub="trial 2: baseline, then stimulation 2 — "
                                                 "median ± interquartile range, n = {n} cells",
        C="Four stimulations, four responses",
        C_sub="median ± interquartile range, n = {n} cells",
        D="Recovery time τ", D_sub="2 trials, one dot per cell and trial\n"
                                  "bar = median ± interquartile range",
        stim_k="stim {k}",
        time="time (s)", y="ΔR/R₀ YFP/CFP  (%)", stim="stim", tau="τ  (s)",
        trial="trial {i}", tick="stim {k}\n{d} s",
        foot=("YFP/CFP ratio (channels C2/C3); ΔR/R₀ relative to the baseline before the first "
              "stimulation; curves smoothed over 5 frames. Cells kept: at least 10 % brighter than the "
              "surrounding neuropil ({n}/{N}). Frames acquired during stimulation are masked "
              "(stimulation light leaks into the YFP channel).\n"
              "A: one image per wavelength, contrast set per image. B–C: trial 1, four 8 to 10 s "
              "stimulations. D: trials 1 and 2 (trial 2: 6, 10, 15 and 20 s); τ fitted with "
              "A·exp(−t/τ) + C on each cell's recovery until the next stimulation, kept if the response "
              "(10 s after − 10 s before) exceeds 3 times the cell's noise, R² ≥ 0.3 and τ is shorter "
              "than the observation window (42 to 53 s)."),
        scale="50 µm"),
}


# ------------------------------------------------------------------ données
def kept_rois():
    """ROIs au moins MIN_CONTRASTE fois plus brillants que l'anneau de neuropile autour."""
    img = np.load(TR_ROOT / "projections" / f"{FOV}.npz")["summed"].astype(float)
    m = np.load(TR_ROOT / "masks" / f"{FOV}_masks.npy")
    others = binary_dilation(m > 0, iterations=2)
    keep = []
    for lab in range(1, m.max() + 1):
        roi = m == lab
        ring = binary_dilation(roi, iterations=12) & ~others
        if ring.any() and img[roi].mean() / img[ring].mean() >= MIN_CONTRASTE:
            keep.append(lab)
    return keep, int(m.max())


def fixed_images():
    """Images fixes de FOV11, recalées sur la référence du FOV."""
    info = json.loads((TR_ROOT / "projections" / f"{FOV}_align.json").read_text())["info"]
    out = []
    for rec, ch in FIXED:
        im = tifffile.imread(config.frame_files(FOV, rec, ch)[0]).squeeze().astype(float)
        dy, dx = info[rec]["fov_shift"]
        out.append(np.nan_to_num(apply_shift(im, dy, dx), nan=0.0))
    return out


def smooth(y, w=SMOOTH):
    """Moyenne glissante centrée ; les images masquées (stim) restent vides."""
    return pd.Series(y).rolling(w, center=True, min_periods=2).mean().to_numpy()


def matrix(g, col):
    w = g.pivot(index="roi", columns="image", values=col).sort_index()
    st = g.drop_duplicates("image").set_index("image").sort_index()
    M = np.where(st["stim"].to_numpy(bool)[None, :], np.nan, w.to_numpy(float))
    return M, st["t_s"].to_numpy(float), st["stim"].to_numpy(bool)


def _model(t, A, tau, C):
    return A * np.exp(-t / tau) + C


def recovery_taus(d):
    """tau de retour de chaque cellule après chaque stimulation des deux trains.

    Ajustement y = base + A·exp(-t/tau) sur le retour, de la fin de la stimulation au
    début de la suivante, la base étant fixée au niveau moyen des 10 s avant la
    stimulation : la réponse revient vers là où elle était. Avec un plateau libre,
    celui-ci absorbe le retour et tau part à la borne pour un tiers des cellules.
    """
    rows = []
    for rec, g in d.groupby("acquisition"):
        stim_t = g[g["stim"]].groupby("pulse")["t_s"].agg(["min", "max"])
        dt = float(np.median(np.diff(np.sort(g["t_s"].unique()))))
        n_stim = g[g["stim"]].groupby("pulse")["image"].nunique()
        for roi, c in g.groupby("roi"):
            for k, (t0, t1) in stim_t.iterrows():
                nxt = stim_t.loc[k + 1, "min"] if k + 1 in stim_t.index else np.inf
                seg = c[~c["stim"] & (c["t_s"] > t1) & (c["t_s"] < nxt)]
                base = c[~c["stim"] & (c["t_s"] < t0) & (c["t_s"] >= t0 - WIN)]["dR_R0_pct"].mean()
                resp = seg[seg["t_s"] <= t1 + WIN]["dR_R0_pct"].mean() - base
                t, y = seg["t_s"].to_numpy() - t1, seg["dR_R0_pct"].to_numpy()
                tau = r2 = np.nan
                try:
                    f = lambda tt, A, tau_: base + A * np.exp(-tt / tau_)
                    p, _ = curve_fit(f, t, y, p0=(min(y[:3].mean() - base, -0.1), 15),
                                     bounds=([-60, 0.5], [0, 300]), maxfev=20000)
                    tau = p[1]
                    r2 = 1 - np.sum((y - f(t, *p)) ** 2) / np.sum((y - y.mean()) ** 2)
                except (RuntimeError, ValueError):
                    pass
                rows.append(dict(acquisition=rec, pulse=int(k), roi=int(roi),
                                 duree_s=int(round((n_stim[k] - 2 * config.STIM_PAD) * dt)),
                                 reponse_pct=resp, fenetre_s=t.max(), tau_s=tau, R2=r2))
    x = pd.DataFrame(rows)
    # Les N_TAU cellules à la réponse médiane la plus forte sur les 8 stimulations.
    best = x.groupby("roi")["reponse_pct"].median().nsmallest(N_TAU).index
    x["meilleure"] = x["roi"].isin(best)
    return x


# ------------------------------------------------------------------ figure
def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_2, length=0, labelsize=11)


def stim_bands(ax, t, stim, label=None):
    idx = np.flatnonzero(stim)
    for blk in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
        ax.axvspan(t[blk[0]], t[blk[-1]], color=STIM, alpha=STIM_ALPHA, lw=0, zorder=0)
        if label:
            ax.text((t[blk[0]] + t[blk[-1]]) / 2, 0.97, label, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=10, color=INK_2)


def heading(ax, letter, title, sub, y=1.16, y_sub=1.03):
    ax.text(-0.02, y, letter, transform=ax.transAxes, fontsize=17, fontweight="bold",
            color=INK, ha="right", va="bottom")
    ax.text(0, y, title, transform=ax.transAxes, fontsize=13.5, color=INK, va="bottom",
            fontweight="bold")
    ax.text(0, y_sub, sub, transform=ax.transAxes, fontsize=10.5, color=INK_2, va="bottom")


def population(ax, M, t):
    S = np.vstack([smooth(r) for r in M])
    with np.errstate(invalid="ignore"):
        q1, med, q3 = (np.nanpercentile(S, q, axis=0) for q in (25, 50, 75))
    ax.fill_between(t, q1, q3, color=ACCENT, alpha=0.22, lw=0, zorder=2)
    ax.plot(t, med, color=ACCENT, linewidth=2.6, zorder=3)
    ax.axhline(0, color=INK_2, linewidth=0.9, linestyle="--", zorder=1)


def zoom_box(masks, roi=ZOOM_ROI, cote=ZOOM_COTE):
    """Carré (y0, y1, x0, x1) de `cote` px centré sur le ROI, borné à l'image."""
    ys, xs = np.nonzero(masks == roi)
    cy, cx = (ys.min() + ys.max()) / 2, (xs.min() + xs.max()) / 2
    n0, n1 = masks.shape
    y0 = int(np.clip(cy - cote / 2, 0, n0 - cote))
    x0 = int(np.clip(cx - cote / 2, 0, n1 - cote))
    return y0, y0 + cote, x0, x0 + cote


def zoom(ax, im, box, color="white", cmap="gray"):
    """Encart en bas à droite agrandissant `box`, son cadre sur l'image, et deux lignes
    reliant les coins haut-droit et bas-gauche du cadre à ceux de l'encart."""
    y0, y1, x0, x1 = box
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor=color,
                               linewidth=1.6, zorder=5))
    pos = [0.55, 0.015, 0.435, 0.435]
    ins = ax.inset_axes(pos)
    crop = im[y0:y1, x0:x1]
    lo, hi = np.percentile(crop, (1, 99.7))
    ins.imshow(crop, cmap=cmap, vmin=lo, vmax=max(hi, lo + 1e-6))
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_edgecolor(color); sp.set_linewidth(1.8)
    # Coins de l'encart en coordonnées image (axe y d'imshow inversé).
    H, W = im.shape
    to_img = lambda fx, fy: (-0.5 + fx * W, -0.5 + (1 - fy) * H)
    for (sx, sy), (fx, fy) in (((x1, y0), (pos[0] + pos[2], pos[1] + pos[3])),
                               ((x0, y1), (pos[0], pos[1]))):
        ix, iy = to_img(fx, fy)
        ax.plot([sx, ix], [sy, iy], color=color, linewidth=1.2, zorder=4)


def build(lang, D, suffixe=""):
    T = TEXT[lang]
    fig = plt.figure(figsize=(16, 9), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, left=0.045, right=0.985, top=0.875, bottom=0.075,
                          wspace=0.14, hspace=0.36, width_ratios=[1.55, 1], height_ratios=[1.55, 1])
    fig.text(0.012, 0.975, T["title"], fontsize=21, fontweight="bold", color=INK, va="top")

    # A — 800 nm et 950 nm côte à côte
    # Placement explicite : la première image part du bord gauche des axes de C, les
    # deux images sont carrées et occupent toute la hauteur de la rangée du haut.
    top = gs[0, 0].get_position(fig)
    left = gs[1, 0].get_position(fig).x0
    fw, fh = fig.get_size_inches()
    side_w = top.height * fh / fw                    # largeur d'un carré, en fraction de figure
    gap = 0.005
    axes_A = []
    for j, (im, lab) in enumerate(zip(D["fixed"], (T["A1"], T["A2"]))):
        ax = fig.add_axes([left + j * (side_w + gap), top.y0, side_w, top.height])
        axes_A.append(ax)
        lo, hi = np.percentile(im, (1, 99.7))
        ax.imshow(im, cmap=CMAPS_A[j], vmin=lo, vmax=hi)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.text(0.03, 0.97, lab, transform=ax.transAxes, color="white", fontsize=11,
                fontweight="bold", va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.25", fc=(0, 0, 0, 0.45), ec="none"))
        if j == 0:
            # Même position relative que les lettres des autres panneaux (voir heading).
            # heading() place la lettre à -0.02 de la largeur des axes de C : même décalage
            # absolu ici, ramené à la largeur de l'image.
            c_w = gs[1, 0].get_position(fig).width
            ax.text(-0.02 * c_w / side_w, 1.02, "A",
                    transform=ax.transAxes, fontsize=17, fontweight="bold",
                    color=INK, ha="right", va="bottom")
            ax.text(0.0, 1.02, T["A"], transform=ax.transAxes, fontsize=13.5,
                    fontweight="bold", color=INK, ha="left", va="bottom")
        else:
            # Barre d'échelle en bas à gauche (le coin droit est pris par l'encart).
            px, x0, y = 50 / config.PIXEL_SIZE_UM, 50, im.shape[0] - 60
            ax.plot([x0, x0 + px], [y, y], color="white", linewidth=4, solid_capstyle="butt")
            ax.text(x0 + px / 2, y - 25, T["scale"], color="white", ha="center", va="bottom",
                    fontsize=11, fontweight="bold")
        zoom(ax, im, D["zoom_box"], cmap=CMAPS_A[j])

    M, t, stim = D["M"], D["t"], D["stim"]
    # B — base, première stimulation, retour (jusqu'avant la deuxième)
    # B s'étend de la fin des images de A (plus la place de l'étiquette y) au bord droit.
    b_left = axes_A[-1].get_position().x1 + 0.06
    b_ref = gs[0, 1].get_position(fig)
    # Le titre de B est à la hauteur de celui de A ; le haut des axes descend pour laisser
    # la place au sous-titre en dessous.
    y_title = top.y1 + 0.02 * top.height
    b_top = y_title - 0.045
    # Tracé recollé : la ligne de base (avant la première stimulation), puis la
    # stimulation B_PULSE et son retour, décalés dans le temps pour suivre directement la
    # ligne de base — comme si cette stimulation commençait à ~30 s. Le sous-titre le dit.
    Mb, tb, sb = D["B_M"], D["B_t"], D["B_stim"]
    lab, _ = ndi_label(sb)
    s_of = lambda k: tb[lab == k]
    seg1 = (tb < s_of(1).min()) & ~sb
    seg2 = (tb >= s_of(B_PULSE).min() - B_AVANT) & (tb < s_of(B_PULSE + 1).min() - 1)
    dt = float(np.median(np.diff(tb)))
    t_b = np.r_[tb[seg1], tb[seg2] - tb[seg2].min() + tb[seg1].max() + dt]
    # Retour décalé vers le haut pour que son maximum (médiane lissée) soit à 0 : le
    # retour, incomplet dans les données, finit ainsi sur la ligne de base.
    R = Mb[:, seg2]
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med_r = np.nanmedian(np.vstack([smooth(r) for r in R]), axis=0)
    D["B_offset"] = float(np.nanmax(med_r))
    M_b = np.hstack([Mb[:, seg1], R - D["B_offset"]])
    s_b = np.r_[sb[seg1], sb[seg2]]
    h_b = b_top - b_ref.y0
    ax = fig.add_axes([b_left, b_ref.y0, b_ref.x1 - b_left, h_b]); style(ax)
    stim_bands(ax, t_b, s_b, T["stim"])
    population(ax, M_b, t_b)
    fin = np.flatnonzero(s_b).max()                  # dernière image de la stimulation
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(M_b, axis=0)
    rec = np.arange(fin + 1, len(t_b))
    rec = rec[np.isfinite(med[rec])]
    x, y = t_b[rec] - t_b[fin], med[rec]
    if B_FIT_BASE is None:
        f = lambda x_, A, tau, C_: C_ + A * np.exp(-x_ / tau)
        p, cv = curve_fit(f, x, y, p0=(y[0] - y[-1], 15, y[-1]), maxfev=20000)
    else:
        f = lambda x_, A, tau: B_FIT_BASE + A * np.exp(-x_ / tau)
        p, cv = curve_fit(f, x, y, p0=(min(y[:3].mean() - B_FIT_BASE, -0.1), 15),
                          bounds=([-60, 0.5], [0, 300]), maxfev=20000)
    xx = np.linspace(x.min(), x.max(), 300)
    ax.plot(t_b[fin] + xx, f(xx, *p), color=ACCENT_2, linewidth=2.2, zorder=4)
    D["B_fit"] = (p[0], p[1], float(np.sqrt(cv[1, 1])))
    ax.set_xlabel(T["time"], fontsize=12, color=INK); ax.set_ylabel(T["y"], fontsize=12, color=INK)
    heading(ax, "B", T["B"], T["B_sub"].format(n=len(Mb)),
            y=(y_title - b_ref.y0) / h_b, y_sub=1 + 0.012 / h_b)

    # C — train complet
    ax = fig.add_subplot(gs[1, 0]); style(ax)
    stim_bands(ax, t, stim, T["stim"])
    population(ax, M, t)
    ax.set_xlabel(T["time"], fontsize=12, color=INK); ax.set_ylabel(T["y"], fontsize=12, color=INK)
    heading(ax, "C", T["C"], T["C_sub"].format(n=len(M)))

    # D — tau : une barre par rang de stimulation, les deux essais réunis
    ax = fig.add_subplot(gs[1, 1]); style(ax)
    r = D["taus"][D["taus"]["meilleure"]]
    pulses = sorted(r["pulse"].unique())
    rng = np.random.default_rng(0)
    stats = []
    for i, k in enumerate(pulses):
        v = r.loc[r["pulse"] == k, "tau_s"].dropna().to_numpy()
        v = v[v <= TAU_MAX]
        q1, med, q3 = np.percentile(v, (25, 50, 75))
        stats.append((k, len(v), med, q1, q3))
        ax.bar(i, med, width=0.62, color=ACCENT, alpha=0.25, edgecolor=ACCENT, linewidth=1.5, zorder=2)
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=26, color=ACCENT, alpha=0.9,
                   edgecolors=SURFACE, linewidths=0.6, zorder=3)
        ax.errorbar(i, med, yerr=[[med - q1], [q3 - med]], color=INK, linewidth=1.7, capsize=8,
                    capthick=1.7, zorder=4)
    ax.set_xticks(range(len(pulses)))
    ax.set_xticklabels([f"{T['stim_k'].format(k=k)}\nn = {n}" for k, n, *_ in stats], fontsize=10.5)
    ax.set_xlim(-0.6, len(pulses) - 0.4)
    ax.set_ylim(0, TAU_MAX * 1.05)
    ax.set_ylabel(T["tau"], fontsize=12, color=INK)
    heading(ax, "D", T["D"], T["D_sub"].format(n=r["roi"].nunique()), y=1.2, y_sub=1.03)
    D["tau_bar"] = stats
    ax_cat = ax                      # axe x catégoriel : pas de virgule décimale

    if lang == "FR":
        comma = FuncFormatter(lambda val, _: f"{val:g}".replace(".", ",").replace("-", "−"))
        for a in fig.axes:
            if not a.get_images():
                a.yaxis.set_major_formatter(comma)
                if a is not ax_cat:
                    a.xaxis.set_major_formatter(comma)
    out = OUT / f"slide_resultats_{lang}{suffixe}.png"
    fig.savefig(out, dpi=160, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.relative_to(HERE)}  (retour de B décalé de {-D['B_offset']:+.2f} pts)")


def main():
    OUT.mkdir(exist_ok=True)
    keep, n_all = kept_rois()
    d = pd.read_csv(TR_ROOT / "results" / "ratio_temps.csv")
    d = d[d["roi"].isin(keep)]
    M, t, stim = matrix(d[d["acquisition"] == TRAINS[0]], "dR_R0_pct")
    B_M, B_t, B_stim = matrix(d[d["acquisition"] == B_TRAIN], "dR_R0_pct")
    taus = recovery_taus(d)
    taus.to_csv(TR_ROOT / "results" / "tau_par_cellule.csv", index=False)
    D = dict(fixed=fixed_images(), M=M, t=t, stim=stim, B_M=B_M, B_t=B_t, B_stim=B_stim, taus=taus, n_keep=len(keep), n_all=n_all,
             zoom_box=zoom_box(np.load(TR_ROOT / "masks" / f"{FOV}_masks.npy")))
    b_ = taus[taus["meilleure"]]
    print(f"{FOV} : {len(keep)}/{n_all} cellules retenues ; tau des {N_TAU} meilleures : ROIs "
          f"{sorted(b_['roi'].unique())}")
    print(b_.groupby("pulse").agg(
          n=("tau_s", "count"), mediane=("tau_s", "median"),
          q1=("tau_s", lambda v: v.quantile(.25)), q3=("tau_s", lambda v: v.quantile(.75)),
          rejetes=("tau_s", lambda v: int((v > TAU_MAX).sum()))).round(1).to_string())
    k_ = b_[b_["tau_s"] <= TAU_MAX].groupby("pulse")["tau_s"]
    print("après rejet des tau > %d s :" % TAU_MAX)
    print(pd.DataFrame({"n": k_.size(), "mediane": k_.median(), "q1": k_.quantile(.25),
                        "q3": k_.quantile(.75)}).round(1).to_string())
    for suffixe in THEMES:
        appliquer_theme(suffixe)
        for lang in TEXT:
            build(lang, D, suffixe)
    A_, tau_, se_ = D["B_fit"]
    print(f"B, remontée ajustée : τ = {tau_:.1f} ± {se_:.1f} s, A = {A_:+.2f} pts "
          f"(base {'libre' if B_FIT_BASE is None else f'fixée à {B_FIT_BASE:g}'})")


if __name__ == "__main__":
    main()
