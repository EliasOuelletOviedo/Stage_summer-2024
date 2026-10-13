"""Slide des résultats préliminaires — tranche du 2026-10-08 (FLIM), sur le modèle de celle
du 2026-09-21 (SCLM+eNpHR_2026-09-21/slide_resultats.py).

Usage :  python3 slide_resultats.py   (de n'importe où, après les étapes 5 à 7 du pipeline)

Produit, dans ce dossier : slide_resultats_FR.png, slide_resultats_EN.png et leurs versions
sur fond noir (_FR_dark.png, _EN_dark.png). `charger()` sert aussi à presentation.py.

Quatre panneaux alignés en 2 x 2, somas Cellpose retenus (étape 3) :
  A  la session entière, ratio YFP/CFP : quatre stimulations, quatre réponses ;
  B  une stimulation, puis le retour : ligne de base, puis stimulation 1 jusqu'à la 2e ;
  C  comme A, pour la durée de vie de la CFP (le donneur) ;
  D  temps de retour τ, une barre par stimulation, un point par cellule.

Choix (rien de tout ceci n'est écrit sur la slide) :
  - A, B, C sont corrigés de la descente graduelle due à l'imagerie, par la méthode de
    l'étape 6 (pente commune, un niveau par intervalle, estimée sur la base et entre les
    stims en écartant les 30 s qui suivent chacune), appliquée à la cellule médiane :
    ratio multiplié par exp(−β·t), τ CFP diminué de β_τ·t. Puis ΔR/R₀ et Δτ par rapport
    à la base avant la 1re stim ;
  - blocs de 3 s (étape 5) ; blocs de stim masqués (la stim ajoute ~4 % de photons dans la
    YFP) ; médiane ± écart interquartile entre cellules, lissés par une moyenne glissante
    triangulaire centrée sur LISSAGE blocs (poids 1/4, 1/2, 1/4 : pas de décalage) ;
  - B : la remontée après la stimulation est ajustée par C + A·exp(−t/τ), t depuis la fin
    de la stim, équilibre C libre (comme à l'étape 7), sur la médiane non lissée des
    cellules ; la courbe ajustée est tracée sur la remontée ;
  - D : τ de l'étape 7 (C + A·exp(−t/τ) sur le retour de chaque cellule, équilibre libre),
    retours mesurables (une baisse, τ hors des bornes), τ ≤ TAU_MAX_S (40 s) ;
    barre = médiane ± écart interquartile.
"""
import json
import sys
import warnings
from pathlib import Path

ICI = Path(__file__).resolve().parent              # <jeu>/présentation/
ROOT = ICI.parent                                  # <jeu>/
sys.path.insert(0, str(ROOT.parent / "sclm_flim"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

from config import BASE_SKIP, EXCLURE_APRES_STIM_S  # noqa: E402


# =====================================================================================
# RÉGLAGES
# =====================================================================================
TEXTES = {
    "FR": dict(
        title="SCLm + eNpHR en tranches organotypiques : résultats préliminaires",
        A="Quatre stimulations, quatre réponses",
        A_sub="médiane ± écart interquartile, n = {n} cellules",
        B="Une stimulation, puis le retour",
        B_sub="ligne de base, puis stimulation 1\nmédiane ± écart interquartile, n = {n} cellules",
        C="Quatre stimulations, durée de vie de la CFP",
        C_sub="médiane ± écart interquartile, n = {n} cellules",
        D="Temps de retour τ",
        D_sub="un point par cellule\nbarre = médiane ± écart interquartile",
        time="temps (s)", y_ratio="ΔR/R₀ YFP/CFP  (%)", y_tau="Δτ CFP  (ps)", stim="stim",
        tau="τ  (s)", tick="stim {k}\nn = {n}"),
    "EN": dict(
        title="SCLm + eNpHR in organotypic slices: preliminary results",
        A="Four stimulations, four responses",
        A_sub="median ± interquartile range, n = {n} cells",
        B="One stimulation, then recovery",
        B_sub="baseline, then stimulation 1\nmedian ± interquartile range, n = {n} cells",
        C="Four stimulations, CFP lifetime",
        C_sub="median ± interquartile range, n = {n} cells",
        D="Recovery time τ",
        D_sub="one dot per cell\nbar = median ± interquartile range",
        time="time (s)", y_ratio="ΔR/R₀ YFP/CFP  (%)", y_tau="Δτ CFP  (ps)", stim="stim",
        tau="τ  (s)", tick="stim {k}\nn = {n}"),
}

# Thèmes : "" = fond clair (celui du 2026-09-21), "_dark" = fond noir. Teintes sombres :
# celles de la palette de référence pour fond sombre, validées sur #000000.
THEMES = {
    "": dict(SURFACE="#fcfcfb", INK="#0b0b0b", INK_2="#52514e", GRID="#e6e6e2",
             ACCENT="#2a78d6", ACCENT_TAU="#1baf7a", FIT="#eb6834", STIM="#f4d35e", STIM_ALPHA=0.45),
    "_dark": dict(SURFACE="#000000", INK="#ffffff", INK_2="#c3c2b7", GRID="#2f2f2c",
                  ACCENT="#3987e5", ACCENT_TAU="#199e70", FIT="#d95926", STIM="#c9a227", STIM_ALPHA=0.35),
}

FIGSIZE = (16, 9)
DPI = 160
# Grille 2 x 2, en fractions de la figure : colonnes de largeurs 1,55 : 1 (comme le
# 2026-09-21), deux rangées de même hauteur. Tous les axes sont placés à la main, pour
# que les bords des quatre panneaux et leurs en-têtes soient alignés.
GAUCHE, DROITE, ECART_COL = 0.055, 0.985, 0.075
BAS, HAUT_LIGNE, ECART_LIGNE = 0.085, 0.30, 0.145
LISSAGE = 3               # blocs de 3 s, nombre impair ; 1 = pas de lissage
TAU_MAX_S = 40            # D : τ plus longs écartés (la fenêtre d'observation fait 108 s)


# =====================================================================================
# Données (produites par le pipeline)
# =====================================================================================
def pente_additive(t, y, stim, fen, exclure=EXCLURE_APRES_STIM_S):
    """Pente commune de y (un niveau par intervalle), hors stims : même règle que l'étape 6."""
    ivs = [(-np.inf, fen[0][0])] + [(b + exclure, fen[i + 1][0] if i + 1 < len(fen) else np.inf)
                                   for i, (_, b) in enumerate(fen)]
    num = den = 0.0
    for a, b in ivs:
        m = ~stim & np.isfinite(y) & (t >= a) & (t < b)
        if m.sum() >= 5:
            x = t[m] - t[m].mean()
            num += (x * (y[m] - y[m].mean())).sum()
            den += (x * x).sum()
    return num / den


def charger():
    d = pd.read_csv(ROOT / "results" / "mesures_blocs.csv")
    d = d[d["famille"] == "cellpose"]
    info = json.loads((ROOT / "projections" / "session.json").read_text())
    fen = [(s["t0"], s["t1"]) for s in info["stims"]]
    blocs = d.drop_duplicates("bloc").set_index("bloc").sort_index()
    t, stim = blocs["t_s"].to_numpy(float), blocs["stim"].to_numpy(bool)
    base = (blocs["phase"].to_numpy() == "avant") & (blocs.index.to_numpy() > BASE_SKIP)
    piv = lambda col: d.pivot(index="roi", columns="bloc", values=col).sort_index(axis=1).to_numpy(float)
    t0 = info["blocs"][0]["t0"]

    # Ratio corrigé de la descente graduelle, puis ΔR/R₀. La pente est celle de la
    # cellule médiane (ln du ratio), estimée comme à l'étape 6 : la courbe tracée est une
    # médiane de cellules, et la cellule médiane dérive moins vite que les somas réunis
    # (photons additionnés, −1,96 %/100 s).
    R = piv("ratio")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        beta = pente_additive(t, np.nanmedian(np.log(R), axis=0), stim, fen)
    R = R * np.exp(-beta * (t - t0))[None, :]
    dR = (R / np.nanmedian(R[:, base], axis=1, keepdims=True) - 1) * 100

    # τ CFP corrigé de la même façon (pente estimée sur la médiane des somas), puis Δτ.
    tau = piv("tau_ch1") * 1000                                      # ps
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        beta_tau = pente_additive(t, np.nanmedian(tau, axis=0), stim, fen)
    tau = tau - beta_tau * (t - t0)[None, :]
    dT = tau - np.nanmedian(tau[:, base], axis=1, keepdims=True)
    dR_tout, dT_tout = dR.copy(), dT.copy()                          # avec les blocs de stim
    for M in (dR, dT):
        M[:, stim] = np.nan                                          # blocs de stim masqués

    taus = pd.read_csv(ROOT / "results" / "temps_de_remontee.csv")
    taus = taus[(taus["niveau"] == "roi") & (taus["stim"] != "moyenne") & taus["mesurable"]
                & (taus["tau_s"] <= TAU_MAX_S)]
    print(f"{ROOT.name} : {len(dR)} somas ; pente du ratio {beta * 1e4:+.2f} %/100 s, "
          f"de τ CFP {beta_tau * 100:+.1f} ps/100 s")
    return dict(t=t, fen=fen, stim=stim, dR=dR, dT=dT, dR_tout=dR_tout, dT_tout=dT_tout,
                n=len(dR), taus=taus)


# =====================================================================================
# Figure
# =====================================================================================
def style(ax, C):
    ax.set_facecolor(C["SURFACE"])
    ax.grid(True, color=C["GRID"], linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C["GRID"])
    ax.tick_params(colors=C["INK_2"], length=0, labelsize=11)


def heading(ax, C, letter, title, sub, y=1.16, y_sub=1.13):
    """Lettre et titre à la même hauteur dans les quatre panneaux ; sous-titre dessous."""
    if letter:
        ax.text(-0.02, y, letter, transform=ax.transAxes, fontsize=17, fontweight="bold",
                color=C["INK"], ha="right", va="bottom")
    ax.text(0, y, title, transform=ax.transAxes, fontsize=13.5, color=C["INK"],
            va="bottom", fontweight="bold")
    ax.text(0, y_sub, sub, transform=ax.transAxes, fontsize=10.5, color=C["INK_2"], va="top")


def stim_bands(ax, C, fen, label):
    for a, b in fen:
        ax.axvspan(a, b, color=C["STIM"], alpha=C["STIM_ALPHA"], lw=0, zorder=0)
        ax.text((a + b) / 2, 0.97, label, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=10, color=C["INK_2"])


def lisser(M):
    """Moyenne glissante triangulaire centrée sur LISSAGE blocs, qui ignore les blocs vides."""
    if LISSAGE <= 1:
        return M.copy()
    k = np.bartlett(LISSAGE + 2)[1:-1]                               # 3 -> 1/2, 1, 1/2
    ok = np.isfinite(M)
    num = ndi.convolve1d(np.where(ok, M, 0.0), k, axis=1, mode="constant")
    den = ndi.convolve1d(ok.astype(float), k, axis=1, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        S = num / den
    S[~ok] = np.nan                                                  # les trous restent vides
    return S


def ajuster_remontee(t, med, t_fin):
    """C + A·exp(−(t − t_fin)/τ) sur la remontée (équilibre C libre) : (A, τ, C), erreurs
    types, et premier et dernier temps ajustés (s depuis t_fin)."""
    ok = np.isfinite(med) & (t > t_fin)
    x, y = t[ok] - t_fin, med[ok]
    f = lambda x_, A, tau, C_: C_ + A * np.exp(-x_ / tau)
    p, cv = curve_fit(f, x, y, p0=[y[0] - y[-1], 15.0, y[-1]],
                      bounds=([-50, 0.3, -20], [50, 300, 20]), maxfev=20000)
    return p, np.sqrt(np.diag(cv)), (x.min(), x.max())


def population(ax, C, M, t, couleur):
    S = lisser(M)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q1, med, q3 = (np.nanpercentile(S, q, axis=0) for q in (25, 50, 75))
    ax.fill_between(t, q1, q3, color=couleur, alpha=0.22, lw=0, zorder=2)
    ax.plot(t, med, color=couleur, linewidth=2.6, zorder=3)
    ax.axhline(0, color=C["INK_2"], linewidth=0.9, linestyle="--", zorder=1)
    ok = np.isfinite(q1) & np.isfinite(q3)
    lo, hi = np.nanmin(q1[ok]), np.nanmax(q3[ok])
    ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.22 * (hi - lo))       # place pour « stim »


def axes_grille(fig):
    lw = (DROITE - GAUCHE - ECART_COL) * 1.55 / 2.55
    rw = DROITE - GAUCHE - ECART_COL - lw
    y_haut = BAS + HAUT_LIGNE + ECART_LIGNE
    rect = {"A": [GAUCHE, y_haut, lw, HAUT_LIGNE], "B": [GAUCHE + lw + ECART_COL, y_haut, rw, HAUT_LIGNE],
            "C": [GAUCHE, BAS, lw, HAUT_LIGNE], "D": [GAUCHE + lw + ECART_COL, BAS, rw, HAUT_LIGNE]}
    return {k: fig.add_axes(v) for k, v in rect.items()}


def build(lang, D, suffixe):
    T, C = TEXTES[lang], THEMES[suffixe]
    fig = plt.figure(figsize=FIGSIZE, facecolor=C["SURFACE"])
    fig.text(0.012, 0.975, T["title"], fontsize=21, fontweight="bold", color=C["INK"], va="top")
    ax = axes_grille(fig)
    t, fen, n = D["t"], D["fen"], D["n"]

    # A — la session entière, ratio
    a = ax["A"]
    style(a, C)
    stim_bands(a, C, fen, T["stim"])
    population(a, C, D["dR"], t, C["ACCENT"])
    a.set_xlabel(T["time"], fontsize=12, color=C["INK"])
    a.set_ylabel(T["y_ratio"], fontsize=12, color=C["INK"])
    heading(a, C, "A", T["A"], T["A_sub"].format(n=n))

    # B — ligne de base, stimulation 1 et son retour, jusqu'à la stimulation 2
    b = ax["B"]
    style(b, C)
    m = t < fen[1][0] if len(fen) > 1 else np.ones_like(t, bool)
    stim_bands(b, C, fen[:1], T["stim"])
    population(b, C, D["dR"][:, m], t[m], C["ACCENT"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(D["dR"][:, m], axis=0)
    (A_, tau_, C_), se, (x0, x1) = ajuster_remontee(t[m], med, fen[0][1])
    xx = np.linspace(x0, x1, 300)                    # du premier au dernier point ajusté
    b.plot(fen[0][1] + xx, C_ + A_ * np.exp(-xx / tau_), color=C["FIT"], linewidth=2.2, zorder=4)
    D["fit_B"] = (A_, tau_, C_, se)
    b.set_xlim(t[m][0] - 3, t[m][-1] + 3)
    b.set_xlabel(T["time"], fontsize=12, color=C["INK"])
    b.set_ylabel(T["y_ratio"], fontsize=12, color=C["INK"])
    heading(b, C, "B", T["B"], T["B_sub"].format(n=n))

    # C — la session entière, durée de vie de la CFP
    c = ax["C"]
    style(c, C)
    stim_bands(c, C, fen, T["stim"])
    population(c, C, D["dT"], t, C["ACCENT_TAU"])
    c.set_xlabel(T["time"], fontsize=12, color=C["INK"])
    c.set_ylabel(T["y_tau"], fontsize=12, color=C["INK"])
    heading(c, C, "C", T["C"], T["C_sub"].format(n=n))
    c.set_xlim(a.get_xlim())

    # D — τ de retour, une barre par stimulation
    dd = ax["D"]
    style(dd, C)
    r = D["taus"]
    stims = sorted(r["stim"].astype(int).unique())
    rng = np.random.default_rng(0)
    ticks = []
    for i, k in enumerate(stims):
        v = r.loc[r["stim"].astype(int) == k, "tau_s"].to_numpy()
        q1, med, q3 = np.percentile(v, (25, 50, 75))
        dd.bar(i, med, width=0.62, color=C["ACCENT"], alpha=0.25, edgecolor=C["ACCENT"],
               linewidth=1.5, zorder=2)
        dd.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=26, color=C["ACCENT"], alpha=0.9,
                   edgecolors=C["SURFACE"], linewidths=0.6, zorder=3)
        dd.errorbar(i, med, yerr=[[med - q1], [q3 - med]], color=C["INK"], linewidth=1.7,
                    capsize=8, capthick=1.7, zorder=4)
        ticks.append(T["tick"].format(k=k, n=len(v)))
    dd.set_xticks(range(len(stims)))
    dd.set_xticklabels(ticks, fontsize=10.5)
    dd.set_xlim(-0.6, len(stims) - 0.4)
    dd.set_ylim(0, r["tau_s"].max() * 1.1)
    dd.set_ylabel(T["tau"], fontsize=12, color=C["INK"])
    heading(dd, C, "D", T["D"], T["D_sub"])

    if lang == "FR":
        virgule = FuncFormatter(lambda val, _: f"{val:g}".replace(".", ",").replace("-", "−"))
        for k_, a_ in ax.items():
            a_.yaxis.set_major_formatter(virgule)
            if k_ != "D":
                a_.xaxis.set_major_formatter(virgule)
    out = ICI / f"slide_resultats_{lang}{suffixe}.png"
    fig.savefig(out, dpi=DPI, facecolor=C["SURFACE"])
    plt.close(fig)
    print(f"  -> {out.name}")


def main():
    D = charger()
    for suffixe in THEMES:
        for lang in TEXTES:
            build(lang, D, suffixe)
    A_, tau_, C_, se = D["fit_B"]
    print(f"  B, remontée après la stimulation 1 : τ = {tau_:.1f} ± {se[1]:.1f} s, "
          f"A = {A_:+.2f} %, équilibre {C_:+.2f} %")


if __name__ == "__main__":
    main()
