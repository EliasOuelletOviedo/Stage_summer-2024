"""Slide des résultats préliminaires — tranche du 2026-10-08 (FLIM), sur le modèle de celle
du 2026-09-21 (SCLM+eNpHR_2026-09-21/slide_resultats.py).

Usage :  python3 slide_resultats.py   (de n'importe où, après les étapes 5 à 7 du pipeline)

Produit, dans ce dossier : slide_resultats_FR.png, slide_resultats_EN.png et leurs versions
sur fond noir (_FR_dark.png, _EN_dark.png). `charger()` et `ajuster_retours()` servent aussi
à presentation.py.

Quatre panneaux alignés en 2 x 2, somas Cellpose retenus (étape 3) ; le ratio en haut, la
durée de vie de la CFP (le donneur) en bas :
  A  la session entière, ratio YFP/CFP : quatre stimulations, quatre réponses, et le
     retour ajusté après chacune ;
  B  temps de retour τ du ratio, une barre par stimulation, un point par cellule ;
  C  comme A, pour la durée de vie de la CFP ;
  D  comme B, pour la durée de vie de la CFP.

Choix (rien de tout ceci n'est écrit sur la slide) :
  - A et C sont corrigés de la descente graduelle due à l'imagerie, par la méthode de
    l'étape 6 (pente commune, un niveau par intervalle, estimée sur la base et entre les
    stims en écartant les 30 s qui suivent chacune), appliquée à la cellule médiane :
    ratio multiplié par exp(−β·t), τ CFP diminué de β_τ·t. Puis ΔR/R₀ et Δτ par rapport
    à la base avant la 1re stim ;
  - blocs de 3 s (étape 5) ; blocs de stim masqués (la stim ajoute ~4 % de photons dans la
    YFP) ; médiane ± écart interquartile entre cellules, lissés par une moyenne glissante
    triangulaire centrée sur LISSAGE blocs (poids 1/4, 1/2, 1/4 : pas de décalage) ;
  - A et C : après chaque stim et jusqu'à la suivante, C_k + A_k·exp(−t/τ), t depuis la
    fin de la stim, sur la médiane non lissée des cellules : un τ commun aux quatre
    retours, une amplitude et un équilibre par stim. Un à un, les retours de la durée de
    vie sont trop bruités pour être ajustés (τ de 29 ± 122 s à la stim 1, à la borne de
    300 s à la stim 4) ; le même modèle pour le ratio garde A et C comparables. Courbes
    tracées du premier au dernier point ajusté ;
  - B : τ de l'étape 7 (C + A·exp(−t/τ) sur le retour de chaque cellule, équilibre libre),
    retours mesurables (une baisse, τ hors des bornes), τ ≤ TAU_MAX_S (40 s) ;
    barre = médiane ± écart interquartile ;
  - D : le même ajustement sur la durée de vie de chaque cellule, telle qu'en C, de la fin
    de chaque stim à la suivante ; retours mesurables (une hausse, τ hors des bornes),
    τ ≤ TAU_MAX_S. Un retour isolé de la durée de vie est deux fois plus bruité que celui
    du ratio (erreur type médiane de τ : 2 τ, contre 1 τ en B) et un retour ajusté sur
    trois va dans le mauvais sens : D montre surtout la dispersion due au bruit ;
  - B et D ont la même échelle de τ.
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
from scipy.optimize import curve_fit, least_squares  # noqa: E402

from config import BASE_SKIP, EXCLURE_APRES_STIM_S  # noqa: E402
from step7_retour import BORNES_TAU, MESURABLE_TAU  # noqa: E402


# =====================================================================================
# RÉGLAGES
# =====================================================================================
TEXTES = {
    "FR": dict(
        title="SCLm + eNpHR en tranches organotypiques : résultats préliminaires",
        A="Quatre stimulations, quatre réponses",
        A_sub="médiane ± écart interquartile, n = {n} cellules",
        B="Temps de retour τ",
        B_sub="un point par cellule\nbarre = médiane ± écart interquartile",
        C="Quatre stimulations, durée de vie de la CFP",
        C_sub="médiane ± écart interquartile, n = {n} cellules",
        D="Temps de retour τ, durée de vie de la CFP",
        D_sub="un point par cellule\nbarre = médiane ± écart interquartile",
        time="temps (s)", y_ratio="ΔR/R₀ YFP/CFP  (%)", y_tau="Δτ CFP  (ps)", stim="stim",
        tau="τ  (s)", tick="stim {k}\nn = {n}"),
    "EN": dict(
        title="SCLm + eNpHR in organotypic slices: preliminary results",
        A="Four stimulations, four responses",
        A_sub="median ± interquartile range, n = {n} cells",
        B="Recovery time τ",
        B_sub="one dot per cell\nbar = median ± interquartile range",
        C="Four stimulations, CFP lifetime",
        C_sub="median ± interquartile range, n = {n} cells",
        D="Recovery time τ, CFP lifetime",
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
TAU_MAX_S = 40            # B, D : τ plus longs écartés (la fenêtre d'observation fait 108 s)


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


def apres_stims(t, fen):
    """Pour chaque stim : sa fin, et les blocs qui la suivent jusqu'à la stim suivante."""
    return [(b, (t > b) & (t < (fen[i + 1][0] if i + 1 < len(fen) else np.inf)))
            for i, (_, b) in enumerate(fen)]


def ajuster_retours(t, med, fen):
    """Après chaque stimulation, C_k + A_k·exp(−(t − fin_k)/τ) jusqu'à la stimulation
    suivante : un τ commun, une amplitude et un équilibre par stim, moindres carrés.
    Une à une, les redescentes de la durée de vie sont trop bruitées pour être ajustées
    (τ de 29 ± 122 s à la stim 1, à la borne de 300 s à la stim 4) ; ensemble, si.
    Rend τ, son erreur type, et [(A_k, C_k, fin_k, premier, dernier temps ajusté)]."""
    segs = []
    for b, m in apres_stims(t, fen):
        m = m & np.isfinite(med)
        segs.append((b, t[m] - b, med[m]))

    def residus(p):
        return np.concatenate([p[2 + 2 * k] + p[1 + 2 * k] * np.exp(-x / p[0]) - y
                               for k, (_, x, y) in enumerate(segs)])

    p0 = [15.0] + [v for _, x, y in segs for v in (y[:3].mean() - y[-10:].mean(), y[-10:].mean())]
    n = len(segs)
    r = least_squares(residus, p0, bounds=([BORNES_TAU[0]] + [-200, -100] * n,
                                           [BORNES_TAU[1]] + [200, 100] * n))
    s2 = (r.fun ** 2).sum() / max(len(r.fun) - len(r.x), 1)
    se = np.sqrt(np.diag(np.linalg.inv(r.jac.T @ r.jac) * s2))
    return r.x[0], se[0], [(r.x[1 + 2 * k], r.x[2 + 2 * k], b, x.min(), x.max())
                           for k, (b, x, _) in enumerate(segs)]


def ajuster_cellules(t, M, fen, signe):
    """Chaque cellule (ligne de M) et chaque stim : C + A·exp(−t/τ) de la fin de la stim à
    la suivante, équilibre libre, comme à l'étape 7. Mesurable : A du signe attendu
    (−1 : une baisse, +1 : une hausse) et τ hors des bornes."""
    f = lambda x, A, tau, C_: C_ + A * np.exp(-x / tau)
    rows = []
    for k, (b, m) in enumerate(apres_stims(t, fen)):
        for i, y in enumerate(M[:, m]):
            ok = np.isfinite(y)
            if ok.sum() < 10:
                continue
            x, y = t[m][ok] - b, y[ok]
            try:
                p, cv = curve_fit(f, x, y, p0=[y[:3].mean() - y[-10:].mean(), 15.0, y[-10:].mean()],
                                  bounds=([-200, BORNES_TAU[0], -200], [200, BORNES_TAU[1], 200]),
                                  maxfev=20000)
            except (RuntimeError, ValueError):
                continue
            rows.append(dict(cellule=i, stim=k + 1, tau_s=p[1], tau_se_s=np.sqrt(cv[1, 1]),
                             mesurable=bool(np.sign(p[0]) == signe
                                            and MESURABLE_TAU[0] < p[1] < MESURABLE_TAU[1])))
    return pd.DataFrame(rows)


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

    # Retours : ajustés sur la médiane des cellules (A, C), puis cellule par cellule (B, D).
    taus = pd.read_csv(ROOT / "results" / "temps_de_remontee.csv")
    taus_R = taus[(taus["niveau"] == "roi") & (taus["stim"] != "moyenne") & taus["mesurable"]
                  & (taus["tau_s"] <= TAU_MAX_S)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit_R = ajuster_retours(t, np.nanmedian(dR, axis=0), fen)
        fit_T = ajuster_retours(t, np.nanmedian(dT, axis=0), fen)
        cel = ajuster_cellules(t, dT, fen, signe=+1)
    taus_T = cel[cel["mesurable"] & (cel["tau_s"] <= TAU_MAX_S)]
    print(f"{ROOT.name} : {len(dR)} somas ; pente du ratio {beta * 1e4:+.2f} %/100 s, "
          f"de τ CFP {beta_tau * 100:+.1f} ps/100 s")
    return dict(t=t, fen=fen, stim=stim, base=base, dR=dR, dT=dT, dR_tout=dR_tout,
                dT_tout=dT_tout, n=len(dR), fit_R=fit_R, fit_T=fit_T, taus_R=taus_R, taus_T=taus_T)


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


def session(ax, C, T, M, t, fen, fit, couleur, ylabel):
    """La session entière, médiane ± écart interquartile, et le retour ajusté après
    chaque stim."""
    style(ax, C)
    stim_bands(ax, C, fen, T["stim"])
    population(ax, C, M, t, couleur)
    tau, _, courbes = fit
    for A_, C_, b, x0, x1 in courbes:
        xx = np.linspace(x0, x1, 300)                # du premier au dernier point ajusté
        ax.plot(b + xx, C_ + A_ * np.exp(-xx / tau), color=C["FIT"], linewidth=2.2, zorder=4)
    ax.set_xlabel(T["time"], fontsize=12, color=C["INK"])
    ax.set_ylabel(ylabel, fontsize=12, color=C["INK"])


def temps_de_retour(ax, C, T, r, couleur, haut):
    """τ de retour, une barre par stimulation (médiane ± écart interquartile), un point
    par cellule."""
    style(ax, C)
    stims = sorted(r["stim"].astype(int).unique())
    rng = np.random.default_rng(0)
    ticks = []
    for i, k in enumerate(stims):
        v = r.loc[r["stim"].astype(int) == k, "tau_s"].to_numpy()
        q1, med, q3 = np.percentile(v, (25, 50, 75))
        ax.bar(i, med, width=0.62, color=couleur, alpha=0.25, edgecolor=couleur,
               linewidth=1.5, zorder=2)
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=26, color=couleur, alpha=0.9,
                   edgecolors=C["SURFACE"], linewidths=0.6, zorder=3)
        ax.errorbar(i, med, yerr=[[med - q1], [q3 - med]], color=C["INK"], linewidth=1.7,
                    capsize=8, capthick=1.7, zorder=4)
        ticks.append(T["tick"].format(k=k, n=len(v)))
    ax.set_xticks(range(len(stims)))
    ax.set_xticklabels(ticks, fontsize=10.5)
    ax.set_xlim(-0.6, len(stims) - 0.4)
    ax.set_ylim(0, haut)
    ax.set_ylabel(T["tau"], fontsize=12, color=C["INK"])


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
    haut = max(D["taus_R"]["tau_s"].max(), D["taus_T"]["tau_s"].max()) * 1.1   # B et D

    # A — la session entière, ratio, et le retour ajusté après chaque stim
    session(ax["A"], C, T, D["dR"], t, fen, D["fit_R"], C["ACCENT"], T["y_ratio"])
    heading(ax["A"], C, "A", T["A"], T["A_sub"].format(n=n))

    # B — τ de retour du ratio, une barre par stimulation
    temps_de_retour(ax["B"], C, T, D["taus_R"], C["ACCENT"], haut)
    heading(ax["B"], C, "B", T["B"], T["B_sub"])

    # C — comme A, durée de vie de la CFP
    session(ax["C"], C, T, D["dT"], t, fen, D["fit_T"], C["ACCENT_TAU"], T["y_tau"])
    heading(ax["C"], C, "C", T["C"], T["C_sub"].format(n=n))
    ax["C"].set_xlim(ax["A"].get_xlim())

    # D — comme B, durée de vie de la CFP
    temps_de_retour(ax["D"], C, T, D["taus_T"], C["ACCENT_TAU"], haut)
    heading(ax["D"], C, "D", T["D"], T["D_sub"])

    if lang == "FR":
        virgule = FuncFormatter(lambda val, _: f"{val:g}".replace(".", ",").replace("-", "−"))
        for k_, a_ in ax.items():
            a_.yaxis.set_major_formatter(virgule)
            if k_ in ("A", "C"):
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
    for nom, (tau, se, courbes) in (("A, ratio (%)", D["fit_R"]), ("C, durée de vie (ps)", D["fit_T"])):
        print(f"  {nom} : τ commun = {tau:.1f} ± {se:.1f} s ; A par stim "
              f"{', '.join(f'{a:+.2f}' for a, *_ in courbes)}")
    for nom, r in (("B, ratio", D["taus_R"]), ("D, durée de vie", D["taus_T"])):
        n_k = r.groupby(r["stim"].astype(int)).size().tolist()
        print(f"  {nom} : n par stim {n_k} ; τ médian {r['tau_s'].median():.1f} s ; "
              f"erreur type médiane {np.median(r['tau_se_s'] / r['tau_s']):.1f} τ")


if __name__ == "__main__":
    main()
