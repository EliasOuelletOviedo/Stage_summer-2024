"""Étape 9 — redescente de la [Cl⁻] après chaque stim, lue par le filtre bayésien (étape 8).

Usage :  python3 step9_redescente.py <dossier> [<dossier> …]   (après les étapes 7 et 8 ; quelques s)

Données : la [Cl⁻] de chaque soma, tranche par tranche (~3 s), lissée et en temps réel
(`results/filtre_par_trame.csv` de l'étape 8). Le filtre n'est pas relancé.

1. **Dérive.** Entre les stims, c monte lentement dans presque tous les somas (le ratio
   descend, étape 6) : autant que la réponse à une stim en une minute environ. Même
   correction qu'à l'étape 6, soma par soma : ln c = a_k + β·t, une pente β par soma, un
   niveau a_k par intervalle sans stim (la base avant la 1re stim, puis de
   EXCLURE_APRES_STIM_S après chaque stim jusqu'à la suivante) ; c corrigé = c·exp(−β(t − t₀)).
   Les niveaux libres laissent les marches des stims hors de la pente.
2. **Alignement.** Δc = c corrigé − sa moyenne des 10 s avant la stim (tranches lues),
   en fonction du temps depuis la fin de la stim, ramené sur une grille commune au pas
   d'une tranche. L'interpolation ne traverse jamais la stim : ses tranches ne sont pas lues.
3. **Redescente.** Δc(x) = C + A·exp(−x/τ), x depuis la fin de la stim, sur la fenêtre de
   l'étape 7 (jusqu'à la stim suivante, commune à toutes). Elle est mesurable si A > 0
   (une hausse, puis une redescente) et si τ ne touche pas ses bornes. Ajustée sur la
   médiane des somas, pour chaque stim et pour la moyenne des stims ; puis pour chaque
   soma, sur la moyenne de ses stims.

À lire avec :
- le c lissé voit les tranches suivantes : il commence à monter avant la fin de la stim
  et il raccourcit le pic ; le c en temps réel est causal, mais il est en retard et bruité ;
- les tranches lissées sont corrélées entre elles : les erreurs types de τ sont optimistes ;
- le τ du ratio (étape 7, tranches d'1 s, ratio corrigé de sa dérive) sert de référence :
  le même retour, vu sans le modèle du capteur ni le lissage ;
- la [Cl⁻] absolue n'est pas calibrée en tranche (étape 8) : seules les variations se lisent.

Sorties : results/filtre_redescente.png, results/filtre_redescente_somas.png,
results/filtre_redescente.csv ; avec plusieurs sessions d'un même jeu,
<jeu>/filtre_redescente_sessions.png.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import OptimizeWarning, curve_fit

import config
from config import EXCLURE_APRES_STIM_S, PRE_WINDOW_S
from figures import GRID, INK, INK_2, PALETTE, STIM_BAND, TRACE, bandes_stim, plt, sauver, style
from step6_ratio_corrige import intervalles
from step7_retour import BORNES_TAU, MESURABLE_TAU, ic95

LECTURES = {"lissé": "c_lisse_mM", "temps réel": "c_temps_reel_mM"}
AVANT_S = 30.0            # montré avant le début de la stim
MIN_TRANCHES = 5          # par intervalle, pour la pente de dérive
MIN_POINTS = 6            # pour ajuster une redescente
BLEU, ORANGE = PALETTE[0], PALETTE[1]
DIVERGENTE = LinearSegmentedColormap.from_list("baisse_hausse", [BLEU, "#e9e9e6", "#d64550"])


# ------------------------------------------------------------------ données
def charger(root):
    D = config.dirs(root)
    info = json.loads((D["projections"] / "session.json").read_text())
    fen = [(s["t0"], s["t1"]) for s in info["stims"]]
    t = pd.read_csv(D["results"] / "filtre_par_trame.csv")
    piv = lambda col: t.pivot(index="t_s", columns="roi", values=col).sort_index()
    lues = piv("lue").to_numpy(bool)
    tc = piv("lue").index.to_numpy(float)
    C = {k: np.where(lues, piv(v).to_numpy(float), np.nan) for k, v in LECTURES.items()}
    return D, fen, tc, C, list(piv("lue").columns)


def pentes(tc, C, fen):
    """β (par s) de ln c = a_k + β·t pour chaque soma, moindres carrés intra-intervalle (étape 6)."""
    iv = intervalles(tc, fen, EXCLURE_APRES_STIM_S)
    L = np.log(C)
    beta = np.zeros(C.shape[1])
    for m in range(C.shape[1]):
        num = den = 0.0
        for a, b in iv:
            k = (tc >= a) & (tc < b) & np.isfinite(L[:, m])
            if k.sum() >= MIN_TRANCHES:
                x = tc[k] - tc[k].mean()
                num += (x * (L[k, m] - L[k, m].mean())).sum()
                den += (x * x).sum()
        beta[m] = num / den if den > 0 else 0.0
    return beta


def grille(tc, fen):
    dt = float(np.median(np.diff(tc)))
    suivants = [fen[i + 1][0] if i + 1 < len(fen) else tc[-1] for i in range(len(fen))]
    fenetre = np.floor(min(s - b for s, (_, b) in zip(suivants, fen)) - 1)     # comme l'étape 7
    duree = max(b - a for a, b in fen)
    xg = dt * np.arange(np.ceil(-(duree + AVANT_S) / dt), np.floor(fenetre / dt) + 1)
    return xg, fenetre, duree, dt


def aligner(tc, X, fen, xg):
    """Δc (M, S, G) sur x = t − fin de la stim ; NaN pendant la stim, et partout sans base."""
    M, S = X.shape[1], len(fen)
    Y = np.full((M, S, len(xg)), np.nan)
    for i, (a, b) in enumerate(fen):
        prec = fen[i - 1][1] if i else -np.inf
        suiv = fen[i + 1][0] if i + 1 < len(fen) else np.inf
        av = (tc >= a - PRE_WINDOW_S) & (tc < a)
        x = tc - b
        for m in range(M):
            y = X[:, m]
            base = np.nanmean(y[av]) if np.isfinite(y[av]).any() else np.nan
            if not np.isfinite(base):
                continue
            for part in ((tc > prec) & (tc < a), (tc > b) & (tc < suiv)):
                k = part & np.isfinite(y)
                if k.sum() < 2:
                    continue
                g = (xg >= x[k].min()) & (xg <= x[k].max())
                Y[m, i, g] = np.interp(xg[g], x[k], y[k] - base)
    return Y


# ------------------------------------------------------------------ redescente
def ajuster(x, y):
    """C + A·exp(−x/τ) sur x > 0 ; mesurable si A > 0 et τ loin des bornes."""
    nan = dict(A=np.nan, tau=np.nan, C=np.nan, se_A=np.nan, se_tau=np.nan, se_C=np.nan, mesurable=False)
    ok = np.isfinite(y) & (x > 0)
    if ok.sum() < MIN_POINTS:
        return nan
    xo, yo = x[ok], y[ok]
    s = max(np.abs(yo).max(), 1e-3)
    f = lambda x, A, tau, C: C + A * np.exp(-x / tau)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            p, cv = curve_fit(f, xo, yo, p0=[yo[0] - yo[-1], 15.0, yo[-1]],
                              bounds=([-10 * s, BORNES_TAU[0], -10 * s], [10 * s, BORNES_TAU[1], 10 * s]),
                              maxfev=20000)
    except (RuntimeError, ValueError):
        return nan
    se = np.sqrt(np.diag(cv))
    return dict(A=p[0], tau=p[1], C=p[2], se_A=se[0], se_tau=se[1], se_C=se[2],
                mesurable=bool(p[0] > 0 and MESURABLE_TAU[0] < p[1] < MESURABLE_TAU[1]))


def courbe(p, x):
    return p["C"] + p["A"] * np.exp(-x / p["tau"])


def x_courbe(xg, y, n=200):
    """De la 1re tranche lue après la stim à la fin de la fenêtre : pas d'extrapolation vers x = 0."""
    ok = np.isfinite(y) & (xg > 0)
    return np.linspace(xg[ok][0] if ok.any() else xg[-1], xg[-1], n)


def tau_ratio(D):
    """τ du ratio, somas réunis (étape 7) : {stim : (τ, se, mesurable)}."""
    f = D["results"] / "temps_de_remontee.csv"
    if not f.exists():
        return {}
    r = pd.read_csv(f)
    r = r[r["niveau"] == "somas réunis"]
    return {str(s): (t, e, bool(m)) for s, t, e, m in zip(r["stim"], r["tau_s"], r["tau_se_s"], r["mesurable"])}


def analyser(root):
    D, fen, tc, C, rois = charger(root)
    xg, fenetre, duree, dt = grille(tc, fen)
    out = dict(D=D, nom=root.name, fen=fen, tc=tc, rois=rois, xg=xg, fenetre=fenetre, duree=duree, dt=dt,
               tau_R=tau_ratio(D))
    rows = []
    for lec, X in C.items():
        beta = pentes(tc, X, fen)
        Xc = X * np.exp(-np.outer(tc - tc[0], beta))
        Y = aligner(tc, Xc, fen, xg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            Ym = np.nanmean(Y, 1)                                   # (M, G) moyenne des stims
            meds = [np.nanmedian(Y[:, i], 0) for i in range(len(fen))] + [np.nanmedian(Ym, 0)]
        fits = [ajuster(xg, m) for m in meds]
        f_roi = [ajuster(xg, y) for y in Ym]
        for s, p in zip([str(i + 1) for i in range(len(fen))] + ["moyenne"], fits):
            rows.append(dict(lecture=lec, niveau="médiane des somas", roi="somas", stim=s, **_ligne(p)))
        for r, p, b in zip(rois, f_roi, beta):
            rows.append(dict(lecture=lec, niveau="soma", roi=r, stim="moyenne", **_ligne(p),
                             derive_pct_100s=100 * 100 * b))
        out[lec] = dict(beta=beta, Xc=Xc, Y=Y, Ym=Ym, meds=meds, fits=fits, f_roi=f_roi)
    out["table"] = pd.DataFrame(rows)
    return out


def _ligne(p):
    return dict(A_mM=p["A"], tau_s=p["tau"], tau_se_s=p["se_tau"], equilibre_mM=p["C"], mesurable=p["mesurable"])


# ------------------------------------------------------------------ figures
def _titre_tau(p):
    if not p["mesurable"]:
        return "pas de redescente mesurable"
    return f"τ = {p['tau']:.0f} ± {p['se_tau']:.0f} s,  A = {p['A']:+.2f} mM"


def _limites(curves, pad=0.15):
    v = np.concatenate([np.ravel(c) for c in curves])
    v = v[np.isfinite(v)]
    lo, hi = (np.percentile(v, (1, 99)) if v.size else (-1, 1))
    m = (hi - lo) * pad or 0.1
    return lo - m, hi + m


def figure_redescente(R, out):
    L, T = R["lissé"], R["temps réel"]
    xg, fen, tc, S = R["xg"], R["fen"], R["tc"], len(R["fen"])
    fig = plt.figure(figsize=(max(13, 3.2 * (S + 1)), 9.4))
    gs = fig.add_gridspec(2, S + 1, height_ratios=[1, 1.15], hspace=0.5, wspace=0.25)

    # A — la session entière : c corrigé de la dérive, relatif à la médiane de chaque soma
    ax = fig.add_subplot(gs[0, :])
    style(ax)
    bandes_stim(ax, fen)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rel = lambda X: X - np.nanmedian(X, 0)
        Rl, Rt = rel(L["Xc"]), rel(T["Xc"])
        med, q1, q3 = np.nanmedian(Rl, 1), *np.nanpercentile(Rl, (25, 75), 1)
        ax.plot(tc, np.nanmedian(Rt, 1), color=ORANGE, lw=0.8, alpha=0.75, label="temps réel (médiane)")
    ax.fill_between(tc, q1, q3, color=BLEU, alpha=0.15, lw=0)
    ax.plot(tc, med, color=BLEU, lw=1.8, label="lissé (médiane et écart interquartile)")
    for i, ((a, b), p) in enumerate(zip(fen, L["fits"])):
        if p["mesurable"]:
            # la médiane des écarts n'est pas l'écart des médianes : niveau calé sur la trace après la stim
            xx = x_courbe(xg, L["meds"][i])
            k = (tc > b) & (tc <= b + xx[-1]) & np.isfinite(med)
            niv = np.mean(med[k] - courbe(p, tc[k] - b)) if k.any() else 0.0
            ax.plot(b + xx, niv + courbe(p, xx), color=INK, lw=2.2, zorder=5)
            ax.text(b + 3, niv + p["C"] + p["A"] * 1.08, f"τ = {p['tau']:.0f} s",
                    color=INK, fontsize=9.5, fontweight="bold", va="bottom")
    ax.set_ylim(*_limites([q1, q3, med]))
    ax.set_xlabel("temps (s)", color=INK)
    ax.set_ylabel("[Cl⁻] corrigée de la dérive\n− médiane du soma (mM)", color=INK)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", ncol=2)
    ax.set_title(f"{R['nom']} — [Cl⁻] des {len(R['rois'])} somas après chaque stim (filtre bayésien, étape 8) ; "
                 "noir : C + A·exp(−t/τ) depuis la fin de la stim", loc="left", color=INK, fontsize=11.5)

    # B — chaque stim, alignée sur sa fin, et la moyenne des stims
    panneaux = [(f"stim {i + 1}", L["Y"][:, i], T["Y"][:, i], L["fits"][i], L["meds"][i]) for i in range(S)]
    panneaux.append((f"moyenne des {S} stims", L["Ym"], T["Ym"], L["fits"][-1], L["meds"][-1]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lims = _limites([np.nanpercentile(p[1], (25, 75), 0) for p in panneaux], pad=0.35)
    for j, (titre, Yl, Yt, p, med) in enumerate(panneaux):
        a2 = fig.add_subplot(gs[1, j])
        style(a2)
        a2.axvspan(-R["duree"], 0, color=STIM_BAND, alpha=0.35, lw=0, zorder=0)
        for y in Yl:
            a2.plot(xg, y, color=TRACE, lw=0.5, alpha=0.45)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            q1, q3 = np.nanpercentile(Yl, (25, 75), 0)
            a2.fill_between(xg, q1, q3, color=BLEU, alpha=0.18, lw=0)
            a2.plot(xg, np.nanmedian(Yt, 0), color=ORANGE, lw=1.0, alpha=0.85)
            a2.plot(xg, np.nanmedian(Yl, 0), color=BLEU, lw=2.0)
        if p["mesurable"]:
            xx = x_courbe(xg, med)
            a2.plot(xx, courbe(p, xx), color=INK, lw=1.6, ls="--")
        a2.axhline(0, color=INK_2, lw=0.8, ls=":")
        a2.set_ylim(*lims)
        a2.set_title(f"{titre}\n{_titre_tau(p)}", color=INK, fontsize=9.5)
        a2.set_xlabel("s depuis la fin de la stim", color=INK, fontsize=9)
        if j == 0:
            a2.set_ylabel("Δ[Cl⁻] (mM)\nrelatif aux 10 s avant", color=INK, fontsize=9.5)
        else:
            a2.tick_params(labelleft=False)
    fig.text(0.01, 0.005, "bandes jaunes : stims (tranches non lues) ; gris : chaque soma (lissé) ; bleu : médiane "
             "et écart interquartile (lissé) ; orange : médiane en temps réel ; tirets : C + A·exp(−t/τ) ; "
             "± : erreur type, optimiste (tranches lissées corrélées)", color=INK_2, fontsize=8)
    sauver(fig, out, dpi=170, bbox_inches="tight")


def figure_somas(R, out):
    L, xg = R["lissé"], R["xg"]
    S = len(R["fen"])
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.2), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))

    # A — chaque soma, moyenne de ses stims, triés par la hausse des 10 s après
    ax = axes[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        h = np.nanmean(L["Ym"][:, (xg > 0) & (xg <= 10)], 1)
    ordre = np.argsort(-np.nan_to_num(h, nan=-np.inf))
    v = np.nanpercentile(np.abs(L["Ym"]), 90) or 1.0           # quelques somas très hauts saturent
    im = ax.imshow(L["Ym"][ordre], aspect="auto", cmap=DIVERGENTE, vmin=-v, vmax=v, interpolation="nearest",
                   extent=(xg[0] - R["dt"] / 2, xg[-1] + R["dt"] / 2, len(ordre) + 0.5, 0.5))
    ax.axvline(-R["duree"], color=INK, lw=0.8, ls=":")
    ax.axvline(0, color=INK, lw=0.8, ls=":")
    ax.set_facecolor(GRID)
    ax.set_xlabel("s depuis la fin de la stim (pointillés : la stim)", color=INK)
    ax.set_ylabel("soma (trié par la hausse des 10 s après)", color=INK)
    ax.set_title(f"Δ[Cl⁻] lissée, moyenne des {S} stims", loc="left", color=INK, fontsize=11)
    ax.tick_params(colors=INK_2, length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Δ[Cl⁻] (mM, saturée au 90e centile)", color=INK)
    cb.outline.set_visible(False)

    # B — τ par stim : la [Cl⁻] (filtre) contre le ratio (étape 7)
    ax = axes[1]
    style(ax)
    noms = [str(i + 1) for i in range(S)] + ["moyenne"]
    for k, (etiq, col, dx, get) in enumerate((
            ("[Cl⁻] lissée (médiane des somas, tranches de 3 s)", BLEU, -0.12,
             lambda i, s: (L["fits"][i]["tau"], L["fits"][i]["se_tau"], L["fits"][i]["mesurable"])),
            ("ratio corrigé (étape 7, somas réunis, tranches d'1 s)", INK_2, 0.12,
             lambda i, s: R["tau_R"].get(s, (np.nan, np.nan, False))))):
        for i, s in enumerate(noms):
            t, se, ok = get(i, s)
            if not np.isfinite(t):
                continue
            lo, hi = ic95(t, se)
            ax.errorbar(i + dx, t, yerr=[[t - lo], [min(hi, 1e3) - t]], fmt="o", ms=8, color=col, lw=1.4,
                        mfc=col if ok else "white", mew=1.6, capsize=0, label=etiq if i == 0 else None)
    ax.set_yscale("log")
    ax.set_ylim(1, 400)
    ax.set_xticks(range(len(noms)), [f"stim {s}" if s != "moyenne" else s for s in noms])
    ax.set_ylabel("τ de la redescente (s), IC à 95 %", color=INK)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_title("Temps de redescente, par stim", loc="left", color=INK, fontsize=11)
    ax.text(0.99, 0.02, "creux : non mesurable", transform=ax.transAxes, ha="right", fontsize=8, color=INK_2)

    # C — distribution entre somas
    ax = axes[2]
    style(ax)
    taus = np.array([p["tau"] for p in L["f_roi"] if p["mesurable"]])
    if taus.size:
        ax.hist(taus, bins=np.geomspace(BORNES_TAU[0], BORNES_TAU[1], 24), color=BLEU, alpha=0.75,
                edgecolor="white", lw=1.5)
    ax.set_xscale("log")
    if taus.size:
        ax.axvline(np.median(taus), color=BLEU, lw=1.6, ls="--",
                   label=f"τ médian des somas : {np.median(taus):.0f} s "
                         f"(quartiles {np.percentile(taus, 25):.0f}–{np.percentile(taus, 75):.0f})")
    tr = R["tau_R"].get("moyenne")
    if tr and tr[2]:
        ax.axvline(tr[0], color=INK_2, lw=1.4, ls=":", label=f"ratio, somas réunis (étape 7) : {tr[0]:.0f} s")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_xlabel("τ du soma (s), moyenne de ses stims", color=INK)
    ax.set_ylabel("somas", color=INK)
    ax.set_title(f"Entre somas : {taus.size} / {len(L['f_roi'])} redescentes mesurables",
                 loc="left", color=INK, fontsize=11)
    fig.suptitle(f"{R['nom']} — redescente de la [Cl⁻] soma par soma (filtre bayésien, c lissé)",
                 x=0.01, ha="left", color=INK, fontsize=12)
    fig.tight_layout()
    sauver(fig, out, dpi=170)


def figure_sessions(Rs, out):
    fig, axes = plt.subplots(1, len(Rs), figsize=(3.6 * len(Rs), 4.6), sharey=True, squeeze=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lims = _limites([np.nanpercentile(R["lissé"]["Ym"], (25, 75), 0) for R in Rs], pad=0.2)
    for ax, R in zip(axes[0], Rs):
        L, xg = R["lissé"], R["xg"]
        style(ax)
        ax.axvspan(-R["duree"], 0, color=STIM_BAND, alpha=0.35, lw=0, zorder=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            q1, q3 = np.nanpercentile(L["Ym"], (25, 75), 0)
            ax.fill_between(xg, q1, q3, color=BLEU, alpha=0.18, lw=0)
            ax.plot(xg, np.nanmedian(L["Ym"], 0), color=BLEU, lw=2.0)
        p = L["fits"][-1]
        if p["mesurable"]:
            xx = x_courbe(xg, L["meds"][-1])
            ax.plot(xx, courbe(p, xx), color=INK, lw=1.5, ls="--")
        ax.axhline(0, color=INK_2, lw=0.8, ls=":")
        ax.set_ylim(*lims)
        ax.set_title(f"{R['nom']}\n{len(R['rois'])} somas, {len(R['fen'])} stims ; {_titre_tau(p)}",
                     color=INK, fontsize=8.8)
        ax.set_xlabel("s depuis la fin de la stim", color=INK, fontsize=9)
    axes[0, 0].set_ylabel("Δ[Cl⁻] lissée (mM)\nmoyenne des stims", color=INK)
    fig.suptitle("Redescente de la [Cl⁻] après les stims, session par session — médiane des somas "
                 "et écart interquartile ; tirets : C + A·exp(−t/τ)", x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    sauver(fig, out, dpi=170)


# ------------------------------------------------------------------ main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")] or [None]
    Rs = []
    for a in args:
        root = config.racine([a] if a else [])
        R = analyser(root)
        D = R["D"]
        print(f"{root.name} — {len(R['rois'])} somas, {len(R['fen'])} stims, fenêtre de {R['fenetre']:.0f} s "
              f"après chaque stim, tranches de {R['dt']:.1f} s")
        R["table"].to_csv(D["results"] / "filtre_redescente.csv", index=False)
        figure_redescente(R, D["results"] / "filtre_redescente.png")
        figure_somas(R, D["results"] / "filtre_redescente_somas.png")
        for lec in LECTURES:
            print(f"  {lec} : dérive médiane {100 * 100 * np.median(R[lec]['beta']):+.2f} %/100 s")
            for s, p in zip([f"stim {i + 1}" for i in range(len(R['fen']))] + ["moyenne"], R[lec]["fits"]):
                print(f"    {s:8s} : A = {p['A']:+.3f} mM, τ = {p['tau']:6.1f} ± {p['se_tau']:5.1f} s, "
                      f"C = {p['C']:+.3f} mM{'' if p['mesurable'] else '  (non mesurable)'}")
            n = sum(p["mesurable"] for p in R[lec]["f_roi"])
            t = [p["tau"] for p in R[lec]["f_roi"] if p["mesurable"]]
            print(f"    somas : {n} / {len(R['rois'])} mesurables"
                  + (f", τ médian {np.median(t):.0f} s (quartiles {np.percentile(t, 25):.0f}–"
                     f"{np.percentile(t, 75):.0f})" if t else ""))
        if R["tau_R"]:
            print("  ratio (étape 7) : " + ", ".join(f"{s} τ = {t:.0f} s{'' if ok else ' (n.m.)'}"
                                                   for s, (t, _, ok) in R["tau_R"].items()))
        Rs.append(R)
    parents = {R["D"]["results"].parent.parent for R in Rs}
    trop_court = [R["nom"] for R in Rs if R["fenetre"] < MIN_POINTS * R["dt"]]
    if trop_court:
        print(f"Hors de la figure des sessions (moins de {MIN_POINTS} tranches entre deux stims) : "
              + ", ".join(trop_court))
    Rs = [R for R in Rs if R["nom"] not in trop_court]
    if len(Rs) > 1 and len(parents) == 1:
        figure_sessions(Rs, parents.pop() / "filtre_redescente_sessions.png")


if __name__ == "__main__":
    main()
