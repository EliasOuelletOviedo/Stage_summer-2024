"""Étape 7 — retour du ratio à l'équilibre après chaque stim : mono-exponentielle et
distribution des temps de remontée.

Usage :  python3 step7_retour.py <dossier>   (après l'étape 6)

Données : le ratio ch2 / ch1 corrigé de sa descente graduelle (étape 6, même méthode,
une pente propre à chaque ROI), en tranches de 1 s comptées depuis la fin de la stim,
rapporté aux 10 s qui la précèdent : y = R / R_avant − 1.

Modèle, sur la même fenêtre après chaque stim (FENETRE jusqu'à la stim suivante) :

    y(x) = C + A · exp(−x / τ)        x = temps depuis la fin de la stim

C, l'équilibre, est libre : la fenêtre (~110 s) fait plusieurs τ, il est bien contraint,
et l'imposer égal à la base d'avant allonge τ quand le retour n'est pas complet
(stim 1 des somas : 12 s libre, 31 s imposé). Un retour est **mesurable** si A < 0 (une
baisse, puis une remontée) et si τ ne touche pas ses bornes (0,5 à 290 s).

Trois niveaux :
- somas réunis (photons de tous les ROIs retenus) : un τ par stim, et un sur la moyenne
  des retours ;
- chaque ROI, sur la moyenne de ses retours : la distribution entre cellules ;
- chaque ROI et chaque stim : dans le CSV et en fond de figure, mais un retour isolé est
  trop bruité pour être lu seul.

Témoin « bruit seul » : pour chaque ROI mesurable, on simule son retour avec son A et
son C, mais le τ commun des somas réunis, plus ses propres résidus décalés
circulairement (même bruit, même autocorrélation), et on le réajuste ; N_SIM fois. Si
la distribution observée n'est pas plus large que celle-ci, les écarts de τ entre
cellules ne se distinguent pas du bruit de mesure.

Sorties : results/retour_equilibre_cellpose.png, results/temps_de_remontee.png,
results/temps_de_remontee.csv
"""
import json
import warnings

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy import stats
from scipy.optimize import OptimizeWarning, curve_fit

import config
from config import PRE_WINDOW_S
from figures import CLAIR, bandes_stim, plt, sauver, style
from step6_ratio_corrige import LISSAGE, ajuster_derive

PAS_S = 1.0               # tranches de 1 s (~17 images)
GARDE_S = 0.0             # rien d'écarté après la fin de la stim : la fuite s'arrête net
BORNES_TAU = (0.3, 300.0)
MESURABLE_TAU = (0.5, 290.0)
N_SIM = 200


# ------------------------------------------------------------------ données
def ratio_corrige(c1, c2, tm, fen):
    """Pente de la descente graduelle (étape 6) estimée sur des tranches de 1 s."""
    k = max(int(round(PAS_S / np.median(np.diff(tm)))), 1)
    m = len(tm) // k * k
    b = lambda v: v[:m].reshape(-1, k).sum(1)
    tt = b(tm) / k
    fit = ajuster_derive(tt, b(c2) / np.maximum(b(c1), 1), fen)
    return fit["beta"]


def retours(c1, c2, tm, fen, beta, fenetre):
    """Pour chaque stim : centres des tranches (s depuis la fin), y = R / R_avant − 1."""
    corr = lambda t: np.exp(-beta * (t - tm[0]))
    out = []
    for a, b in fen:
        av = (tm >= a - PRE_WINDOW_S) & (tm < a)
        r_av = c2[av].sum() / c1[av].sum() * corr(tm[av].mean())
        x = tm - b
        m = (x > GARDE_S) & (x <= fenetre)
        i = np.floor((x[m] - GARDE_S) / PAS_S).astype(int)
        n1 = np.bincount(i, weights=c1[m])
        n2 = np.bincount(i, weights=c2[m])
        tc = np.bincount(i, weights=tm[m]) / np.bincount(i)
        xc = GARDE_S + (np.arange(len(n1)) + 0.5) * PAS_S
        out.append((xc, n2 / np.maximum(n1, 1) * corr(tc) / r_av - 1, r_av))
    return out


def ajuster(x, y):
    """C + A·exp(−x/τ) par moindres carrés ; dict avec erreurs types et résidus."""
    nan = dict(A=np.nan, tau=np.nan, C=np.nan, se_A=np.nan, se_tau=np.nan, se_C=np.nan,
               mesurable=False, residus=None)
    ok = np.isfinite(y)
    if ok.sum() < 10:
        return nan
    f = lambda x, A, tau, C: C + A * np.exp(-x / tau)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            p, cv = curve_fit(f, x[ok], y[ok], p0=[-0.02, 15.0, 0.0],
                              bounds=([-1, BORNES_TAU[0], -0.2], [1, BORNES_TAU[1], 0.2]),
                              maxfev=20000)
    except (RuntimeError, ValueError):
        return nan
    se = np.sqrt(np.diag(cv))
    res = np.full(len(y), np.nan)
    res[ok] = y[ok] - f(x[ok], *p)
    return dict(A=p[0], tau=p[1], C=p[2], se_A=se[0], se_tau=se[1], se_C=se[2],
                mesurable=bool(p[0] < 0 and MESURABLE_TAU[0] < p[1] < MESURABLE_TAU[1]),
                residus=res)


def ic95(t, se):
    """Intervalle à 95 % de τ, en log-normal (τ est positif et asymétrique)."""
    s = 1.96 * se / t
    return t * np.exp(-s), t * np.exp(s)


def temoin_bruit(x, fits_roi, tau_commun, rng):
    """τ réajustés quand chaque ROI mesurable suit le τ commun, avec son propre bruit."""
    sims = []
    for r in fits_roi:
        if not r["mesurable"]:
            continue
        res = r["residus"]
        ok = np.isfinite(res)
        for _ in range(N_SIM):
            y = r["C"] + r["A"] * np.exp(-x / tau_commun)
            y[ok] += np.roll(res[ok], rng.integers(ok.sum()))
            y[~ok] = np.nan
            s = ajuster(x, y)
            if s["mesurable"]:
                sims.append(s["tau"])
    return np.array(sims)


# ------------------------------------------------------------------ figures
def figure_somas(tm, c1, c2, fen, beta, rets, fits, moyen, etiquette, out, T=CLAIR):
    pal = T["palette"]
    fig = plt.figure(figsize=(13, 8.4), facecolor=T["surface"])
    gs = fig.add_gridspec(2, len(fen) + 1, height_ratios=[1.25, 1], hspace=0.42, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    style(ax, T)
    bandes_stim(ax, fen, T)
    lis = lambda v: ndi.uniform_filter1d(v.astype(float), LISSAGE)
    corr = np.exp(-beta * (tm - tm[0]))
    rc = c2 / np.maximum(c1, 1) * corr
    ax.plot(tm, rc, color=T["ink_2"], lw=0.3, alpha=0.3)
    ax.plot(tm, lis(c2) / lis(c1) * corr, color=T["ink"], lw=1.2)
    for i, ((a, b), (xc, y, r_av), p) in enumerate(zip(fen, rets, fits)):
        xx = np.linspace(GARDE_S, xc[-1], 300)
        ax.plot(b + xx, r_av * (1 + p["C"] + p["A"] * np.exp(-xx / p["tau"])),
                color=pal[i % len(pal)], lw=2.6, zorder=5)
        ax.text(b + 6, r_av * (1 + p["A"] * 1.35), f"τ = {p['tau']:.1f} s",
                color=pal[i % len(pal)], fontsize=10, fontweight="bold", va="top")
    v = rc[np.isfinite(rc)]
    lo, hi = np.percentile(v, (0.5, 99.5))
    ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.1 * (hi - lo))
    ax.set_ylabel("ratio ch2 / ch1\n(corrigé)", color=T["ink"], fontsize=10)
    ax.set_xlabel("temps (s)", color=T["ink"], fontsize=10)
    ax.set_title(f"{etiquette[0].upper() + etiquette[1:]} — retour à l'équilibre après chaque stim : "
                 "C + A·exp(−t/τ), t depuis la fin de la stim", color=T["ink"], fontsize=12, loc="left")

    panneaux = [(f"stim {i + 1}", r[0], r[1], p, pal[i % len(pal)])
                for i, (r, p) in enumerate(zip(rets, fits))]
    panneaux.append((f"moyenne des {len(rets)}", rets[0][0], np.nanmean([r[1] for r in rets], axis=0), moyen,
                     T["ink"]))
    lims = np.nanpercentile(np.concatenate([p[2] for p in panneaux]) * 100, (1, 99))
    for j, (titre, xc, y, p, c) in enumerate(panneaux):
        a2 = fig.add_subplot(gs[1, j])
        style(a2, T)
        a2.scatter(xc, y * 100, s=7, color=c, alpha=T["alpha_points"], lw=0)
        xx = np.linspace(GARDE_S, xc[-1], 300)
        a2.plot(xx, (p["C"] + p["A"] * np.exp(-xx / p["tau"])) * 100, color=T["ink"], lw=1.8)
        a2.axhline(0, color=T["ink_2"], lw=0.8, ls="--")
        a2.axhline(p["C"] * 100, color=c, lw=0.8, ls=":")
        a2.set_ylim(lims[0] - 0.4, lims[1] + 0.4)
        a2.set_title(f"{titre}\nτ = {p['tau']:.1f} ± {p['se_tau']:.1f} s", color=T["ink"], fontsize=10)
        a2.text(0.97, 0.05, f"A = {p['A'] * 100:+.2f} %\néquilibre {p['C'] * 100:+.2f} %",
                transform=a2.transAxes, ha="right", va="bottom", fontsize=8, color=T["ink_2"])
        a2.set_xlabel("s depuis la fin", color=T["ink"], fontsize=9)
        if j == 0:
            a2.set_ylabel("R / R$_{avant}$ − 1 (%)", color=T["ink"], fontsize=9.5)
    fig.text(0.01, 0.005, "points : tranches de 1 s ; tirets : niveau des 10 s avant la stim ; "
             "pointillés : équilibre ajusté (C) ; ± : erreur type", color=T["ink_2"], fontsize=8.5)
    sauver(fig, out, dpi=170, T=T, bbox_inches="tight")


def figure_distribution(df, sims, tau_commun, se_commun, ks, out, T=CLAIR):
    pal = T["palette"]
    roi = df[(df["niveau"] == "roi") & (df["stim"] == "moyenne")]
    mes = roi[roi["mesurable"]].sort_values("tau_s")
    som = df[df["niveau"] == "somas réunis"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.4), facecolor=T["surface"],
                             gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
    lo_c, hi_c = ic95(tau_commun, se_commun)
    xl = (2, 300)

    # A — chaque soma, moyenne de ses retours, avec son intervalle
    ax = axes[0]
    style(ax, T)
    ax.axvspan(lo_c, hi_c, color=T["ink"], alpha=0.10, lw=0)
    ax.axvline(tau_commun, color=T["ink"], lw=1.4)
    lo, hi = ic95(mes["tau_s"].to_numpy(), mes["tau_se_s"].to_numpy())
    yy = np.arange(len(mes))
    ax.hlines(yy, np.clip(lo, *xl), np.clip(hi, *xl), color=pal[0], lw=1.6, alpha=0.6)
    ax.scatter(mes["tau_s"], yy, color=pal[0], s=26, zorder=3)
    ax.set_yticks(yy)
    ax.set_yticklabels(mes["roi"], fontsize=7.5, color=T["ink_2"])
    ax.set_xscale("log")
    ax.set_xlim(*xl)
    ax.set_xlabel("τ de remontée (s, échelle log)", color=T["ink"])
    ax.set_ylabel("ROI", color=T["ink"])
    n_non = len(roi) - len(mes)
    n_st = df.loc[df["stim"] != "moyenne", "stim"].nunique()
    ax.set_title(f"A · Chaque soma, moyenne de ses {n_st} retours", color=T["ink"], fontsize=11, loc="left")
    ax.text(0.98, 0.02, f"{len(mes)} mesurables sur {len(roi)}\n({n_non} sans baisse nette ou\nτ hors bornes)\n"
            "trait : IC 95 %", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=T["ink_2"])

    # B — observé contre bruit seul
    ax = axes[1]
    style(ax, T)
    obs = np.sort(mes["tau_s"].to_numpy())
    sim = np.sort(sims)
    ax.step(sim, np.arange(1, len(sim) + 1) / len(sim), where="post", color=T["neutre"], lw=2.2,
            label=f"bruit seul : tous à {tau_commun:.1f} s\n({len(sim):,} retours simulés)")
    ax.step(obs, np.arange(1, len(obs) + 1) / len(obs), where="post", color=pal[0], lw=2.2,
            label=f"observé ({len(obs)} somas)")
    for q, ls in ((0.25, ":"), (0.5, "--"), (0.75, ":")):
        ax.axhline(q, color=T["repere"], lw=0.8, ls=ls)
    ax.set_xscale("log")
    ax.set_xlim(*xl)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("τ de remontée (s, échelle log)", color=T["ink"])
    ax.set_ylabel("fraction des somas (cumulée)", color=T["ink"])
    q_o, q_s = np.percentile(obs, [25, 50, 75]), np.percentile(sim, [25, 50, 75])
    ax.set_title("B · Plus étalé que le bruit de mesure ?", color=T["ink"], fontsize=11, loc="left")
    ax.text(0.98, 0.02, f"médiane [Q1–Q3]\nobservé : {q_o[1]:.1f} [{q_o[0]:.1f}–{q_o[2]:.1f}] s\n"
            f"bruit seul : {q_s[1]:.1f} [{q_s[0]:.1f}–{q_s[2]:.1f}] s\nKolmogorov-Smirnov p = {ks:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=T["ink_2"])
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", labelcolor=T["ink"])

    # C — stim par stim
    ax = axes[2]
    style(ax, T)
    rs = df[(df["niveau"] == "roi") & (df["stim"] != "moyenne") & df["mesurable"]]
    rng = np.random.default_rng(1)
    stims = sorted(rs["stim"].unique(), key=int)
    etiquettes = []
    for i in map(int, stims):
        v = rs.loc[rs["stim"] == str(i), "tau_s"].to_numpy()
        c = pal[(i - 1) % len(pal)]
        if len(v) >= 3:
            ax.boxplot([v], positions=[i - 0.08], widths=0.32, showfliers=False, whis=0,
                       showcaps=False, patch_artist=True, manage_ticks=False,
                       boxprops=dict(facecolor=c, alpha=0.22, edgecolor=c),
                       medianprops=dict(color=c, lw=2.2), whiskerprops=dict(lw=0))
        ax.scatter(np.full(len(v), i - 0.08) + rng.uniform(-0.13, 0.13, len(v)), v, s=14,
                   color=c, alpha=T["alpha_points"], lw=0, zorder=3)
        p = som[som["stim"] == str(i)].iloc[0]
        lo, hi = ic95(p["tau_s"], p["tau_se_s"])
        ax.errorbar(i + 0.3, p["tau_s"], yerr=[[p["tau_s"] - lo], [hi - p["tau_s"]]], fmt="D",
                    color=T["ink"], ms=6, capsize=3, zorder=5)
        n_tot = ((df["niveau"] == "roi") & (df["stim"] == str(i))).sum()
        etiquettes.append(f"stim {i}\n{len(v)}/{n_tot}")
    ax.set_yscale("log")
    ax.set_ylim(*xl)
    ax.set_xlim(0.5, len(stims) + 0.6)
    ax.set_xticks([int(i) for i in stims])
    ax.set_xticklabels(etiquettes, color=T["ink"], fontsize=9)
    ax.set_ylabel("τ de remontée (s, échelle log)", color=T["ink"])
    ax.set_title("C · Stim par stim", color=T["ink"], fontsize=11, loc="left")
    ax.text(0.98, 0.98, "points et boîte : chaque soma\n(retour isolé, très bruité)\n"
            "◆ : somas réunis ± IC 95 %\nsous l'axe : retours mesurables", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.5, color=T["ink_2"])
    fig.suptitle("Temps de remontée du ratio après les stims — somas Cellpose ; "
                 f"trait vertical et bande (A) : somas réunis, {tau_commun:.1f} ± {se_commun:.1f} s "
                 "(IC 95 %)", color=T["ink"], fontsize=12.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    sauver(fig, out, dpi=170, T=T)


# ------------------------------------------------------------------ calcul
def calculer(root):
    """Tous les ajustements et le témoin « bruit seul ». Rend un dict pour les figures."""
    D = config.dirs(root)
    info = json.loads((D["projections"] / "session.json").read_text())
    f = np.load(D["projections"] / "flim.npz")
    d = np.load(D["mesures"] / "declins.npz")
    kept = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
    t = f["trames_t"]
    n = len(t) - 1                      # dernière image coupée par l'arrêt : écartée
    tm = (t + np.median(np.diff(t)) / 2)[:n]
    fen = [(s["t0"], s["t1"]) for s in info["stims"]]
    if not fen:
        raise SystemExit("aucune stim dans session.json")
    if fen[0][0] < PRE_WINDOW_S:
        raise SystemExit(f"pas de ligne de base (1re stim à {fen[0][0]:.1f} s) : étape 7 sautée")
    suivants = [fen[i + 1][0] if i + 1 < len(fen) else tm[-1] for i in range(len(fen))]
    fenetre = np.floor(min(s - b for s, (_, b) in zip(suivants, fen)) - 1)   # commune à toutes
    P = d["par_image_cellpose"][:, :n].astype(float)
    rows = []

    # --- somas réunis
    pc = P[kept].sum(0)
    beta = ratio_corrige(pc[:, 0], pc[:, 1], tm, fen)
    rets = retours(pc[:, 0], pc[:, 1], tm, fen, beta, fenetre)
    fits = [ajuster(x, y) for x, y, _ in rets]
    moyen = ajuster(rets[0][0], np.nanmean([r[1] for r in rets], axis=0))
    for i, p in enumerate(fits + [moyen]):
        rows.append(dict(niveau="somas réunis", roi="somas", stim=str(i + 1) if i < len(fits) else "moyenne",
                         beta_pct_100s=beta * 1e4, A_pct=p["A"] * 100, tau_s=p["tau"], tau_se_s=p["se_tau"],
                         equilibre_pct=p["C"] * 100, mesurable=p["mesurable"]))

    # --- chaque ROI
    fits_roi = []
    for l in kept:
        b_l = ratio_corrige(P[l, :, 0], P[l, :, 1], tm, fen)
        r_l = retours(P[l, :, 0], P[l, :, 1], tm, fen, b_l, fenetre)
        for i, (x, y, _) in enumerate(r_l):
            p = ajuster(x, y)
            rows.append(dict(niveau="roi", roi=str(l), stim=str(i + 1), beta_pct_100s=b_l * 1e4,
                             A_pct=p["A"] * 100, tau_s=p["tau"], tau_se_s=p["se_tau"],
                             equilibre_pct=p["C"] * 100, mesurable=p["mesurable"]))
        p = ajuster(r_l[0][0], np.nanmean([r[1] for r in r_l], axis=0))
        fits_roi.append(p)
        rows.append(dict(niveau="roi", roi=str(l), stim="moyenne", beta_pct_100s=b_l * 1e4,
                         A_pct=p["A"] * 100, tau_s=p["tau"], tau_se_s=p["se_tau"],
                         equilibre_pct=p["C"] * 100, mesurable=p["mesurable"]))
    df = pd.DataFrame(rows)
    mes = df[(df["niveau"] == "roi") & (df["stim"] == "moyenne") & df["mesurable"]]
    sims = temoin_bruit(rets[0][0], fits_roi, moyen["tau"], np.random.default_rng(0))
    ks = stats.ks_2samp(mes["tau_s"], sims).pvalue
    return dict(D=D, kept=kept, tm=tm, fen=fen, fenetre=fenetre, pc=pc, beta=beta, rets=rets,
                fits=fits, moyen=moyen, df=df, sims=sims, ks=ks,
                etiquette=f"{len(kept)} somas Cellpose réunis")


def dessiner(R, dossier, T=CLAIR):
    """Les deux figures de l'étape, dans `dossier`, avec le thème `T`."""
    figure_somas(R["tm"], R["pc"][:, 0], R["pc"][:, 1], R["fen"], R["beta"], R["rets"], R["fits"],
                 R["moyen"], R["etiquette"], dossier / "retour_equilibre_cellpose.png", T)
    figure_distribution(R["df"], R["sims"], R["moyen"]["tau"], R["moyen"]["se_tau"], R["ks"],
                        dossier / "temps_de_remontee.png", T)


# ------------------------------------------------------------------ main
def main():
    root = config.racine()
    R = calculer(root)
    D, df, moyen = R["D"], R["df"], R["moyen"]
    print(f"{root.name} — fenêtre de {R['fenetre']:.0f} s après chaque stim, tranches de {PAS_S:.0f} s")
    print(f"  somas réunis (β = {R['beta'] * 1e4:+.2f} %/100 s) :")
    for i, p in enumerate(R["fits"] + [moyen]):
        nom = f"stim {i + 1}" if i < len(R["fits"]) else f"moyenne des {len(R['fits'])}"
        print(f"    {nom:14s} τ = {p['tau']:5.1f} ± {p['se_tau']:4.1f} s   A = {p['A'] * 100:+.2f} %   "
              f"équilibre {p['C'] * 100:+.2f} %")
    df.to_csv(D["results"] / "temps_de_remontee.csv", index=False)
    mes = df[(df["niveau"] == "roi") & (df["stim"] == "moyenne") & df["mesurable"]]
    q = np.percentile(mes["tau_s"], [25, 50, 75])
    qs = np.percentile(R["sims"], [25, 50, 75])
    print(f"\n  chaque soma (moyenne de ses {len(R['fits'])} retours) : {len(mes)}/{len(R['kept'])} mesurables, "
          f"τ médian {q[1]:.1f} s [{q[0]:.1f}–{q[2]:.1f}]")
    print(f"  bruit seul (tous à {moyen['tau']:.1f} s) : τ médian {qs[1]:.1f} s [{qs[0]:.1f}–{qs[2]:.1f}], "
          f"Kolmogorov-Smirnov p = {R['ks']:.2f}")
    rs = df[(df["niveau"] == "roi") & (df["stim"] != "moyenne")]
    print(f"  retours isolés (ROI x stim) : {int(rs['mesurable'].sum())}/{len(rs)} mesurables, "
          f"erreur type médiane {np.nanmedian(rs.loc[rs['mesurable'], 'tau_se_s'] / rs.loc[rs['mesurable'], 'tau_s']) * 100:.0f} % de τ")
    dessiner(R, D["results"])
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
