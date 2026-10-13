"""Étape 5 — ratio FRET et durées de vie dans le temps, réponse à chaque stimulation.

Usage :  python3 step5_analyse.py <dossier>

Pour chaque ROI et chaque tranche d'environ ANALYSE_BLOC_S (3 s : les blocs FLIMApp
regroupés — tels quels à 50 images, par 3 à 17 images) :
- photons ch1 (CFP) et ch2 (YFP), ratio = ch2 / ch1, valeurs brutes ;
- τ ch1 et τ ch2 : mono-exponentielle reconvoluée (flim_fit.Modele), décalage de l'IRF
  fixé par ROI ; τm, le temps moyen sans ajustement.

Familles de ROIs : `cellpose` (triés à l'étape 3), `flimapp` (dessinés pendant
l'expérience, F1…F8), `neuropile` (hors somas) et `champ` (toute l'image).

R0, τ0 : médiane avant la première stim (le 1er bloc écarté). Les blocs qui touchent
une stim sont gardés dans le CSV (`stim = True`) mais exclus des figures et des
mesures : la stim y ajoute ~4 % de photons dans ch2.

Réponse à une stim : moyenne des blocs centrés dans les 10 s après la fin, moins celle
des blocs centrés dans les 10 s avant le début. Le **témoin** applique exactement la
même mesure à une fausse stim de même durée, placée au milieu de chaque intervalle
sans stim : il dit ce que la dérive et le bruit produisent seuls.

Sorties dans results/ : mesures_blocs.csv, resume_par_roi.csv, avant_apres_stim.csv,
resume_stims.csv, images_fixes_C2.csv et les figures.
"""
import json
import warnings

import numpy as np
import pandas as pd
from scipy import stats

import config
from config import (ANALYSE_BLOC_S, BASE_SKIP, CANAUX, MIN_VALID_FRAC, NOMS_CANAUX,
                    PRE_WINDOW_S, RESPONSE_WINDOW_S)
from figures import (INK, INK_2, PALETTE, SURFACE, TRACE, bandes_stim, plt, sauver, style,
                     ylim_robuste)
from flim_fit import Modele

COL_TEMOIN = "#8d8d86"


# ------------------------------------------------------------------ mesures
def ajuster_roi(H, modeles):
    """H [nb, 2, 256] -> dict de séries [nb] par canal : τ, τm, χ² réduit, photons."""
    out = {}
    for k in CANAUX:
        M = modeles[k]
        h = H[:, k - 1].astype(float)
        tot = h.sum(0)
        tau_g, s, _, chi_g = M.ajuster(tot)
        r = [M.ajuster(h[b], s=s, tau0=tau_g) for b in range(len(h))]
        out[f"photons_ch{k}"] = h.sum(1)
        out[f"tau_ch{k}"] = np.array([x[0] for x in r])
        out[f"fond_ch{k}"] = np.array([x[2] for x in r])
        out[f"chi2_ch{k}"] = np.array([x[3] for x in r])
        out[f"tm_ch{k}"] = np.array([M.temps_moyen(h[b]) for b in range(len(h))])
        out[f"_global_ch{k}"] = (tau_g, s, chi_g)
    return out


def tranches(info):
    """Regroupement des blocs FLIMApp en tranches d'environ ANALYSE_BLOC_S : (k, n, t0, t1,
    stim) — k blocs par tranche, n tranches (les blocs en trop à la fin sont laissés)."""
    b0 = np.array([b["t0"] for b in info["blocs"]])
    b1 = np.array([b["t1"] for b in info["blocs"]])
    st = np.array([b["stim"] for b in info["blocs"]])
    k = max(1, int(round(ANALYSE_BLOC_S / np.median(b1 - b0))))
    n = len(b0) // k
    return k, n, b0[:n * k:k], b1[k - 1:n * k:k], st[:n * k].reshape(n, k).any(1)


def regrouper(x, k, n, f=np.sum):
    """Somme (ou `f`) de k blocs consécutifs le long du premier axe."""
    x = np.asarray(x)
    return f(x[:n * k].reshape(n, k, *x.shape[1:]), axis=1)


def phases(t0, t1, stim, fen):
    pulse = np.array([sum(e <= a for _, e in fen) for a in t0])
    premier = fen[0][0] if fen else np.inf
    phase = np.where(stim, "pendant", np.where(t1 <= premier, "avant", "apres"))
    return t0, t1, (t0 + t1) / 2, stim, phase, pulse


def temoins(fen, t_fin):
    """Fausses stims de même durée, au milieu de chaque intervalle sans stim, gardées
    si leurs fenêtres avant / après ne touchent aucune vraie stim."""
    out = []
    for i, (a, b) in enumerate(fen):
        suivant = fen[i + 1][0] if i + 1 < len(fen) else t_fin
        mid = (b + suivant) / 2
        ta, tb = mid - (b - a) / 2, mid + (b - a) / 2
        if ta - PRE_WINDOW_S >= b and tb + RESPONSE_WINDOW_S <= suivant:
            out.append((ta, tb))
    return out


def tableau(root, D, info):
    S = config.session(root)
    I, dt = config.irf(root)
    modeles = {k: Modele(I[k], dt) for k in CANAUX}
    d = np.load(D["mesures"] / "declins.npz")
    f = np.load(D["projections"] / "flim.npz")
    sel = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
    fen = [(s["t0"], s["t1"]) for s in info["stims"]]
    k, nb, b0, b1, bst = tranches(info)
    t0, t1, tc, stim, phase, pulse = phases(b0, b1, bst, fen)
    if k > 1:
        print(f"  blocs FLIMApp de {np.median(b1 - b0) / k:.2f} s regroupés par {k} : "
              f"{nb} tranches de {np.median(b1 - b0):.1f} s")
    G = lambda x: regrouper(x, k, nb)

    rois = [("cellpose", str(l), G(d["cellpose"][l]), regrouper(d["frac_cellpose"][:, l], k, nb, np.min))
            for l in sel]
    rois += [("flimapp", f"F{i + 1}", G(d["flimapp"][i]), np.ones(nb)) for i in range(len(d["flimapp"]))]
    rois += [("neuropile", "neuropile", G(d["neuropile"]), np.ones(nb)),
             ("champ", "champ", G(np.stack([f["declin_bloc_ch1"], f["declin_bloc_ch2"]], 1)), np.ones(nb))]
    rows, glob = [], []
    for fam, roi, H, frac in rois:
        r = ajuster_roi(H, modeles)
        ok = frac >= MIN_VALID_FRAC
        for k in CANAUX:
            tg, s, chi = r.pop(f"_global_ch{k}")
            glob.append(dict(famille=fam, roi=roi, canal=k, tau_global_ns=tg, decalage_irf=s,
                             chi2_global=chi))
        df = pd.DataFrame(r)
        df.loc[~ok, :] = np.nan
        with np.errstate(invalid="ignore", divide="ignore"):
            df["ratio"] = df["photons_ch2"] / df["photons_ch1"]
        df.insert(0, "famille", fam)
        df.insert(1, "roi", roi)
        df.insert(2, "bloc", np.arange(1, nb + 1))
        df.insert(3, "t_s", np.round(tc, 3))
        df.insert(4, "stim", stim)
        df.insert(5, "phase", phase)
        df.insert(6, "pulse", pulse)
        df["frac_valide"] = frac
        rows.append(df)
        if fam != "cellpose":
            print(f"    {fam:9s} {roi:9s} τ ch1 {np.nanmedian(r['tau_ch1']):.3f} ns, "
                  f"τ ch2 {np.nanmedian(r['tau_ch2']):.3f} ns, ratio {np.nanmedian(df['ratio']):.2f}")
    df = pd.concat(rows, ignore_index=True)

    # --- normalisation par la base avant la 1re stim
    base = (df["phase"] == "avant") & (df["bloc"] > BASE_SKIP)
    g = [df["famille"], df["roi"]]
    for col, ref, out, mult in [("ratio", "R0", "dR_R0_pct", None),
                                ("tau_ch1", "tau0_ch1", "dtau_ch1_ps", 1000),
                                ("tau_ch2", "tau0_ch2", "dtau_ch2_ps", 1000),
                                ("tm_ch1", "tm0_ch1", "dtm_ch1_ps", 1000),
                                ("photons_ch1", "N0_ch1", "dN_ch1_pct", None),
                                ("photons_ch2", "N0_ch2", "dN_ch2_pct", None)]:
        df[ref] = df[col].where(base).groupby(g).transform("median")
        df[out] = (df[col] - df[ref]) * mult if mult else (df[col] / df[ref] - 1) * 100
    return df, pd.DataFrame(glob), fen, S


MESURES = {"dR_R0_pct": ("ΔR / R$_0$", "%"), "dtau_ch1_ps": ("Δτ CFP (ch1)", "ps"),
           "dtau_ch2_ps": ("Δτ YFP (ch2)", "ps"), "dtm_ch1_ps": ("Δτm CFP (ch1)", "ps"),
           "dN_ch1_pct": ("Δ photons CFP (ch1)", "%"), "dN_ch2_pct": ("Δ photons YFP (ch2)", "%")}


def avant_apres(df, fen, tem):
    """Une ligne par ROI x stim (vraie ou témoin) : moyennes avant / après et leur écart."""
    rows = []
    ok = ~df["stim"]
    for type_, wins in (("stim", fen), ("témoin", tem)):
        for i, (a, b) in enumerate(wins):
            av = df[ok & (df["t_s"] >= a - PRE_WINDOW_S) & (df["t_s"] < a)]
            ap = df[ok & (df["t_s"] > b) & (df["t_s"] <= b + RESPONSE_WINDOW_S)]
            ma = av.groupby(["famille", "roi"])[list(MESURES)].mean()
            mp = ap.groupby(["famille", "roi"])[list(MESURES)].mean()
            w = pd.concat({"avant": ma, "apres": mp}, axis=1)
            w.columns = [f"{m}_{p}" for p, m in w.columns]
            for m in MESURES:
                w[f"{m}_ecart"] = w[f"{m}_apres"] - w[f"{m}_avant"]
            w = w.reset_index()
            w.insert(0, "type", type_)
            w.insert(1, "numero", i + 1)
            w.insert(2, "t0_s", round(a, 2))
            w.insert(3, "t1_s", round(b, 2))
            w["n_blocs_avant"], w["n_blocs_apres"] = av["bloc"].nunique(), ap["bloc"].nunique()
            rows.append(w)
    return pd.concat(rows, ignore_index=True)


def recuperation(df, fen, t_fin):
    """Part de la réponse revenue juste avant la stim suivante (10 s avant), par ROI."""
    rows = []
    ok = ~df["stim"]
    for i, (a, b) in enumerate(fen):
        suivant = fen[i + 1][0] if i + 1 < len(fen) else t_fin
        av = df[ok & (df["t_s"] >= a - PRE_WINDOW_S) & (df["t_s"] < a)]
        ap = df[ok & (df["t_s"] > b) & (df["t_s"] <= b + RESPONSE_WINDOW_S)]
        fin = df[ok & (df["t_s"] >= suivant - PRE_WINDOW_S) & (df["t_s"] < suivant)]
        for m in ("dR_R0_pct", "dtau_ch1_ps"):
            x = [d.groupby(["famille", "roi"])[m].mean() for d in (av, ap, fin)]
            w = pd.DataFrame({"avant": x[0], "apres": x[1], "fin": x[2]}).reset_index()
            w["mesure"], w["stim"] = m, i + 1
            w["delai_fin_s"] = round(suivant - b, 1)
            rows.append(w)
    return pd.concat(rows, ignore_index=True)


def resume(pp):
    rows = []
    for (typ, fam, num), g in pp.groupby(["type", "famille", "numero"], sort=False):
        row = {"type": typ, "famille": fam, "numero": num, "n": len(g)}
        for m in MESURES:
            e = g[f"{m}_ecart"].dropna()
            row[f"{m}_median"] = e.median()
            row[f"{m}_n_baisse"] = int((e < 0).sum())
            row[f"{m}_p_wilcoxon"] = (stats.wilcoxon(e).pvalue if len(e) >= 6 and e.abs().sum() > 0
                                      else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ figures
def _matrice(sub, col):
    w = sub.pivot(index="roi", columns="bloc", values=col)
    w = w.loc[sorted(w.index, key=lambda r: (len(r), r))]
    t = sub.drop_duplicates("bloc").set_index("bloc")["t_s"].sort_index()
    st = sub.drop_duplicates("bloc").set_index("bloc")["stim"].sort_index().to_numpy(bool)
    M = w.to_numpy(float, copy=True)
    M[:, st] = np.nan
    return M, t.to_numpy(float), w.index.to_numpy()


def _traits(ax, sub, col, color):
    M, t, _ = _matrice(sub, col)
    for row in M:
        ax.plot(t, row, color=TRACE, linewidth=0.6, alpha=0.6, zorder=1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q1, med, q3 = (np.nanpercentile(M, q, axis=0) for q in (25, 50, 75))
    ax.fill_between(t, q1, q3, color=color, alpha=0.22, lw=0, zorder=2)
    ax.plot(t, med, color=color, linewidth=2.2, zorder=3)
    ylim_robuste(ax, M)


def figure_traces(df, fen, colonnes, titre, note, out):
    sub = df[df["famille"] == "cellpose"]
    fig, axes = plt.subplots(len(colonnes), 1, figsize=(10, 2.9 * len(colonnes) + 0.6),
                             sharex=True, facecolor=SURFACE)
    for ax, (col, lab, color) in zip(axes, colonnes):
        style(ax)
        bandes_stim(ax, fen)
        _traits(ax, sub, col, color)
        if col.startswith("d"):
            ax.axhline(0, color=INK_2, linewidth=0.9, linestyle="--", zorder=2)
        ax.set_ylabel(lab, color=INK, fontsize=10.5)
    axes[-1].set_xlabel("temps (s)", color=INK, fontsize=11)
    axes[0].set_title(titre, color=INK, fontsize=13, pad=34, loc="left")
    axes[0].text(0, 1.02, note, transform=axes[0].transAxes, color=INK_2, fontsize=9.5,
                 va="bottom", linespacing=1.5)
    fig.tight_layout()
    sauver(fig, out, dpi=190)


def figure_decomposition(df, fen, tem, out):
    """Ce qui bouge : ratio, τ CFP, photons de chaque canal — somas, neuropile, champ."""
    cols = ["dR_R0_pct", "dtau_ch1_ps", "dN_ch1_pct", "dN_ch2_pct"]
    fams = [("cellpose", "somas Cellpose (médiane, écart interquartile)", PALETTE[0]),
            ("neuropile", "neuropile", PALETTE[1]), ("champ", "champ entier", INK_2)]
    fig, axes = plt.subplots(len(cols), 1, figsize=(10, 10.5), sharex=True, facecolor=SURFACE)
    for ax, col in zip(axes, cols):
        style(ax)
        bandes_stim(ax, fen)
        for a, b in tem:
            ax.axvspan(a, b, color=COL_TEMOIN, alpha=0.10, lw=0, zorder=0)
        allv = []
        for fam, lab, c in fams:
            sub = df[df["famille"] == fam]
            M, t, _ = _matrice(sub, col)
            allv.append(M)
            with np.errstate(invalid="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med = np.nanmedian(M, axis=0)
                if len(M) > 1:
                    q1, q3 = np.nanpercentile(M, 25, axis=0), np.nanpercentile(M, 75, axis=0)
                    ax.fill_between(t, q1, q3, color=c, alpha=0.18, lw=0)
            ax.plot(t, med, color=c, lw=1.8, label=lab)
        ax.axhline(0, color=INK_2, lw=0.9, ls="--")
        nom, unite = MESURES[col]
        ax.set_ylabel(f"{nom} ({unite})", color=INK, fontsize=10)
        ylim_robuste(ax, np.concatenate([np.ravel(np.nanmedian(m, 0)) for m in allv]), pad=0.25)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left", ncol=3)
    axes[-1].set_xlabel("temps (s)", color=INK, fontsize=11)
    axes[0].set_title("Ce qui change à chaque stimulation — par rapport à la base avant la "
                      "1re stim", color=INK, fontsize=12.5, loc="left", pad=26)
    axes[0].text(0, 1.02, "bandes jaunes : stimulations (blocs masqués) ; bandes grises : "
                 "témoins sans stim, mêmes durées", transform=axes[0].transAxes, color=INK_2,
                 fontsize=9.5, va="bottom")
    fig.tight_layout()
    sauver(fig, out, dpi=180)


def figure_alignee(df, fen, tem, out):
    """Chaque stim (et chaque témoin) alignée sur son début, rapportée aux 10 s d'avant."""
    sub = df[df["famille"] == "cellpose"]
    cols = [("dR_R0_pct", "ΔR / R$_0$ (points)"), ("dtau_ch1_ps", "Δτ CFP (ps)")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURFACE)
    for ax, (col, lab) in zip(axes, cols):
        style(ax)
        M, t, _ = _matrice(sub, col)
        dur = fen[0][1] - fen[0][0] if fen else 10
        ax.axvspan(0, dur, color="#f4d35e", alpha=0.35, lw=0)
        for typ, wins in (("témoin", tem), ("stim", fen)):
            for i, (a, b) in enumerate(wins):
                rel = t - a
                av = (rel >= -PRE_WINDOW_S) & (rel < 0)
                m = (rel >= -30) & (rel <= 105)
                with np.errstate(invalid="ignore"), warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    R = M - np.nanmean(M[:, av], axis=1, keepdims=True)
                    med = np.nanmedian(R, axis=0)
                if typ == "stim":
                    ax.plot(rel[m], med[m], color=PALETTE[i % len(PALETTE)], lw=2.0,
                            marker="o", ms=3, label=f"stim {i + 1}")
                else:
                    ax.plot(rel[m], med[m], color=COL_TEMOIN, lw=1.1, alpha=0.8,
                            label="témoins" if i == 0 else None)
        ax.axhline(0, color=INK_2, lw=0.9, ls="--")
        ax.set_xlabel("temps depuis le début de la stim (s)", color=INK)
        ax.set_ylabel(lab, color=INK)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"Réponse de chaque stimulation, somas Cellpose (médiane, n = "
                 f"{sub['roi'].nunique()}) — chaque courbe part des 10 s qui précèdent sa stim",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    sauver(fig, out, dpi=180)


def figure_avant_apres(pp, out):
    sub = pp[pp["famille"] == "cellpose"]
    cols = [("dR_R0_pct", "ΔR / R$_0$ (%)"), ("dtau_ch1_ps", "Δτ CFP (ps)")]
    groupes = [(f"stim {n}", sub[(sub["type"] == "stim") & (sub["numero"] == n)])
               for n in sorted(sub.loc[sub["type"] == "stim", "numero"].unique())]
    tem = sub[sub["type"] == "témoin"]
    if len(tem):
        groupes.append(("témoins\n(réunis)", tem))
    fig, axes = plt.subplots(2, len(groupes), figsize=(2.5 * len(groupes) + 1.2, 8.4),
                             sharey="row", facecolor=SURFACE, squeeze=False)
    for r, (col, lab) in enumerate(cols):
        allv = []
        for c, (nom, g) in enumerate(groupes):
            ax = axes[r][c]
            style(ax)
            w = g[[f"{col}_avant", f"{col}_apres"]].dropna().to_numpy()
            allv.append(w)
            alpha = min(0.6, max(0.08, 25 / max(len(w), 1)))
            for row in w:
                ax.plot([0, 1], row, color=TRACE, lw=0.7, alpha=alpha)
            couleur = COL_TEMOIN if nom.startswith("témoin") else PALETTE[(c) % len(PALETTE)]
            for j in (0, 1):
                ax.scatter(np.full(len(w), j), w[:, j], s=11, color=couleur, alpha=0.8,
                           edgecolors=SURFACE, linewidths=0.5, zorder=3)
            med = np.median(w, axis=0)
            ax.plot([0, 1], med, color=INK, lw=2.2, marker="o", ms=6, zorder=4)
            e = w[:, 1] - w[:, 0]
            p = stats.wilcoxon(e).pvalue if len(e) >= 6 else np.nan
            ax.set_title(f"{nom}\nécart méd. {np.median(e):+.2f} ; {int((e < 0).sum())}/{len(e)} ↓\n"
                         f"p = {p:.1g}", color=INK, fontsize=9.5)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"{PRE_WINDOW_S:.0f} s\navant", f"{RESPONSE_WINDOW_S:.0f} s\naprès"],
                               color=INK, fontsize=9)
            ax.set_xlim(-0.35, 1.35)
            ax.axhline(0, color=INK_2, lw=0.8, ls="--")
        axes[r][0].set_ylabel(lab, color=INK, fontsize=11)
        ylim_robuste(axes[r][0], np.concatenate([a.ravel() for a in allv]))
    fig.suptitle("Avant / après chaque stimulation, somas Cellpose — un trait par ROI, médiane "
                 "en noir ; p : Wilcoxon sur l'écart", color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    sauver(fig, out, dpi=180)


def figure_heatmap(df, pp, fen, out):
    sub = df[df["famille"] == "cellpose"]
    st = pp[(pp["famille"] == "cellpose") & (pp["type"] == "stim")]
    score = st.groupby("roi")["dtau_ch1_ps_ecart"].mean()
    fig, axes = plt.subplots(1, 2, figsize=(15, 0.16 * sub["roi"].nunique() + 2.2), facecolor=SURFACE)
    for ax, col, cmap_lab in zip(axes, ["dR_R0_pct", "dtau_ch1_ps"], ["ΔR/R$_0$ (%)", "Δτ CFP (ps)"]):
        M, t, rois = _matrice(sub, col)
        o = np.argsort(-score.reindex(rois).fillna(-1e9).to_numpy())
        M, r = M[o], rois[o]
        lim = np.nanpercentile(np.abs(M), 97)
        cm = "RdBu_r" if col.startswith("dtau") else "RdBu"
        h = ax.imshow(M, aspect="auto", cmap=cm, vmin=-lim, vmax=lim, interpolation="nearest",
                      extent=(t[0], t[-1], len(r) - 0.5, -0.5))
        for a, b in fen:
            ax.axvline(a, color=INK, lw=0.9, ls="--")
            ax.axvline(b, color=INK, lw=0.9, ls="--")
        ax.set_yticks(range(len(r)))
        ax.set_yticklabels(r, fontsize=6.5, color=INK_2)
        ax.set_xlabel("temps (s)", color=INK)
        ax.tick_params(colors=INK_2, length=0)
        cb = fig.colorbar(h, ax=ax, fraction=0.03, pad=0.01)
        cb.set_label(cmap_lab, color=INK_2)
        ax.set_title(cmap_lab + " par ROI", color=INK, fontsize=11, loc="left")
    axes[0].set_ylabel("ROI (trié par Δτ CFP moyen après les stims)", color=INK, fontsize=9.5)
    fig.suptitle(f"Somas Cellpose — {len(fen)} stimulations (traits tiretés) ; rouge = moins de "
                 "FRET dans les deux panneaux (ratio en baisse, τ CFP en hausse) ; blanc = blocs "
                 "de stim ou pas de changement", color=INK, fontsize=11.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    sauver(fig, out, dpi=170)


def figure_flimapp(df, fr, fen, out):
    sub = df[df["famille"] == "flimapp"]
    rois = sorted(sub["roi"].unique(), key=lambda r: int(r[1:]))
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), sharex=True, facecolor=SURFACE)
    for ax, (col, lab) in zip(axes, [("dR_R0_pct", "ΔR / R$_0$ (%)"), ("dtau_ch1_ps", "Δτ CFP (ps)")]):
        style(ax)
        bandes_stim(ax, fen)
        for i, r in enumerate(rois):
            g = sub[sub["roi"] == r].sort_values("bloc")
            y = g[col].where(~g["stim"])
            ax.plot(g["t_s"], y, color=plt.cm.tab10(i), lw=1.1, label=r)
        ax.axhline(0, color=INK_2, lw=0.9, ls="--")
        ax.set_ylabel(lab, color=INK)
        ylim_robuste(ax, sub[col].where(~sub["stim"]))
    axes[0].legend(frameon=False, fontsize=8.5, ncol=8, loc="lower left")
    axes[1].set_xlabel("temps (s)", color=INK)
    axes[0].set_title("ROIs dessinées dans FLIMApp pendant l'expérience (F1…F8)", color=INK,
                      fontsize=12.5, loc="left")
    fig.tight_layout()
    sauver(fig, out, dpi=180)


def figure_validation(df, fr, out):
    sub = df[df["famille"] == "flimapp"].copy()
    sub["cle"] = sub["roi"].str[1:].astype(int) * 10000 + sub["bloc"]
    ref = pd.DataFrame({"cle": [int(r["roi_index"]) * 10000 + int(r["pass"]) for r in fr],
                        **{f"flimapp_ch{k}": [float(r[f"lifetime_ch{k}_ns"]) for r in fr] for k in CANAUX}})
    m = sub.merge(ref, on="cle")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), facecolor=SURFACE)
    for ax, k in zip(axes, CANAUX):
        style(ax)
        x, y = m[f"flimapp_ch{k}"], m[f"tau_ch{k}"]
        ax.scatter(x, y, s=5, alpha=0.4, color=PALETTE[k - 1])
        lo, hi = np.nanpercentile(np.r_[x, y], (0.5, 99.5))
        ax.plot([lo, hi], [lo, hi], color=INK_2, lw=0.9, ls="--")
        r = np.corrcoef(x, y)[0, 1]
        ax.set_title(f"ch{k} ({NOMS_CANAUX[k]}) : r = {r:.3f}, écart médian "
                     f"{np.nanmedian(y - x) * 1000:+.0f} ps", color=INK, fontsize=10.5)
        ax.set_xlabel("τ FLIMApp, en direct (ns)", color=INK)
        ax.set_ylabel("τ de ce pipeline (ns)", color=INK)
    fig.suptitle("Durées de vie bloc par bloc des 8 ROIs FLIMApp : ce pipeline contre le direct",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    sauver(fig, out, dpi=160)


# ------------------------------------------------------------------ images fixes
def images_fixes(root, D, info, sel, rep):
    """C2 moyen par ROI Cellpose sur chaque image fixe Bliq (repère FLIM), et indice
    d'expression de eNpHR-YFP = C2 à 950 nm / C2 à 800 nm, aux images du plan de la
    session seulement."""
    p = D["projections"] / "fovs.npz"
    if not p.exists():
        return pd.DataFrame()
    fovs = np.load(p)
    masks = np.load(D["masks"] / "cellpose_masks.npy")
    rows = []
    for nom, t in info["fovs"].items():
        im = fovs[f"{nom} C2".replace(" ", "_")]
        ok = np.isfinite(im)
        lab = np.where(ok, masks, 0)
        s = np.bincount(lab.ravel(), weights=np.nan_to_num(im).ravel(), minlength=masks.max() + 1)
        n = np.bincount(lab.ravel(), minlength=masks.max() + 1)
        npx = np.bincount(masks.ravel(), minlength=masks.max() + 1)
        for l in sel:
            rows.append(dict(image=nom, nm=t["nm"], moment=t["moment"], meme_plan=t["meme_plan"],
                             roi=str(l), C2_moyen=s[l] / n[l] if n[l] else np.nan,
                             frac_valide=n[l] / npx[l] if npx[l] else 0))
    fx = pd.DataFrame(rows)
    fx.to_csv(D["results"] / "images_fixes_C2.csv", index=False)
    return fx


def figure_expression(fx, rep, out):
    pl = fx[fx["meme_plan"]]
    w = pl.pivot_table(index="roi", columns="nm", values="C2_moyen")
    if not {800, 950} <= set(w.columns):
        return None
    idx = (w[950] / w[800]).rename("indice_950_800")
    m = rep.join(idx, how="inner")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), facecolor=SURFACE)
    style(axes[0])
    for row in w[[800, 950]].to_numpy():
        axes[0].plot([0, 1], row, color=TRACE, lw=0.7, alpha=0.6)
    axes[0].plot([0, 1], w[[800, 950]].median().to_numpy(), color=INK, lw=2.2, marker="o")
    axes[0].set_xticks([0, 1])
    axes[0].set_xticklabels(["800 nm", "950 nm"], color=INK)
    axes[0].set_xlim(-0.3, 1.3)
    axes[0].set_ylabel("C2 moyen par ROI (counts)", color=INK)
    axes[0].set_title(f"Images fixes « après », n = {len(w)}", color=INK, fontsize=10.5)
    out_stats = {}
    for ax, col, lab in [(axes[1], "dR_R0_pct_ecart_stims", "réponse ΔR/R$_0$ (points)"),
                         (axes[2], "dtau_ch1_ps_ecart_stims", "réponse Δτ CFP (ps)")]:
        style(ax)
        ax.scatter(m["indice_950_800"], m[col], s=22, color=PALETTE[0], alpha=0.85)
        rho, p = stats.spearmanr(m["indice_950_800"], m[col], nan_policy="omit")
        out_stats[col] = (rho, p)
        ax.axhline(0, color=INK_2, lw=0.8, ls="--")
        ax.set_xlabel("indice eNpHR : C2 950 nm / C2 800 nm", color=INK)
        ax.set_ylabel(lab + ", moyenne des stims", color=INK, fontsize=9.5)
        ax.set_title(f"Spearman ρ = {rho:+.2f}, p = {p:.2g}", color=INK, fontsize=10.5)
    fig.suptitle("La réponse dépend-elle de l'expression de eNpHR-YFP ? (images au plan de la "
                 "session seulement)", color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    sauver(fig, out, dpi=170)
    return idx, out_stats


# ------------------------------------------------------------------ main
def main():
    root = config.racine()
    D = config.dirs(root)
    info = json.loads((D["projections"] / "session.json").read_text())
    print(f"{root.name} — ajustements…")
    df, glob, fen, S = tableau(root, D, info)
    if not ((df["famille"] == "cellpose") & (df["phase"] == "avant") & (df["bloc"] > BASE_SKIP)).any():
        # Pas de ligne de base avant la 1re stim : ni R0, ni réponse avant / après, ni témoin.
        df.to_csv(D["results"] / "mesures_blocs.csv", index=False)
        glob.to_csv(D["qc"] / "ajustement_global.csv", index=False)
        sel = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
        images_fixes(root, D, info, sel, None)
        figure_traces(df, fen, [("ratio", "ratio ch2 / ch1 (YFP / CFP)", PALETTE[0]),
                                ("tau_ch1", "τ CFP, ch1 (ns)", PALETTE[2])],
                      "Ratio et durée de vie du donneur, valeurs brutes — somas Cellpose",
                      f"{df.loc[df['famille'] == 'cellpose', 'roi'].nunique()} ROIs · pas de ligne de "
                      f"base : la 1re stim commence à {fen[0][0]:.1f} s\nbandes jaunes : stimulations",
                      D["results"] / "trace_brute.png")
        print(f"\n  Pas de ligne de base (1re stim à {fen[0][0]:.1f} s) : seules les valeurs brutes "
              f"sont écrites (mesures_blocs.csv, trace_brute.png).\nRésultats -> {D['results']}")
        return
    t_fin = info["blocs"][-1]["t1"]
    tem = temoins(fen, t_fin)
    df.to_csv(D["results"] / "mesures_blocs.csv", index=False)
    glob.to_csv(D["qc"] / "ajustement_global.csv", index=False)

    pp = avant_apres(df, fen, tem)
    pp.to_csv(D["results"] / "avant_apres_stim.csv", index=False)
    rs = resume(pp)
    rs.to_csv(D["results"] / "resume_stims.csv", index=False)
    rec = recuperation(df, fen, t_fin)
    rec.to_csv(D["results"] / "recuperation.csv", index=False)

    # --- résumé par ROI : base et réponse moyenne des stims / des témoins
    base = df[(df["phase"] == "avant") & (df["bloc"] > BASE_SKIP)].groupby(["famille", "roi"])
    rep = pd.DataFrame({"R0": base["ratio"].median(), "tau0_ch1_ns": base["tau_ch1"].median(),
                        "tau0_ch2_ns": base["tau_ch2"].median(),
                        "photons_ch1_par_bloc": base["photons_ch1"].median(),
                        "photons_ch2_par_bloc": base["photons_ch2"].median()})
    for typ in ("stim", "témoin"):
        g = pp[pp["type"] == typ].groupby(["famille", "roi"])
        for m in ("dR_R0_pct", "dtau_ch1_ps", "dtau_ch2_ps", "dN_ch1_pct", "dN_ch2_pct"):
            rep[f"{m}_ecart_{'stims' if typ == 'stim' else 'temoins'}"] = g[f"{m}_ecart"].mean()
    rep = rep.reset_index()

    print(f"\n  {len(fen)} stims, {len(tem)} témoins ; somas : {df.loc[df['famille'] == 'cellpose', 'roi'].nunique()}")
    cols = ["dR_R0_pct", "dtau_ch1_ps", "dtau_ch2_ps", "dN_ch1_pct", "dN_ch2_pct"]
    for fam in ("cellpose", "flimapp", "neuropile", "champ"):
        print(f"\n  {fam} — écart après − avant (médiane sur les ROIs)")
        print("    " + " ".join(f"{MESURES[c][0].replace('$', '').replace('_0', '0'):>22s}" for c in cols))
        for _, r in rs[rs["famille"] == fam].iterrows():
            lab = f"{r['type']} {r['numero']}"
            print(f"    {lab:9s}" + " ".join(f"{r[c + '_median']:+18.2f} {MESURES[c][1]:3s}" for c in cols)
                  + (f"   p(Δτ1) = {r['dtau_ch1_ps_p_wilcoxon']:.1g}, p(ΔR) = {r['dR_R0_pct_p_wilcoxon']:.1g}"
                     if fam == "cellpose" else ""))

    cp = rep[rep["famille"] == "cellpose"].set_index("roi")
    for m in ("dR_R0_pct", "dtau_ch1_ps"):
        a, b = cp[f"{m}_ecart_stims"], cp[f"{m}_ecart_temoins"]
        p = stats.wilcoxon(a - b, nan_policy="omit").pvalue
        print(f"\n  {m} : stims {a.median():+.2f} contre témoins {b.median():+.2f} "
              f"(même ROI, Wilcoxon p = {p:.2g}, {int((a < b).sum())}/{len(a)} ROIs plus bas après stim)")
    rc = rec[rec["famille"] == "cellpose"]
    for (m, k), g in rc.groupby(["mesure", "stim"]):
        rep_ = (g["apres"] - g["avant"]).median()
        rest = (g["fin"] - g["apres"]).median()
        print(f"  {m} stim {k} : réponse {rep_:+.2f}, puis {rest:+.2f} de plus jusqu'aux 10 s "
              f"qui précèdent la stim suivante ({g['delai_fin_s'].iloc[0]:.0f} s après la fin)")

    # --- figures
    sel = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
    n = len(sel)
    figure_traces(df, fen, [("ratio", "ratio ch2 / ch1 (YFP / CFP)", PALETTE[0]),
                            ("dR_R0_pct", "ΔR / R$_0$ (%)", PALETTE[0])],
                  "Ratio FRET dans le temps, somas Cellpose",
                  f"{n} ROIs, blocs de 3 s, valeurs brutes · médiane et écart interquartile\n"
                  "bandes jaunes : stimulations (blocs masqués : la stim ajoute ~4 % de photons dans ch2)",
                  D["results"] / "trace_ratio.png")
    figure_traces(df, fen, [("tau_ch1", "τ CFP, ch1 (ns)", PALETTE[2]),
                            ("dtau_ch1_ps", "Δτ CFP (ps)", PALETTE[2])],
                  "Durée de vie du donneur (CFP) dans le temps, somas Cellpose",
                  f"{n} ROIs, mono-exponentielle reconvoluée, blocs de 3 s · médiane et écart "
                  "interquartile\nune hausse de τ CFP = moins de FRET (entrée de Cl⁻ ou perte de "
                  "l'accepteur)", D["results"] / "trace_duree_de_vie.png")
    figure_decomposition(df, fen, tem, D["results"] / "decomposition.png")
    figure_alignee(df, fen, tem, D["results"] / "reponse_alignee.png")
    figure_avant_apres(pp, D["results"] / "avant_apres_stim.png")
    figure_heatmap(df, pp, fen, D["results"] / "heatmap.png")
    fr = config.frames_flimapp(root)
    if (df["famille"] == "flimapp").any():
        figure_flimapp(df, fr, fen, D["results"] / "trace_flimapp.png")
        figure_validation(df, fr, D["qc"] / "validation_duree_de_vie.png")

    fx = images_fixes(root, D, info, sel, rep)
    if len(fx):
        r = figure_expression(fx, cp, D["results"] / "expression_eNpHR.png")
        if r:
            idx, st = r
            rep = rep.merge(idx.reset_index().assign(famille="cellpose"), on=["famille", "roi"], how="left")
            print("\n  indice eNpHR (C2 950 / 800 nm, plan de la session) contre la réponse : "
                  + ", ".join(f"{k} ρ = {v[0]:+.2f} (p = {v[1]:.2g})" for k, v in st.items()))
    rep.to_csv(D["results"] / "resume_par_roi.csv", index=False)
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
