"""Étape 6 — la réponse à la stim dépend-elle du ratio basal ?

Question : les cellules à R basal élevé (peu de Cl-) répondent-elles plus fort ?

Pour chaque ROI et chaque stim étiquetée (VIS / IR) :
  R basal  = ratio C2/C3 moyen au début de la ligne de base, de l'image 3 jusqu'à 10 s
             avant la stim — une fenêtre distincte de celle du ΔR, pour que le bruit
             de mesure ne crée pas à lui seul une corrélation (régression vers la moyenne) ;
  ΔR       = ratio moyen des 10 s après la stim - celui des 10 s avant ;
  ΔR/R0    = ΔR rapporté au ratio des 10 s avant, en %.

La figure ne montre que ΔR absolu contre R basal après toutes les stims VIS (1, 10,
50, 100 % et FOV2/stim_VIS), un point par ROI et par stim — un même ROI y figure donc
une fois par stim VIS de sa tranche.
Le script imprime en plus les tendances qui permettent de l'interpréter :
  A  ΔR absolu contre R basal, VIS 50 % ;
  B  ΔR/R0 contre R basal, VIS 50 % et conditions témoins (IR, VIS 1 %) ;
  C  réponse propre à la stim : ΔR/R0 à VIS 50 % moins ΔR/R0 au témoin de la même
     cellule, contre R basal.

Les tendances sont estimées **à l'intérieur de chaque FOV** (valeurs centrées sur la
moyenne de la tranche, pente de Theil-Sen, rho de Spearman sur les rangs intra-FOV) :
les tranches diffèrent à la fois par leur R basal et par leur amplitude de réponse, et
une corrélation calculée tous FOVs confondus mêlerait les deux.

Usage :  python3 step6_deltaR_vs_Rbasal.py   (après step4_ratio.py)
"""
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, theilslopes

import config
from config import dirs
from step4_ratio import BASE_SKIP, INK, INK_2, SURFACE, _style

WIN = 10.0            # s avant / après la stim pour le ΔR
POINT = "#3a3a36"


# ------------------------------------------------------------------ mesures
def measure(df):
    rows = []
    for (fov, rec), g in df.groupby(["fov", "acquisition"], sort=False):
        typ = config.stim_type(rec)
        if typ is None or not g["stim"].any():
            continue
        t0, t1 = g.loc[g["stim"], "t_s"].min(), g.loc[g["stim"], "t_s"].max()
        bas = g[(g["phase"] == "avant") & (g["t_s"] < t0 - WIN) & (g["image"] > BASE_SKIP)]
        av = g[(g["phase"] == "avant") & (g["t_s"] >= t0 - WIN)]
        ap = g[(g["phase"] == "apres") & (g["t_s"] <= t1 + WIN)]
        w = pd.DataFrame({"R_basal": bas.groupby("roi")["ratio"].mean(),
                          "R_avant": av.groupby("roi")["ratio"].mean(),
                          "R_apres": ap.groupby("roi")["ratio"].mean()})
        w["dR"] = w["R_apres"] - w["R_avant"]
        w["dR_R0_pct"] = 100 * w["dR"] / w["R_avant"]
        pw = rec.split("VIS_")[1].split("-")[0] if "VIS_" in rec else None
        cond = "IR" if typ == "IR" else (f"VIS {pw} %" if pw else "VIS ?")
        rows.append(w.reset_index().assign(fov=fov, acquisition=rec, condition=cond))
    out = pd.concat(rows, ignore_index=True).dropna(subset=["R_basal", "dR"])
    return out[["fov", "acquisition", "condition", "roi", "R_basal", "R_avant", "R_apres",
                "dR", "dR_R0_pct"]]


def within(sub, y, x="R_basal", group="fov"):
    """Tendance intra-FOV : rho de Spearman sur les rangs et pente de Theil-Sen [IC95]."""
    rx = sub.groupby(group)[x].rank(pct=True)
    ry = sub.groupby(group)[y].rank(pct=True)
    rho, p = spearmanr(rx, ry)
    xc = sub[x] - sub.groupby(group)[x].transform("mean")
    yc = sub[y] - sub.groupby(group)[y].transform("mean")
    s, _, lo, hi = theilslopes(yc, xc)
    return dict(rho=rho, p=p, pente=s, lo=lo, hi=hi, n=len(sub))


# ------------------------------------------------------------------ figure
def _pfmt(p):
    return f"p = {p:.1g}" if p >= 1e-3 else "p < 0,001"


def figure(vis, st, D, ylabel="ΔR après stim VIS  (ratio C2 / C3)", out=None):
    """ΔR absolu contre R basal après toutes les stims VIS : un point par ROI et par stim,
    et la droite de tendance intra-acquisition (pente de Theil-Sen passant par le centre
    du nuage), sans autre élément."""
    fig, ax = plt.subplots(figsize=(6.6, 5.2), facecolor=SURFACE)
    _style(ax)
    ax.axhline(0, color=INK_2, linewidth=0.9, linestyle="--", zorder=1)
    ax.scatter(vis["R_basal"], vis["dR"], s=16, color=POINT, alpha=0.7,
               edgecolors=SURFACE, linewidths=0.5, zorder=3)
    # Pente estimée à l'intérieur de chaque FOV x stim : les tranches diffèrent à la fois
    # par leur R basal et par leur amplitude, une droite toutes tranches confondues les
    # mêlerait.
    xs = np.array([vis["R_basal"].quantile(0.01), vis["R_basal"].quantile(0.99)])
    ax.plot(xs, vis["dR"].mean() + st["pente"] * (xs - vis["R_basal"].mean()),
            color=INK, linewidth=2.0, zorder=4)
    ax.set_xlabel("R basal (ratio C2 / C3)", color=INK, fontsize=11)
    ax.set_ylabel(ylabel, color=INK, fontsize=11)
    fig.tight_layout()
    out = out or D["results"] / "deltaR_vs_Rbasal.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


# ------------------------------------------------------------------ main
def main():
    D = dirs()
    path = D["results"] / "ratio_temps.csv"
    if not path.exists():
        sys.exit(f"{path} absent — lance d'abord step4_ratio.py")
    df = pd.read_csv(path, usecols=["fov", "acquisition", "roi", "image", "t_s",
                                    "phase", "stim", "ratio"])
    x = measure(df)
    x.to_csv(D["results"] / "deltaR_vs_Rbasal.csv", index=False)

    vis = x[x["condition"] == "VIS 50 %"]
    ctl = x[x["condition"].isin(["IR", "VIS 1 %"])]
    # Même cellule à VIS 50 % et au témoin de sa tranche (IR ou VIS 1 %).
    paired = vis.merge(ctl[["fov", "roi", "R_basal", "dR_R0_pct"]], on=["fov", "roi"],
                       suffixes=("", "_temoin"))
    paired["R_basal"] = (paired["R_basal"] + paired["R_basal_temoin"]) / 2
    paired["specifique"] = paired["dR_R0_pct"] - paired["dR_R0_pct_temoin"]

    vis_all = x[x["condition"].str.startswith("VIS")]
    stats = {"A_tout": within(vis_all, "dR", group=["fov", "acquisition"]),
             "A": within(vis, "dR"),
             "B_vis": within(vis, "dR_R0_pct"),
             "B_ctl": within(ctl, "dR_R0_pct", group=["fov", "acquisition"]),
             "C": within(paired, "specifique")}
    labels = {"A_tout": "ΔR absolu, toutes stims VIS", "A": "ΔR absolu, VIS 50 %", "B_vis": "ΔR/R0, VIS 50 %",
              "B_ctl": "ΔR/R0, témoins (IR, VIS 1 %)", "C": "VIS 50 % - témoin (même cellule)"}
    print("tendance intra-FOV en fonction du R basal :")
    for k, s in stats.items():
        unit = "ratio" if k.startswith("A") else "%"
        print(f"  {labels[k]:34s} rho = {s['rho']:+.2f} ({_pfmt(s['p'])}, n = {s['n']})   "
              f"pente {s['pente']:+.3g} {unit}/unité de R [IC95 {s['lo']:+.3g}, {s['hi']:+.3g}]")
    figure(vis_all, stats["A_tout"], D)
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
