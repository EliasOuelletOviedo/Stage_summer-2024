"""Étape 11 — l'écart entre champs vient-il d'un fond non lié au capteur ?

Usage :  python3 step11_fond.py <dossier>          (quelques secondes ; après l'étape 9)

Les validations croisées (étape 10) montrent que l'erreur d'une calibration vient surtout de l'écart entre champs :
à 20 mM, les FOV pâles (fov3, fov5 le 2026-11-09) sont lus trop haut, surtout par les lectures fondées sur le
rapport. Hypothèse : un fond (bruit des détecteurs, lumière parasite, autofluorescence à déclin lent) pèse davantage
dans les champs pâles et abaisse le rapport.

L'étape 9 (modèle à deux populations) ajuste, pour chaque soma et chaque canal, un fond uniforme sur la période,
part φ des photons de la fenêtre. On compare ce fond d'un FOV à l'autre (en part, et en photons par pixel et par
cycle), puis on corrige le rapport : R_f = Y·(1 − φ_Y) / (B·(1 − φ_B)) (photons de la fenêtre). Test : la
validation croisée par FOV de l'étape 10 (isotherme sur les médianes, le k-ième FOV de chaque [Cl⁻] imposée retiré
au pli k), sur le rapport brut de la fenêtre et sur R_f.

Limite : un fond uniforme ne capte que ce qui est plat dans le temps ; une autofluorescence à déclin court se
confond avec le déclin du capteur et reste dans R_f.

Sorties : results/fond_fovs.csv, results/fond_validation.csv, results/fond.png
"""
import re
import sys

import numpy as np
import pandas as pd
import scipy.stats as st

import config as C
import figures as G
from figures import plt
from step7_calibration import ajuster, inverser

CONDITIONS = ("0 mM", "20 mM", "131 mM")


def num_fov(nom):
    return int(re.search(r"fov(\d+)$", nom).group(1))


def validation_par_fov(T, col):
    """Erreur quadratique moyenne des médianes de FOV lues hors calibration, par condition ; basal par pli."""
    conc = [C.CL_MM[c] for c in CONDITIONS]
    imp = T[T.condition.isin(CONDITIONS)]
    err = {c: [] for c in CONDITIONS}
    basal = []
    for k in sorted({num_fov(a) for a in imp.acquisition}):
        tr = imp[imp.acquisition.map(num_fov) != k]
        f = ajuster(conc, [tr[tr.condition == c][col].median() for c in CONDITIONS])
        if f is None:
            return {c: np.nan for c in CONDITIONS}, np.nan
        for c in CONDITIONS:
            te = imp[(imp.condition == c) & (imp.acquisition.map(num_fov) == k)]
            if len(te):
                err[c].append(float(np.median(inverser(te[col].values, *f))) - C.CL_MM[c])
        basal.append(float(np.median(inverser(T[T.condition == "basal"][col].values, *f))))
    rmse = {c: (float(np.sqrt(np.mean(np.square(e)))) if np.all(np.isfinite(e)) else np.inf) for c, e in err.items()}
    return rmse, (float(np.median(basal)), float(np.min(basal)), float(np.max(basal)))


def main():
    root = C.racine()
    D = C.dirs(root)
    p = D["results"] / "deux_canaux_somas.csv"
    if not p.exists():
        sys.exit("results/deux_canaux_somas.csv absent : lance d'abord step9_flim_deux_canaux.py")
    T = pd.read_csv(p)
    s6 = pd.read_csv(D["results"] / "somas.csv")
    s6 = s6[s6.type == "soma"][["acquisition", "roi", "aire_px", "cycles"]]
    T = T.merge(s6, on=["acquisition", "roi"], how="left")
    T["R_fenetre"] = T.photons_yfp / T.photons_cfp
    T["R_fond"] = T.photons_yfp * (1 - T.fond_jaune) / (T.photons_cfp * (1 - T.fond_bleu))
    for k, canal in (("bleu", "cfp"), ("jaune", "yfp")):
        T[f"fond_{k}_px_cycle"] = T[f"fond_{k}"] * T[f"photons_{canal}"] / (T.aire_px * T.cycles)
    T["intensite"] = (T.photons_cfp + T.photons_yfp) / (T.aire_px * T.cycles)

    F = T.groupby(["acquisition", "condition"]).agg(
        somas=("roi", "size"), intensite=("intensite", "median"),
        fond_bleu=("fond_bleu", "median"), fond_jaune=("fond_jaune", "median"),
        fond_bleu_px_cycle=("fond_bleu_px_cycle", "median"), fond_jaune_px_cycle=("fond_jaune_px_cycle", "median"),
        R_fenetre=("R_fenetre", "median"), R_fond=("R_fond", "median")).reset_index()
    F.to_csv(D["results"] / "fond_fovs.csv", index=False)
    pd.set_option("display.width", 200)
    print(F.round(4).to_string(index=False))

    # écart relatif de chaque FOV à la médiane de sa condition, contre son intensité
    for col in ("R_fenetre", "R_fond"):
        F[f"ecart_{col}"] = F.groupby("condition")[col].transform(lambda v: v / v.median() - 1)
    imp = F[F.condition.isin(CONDITIONS)]
    for col in ("R_fenetre", "R_fond"):
        r = st.spearmanr(imp.intensite, imp[f"ecart_{col}"])[0]
        cv = imp.groupby("condition")[col].agg(lambda v: v.std() / v.mean())
        print(f"  {col} : ρ(intensité du FOV, écart à sa condition) = {r:+.2f} ; dispersion des médianes de FOV "
              f"(CV) : " + ", ".join(f"{c} {100 * cv[c]:.0f} %" for c in CONDITIONS))
    for k in ("bleu", "jaune"):
        r = st.spearmanr(F.intensite, F[f"fond_{k}"])[0]
        print(f"  fond {k} (part) contre l'intensité du FOV : ρ = {r:+.2f} ; en photons/px/cycle : médiane "
              f"{F[f'fond_{k}_px_cycle'].median():.2f}, de {F[f'fond_{k}_px_cycle'].min():.2f} à {F[f'fond_{k}_px_cycle'].max():.2f}")

    rows = []
    for col, nom in (("R_fenetre", "rapport brut (fenêtre)"), ("R_fond", "rapport, fond retiré")):
        rmse, b = validation_par_fov(T, col)
        rows.append(dict(lecture=nom, **{f"rmse_{c.replace(' ', '')}": v for c, v in rmse.items()},
                         basal_mM=b[0], basal_min=b[1], basal_max=b[2]))
    V = pd.DataFrame(rows)
    V.to_csv(D["results"] / "fond_validation.csv", index=False)
    print("\nValidation croisée par FOV (erreur quadratique moyenne des médianes de FOV, mM) :")
    print(V.round(2).to_string(index=False))
    figure(D, T, F, V)


def figure(D, T, F, V):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax in axes:
        G.style(ax)
    for cond in C.CONDITIONS_FIXES:
        f = F[F.condition == cond]
        col = C.COULEUR_CONDITION[cond]
        axes[0].plot(f.intensite, 100 * f.fond_bleu, "o", color=col, ms=7, label=f"{cond} — CFP")
        axes[0].plot(f.intensite, 100 * f.fond_jaune, "s", color=col, ms=6, mfc="none", label=f"{cond} — YFP")
        for _, r in f.iterrows():
            axes[1].annotate("", xy=(r.intensite, r.R_fond), xytext=(r.intensite, r.R_fenetre),
                             arrowprops=dict(arrowstyle="->", color=col, lw=1.2))
            axes[1].text(r.intensite, r.R_fond, f" {r.acquisition.split('_')[1]}", fontsize=7, color=col, va="center")
        axes[1].plot(f.intensite, f.R_fenetre, "o", color=col, ms=5, alpha=0.6, label=cond)
    axes[0].set_xscale("log"); axes[1].set_xscale("log")
    axes[0].set_xlabel("intensité du FOV (photons / px / cycle, médiane des somas)", color=G.INK)
    axes[0].set_ylabel("part de fond uniforme (%)", color=G.INK)
    axes[0].legend(frameon=False, fontsize=7, ncol=2)
    G.titre(axes[0], "Fond ajusté par l'étape 9, médiane des somas de chaque FOV", "ronds : CFP ; carrés : YFP")
    axes[1].set_xlabel("intensité du FOV (photons / px / cycle)", color=G.INK)
    axes[1].set_ylabel("rapport YFP/CFP, médiane des somas", color=G.INK)
    axes[1].legend(frameon=False, fontsize=8)
    G.titre(axes[1], "Rapport brut → rapport fond retiré, par FOV", "flèche : effet de la correction")
    ax = axes[2]
    x = np.arange(3)
    for i, (_, r) in enumerate(V.iterrows()):
        vals = [r["rmse_0mM"], r["rmse_20mM"], r["rmse_131mM"]]
        ax.bar(x + (i - 0.5) * 0.38, np.nan_to_num(vals, posinf=0), 0.36, color=G.PALETTE[i], label=r.lecture)
        for xi, v in zip(x, vals):
            ax.text(xi + (i - 0.5) * 0.38, (0 if not np.isfinite(v) else v), "∞" if not np.isfinite(v) else f"{v:.1f}",
                    ha="center", va="bottom", fontsize=8, color=G.INK)
    ax.set_xticks(x); ax.set_xticklabels(CONDITIONS)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylabel("erreur quadratique moyenne des médianes de FOV (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=8)
    G.titre(ax, "Validation croisée par FOV (isotherme)", "le k-ième FOV de chaque condition lu hors calibration")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "fond.png")


if __name__ == "__main__":
    main()
