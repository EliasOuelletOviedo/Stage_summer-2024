"""Étape 15 — calibration avec la fraction FRET x du déclin bleu, au lieu du temps moyen.

Usage :  python3 step15_fraction_fret.py <dossier> [--bootstrap 2000]

Déclin CFP à deux populations :

    d_D(t) = x · exp(−t/τ_DA)  +  (1 − x) · exp(−t/τ_D)
             (FRET)

Plus un terme de diffusion en forme d'IRF (exponentielle de τ_diff ≤ quelques dizaines de ps,
ajustée), exclu de x : x = a_DA / (a_DA + a_D).

1. τ_DA, τ_D et τ_diff : communes à tous les somas gardés, ajustées (Poisson, IRF de l'étape 5, décalage
   médian des FOV, fond libre) sur la somme des déclins CFP de tous ces somas, conditions réunies.
2. x : par soma, τ fixées, amplitudes pré-exponentielles et fond libres (`Modele.amplitudes`).
   x est la part d'amplitude (pas d'intensité) de la population FRET, comme dans l'équation.
3. Isotherme et bootstrap hiérarchique de l'étape 7, pour R, τ mono et x.

Sorties : results/fraction_fret_somas.csv, results/fraction_fret.json,
results/calibration_fraction_fret.png (même figure que calibration.png, 3ᵉ panneau = x).
"""
import argparse
import json

import numpy as np
import pandas as pd

import config as C
import figures as G
from figures import plt
from flim import Modele
from step7_calibration import ajuster, inverser, isotherme, tirage

LECT = [("R", "rapport YFP/CFP brut", ""), ("tau_mono", "τ CFP, mono-exponentielle", "ns"),
        ("x", "x, fraction FRET du déclin CFP", "")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--bootstrap", type=int, default=2000)
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    irf = json.loads((D["qc"] / "irf.json").read_text())
    M = Modele(C.charger_irf(irf["voie_cfp"], root))
    T = pd.read_csv(D["results"] / "somas.csv")

    # déclins CFP des somas gardés (cycles gardés réunis, comme à l'étape 6)
    declins, cles = {}, []
    for a in C.acquisitions(root):
        if a["type"] != "fov" or not (D["mesures"] / f"{a['nom']}.npz").exists():
            continue
        m = np.load(D["mesures"] / f"{a['nom']}.npz")
        tot = m["declins"][m["garde"]].sum(0)                                # (n+2, 2, 256)
        s = T[(T.acquisition == a["nom"]) & (T.type == "soma") & T.garde]
        for roi in s.roi.astype(int):
            declins[(a["nom"], roi)] = (tot[roi, 0].astype(float), irf["fovs"][a["nom"]]["s"])

    # 1. τ communes
    somme = np.sum([y for y, _ in declins.values()], axis=0)
    g = M.ajuster_multi(somme, 3, s=irf["s_median"], depart=[2.0, 0.5, 0.02])
    tau_D, tau_DA, tau_diff = map(float, g["taus"])
    g["x"] = float(g["amps"][1] / (g["amps"][0] + g["amps"][1]))
    print(f"  τ_D = {tau_D:.3f} ns, τ_DA = {tau_DA:.3f} ns, diffusion τ = {1e3*tau_diff:.0f} ps (amplitude {g['amps'][2]:.3f}), x global = {g['x']:.3f}, "
          f"fond {g['b']:.3f}, déviance/ddl {g['dev_ddl']:.2f} ({len(declins)} somas, {somme.sum():.3g} photons)")

    # 2. x par soma
    lignes = []
    for (nom, roi), (y, s_fov) in declins.items():
        a, b, dd = M.amplitudes(y, [tau_D, tau_DA, tau_diff], s_fov)
        lignes.append(dict(acquisition=nom, roi=roi, x=float(a[1] / (a[0] + a[1])), diffusion=float(a[2]), b_x=b, dev_ddl_x=dd))
    X = pd.DataFrame(lignes)
    X.to_csv(D["results"] / "fraction_fret_somas.csv", index=False)
    S = T[(T.type == "soma") & T.garde].merge(X, on=["acquisition", "roi"])
    conds = [c for c in C.CONDITIONS_FIXES if (S.condition == c).any()]
    par = {c: S[S.condition == c] for c in conds}
    imposees = [c for c in conds if c != "basal"]
    cs = [C.CL_MM[c] for c in imposees]

    # 3. isothermes et bootstrap (étape 7)
    rng = np.random.default_rng(1)
    CAL, boot = {}, {}
    for col, nom, _ in LECT:
        med = [par[c][col].median() for c in imposees]
        fit = ajuster(cs, med)
        B = []
        for _ in range(args.bootstrap):
            f = ajuster(cs, [tirage(par[c], rng)[col].median() for c in imposees])
            xb = tirage(par["basal"], rng)[col].median() if "basal" in par else np.nan
            B.append((np.nan,) * 4 if f is None else (*f, float(inverser(xb, *f))))
        B = np.array(B); boot[col] = B
        ok = np.isfinite(B[:, 2])
        row = dict(medianes=med)
        if fit is not None:
            row.update(X0=fit[0], Xinf=fit[1], K_mM=fit[2])
            for j, k in enumerate(("X0", "Xinf", "K_mM", "basal_cl_mM")):
                v = B[ok, j]; v = v[np.isfinite(v)]
                row[k + "_ic"] = np.percentile(v, [2.5, 97.5]).tolist()
            if "basal" in par:
                row["basal_X"] = float(par["basal"][col].median())
                row["basal_cl_mM"] = float(inverser(row["basal_X"], *fit))
        CAL[col] = row
        print(f"  {nom:>32} : " + (f"K {fit[2]:.1f} mM [{row['K_mM_ic'][0]:.1f} ; {row['K_mM_ic'][1]:.1f}], "
                                   f"basal {row.get('basal_cl_mM', np.nan):.1f} mM" if fit else "pas d'isotherme"))
    json.dump(dict(tau_D=tau_D, tau_DA=tau_DA, tau_diffusion=tau_diff, ajustement_global=dict(x=g["x"], amplitudes=g["amps"].tolist(), b=g["b"], dev_ddl=g["dev_ddl"]),
                   calibration=CAL), open(D["results"] / "fraction_fret.json", "w"), indent=2, default=float)

    # figure : comme calibration.png de l'étape 7
    rng = np.random.default_rng(0)
    cc = np.linspace(0, 140, 300)
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.8))
    for ax, (col, nom, unite) in zip(axes, LECT):
        G.style(ax)
        row = CAL[col]
        for cond in imposees:
            s, x = par[cond], C.CL_MM[cond]
            ax.plot(x + rng.normal(0, 1.6, len(s)), s[col], "o", ms=3, alpha=0.45, color=C.COULEUR_CONDITION[cond], mec="none")
            fm = s.groupby("acquisition")[col].median()
            ax.plot(np.full(len(fm), x) + 5, fm, "_", ms=10, mew=1.6, color=C.COULEUR_CONDITION[cond])
            ax.plot(x, s[col].median(), "D", ms=7, color=G.INK, mfc=C.COULEUR_CONDITION[cond], zorder=5)
        if "K_mM" in row:
            B = boot[col]; ok = np.isfinite(B[:, 2])
            courbes = np.array([isotherme(cc, *b[:3]) for b in B[ok][:400]])
            ax.fill_between(cc, *np.percentile(courbes, [2.5, 97.5], axis=0), color=G.TRACE, alpha=0.35, lw=0)
            ax.plot(cc, isotherme(cc, row["X0"], row["Xinf"], row["K_mM"]), color=G.INK, lw=1.4)
            if "basal_cl_mM" in row and np.isfinite(row["basal_cl_mM"]):
                cb, (lo, hi) = row["basal_cl_mM"], row["basal_cl_mM_ic"]
                ax.errorbar(cb, row["basal_X"], xerr=[[max(cb - lo, 0)], [max(hi - cb, 0)]], fmt="s", ms=8,
                            color=C.COULEUR_CONDITION["basal"], capsize=3, zorder=6,
                            label=f"basal : {cb:.1f} mM [{lo:.1f} ; {hi:.1f}]")
                ax.legend(frameon=False, fontsize=8.5, loc="best")
            sous = (f"K = {row['K_mM']:.1f} mM [{row['K_mM_ic'][0]:.1f} ; {row['K_mM_ic'][1]:.1f}], "
                    f"X₀ = {row['X0']:.3g}, X∞ = {row['Xinf']:.3g}")
        else:
            sous = "pas d'isotherme : médianes non monotones"
        if col == "x":
            ax.text(0.98, 0.98, f"τ_DA = {tau_DA:.3f} ns\nτ_D = {tau_D:.3f} ns\ndiffusion exclue (τ = {1e3*tau_diff:.0f} ps)\n(communes aux somas)", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8.5, color=G.INK_2)
        ax.set_xlabel("[Cl⁻] imposée (mM)", color=G.INK)
        ax.set_ylabel(nom + (f" ({unite})" if unite else ""), color=G.INK)
        G.titre(ax, nom if nom[0] in "τx" else nom[0].upper() + nom[1:], sous, fontsize=11)
    fig.text(0.01, 0.005, "points : somas ; traits : médiane de chaque FOV ; losanges : médiane de la "
             "condition ; bande : IC 95 % (bootstrap hiérarchique) ; carré : médiane du basal, inversée ; "
             "x : d_D = x·e^(−t/τ_DA) + (1 − x)·e^(−t/τ_D), diffusion (forme d'IRF) ajustée à part", color=G.INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    G.sauver(fig, D["results"] / "calibration_fraction_fret.png")
    print("  -> results/calibration_fraction_fret.png, fraction_fret_somas.csv, fraction_fret.json")


if __name__ == "__main__":
    main()
