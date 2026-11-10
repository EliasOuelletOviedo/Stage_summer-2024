"""Étape 7 — calibration en tissu et [Cl⁻]i basale.

Usage :  python3 step7_calibration.py <dossier> [--bootstrap 2000]

Pour chaque lecture X (rapport YFP/CFP brut ; corrigé du temps mort ; τ CFP mono-exponentielle ;
temps moyen τm ; τ moyen des trois populations), isotherme d'un capteur à un site :

    X(c) = X∞ + (X₀ − X∞) · K / (K + c)          c = [Cl⁻] imposée (0, 20, 131 mM)

ajustée sur les médianes des somas de chaque condition. Trois concentrations pour trois
paramètres : la courbe passe exactement par les trois médianes ; son incertitude vient du
**bootstrap hiérarchique** (on retire les FOV de chaque condition, puis les somas de chaque
FOV tiré), qui respecte le fait que les somas d'un même champ partagent profondeur,
puissance et état de la tranche. Inversion : c = K · (X₀ − X) / (X − X∞).

[Cl⁻]i basale : médiane des somas du basal, inversée (avec son bootstrap) ; et cellule par
cellule, avec la courbe ajustée. **Auto-cohérence** : les somas des conditions imposées,
relus un par un par la même courbe, montrent la précision d'une lecture sur une seule
cellule (des cellules toutes à 20 mM ne se relisent pas toutes à 20 mM).

Sorties : results/calibration.csv, results/basal_cellules.csv, results/conditions.csv,
results/autocoherence.csv, results/calibration.png, results/basal.png,
results/rapport_vs_tau.png, results/temps_mort.png
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import config as C
import figures as G
from figures import plt

LECTURES = [
    ("R", "rapport YFP/CFP brut", ""),
    ("R_corrige", "rapport YFP/CFP corrigé du temps mort", ""),
    ("tau_mono", "τ CFP, mono-exponentielle", "ns"),
    ("tau_m", "τm CFP, temps moyen", "ns"),
    ("tau_3pop", "τ CFP, 3 populations", "ns"),
]
CL_IMPOSEES = (0.0, 20.0, 131.0)


def isotherme(c, X0, Xinf, K):
    c = np.asarray(c, float)
    return Xinf + (X0 - Xinf) * K / (K + c)


def ajuster(cs, xs):
    """(X0, X∞, K) passant par les points (cs, xs) ; None si la courbe n'existe pas (non
    monotone, ou une pente qui ne s'aplatit pas : K hors de [0,1 ; 10⁴] mM)."""
    cs, xs = np.asarray(cs, float), np.asarray(xs, float)
    if not np.all(np.isfinite(xs)):
        return None
    d = np.diff(xs)
    if not (np.all(d < 0) or np.all(d > 0)):
        return None

    def res(p):
        return isotherme(cs, p[0], p[1], np.exp(p[2])) - xs
    x0 = [xs[0], xs[-1] - 0.3 * (xs[0] - xs[-1]), np.log(30.0)]
    r = least_squares(res, x0, bounds=([-np.inf, -np.inf, np.log(0.1)], [np.inf, np.inf, np.log(1e4)]))
    if np.max(np.abs(r.fun)) > 1e-6 * max(np.abs(xs).max(), 1):
        return None
    K = float(np.exp(r.x[2]))
    if K <= 0.1 * 1.0001 or K >= 1e4 * 0.9999:
        return None
    return float(r.x[0]), float(r.x[1]), K


def inverser(x, X0, Xinf, K):
    x = np.asarray(x, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        c = K * (X0 - x) / (x - Xinf)
    # au-delà de X0 (côté 0 mM) : c < 0 ; au-delà de X∞ : pas de solution
    return np.where((x - Xinf) * (X0 - Xinf) > 0, c, np.inf)


def tirage(df, rng):
    """Un tirage hiérarchique : FOV avec remise, puis somas avec remise dans chaque FOV."""
    fovs = df.acquisition.unique()
    pris = rng.choice(fovs, len(fovs), replace=True)
    parts = []
    for f in pris:
        s = df[df.acquisition == f]
        parts.append(s.iloc[rng.integers(0, len(s), len(s))])
    return pd.concat(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--bootstrap", type=int, default=2000)
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    p = D["results"] / "somas.csv"
    if not p.exists():
        sys.exit("results/somas.csv absent : lance d'abord step6_somas.py")
    T = pd.read_csv(p)
    S = T[(T.type == "soma") & T.garde].copy()
    conds = [c for c in C.CONDITIONS_FIXES if (S.condition == c).any()]
    rng = np.random.default_rng(1)

    # --- résumé par condition
    rows = []
    for cond in conds:
        s = S[S.condition == cond]
        r = dict(condition=cond, cl_mM=C.CL_MM[cond], fovs=s.acquisition.nunique(), somas=len(s),
                 photons_cfp_median=s.photons_cfp.median(), mu_yfp_median=s.mu_yfp.median(),
                 mu_cfp_median=s.mu_cfp.median())
        for col, _, _ in LECTURES:
            q = s[col].quantile([0.25, 0.5, 0.75]).values
            r.update({f"{col}_q1": q[0], f"{col}_median": q[1], f"{col}_q3": q[2],
                      f"{col}_ecart_fov": s.groupby("acquisition")[col].median().std()})
        rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(D["results"] / "conditions.csv", index=False)

    # --- calibration et bootstrap
    cal_rows, basal_rows, auto_rows, boot = [], [], [], {}
    par_cond = {c: S[S.condition == c] for c in conds}
    imposees = [c for c in conds if c != "basal"]
    cs = [C.CL_MM[c] for c in imposees]
    for col, nom, unite in LECTURES:
        med = [par_cond[c][col].median() for c in imposees]
        fit = ajuster(cs, med)
        B = []
        for _ in range(args.bootstrap):
            mb = [tirage(par_cond[c], rng)[col].median() for c in imposees]
            f = ajuster(cs, mb)
            xb = tirage(par_cond["basal"], rng)[col].median() if "basal" in par_cond else np.nan
            B.append((np.nan,) * 4 if f is None else (*f, float(inverser(xb, *f))))
        B = np.array(B)
        boot[col] = B
        ok = np.isfinite(B[:, 2])
        q = lambda v: np.nanpercentile(v[ok & np.isfinite(v)], [2.5, 50, 97.5]) if (ok & np.isfinite(v)).any() else [np.nan] * 3
        row = dict(lecture=col, nom=nom, unite=unite, X_0mM=med[0], X_20mM=med[1], X_131mM=med[2],
                   bootstrap=args.bootstrap, echecs=int((~ok).sum()))
        if fit is not None:
            X0, Xinf, K = fit
            row.update(X0=X0, Xinf=Xinf, K_mM=K, dynamique=X0 / Xinf if Xinf else np.nan)
            for j, k in enumerate(("X0", "Xinf", "K_mM")):
                lo, mid, hi = q(B[:, j])
                row.update({f"{k}_ic_bas": lo, f"{k}_ic_haut": hi})
            if "basal" in par_cond:
                xb = par_cond["basal"][col].median()
                cb = float(inverser(xb, X0, Xinf, K))
                lo, mid, hi = q(B[:, 3])
                row.update(basal_X=xb, basal_cl_mM=cb, basal_cl_ic_bas=lo, basal_cl_ic_haut=hi)
                cel = inverser(par_cond["basal"][col].values, X0, Xinf, K)
                for (_, sr), v in zip(par_cond["basal"].iterrows(), cel):
                    basal_rows.append(dict(lecture=col, acquisition=sr.acquisition, roi=sr.roi, X=sr[col], cl_mM=v))
            # pente locale : combien de mM vaut un écart de X autour de chaque concentration
            for c in imposees:
                v = inverser(par_cond[c][col].values, X0, Xinf, K)
                for (_, sr), vv in zip(par_cond[c].iterrows(), v):
                    auto_rows.append(dict(lecture=col, condition=c, cl_mM=C.CL_MM[c], acquisition=sr.acquisition,
                                          roi=sr.roi, X=sr[col], cl_relu_mM=vv))
        cal_rows.append(row)
        print(f"  {nom:>40} : " + (f"X0 {fit[0]:.4g}, X∞ {fit[1]:.4g}, K {fit[2]:.1f} mM "
                                   f"[{row.get('K_mM_ic_bas', np.nan):.1f} ; {row.get('K_mM_ic_haut', np.nan):.1f}]"
                                   f" ; basal {row.get('basal_cl_mM', np.nan):.1f} mM "
                                   f"[{row.get('basal_cl_ic_bas', np.nan):.1f} ; {row.get('basal_cl_ic_haut', np.nan):.1f}]"
                                   f" ; échecs du bootstrap {row['echecs']}" if fit else "pas d'isotherme (non monotone)"))
    CAL = pd.DataFrame(cal_rows)
    CAL.to_csv(D["results"] / "calibration.csv", index=False)
    BAS = pd.DataFrame(basal_rows)
    BAS.to_csv(D["results"] / "basal_cellules.csv", index=False)
    AUTO = pd.DataFrame(auto_rows)
    AUTO.to_csv(D["results"] / "autocoherence.csv", index=False)
    np.savez_compressed(D["mesures"] / "bootstrap_calibration.npz", **boot)
    print("  -> results/calibration.csv, basal_cellules.csv, conditions.csv, autocoherence.csv")
    figures(D, S, CAL, BAS, AUTO, boot, conds)


def figures(D, S, CAL, BAS, AUTO, boot, conds):
    rng = np.random.default_rng(0)
    cc = np.linspace(0, 140, 300)
    # --- calibration : 3 lectures principales
    lect = [l for l in LECTURES if l[0] in ("R", "tau_mono", "tau_m")]
    fig, axes = plt.subplots(1, len(lect), figsize=(5.4 * len(lect), 4.8))
    for ax, (col, nom, unite) in zip(axes, lect):
        G.style(ax)
        row = CAL[CAL.lecture == col].iloc[0]
        for cond in conds:
            s = S[S.condition == cond]
            x = C.CL_MM[cond] if cond != "basal" else np.nan
            if cond == "basal":
                continue
            j = rng.normal(0, 1.6, len(s))
            ax.plot(x + j, s[col], "o", ms=3, alpha=0.45, color=C.COULEUR_CONDITION[cond], mec="none")
            fm = s.groupby("acquisition")[col].median()
            ax.plot(np.full(len(fm), x) + 5, fm, "_", ms=10, mew=1.6, color=C.COULEUR_CONDITION[cond])
            ax.plot(x, s[col].median(), "D", ms=7, color=G.INK, mfc=C.COULEUR_CONDITION[cond], zorder=5)
        if np.isfinite(row.get("K_mM", np.nan)):
            B = boot[col]
            ok = np.isfinite(B[:, 2])
            courbes = np.array([isotherme(cc, *b[:3]) for b in B[ok][:400]])
            ax.fill_between(cc, *np.percentile(courbes, [2.5, 97.5], axis=0), color=G.TRACE, alpha=0.35, lw=0)
            ax.plot(cc, isotherme(cc, row.X0, row.Xinf, row.K_mM), color=G.INK, lw=1.4)
            if "basal" in conds and np.isfinite(row.get("basal_cl_mM", np.nan)):
                xb, cb = row.basal_X, row.basal_cl_mM
                ax.errorbar(cb, xb, xerr=[[max(cb - row.basal_cl_ic_bas, 0)], [max(row.basal_cl_ic_haut - cb, 0)]],
                            fmt="s", ms=8, color=C.COULEUR_CONDITION["basal"], capsize=3, zorder=6,
                            label=f"basal : {cb:.1f} mM [{row.basal_cl_ic_bas:.1f} ; {row.basal_cl_ic_haut:.1f}]")
                ax.legend(frameon=False, fontsize=8.5, loc="best")
            sous = (f"K = {row.K_mM:.1f} mM [{row.K_mM_ic_bas:.1f} ; {row.K_mM_ic_haut:.1f}], "
                    f"X₀ = {row.X0:.3g}, X∞ = {row.Xinf:.3g}")
        else:
            sous = "pas d'isotherme : médianes non monotones"
        ax.set_xlabel("[Cl⁻] imposée (mM)", color=G.INK)
        ax.set_ylabel(f"{nom}" + (f" ({unite})" if unite else ""), color=G.INK)
        G.titre(ax, nom[0].upper() + nom[1:], sous, fontsize=11)
    fig.text(0.01, 0.005, "points : somas ; traits : médiane de chaque FOV ; losanges : médiane de la "
             "condition ; bande : IC 95 % (bootstrap hiérarchique) ; carré : médiane du basal, inversée",
             color=G.INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    G.sauver(fig, D["results"] / "calibration.png")

    # --- basal cellule par cellule et auto-cohérence
    lect = [l for l in LECTURES if l[0] in ("R", "tau_mono")]
    fig, axes = plt.subplots(1, len(lect), figsize=(6.4 * len(lect), 4.6))
    bins = np.linspace(-20, 160, 73)
    for ax, (col, nom, _) in zip(axes, lect):
        G.style(ax)
        for cond in [c for c in conds if c != "basal"]:
            v = AUTO[(AUTO.lecture == col) & (AUTO.condition == cond)].cl_relu_mM
            v = np.clip(v[np.isfinite(v)], bins[0], bins[-1])
            ax.hist(v, bins=bins, histtype="step", lw=1.3, color=C.COULEUR_CONDITION[cond],
                    label=f"{cond} relus : médiane {np.median(v):.1f}, IQR {np.subtract(*np.percentile(v, [75, 25])):.1f} mM")
        v = BAS[BAS.lecture == col].cl_mM
        v = np.clip(v[np.isfinite(v)], bins[0], bins[-1])
        ax.hist(v, bins=bins, color=C.COULEUR_CONDITION["basal"], alpha=0.55,
                label=f"basal : médiane {np.median(v):.1f}, IQR {np.subtract(*np.percentile(v, [75, 25])):.1f} mM")
        ax.set_xlabel("[Cl⁻] lue par la calibration (mM ; bornée à [−20, 160])", color=G.INK)
        ax.set_ylabel("somas", color=G.INK)
        ax.legend(frameon=False, fontsize=8)
        G.titre(ax, f"Lecture cellule par cellule — {nom}",
                "contours : cellules à [Cl⁻] imposée, relues par la courbe (précision d'une cellule seule)",
                fontsize=11)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "basal.png")

    # --- rapport contre τ, soma par soma
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    G.style(ax)
    for cond in conds:
        s = S[S.condition == cond]
        ax.plot(s.tau_mono, s.R, "o", ms=3.5, alpha=0.6, color=C.COULEUR_CONDITION[cond], mec="none",
                label=f"{cond} (n = {len(s)})")
    ax.set_xlabel("τ CFP, mono-exponentielle (ns)", color=G.INK)
    ax.set_ylabel("rapport YFP/CFP brut", color=G.INK)
    ax.legend(frameon=False, fontsize=9)
    G.titre(ax, "Les deux lectures du FRET, soma par soma", "moins de FRET : rapport plus bas, τ CFP plus long")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "rapport_vs_tau.png")

    # --- temps mort
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax in axes:
        G.style(ax)
    for cond in conds:
        s = S[S.condition == cond]
        axes[0].plot(s.mu_yfp, (s.R_corrige / s.R - 1) * 100, "o", ms=3, alpha=0.6, color=C.COULEUR_CONDITION[cond],
                     mec="none", label=cond)
        axes[1].plot(s.mu_cfp, s.mu_yfp, "o", ms=3, alpha=0.6, color=C.COULEUR_CONDITION[cond], mec="none", label=cond)
    axes[0].set_xlabel("μ YFP : photons enregistrés par impulsion (moyenne pondérée du soma)", color=G.INK)
    axes[0].set_ylabel("rapport corrigé / brut − 1 (%)", color=G.INK)
    G.titre(axes[0], "Effet du temps mort sur le rapport", "au plus un photon par impulsion et par entrée")
    axes[1].set_xlabel("μ CFP", color=G.INK)
    axes[1].set_ylabel("μ YFP", color=G.INK)
    G.titre(axes[1], "Taux par impulsion des deux canaux", "pileup/ : τ biaisé de −4 à −8 % par 0,1 photon/impulsion")
    axes[1].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "temps_mort.png")


if __name__ == "__main__":
    main()
