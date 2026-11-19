"""Étape 5 — peut-on lire une [Cl⁻] absolue avec la calibration en tranche du 2026-11-09 ?

Usage :  python3 step5_absolu.py <dossier>          (~30 s ; après l'étape 3)

La seule calibration en tissu disponible est celle du 2026-11-09 (ionophores, imagerie au
polygone, SPCM ; ../SCLM_calibration_2026-11-09). Trois lectures, au même capteur :
  - le rapport YFP/CFP par son isotherme (X0, X∞, K) — son échelle dépend de η, donc du
    trajet optique et d'une éventuelle YFP sans FRET (eNpHR-YFP) ;
  - τ CFP mono par son isotherme — même définition (IRF décalée d'une valeur ajustée sur
    toute la session puis fixée, fond libre, fenêtre ~0,83–12,4 ns) ;
  - la forme du déclin CFP comparée aux motifs du capteur seul à 0 et 131 mM du 2026-11-09
    (étape 13 là-bas : composante prompte du tissu retirée), ramenés sur les canaux de
    FLIMApp par la fonction cumulée, alignés en temps (décalage ajusté), avec un fond plat
    libre : part f du motif « 131 mM », pour trois débuts de fenêtre.
Les concentrations négatives sont gardées.

Sorties : results/absolu.csv, results/absolu_temps.csv, results/absolu.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import config as C
import figures as G
from figures import plt
from step3_series import FIN_S, blocs_stim, moyennes

CAL = Path(__file__).parent / "../SCLM_calibration_2026-11-09"
P_0911 = 1e9 / 80.2475e6                   # le 2026-11-09 : 256 canaux sur la période, fenêtre [17, 254)
FEN_0911 = (17, 254)
FENETRES = ((36, 240), (40, 240), (45, 240))


def inverser(x, X0, Xi, K):
    """[Cl⁻] d'une lecture par l'isotherme X = X∞ + (X0 − X∞)·K/(K + c), négatifs gardés."""
    return K * (X0 - x) / (x - Xi)


def motif(z, cond, ch="sB"):
    return np.mean([z[f"{cond}_fov{k}|{ch}"] / z[f"{cond}_fov{k}|{ch}"].sum() for k in range(1, 6)], 0)


def vers_grille(m, decalage):
    """Motif du 2026-11-09 (canaux [17, 254) de P/256) → canaux de FLIMApp, décalé de `decalage` ns."""
    bords = np.arange(FEN_0911[0], FEN_0911[1] + 1) * P_0911 / 256 + decalage
    cum = np.r_[0, np.cumsum(m)]
    return np.diff(np.interp(np.arange(C.N_CANAUX + 1) * C.DT_NS, bords, cum, left=0, right=cum[-1]))


def forme(y, m0, m1, lo, hi):
    """f (part du motif 131 mM), fond plat, χ²/ddl et décalage, pour une fenêtre."""
    def fit(dec, tout=False):
        p0, p1 = vers_grille(m0, dec)[lo:hi], vers_grille(m1, dec)[lo:hi]
        p0, p1 = p0 / p0.sum(), p1 / p1.sum()
        b = y[lo:hi] / y[lo:hi].sum()
        X = np.c_[p1 - p0, np.ones(hi - lo) / (hi - lo)]
        w = 1 / np.maximum(b, 1e-12)
        coef = np.linalg.solve((X * w[:, None]).T @ X, (X * w[:, None]).T @ (b - p0))
        r = b - p0 - X @ coef
        chi = float((r ** 2 * w).sum() * y[lo:hi].sum())
        return (chi, coef, np.sqrt(np.mean((r / b) ** 2))) if tout else chi
    o = minimize_scalar(fit, bounds=(-0.6, 0.3), method="bounded")
    chi, coef, rms = fit(o.x, True)
    return dict(f=float(coef[0]), fond=float(coef[1]), chi2_ddl=chi / (hi - lo - 3), rms_pct=100 * rms, decalage_ps=1000 * o.x)


def main():
    root = C.racine()
    D = C.dirs(root)
    T = pd.read_csv(D["results"] / "series_passes.csv")
    cal = pd.read_csv(CAL / "results/calibration.csv").set_index("lecture")
    P = C.charger_passes(D)
    t = T.t_s.values
    blocs = blocs_stim(T.stim.values.astype(bool))
    lignes = []
    base = t < t[blocs[0][0]]
    fin = (t >= FIN_S - 90) & (t < FIN_S) & ~T.stim.values.astype(bool)
    for cle, lec, col in (("rapport", "R", "R"), ("τ CFP mono", "tau_mono", "tau_cfp")):
        X0, Xi, K = cal.loc[lec, ["X0", "Xinf", "K_mM"]]
        for moment, m in (("base (avant la 1re stim)", base), (f"{FIN_S - 90:.0f}–{FIN_S:.0f} s, hors stims", fin)):
            x = T[col][m].mean()
            lignes.append(dict(lecture=cle, moment=moment, valeur=x, X0_0911=X0, Xinf_0911=Xi, K_0911=K,
                               c_mM=inverser(x, X0, Xi, K)))
    # forme du déclin CFP contre les motifs du 2026-11-09
    z = np.load(CAL / "mesures/additif.npz")
    m0, m1 = motif(z, "0mM"), motif(z, "131mM")
    y = P["H_cfp"][:, 0][base].sum(0).astype(float)
    for lo, hi in FENETRES:
        r = forme(y, m0, m1, lo, hi)
        lignes.append(dict(lecture="forme CFP (motifs du capteur seul)", moment=f"base, fenêtre [{lo}, {hi})",
                           valeur=r["f"], c_mM=np.nan, **{k: v for k, v in r.items() if k != "f"}))
    A = pd.DataFrame(lignes)
    A.to_csv(D["results"] / "absolu.csv", index=False)
    pd.set_option("display.width", 230)
    print(A.round(4).to_string(index=False))

    # lectures au fil de la session (moyennes de 5 s)
    tt, R5 = moyennes(t, T.R.values)
    _, tau5 = moyennes(t, T.tau_cfp.values)
    cR = inverser(R5, *cal.loc["R", ["X0", "Xinf", "K_mM"]])
    cT = inverser(tau5, *cal.loc["tau_mono", ["X0", "Xinf", "K_mM"]])
    pd.DataFrame(dict(t_s=tt, R=R5, tau_cfp=tau5, c_rapport_mM=cR, c_tau_mM=cT)).to_csv(D["results"] / "absolu_temps.csv", index=False)
    ok = tt < FIN_S
    for nom, c in (("rapport", cR), ("τ CFP", cT)):
        p = np.polyfit(tt[ok], c[ok], 1)
        print(f"  {nom} lu avec la calibration du 2026-11-09 : {c[ok][:12].mean():+.2f} mM au début, "
              f"pente {100 * p[0]:+.2f} mM/100 s avant {FIN_S:.0f} s")
    figure(D, tt, cR, cT, blocs, t)


def figure(D, tt, cR, cT, blocs, t):
    fig, ax = plt.subplots(figsize=(12, 4.6))
    G.style(ax)
    G.bandes(ax, [(t[a], t[b - 1]) for a, b in blocs])
    ax.plot(tt, cR, color=G.INK, lw=1.4, label="rapport YFP/CFP")
    ax.plot(tt, cT, color=G.COULEUR_CANAL["cfp"], lw=1.4, label="τ CFP mono")
    ax.axhline(0, color=G.TRACE, lw=0.8)
    ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.set_ylabel("[Cl⁻] lue (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=9)
    G.titre(ax, "Lectures absolues avec la calibration en tranche du 2026-11-09 (moyennes de 5 s)",
            "deux lectures du même capteur qui ne s'accordent ni sur la valeur ni sur la dérive : pas de [Cl⁻] absolue ici",
            fontsize=11)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "absolu.png")


if __name__ == "__main__":
    main()
