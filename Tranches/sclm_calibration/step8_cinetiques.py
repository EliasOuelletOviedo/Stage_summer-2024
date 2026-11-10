"""Étape 8 — cinétiques des changements de solution (time-lapses).

Usage :  python3 step8_cinetiques.py <dossier>

Chaque transition est un time-lapse lancé au changement de solution (t = 0 au début du
premier cycle, d'après l'expérimentateur) : un cycle de 10 s toutes les ~41 s, sur un champ
qui n'est aucun des FOV des conditions fixes. Pour chaque cycle gardé (étape 1 : durée,
plan focal) :

  - somas réunis (photons de tous les somas gardés) : rapport YFP/CFP brut et corrigé du
    temps mort, τ CFP mono-exponentielle (décalage fixé sur le déclin de toute la
    transition), τm, τ des trois populations (τ communes de l'étape 5) ;
  - chaque soma : rapport et τ CFP (décalage fixé sur son propre déclin total) ;
  - à la seconde près, dans chaque cycle : rapport des somas réunis.

Ajustement : X(t) = X∞ + (X₀ − X∞)·exp(−t/τ_éq) sur les somas réunis et sur chaque soma,
au moindre carré pondéré par l'incertitude de Poisson de chaque point. Pour basal→0,
qui ne suit pas une exponentielle (le rapport descend d'abord), une somme de deux
exponentielles est aussi ajustée. Les FOV de la condition d'arrivée ont été pris plus tard :
la fraction de la transition accomplie à leur début est reportée.

Sorties : results/transitions_cycles.csv, results/transitions_somas.csv,
results/transitions_ajustements.csv, results/transition_<nom>.png, results/transitions.png
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import config as C
import figures as G
from figures import plt
from flim import Modele
from mesures import comptes_cycles, selection


def expo(t, X0, Xinf, tau):
    return Xinf + (X0 - Xinf) * np.exp(-t / tau)


def expo2(t, Xinf, A1, tau1, A2, tau2):
    return Xinf + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)


def ajuster_expo(t, x, sx):
    ok = np.isfinite(x) & np.isfinite(sx) & (sx > 0)
    t, x, sx = t[ok], x[ok], sx[ok]
    if len(t) < 5:
        return None
    best = None
    for tau0 in (0.5, 2, 5, 15, 40):
        try:
            p, cov = curve_fit(expo, t, x, p0=[x[0], x[-1], tau0], sigma=sx, absolute_sigma=True,
                               bounds=([-np.inf, -np.inf, 0.02], [np.inf, np.inf, 500]), maxfev=20000)
        except Exception:
            continue
        chi2 = float((((x - expo(t, *p)) / sx) ** 2).sum())
        if best is None or chi2 < best[2]:
            best = (p, cov, chi2)
    if best is None:
        return None
    p, cov, chi2 = best
    if p[2] >= 499:                       # τ à la borne : pas de plateau dans la fenêtre, échec
        return dict(echec="pas de plateau (τ_éq à la borne de 500 min)", n=len(t))
    e = np.sqrt(np.clip(np.diag(cov), 0, None)) * np.sqrt(max(chi2 / max(len(t) - 3, 1), 1))
    return dict(X0=p[0], Xinf=p[1], tau_min=p[2], tau_err=e[2], Xinf_err=e[1], chi2_ddl=chi2 / max(len(t) - 3, 1), n=len(t))


def ajuster_expo2(t, x, sx):
    ok = np.isfinite(x) & np.isfinite(sx) & (sx > 0)
    t, x, sx = t[ok], x[ok], sx[ok]
    best = None
    for t1, t2 in ((1, 15), (2, 25), (3, 40), (5, 60)):
        for sgn in (1, -1):
            try:
                p, cov = curve_fit(expo2, t, x, p0=[x[-1], sgn * 0.3 * abs(x[0] - x[-1]) + 1e-3, t1,
                                                    (x[0] - x[-1]), t2], sigma=sx, absolute_sigma=True,
                                   bounds=([-np.inf, -np.inf, 0.05, -np.inf, 0.05], [np.inf] * 2 + [300, np.inf, 600]),
                                   maxfev=40000)
            except Exception:
                continue
            chi2 = float((((x - expo2(t, *p)) / sx) ** 2).sum())
            if best is None or chi2 < best[1]:
                best = (p, chi2)
    if best is None:
        return None
    p, chi2 = best
    if p[2] > p[4]:                                    # τ1 la plus courte
        p = np.array([p[0], p[3], p[4], p[1], p[2]])
    return dict(Xinf=p[0], A1=p[1], tau1_min=p[2], A2=p[3], tau2_min=p[4], chi2_ddl=chi2 / max(len(t) - 5, 1), n=len(t))


def sigma_ratio(nc, ny):
    """Écart-type de Poisson du rapport ny/nc."""
    r = ny / np.maximum(nc, 1)
    return r * np.sqrt(1 / np.maximum(ny, 1) + 1 / np.maximum(nc, 1))


def main():
    root = C.racine()
    D = C.dirs(root)
    irf = json.loads((D["qc"] / "irf.json").read_text())
    M = Modele(C.charger_irf(irf["voie_cfp"], root))
    taus = np.array(irf["taus_cfp"])
    s_glob = irf["s_median"]
    T = pd.read_csv(D["results"] / "somas.csv") if (D["results"] / "somas.csv").exists() else None
    A_all = C.acquisitions(root)
    A = [a for a in A_all if a["type"] == "transition" and (D["mesures"] / f"{a['nom']}.npz").exists()]
    if not A:
        sys.exit("Aucune transition mesurée — lance d'abord step4_declins.py")
    cyc_rows, soma_rows, fit_rows, secs = [], [], [], {}
    for a in A:
        z = np.load(D["projections"] / f"{a['nom']}.npz")
        m = np.load(D["mesures"] / f"{a['nom']}.npz")
        n = int(m["n_rois"])
        k_somas = selection(D, a["nom"], n)
        g = m["garde"]
        t0 = m["debut"][0]
        t_mid = (m["debut"] + m["duree_s"] / 2 - t0) / 60                  # minutes
        dec = m["declins"].astype(float)                                   # (cy, n+2, 2, 256)
        cc = comptes_cycles(z, m["etiquettes"], n + 1, m["decalages"])     # (cy, n+2, 2, 3)
        reunis = dec[:, k_somas].sum(1)                                    # (cy, 2, 256)
        cc_r = cc[:, k_somas].sum(1)                                       # (cy, 2, 3)
        tot_r = reunis[g].sum(0)
        _, s_r, _, _ = M.ajuster_mono(tot_r[0])
        for k in range(len(dec)):
            y = reunis[k, 0]
            tau, _, b, dd = M.ajuster_mono(y, s=s_r)
            amp, b3, dd3 = M.amplitudes(y, taus, s_glob)
            nc, ny = y.sum(), reunis[k, 1].sum()
            cyc_rows.append(dict(transition=a["nom"], cycle=k, t_min=t_mid[k], garde=bool(g[k]),
                                 photons_cfp=nc, photons_yfp=ny, R=ny / max(nc, 1), R_sigma=sigma_ratio(nc, ny),
                                 R_corrige=cc_r[k, 1, 1] / max(cc_r[k, 0, 1], 1e-9),
                                 mu_yfp=cc_r[k, 1, 2] / max(cc_r[k, 1, 0], 1e-9),
                                 tau_mono=tau, dev_ddl_mono=dd, tau_m=M.temps_moyen(y, s_glob),
                                 a1=amp[0], a2=amp[1], a3=amp[2],
                                 tau_3pop=float((amp * taus ** 2).sum() / (amp * taus).sum()) if np.all(np.isfinite(amp)) else np.nan))
        # chaque soma
        for lab in k_somas:
            tot = dec[g, lab].sum(0)
            if tot[0].sum() < C.MIN_PHOTONS_FIT:
                continue
            _, s_l, _, _ = M.ajuster_mono(tot[0])
            for k in range(len(dec)):
                nc, ny = dec[k, lab, 0].sum(), dec[k, lab, 1].sum()
                tau = M.ajuster_mono(dec[k, lab, 0], s=s_l)[0] if g[k] else np.nan
                soma_rows.append(dict(transition=a["nom"], roi=int(lab), cycle=k, t_min=t_mid[k], garde=bool(g[k]),
                                      photons_cfp=nc, photons_yfp=ny, R=ny / max(nc, 1), tau_mono=tau))
        # à la seconde
        sec = m["secondes"][:, k_somas].sum(1).astype(float)              # (cy, 2, 11)
        ts = (m["debut"][:, None] + np.arange(sec.shape[2])[None, :] + 0.5 - t0) / 60
        plein = np.arange(sec.shape[2])[None, :] < np.floor(m["duree_s"])[:, None]
        secs[a["nom"]] = (ts, sec, plein & g[:, None])
        # ajustements
        Cc = pd.DataFrame([r for r in cyc_rows if r["transition"] == a["nom"]])
        Ck = Cc[Cc.garde]
        for col in ("R", "R_corrige", "tau_mono", "tau_m", "tau_3pop"):
            if col.startswith("R"):
                sx = Ck.R_sigma.values * (Ck[col].values / Ck.R.values)
            else:
                # incertitude de τ par cycle : dispersion entre cycles voisins (différences successives)
                v = Ck[col].values
                sx = np.full(len(v), np.nanstd(np.diff(v)) / np.sqrt(2) if len(v) > 3 else np.nan)
            f = ajuster_expo(Ck.t_min.values, Ck[col].values, sx)
            row = dict(transition=a["nom"], lecture=col, modele="1 exponentielle", **(f or {}))
            fit_rows.append(row)
            if a["nom"].startswith("basal") and col in ("R", "tau_mono"):
                f2 = ajuster_expo2(Ck.t_min.values, Ck[col].values, sx)
                if f2:
                    fit_rows.append(dict(transition=a["nom"], lecture=col, modele="2 exponentielles", **f2))
        # somas un par un (rapport)
        Sk = pd.DataFrame([r for r in soma_rows if r["transition"] == a["nom"] and r["garde"]])
        for lab, s in Sk.groupby("roi"):
            f = ajuster_expo(s.t_min.values, s.R.values, sigma_ratio(s.photons_cfp.values, s.photons_yfp.values))
            if f:
                fit_rows.append(dict(transition=a["nom"], lecture="R", modele="1 exponentielle, soma", roi=int(lab), **f))
        # fraction accomplie au début des FOV de la condition d'arrivée
        fovs = [b for b in A_all if b["type"] == "fov" and b["condition"] == a["apres"]]
        if fovs:
            zf = [np.load(D["projections"] / f"{b['nom']}.npz")["debut"][0] for b in fovs]
            for r in fit_rows:
                if r["transition"] == a["nom"] and r.get("modele") == "1 exponentielle" and "tau_min" in r:
                    r["fov_debut_min"] = (min(zf) - t0) / 60
                    r["fov_fin_min"] = (max(zf) - t0) / 60
                    r["fraction_au_premier_fov"] = 1 - np.exp(-r["fov_debut_min"] / r["tau_min"])
        fr = [r for r in fit_rows if r["transition"] == a["nom"] and r.get("lecture") == "R" and r.get("modele") == "1 exponentielle"]
        if fr and "tau_min" in fr[0] and np.isfinite(fr[0]["tau_min"]):
            print(f"  {a['nom']:>10} : {len(k_somas)} somas, {g.sum()}/{len(g)} cycles ; rapport "
                  f"{Ck.R.iloc[0]:.3f} -> {Ck.R.iloc[-1]:.3f}, τ_éq {fr[0]['tau_min']:.2f} ± {fr[0]['tau_err']:.2f} min, "
                  f"X∞ {fr[0]['Xinf']:.3f}")
    CY = pd.DataFrame(cyc_rows)
    SO = pd.DataFrame(soma_rows)
    FI = pd.DataFrame(fit_rows)
    CY.to_csv(D["results"] / "transitions_cycles.csv", index=False)
    SO.to_csv(D["results"] / "transitions_somas.csv", index=False)
    FI.to_csv(D["results"] / "transitions_ajustements.csv", index=False)
    print("  -> results/transitions_cycles.csv, transitions_somas.csv, transitions_ajustements.csv")
    figures(D, A, CY, SO, FI, secs, T)


def figures(D, A, CY, SO, FI, secs, T):
    for a in A:
        Cc = CY[CY.transition == a["nom"]]
        S = SO[SO.transition == a["nom"]]
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        for ax in axes:
            G.style(ax)
        g = Cc.garde.values
        for ax, col, nom in ((axes[0], "R", "rapport YFP/CFP"), (axes[1], "tau_mono", "τ CFP mono-exp. (ns)")):
            for lab, s in S.groupby("roi"):
                ax.plot(s.t_min[s.garde], s[col][s.garde], color=G.TRACE, lw=0.5, alpha=0.6)
            ax.plot(Cc.t_min[g], Cc[col][g], "o-", color=G.INK, ms=4, lw=1.2, label="somas réunis")
            ax.plot(Cc.t_min[~g], Cc[col][~g], "o", mfc="none", color=G.PALETTE[4], ms=6, label="cycle écarté (plan focal)")
            med = S[S.garde].groupby("cycle")[col].median()
            ax.plot(Cc.t_min.values[med.index], med.values, "-", color=G.PALETTE[0], lw=1.2, label="médiane des somas")
            f = FI[(FI.transition == a["nom"]) & (FI.lecture == col) & (FI.modele == "1 exponentielle")]
            tt = np.linspace(0, Cc.t_min.max(), 300)
            if len(f) and "tau_min" in f.columns and np.isfinite(f.iloc[0]["tau_min"]):
                r = f.iloc[0]
                ax.plot(tt, expo(tt, r.X0, r.Xinf, r.tau_min), "--", color=G.PALETTE[1], lw=1.3,
                        label=f"1 exp. : τ_éq {r.tau_min:.1f} ± {r.tau_err:.1f} min")
            f2 = FI[(FI.transition == a["nom"]) & (FI.lecture == col) & (FI.modele == "2 exponentielles")]
            if len(f2):
                r = f2.iloc[0]
                ax.plot(tt, expo2(tt, r.Xinf, r.A1, r.tau1_min, r.A2, r.tau2_min), ":", color=G.PALETTE[2], lw=1.5,
                        label=f"2 exp. : {r.tau1_min:.1f} et {r.tau2_min:.1f} min")
            if T is not None:
                for cond, ls in ((a["avant"], (0, (4, 3))), (a["apres"], "-")):
                    s = T[(T.condition == cond) & (T.type == "soma") & T.garde]
                    if len(s):
                        ax.axhline(s[col].median(), color=C.COULEUR_CONDITION.get(cond, G.INK_2), lw=1, ls=ls, alpha=0.8)
                        ax.text(Cc.t_min.max(), s[col].median(), f" FOV {cond}", color=C.COULEUR_CONDITION.get(cond, G.INK_2),
                                fontsize=8, va="center")
            ax.set_ylabel(nom, color=G.INK)
            ax.legend(frameon=False, fontsize=8, loc="best")
        ts, sec, ok = secs[a["nom"]]
        r = sec[:, 1] / np.maximum(sec[:, 0], 1)
        axes[0].plot(ts[ok], r[ok], ".", color=G.PALETTE[0], ms=2, alpha=0.5)
        axes[2].plot(Cc.t_min, Cc.photons_cfp / 1e6, "o-", ms=3, color=G.COULEUR_CANAL["cfp"], label="CFP")
        axes[2].plot(Cc.t_min, Cc.photons_yfp / 1e6, "o-", ms=3, color=G.COULEUR_CANAL["yfp"], label="YFP")
        axes[2].set_ylabel("photons par cycle,\nsomas réunis (M)", color=G.INK)
        axes[2].legend(frameon=False, fontsize=9)
        axes[2].set_xlabel("minutes depuis le changement de solution (début du time-lapse)", color=G.INK)
        G.titre(axes[0], f"Transition {a['condition']} mM — {S.roi.nunique()} somas",
                "points fins : rapport des somas réunis seconde par seconde ; traits gris : chaque soma ; "
                "lignes horizontales : médiane des FOV de la condition (autres cellules)")
        fig.tight_layout()
        G.sauver(fig, D["results"] / f"transition_{a['nom']}.png")


if __name__ == "__main__":
    main()
