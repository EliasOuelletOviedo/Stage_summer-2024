"""Étape 6 — variation du rapport à chaque stim, par ROI : latence, τ d'activation, τ de retour.

Usage :  python3 step6_cinetique.py <dossier>          (~1 min ; après l'étape 3)

Une série par ROI : les passes de même code de routage (étape 1). En mode ROI (spirales des
galvos), le code est l'indice de la ROI visitée. Sans mode ROI, toutes les passes portent le
code 1 et la série est celle de tout ce que le laser a balayé : le 2026-11-11, le champ
entier du polygone (FLIMApp n'enregistre pas les horloges de ligne et de trame en mode
Realtime par passes, il n'y a donc pas d'image à découper en ROIs).

Pour chaque ROI : ln R passe par passe, redressé par la base quadratique ajustée avec une
réponse M1 (étape 4), en %, rapporté aux 10 s avant chaque stim. Puis, sur chaque stim et
sur la moyenne des stims, un modèle marche-retour avec latence :

    x(t) = 0                                        t < δ
         = A·(1 − exp(−(t − δ)/τ_act))              δ ≤ t < δ + T      (activation)
         = A·(1 − exp(−T/τ_act))·exp(−(t − δ − T)/τ_ret)    t ≥ δ + T  (retour)

t compté depuis le début de la stim, T sa durée (~10 s), plus un décalage constant libre.
δ : latence ; τ_act : constante de temps de la descente pendant la stim ; τ_ret : celle de
la remontée après. Moindres carrés pondérés par le bruit de photons (étape 3), plusieurs
départs ; erreurs types par le jacobien (χ²/ddl compris) ; IC à 95 % des τ en log-normal.
Ajusté sur chaque stim (sa fenêtre : de 20 s avant à 110 s après son début, sans entrer
dans les 10 s qui précèdent la stim suivante ni dépasser FIN_S), puis sur toutes les stims
ensemble : δ, τ_act et τ_ret communs, amplitude et décalage propres à chaque stim (ligne
« stim 0 » du CSV), IC à 95 % par rééchantillonnage des stims (bootstrap). Une constante est
« mesurable » si A < 0 (une baisse du rapport) et si son IC à 95 % tient dans [0,5 ; 290] s.
Avec une stim de 10 s, τ_act n'est mesurable que s'il est nettement plus court que la stim :
plus long, la descente est une rampe dont seule la pente (A / τ_act) se lit.

Sorties : results/cinetique.csv, results/cinetique.png
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import config as C
import figures as G
import spc as S
from figures import plt
from step3_series import FIN_S, blocs_stim
from step4_reponses import ajuster

BORNES = dict(delta=(0.0, 15.0), tau=(0.3, 300.0))
MESURABLE = (0.5, 290.0)
AVANT, APRES = 20.0, 110.0          # fenêtre autour de chaque stim (s)
N_BOOT = 500
GRAINE = 20261111


def modele(p, tp, T):
    c0, A, d, ta, tr = p
    x = np.zeros_like(tp)
    m1 = (tp >= d) & (tp < d + T)
    x[m1] = A * (1 - np.exp(-(tp[m1] - d) / ta))
    m2 = tp >= d + T
    x[m2] = A * (1 - np.exp(-T / ta)) * np.exp(-(tp[m2] - d - T) / tr)
    return c0 + x


def ajuster_marche(tp, y, s, T):
    lb = [-np.inf, -np.inf, BORNES["delta"][0], BORNES["tau"][0], BORNES["tau"][0]]
    ub = [np.inf, np.inf, BORNES["delta"][1], BORNES["tau"][1], BORNES["tau"][1]]
    best = None
    for d0 in (0.0, 2.0, 5.0):
        for ta0 in (2.0, 10.0, 50.0):
            for tr0 in (5.0, 20.0, 60.0):
                r = least_squares(lambda p: (y - modele(p, tp, T)) / s, [0.0, -1.0, d0, ta0, tr0],
                                  bounds=(lb, ub), x_scale="jac")
                if best is None or r.cost < best.cost:
                    best = r
    n, k = len(y), len(best.x)
    chi2 = float((best.fun ** 2).sum())
    try:
        cov = np.linalg.pinv(best.jac.T @ best.jac) * chi2 / (n - k)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        se = np.full(k, np.nan)
    c0, A, d, ta, tr = best.x
    amp = A * (1 - np.exp(-T / ta))
    ic = lambda tau, e: (tau * np.exp(-1.96 * e / tau), tau * np.exp(1.96 * e / tau)) if np.isfinite(e) and tau > 0 else (np.nan, np.nan)
    mes = lambda tau, e: bool(A < 0 and MESURABLE[0] < ic(tau, e)[0] and ic(tau, e)[1] < MESURABLE[1])
    return dict(decalage=c0, A=A, latence_s=d, se_latence=se[2], tau_act_s=ta, se_tau_act=se[3],
                tau_act_ic95=ic(ta, se[3]), tau_ret_s=tr, se_tau_ret=se[4], tau_ret_ic95=ic(tr, se[4]),
                amplitude_pct=amp, pente_initiale_pct_s=A / ta, chi2_ddl=chi2 / (n - k),
                act_mesurable=mes(ta, se[3]), ret_mesurable=mes(tr, se[4]), p=best.x)


def ajuster_conjoint(segs, depart=None):
    """δ, τ_act, τ_ret communs ; A_k et c_k propres à chaque stim. segs : [(tp, yy, sp, T)]."""
    K = len(segs)
    lb = [BORNES["delta"][0], BORNES["tau"][0], BORNES["tau"][0]] + [-np.inf] * (2 * K)
    ub = [BORNES["delta"][1], BORNES["tau"][1], BORNES["tau"][1]] + [np.inf] * (2 * K)

    def res(q):
        d, ta, tr = q[:3]
        return np.concatenate([(yy - modele([q[3 + K + k], q[3 + k], d, ta, tr], tp, T)) / sp
                               for k, (tp, yy, sp, T) in enumerate(segs)])

    departs = [depart] if depart is not None else [
        [d0, ta0, tr0] + [-1.0] * K + [0.0] * K for d0 in (0.0, 2.0, 5.0) for ta0 in (2.0, 10.0, 50.0)
        for tr0 in (5.0, 20.0, 60.0)]
    best = None
    for q0 in departs:
        r = least_squares(res, q0, bounds=(lb, ub), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    n, k = len(best.fun), len(best.x)
    chi2 = float((best.fun ** 2).sum())
    try:
        se = np.sqrt(np.diag(np.linalg.pinv(best.jac.T @ best.jac) * chi2 / (n - k)))
    except np.linalg.LinAlgError:
        se = np.full(k, np.nan)
    d, ta, tr = best.x[:3]
    A = best.x[3:3 + K]
    amps = A * (1 - np.exp(-np.array([sg[3] for sg in segs]) / ta))
    return dict(q=best.x, latence_s=d, se_latence=se[0], tau_act_s=ta, se_tau_act=se[1], tau_ret_s=tr,
                se_tau_ret=se[2], amplitude_pct=float(amps.mean()), amplitudes=amps, A=float(A.mean()),
                pente_initiale_pct_s=float(A.mean() / ta), chi2_ddl=chi2 / (n - k))


def main():
    root = C.racine()
    D = C.dirs(root)
    P = C.charger_passes(D)
    B = pd.read_csv(D["results"] / "bruit.csv")
    rng = np.random.default_rng(GRAINE)
    t_all, stim_all, roi_all = P["t_s"], P["stim"].astype(bool), P["roi_passe"]
    nY_all = P["H_yfp"][:, S.ACQ].sum(1).astype(float)
    nC_all = P["H_cfp"][:, S.ACQ].sum(1).astype(float)
    lignes, traces = [], {}
    for roi in [r for r in np.unique(roi_all) if r > 0]:
        sel = (roi_all == roi) & (t_all < FIN_S)
        t, stim = t_all[sel], stim_all[sel]
        nY, nC = nY_all[sel], nC_all[sel]
        x = np.log(nY / nC)
        s = np.sqrt(1 / nY + 1 / nC)
        dt = float(np.median(np.diff(t)))
        r = ajuster(t, x, s, stim.astype(float), "M1", dt)
        y = 100 * (x - (r["pred"] - r["reponse"]))                  # %, base retirée
        sp = 100 * s
        blocs = blocs_stim(stim)
        debuts = [t[a] for a, _ in blocs]
        alignes = []
        for k, (a, b) in enumerate(blocs, start=1):
            t0, T = t[a], t[b - 1] - t[a] + dt
            fin_fen = debuts[k] - C.PRE_S if k < len(blocs) else FIN_S
            m = (t >= t0 - AVANT) & (t < min(t0 + APRES, fin_fen))
            tp = t[m] - t0
            yy = y[m] - y[(t >= t0 - C.PRE_S) & (t < t0)].mean()
            f = ajuster_marche(tp, yy, sp[m], T)
            lignes.append(dict(roi=int(roi), stim=k, **{kk: v for kk, v in f.items() if kk != "p"}))
            alignes.append((tp, yy, sp[m], f["p"], T))
        # toutes les stims ensemble : δ, τ_act, τ_ret communs
        segs = [(tp, yy, spm, T) for tp, yy, spm, _, T in alignes]
        f = ajuster_conjoint(segs)
        boot = []
        for _ in range(N_BOOT):
            ii = rng.integers(len(segs), size=len(segs))
            q0 = list(f["q"][:3]) + list(f["q"][3:3 + len(segs)][ii]) + list(f["q"][3 + len(segs):][ii])
            fb = ajuster_conjoint([segs[i] for i in ii], depart=q0)
            boot.append([fb["latence_s"], fb["tau_act_s"], fb["tau_ret_s"], fb["amplitude_pct"]])
        q = np.percentile(np.array(boot), [2.5, 97.5], axis=0)
        ic = lambda tau, e: (tau * np.exp(-1.96 * e / tau), tau * np.exp(1.96 * e / tau)) if np.isfinite(e) else (np.nan, np.nan)
        mes = lambda lo_hi: bool(f["A"] < 0 and MESURABLE[0] < lo_hi[0] and lo_hi[1] < MESURABLE[1])
        lignes.append(dict(roi=int(roi), stim=0, decalage=np.nan, A=f["A"], latence_s=f["latence_s"],
                           se_latence=f["se_latence"], tau_act_s=f["tau_act_s"], se_tau_act=f["se_tau_act"],
                           tau_act_ic95=ic(f["tau_act_s"], f["se_tau_act"]), tau_ret_s=f["tau_ret_s"],
                           se_tau_ret=f["se_tau_ret"], tau_ret_ic95=ic(f["tau_ret_s"], f["se_tau_ret"]),
                           amplitude_pct=f["amplitude_pct"], pente_initiale_pct_s=f["pente_initiale_pct_s"],
                           chi2_ddl=f["chi2_ddl"], act_mesurable=mes(tuple(q[:, 1])), ret_mesurable=mes(tuple(q[:, 2])),
                           latence_boot_ic95=tuple(q[:, 0]), tau_act_boot_ic95=tuple(q[:, 1]),
                           tau_ret_boot_ic95=tuple(q[:, 2]), amplitude_boot_ic95=tuple(q[:, 3]),
                           amplitudes_par_stim=tuple(np.round(f["amplitudes"], 3))))
        # moyenne des stims alignées, pour la figure (instants communs à au moins deux stims)
        g = np.round(np.arange(-AVANT, APRES, dt) / dt).astype(int)
        M = np.full((len(alignes), len(g)), np.nan)
        for i, (tp, yy, _, _, _) in enumerate(alignes):
            kk = np.round(tp / dt).astype(int)
            j = np.searchsorted(g, kk)
            ok = (j < len(g)) & (g[np.minimum(j, len(g) - 1)] == kk)
            M[i, j[ok]] = yy[ok]
        garde = (~np.isnan(M)).sum(0) >= 2
        T_moy = float(np.mean([a[4] for a in alignes]))
        p_moy = [0.0, f["A"], f["latence_s"], f["tau_act_s"], f["tau_ret_s"]]
        traces[int(roi)] = dict(alignes=alignes, tg=g[garde] * dt, moy=np.nanmean(M[:, garde], 0), p=p_moy, T=T_moy)
        print(f"  ROI {roi} : {len(blocs)} stims, base M1 retirée (χ²/ddl {r['chi2_ddl']:.2f})")
    L = pd.DataFrame(lignes)
    sortie = L.copy()
    for col in [c for c in L.columns if c.endswith("ic95")]:                # intervalles : deux colonnes
        sortie[f"{col}_bas"] = [float(v[0]) if isinstance(v, tuple) else np.nan for v in L[col]]
        sortie[f"{col}_haut"] = [float(v[1]) if isinstance(v, tuple) else np.nan for v in L[col]]
        sortie = sortie.drop(columns=col)
    if "amplitudes_par_stim" in sortie:
        sortie["amplitudes_par_stim"] = [" / ".join(f"{float(x):.3f}" for x in v) if isinstance(v, tuple) else ""
                                         for v in L["amplitudes_par_stim"]]
    sortie.to_csv(D["results"] / "cinetique.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["roi", "stim", "amplitude_pct", "latence_s", "se_latence", "tau_act_s", "se_tau_act", "tau_ret_s",
            "se_tau_ret", "pente_initiale_pct_s", "chi2_ddl", "act_mesurable", "ret_mesurable"]
    print(L[cols].round(3).to_string(index=False))
    m = L[L.stim == 0].iloc[0]
    print(f"  toutes les stims ensemble, IC à 95 % par bootstrap des stims : latence {np.round(m.latence_boot_ic95, 1)} s, "
          f"τ_act {np.round(m.tau_act_boot_ic95, 1)} s, τ_ret {np.round(m.tau_ret_boot_ic95, 1)} s, "
          f"amplitude {np.round(m.amplitude_boot_ic95, 2)} %")
    figure(D, traces, L)


def figure(D, traces, L):
    for roi, tr in traces.items():
        n = len(tr["alignes"])
        fig, axes = plt.subplots(1, n + 2, figsize=(4.2 * (n + 2), 4.4), sharey=False)
        for ax in axes:
            G.style(ax)
        for i, (tp, yy, sp, p, T) in enumerate(tr["alignes"]):
            ax = axes[i]
            k = max(int(round(1.0 / np.median(np.diff(tp)))), 1)
            nb = len(tp) // k * k
            ax.plot(tp, yy, color=G.TRACE, lw=0.5)
            ax.plot(tp[:nb].reshape(-1, k).mean(1), yy[:nb].reshape(-1, k).mean(1), color=G.INK, lw=1.1)
            g = np.linspace(tp.min(), tp.max(), 800)
            ax.plot(g, modele(p, g, T), color=G.PALETTE[1], lw=1.6)
            ax.axvspan(0, T, color=G.STIM, alpha=0.45, lw=0, zorder=0)
            ax.axhline(0, color=G.TRACE, lw=0.8)
            l = L[(L.roi == roi) & (L.stim == i + 1)].iloc[0]
            ax.set_xlabel("temps depuis le début de la stim (s)", color=G.INK)
            if i == 0:
                ax.set_ylabel("Δ ln R (%)", color=G.INK)
            G.titre(ax, f"stim {i + 1}", f"δ {l.latence_s:.1f} s · τ_act {l.tau_act_s:.1f} s · τ_ret {l.tau_ret_s:.1f} s",
                    fontsize=10.5)
        ax = axes[n]
        ax.plot(tr["tg"], tr["moy"], color=G.TRACE, lw=0.6)
        k = max(int(round(1.0 / np.median(np.diff(tr["tg"])))), 1)
        nb = len(tr["tg"]) // k * k
        ax.plot(tr["tg"][:nb].reshape(-1, k).mean(1), tr["moy"][:nb].reshape(-1, k).mean(1), color=G.INK, lw=1.2)
        g = np.linspace(tr["tg"].min(), tr["tg"].max(), 800)
        ax.plot(g, modele(tr["p"], g, tr["T"]), color=G.PALETTE[1], lw=1.8)
        ax.axvspan(0, tr["T"], color=G.STIM, alpha=0.45, lw=0, zorder=0)
        ax.axhline(0, color=G.TRACE, lw=0.8)
        l = L[(L.roi == roi) & (L.stim == 0)].iloc[0]
        ax.set_xlabel("temps depuis le début de la stim (s)", color=G.INK)
        G.titre(ax, "moyenne des stims ; modèle commun", f"δ {l.latence_s:.1f} s · τ_act {l.tau_act_s:.1f} s · "
                f"τ_ret {l.tau_ret_s:.1f} s · {l.amplitude_pct:+.2f} %", fontsize=10.5)
        ax = axes[n + 1]
        l = L[L.roi == roi]
        x = np.arange(len(l))
        lab = [f"stim {s}" if s else "toutes" for s in l.stim]
        for j, (col, se, c, nom, mes) in enumerate((("tau_act_s", "se_tau_act", G.PALETTE[0], "τ d'activation", "act_mesurable"),
                                                    ("tau_ret_s", "se_tau_ret", G.PALETTE[1], "τ de retour", "ret_mesurable"),
                                                    ("latence_s", "se_latence", G.PALETTE[2], "latence δ", None))):
            ok = l[mes].values.astype(bool) if mes else np.ones(len(l), bool)
            xx = x + (j - 1) * 0.22
            ax.errorbar(xx[ok], l[col][ok], yerr=1.96 * l[se][ok], fmt="o", color=c, ms=5, capsize=2, lw=1, label=nom)
            if (~ok).any():
                ax.plot(xx[~ok], np.minimum(l[col][~ok], 58), "o", mfc="none", color=c, ms=6,
                        label=f"{nom} non mesurable" if j == 0 else None)
        ax.set_xticks(x)
        ax.set_xticklabels(lab, fontsize=8.5)
        ax.set_ylabel("s", color=G.INK)
        ax.set_ylim(0, 60)
        ax.legend(frameon=False, fontsize=8)
        G.titre(ax, "Constantes de temps", "± 1,96 erreur type ; cercles vides : non mesurable (au plus 58 s affiché)", fontsize=10.5)
        fig.tight_layout()
        G.sauver(fig, D["results"] / ("cinetique.png" if len(traces) == 1 else f"cinetique_roi{roi}.png"))


if __name__ == "__main__":
    main()
