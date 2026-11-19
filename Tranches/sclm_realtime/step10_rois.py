"""Étape 10 — par ROI : photons de chaque passe, [Cl⁻] du rapport ancré sur τ, avant / après, cinétique.

Usage :  python3 step10_rois.py <dossier>          (~3 min ; après les étapes 3, 7, 8 et 9)

1. Photons d'acquisition (code ≠ 0, dans une passe) rangés par passe, par ROI (masques de
   l'étape 9 ; 0 = hors des ROIs), par canal et par canal de temps. Position de chaque photon
   comme à l'étape 8 (phase des lignes suivie trame par trame), ramenée sur l'image de
   référence par le décalage de son bloc (interpolé entre les centres des blocs de 10 s).
   Contrôle : la somme des ROIs (0 compris) redonne les photons de chaque passe (étape 1).
2. [Cl⁻] de chaque ROI comme à l'étape 7, ancrage unique : déclin CFP réuni des 60 s avant la
   1re stim → τ mono (décalage de l'IRF de la session, étape 3) → c_τ par l'isotherme du
   2026-11-09 ; η du ROI = R de cette fenêtre / R_iso(c_τ) ; chaque fenêtre de temps :
   R_iso⁻¹(R / η), R étant le rapport des photons sommés sur la fenêtre. Pas la moyenne des
   lectures passe par passe : à 100–400 photons CFP par passe, l'inversion (convexe) de
   l'isotherme la biaise vers le haut (+0,06 à +0,47 mM sur la fenêtre de l'ancre).
   Un ROI dont l'ancre dépasse C_MAX_LISIBLE (τ trop près du plateau de l'isotherme, où la
   conversion amplifie le bruit) est écarté des lectures en mM, gardé pour le rapport.
3. Avant / après chaque stim, par ROI : Δ[Cl⁻] (par rapport à l'ancre du ROI), ΔR/R₀ (R₀ : la
   fenêtre de l'ancre) et Δτ CFP (τ mono des déclins réunis de chaque fenêtre de 10 s), valeurs
   brutes ; Wilcoxon entre ROIs. Écarts redressés (colonnes *_redresse) : 10 s après − 10 s
   avant sur ln R redressé par la base de M1 du ROI (comme l'étape 7), convertis en mM au
   niveau des 10 s avant.
4. Cinétique par ROI : modèle marche-retour de l'étape 6, toutes les stims ensemble (latence,
   τ d'activation, τ de retour communs ; ln R redressé par la base de M1 du ROI).
5. Ancrage unique par ROI (comme hybride_ancrage_unique.png de l'étape 7) : [Cl⁻] de chaque
   soma dans le temps, la moyenne des ROIs lisibles en haut.

Sorties : mesures/rois_passes.npz, results/rois_cl.csv, results/rois_avant_apres.csv,
results/rois_cinetique.csv, results/rois_series.png, results/avant_apres_rois.png,
results/avant_apres_rois_reunies.png, results/cinetique_rois.png,
results/hybride_ancrage_unique_rois.png
"""
import time

import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import center_of_mass
from matplotlib import patheffects as pe

import config as C
import figures as G
import flim
import spc as S
from figures import plt
from step3_series import FIN_S, blocs_stim, moyennes
from step4_reponses import ajuster
from step5_absolu import CAL, inverser
from step6_cinetique import ajuster_conjoint, modele
from step7_hybride import dlnR_dc, iso

NCAN = 256
C_MAX_LISIBLE = 40.0      # au-delà, τ de l'ancre trop près du plateau de l'isotherme : [Cl⁻] du ROI non lisible


def ranger(root, D, P, masks):
    """Photons d'acquisition par passe, ROI, canal (0 CFP, 1 YFP) et canal de temps."""
    ses = C.session(root)
    z = np.load(D["mesures"] / "images.npz")
    T0, NB, A0, L0, roul, lignes = float(z["T0"]), int(z["NB"]), float(z["A0"]), int(z["L0"]), int(z["roulement"]), \
        int(z["lignes"])
    tc, shu, tb, dep = z["tc"], z["decalages"], z["t_bloc"], z["dep_bloc"]
    deb, fin = P["t_debut_s"], P["t_fin_s"]
    n_p, n_r = len(deb), int(masks.max()) + 1
    H = np.zeros(n_p * n_r * 2 * NCAN, np.int64)
    for k, nom in enumerate(("cfp", "yfp")):
        for t, ph, _, champ, adc, _ in S._blocs(ses["spc"][C.CANAL[nom]]["spc"]):
            t = t[ph] * ses["tic_s"]
            code, can = champ[ph], (4095 - adc[ph]) // 16
            i = np.searchsorted(deb, t, side="right") - 1
            ok = (i >= 0) & (code != 0)
            ok[ok] &= t[ok] < fin[i[ok]]
            t, can, i = t[ok], can[ok], i[ok]
            u = (t - A0) / T0 * NB - np.interp(t, tc, shu)
            L = np.floor(u / NB).astype(np.int64)
            x = u - L * NB
            y = ((L - L0) % lignes + roul) % lignes
            yr = np.rint(y - np.interp(t, tb, dep[:, 0])).astype(np.int64) % lignes
            xr = np.rint(x - np.interp(t, tb, dep[:, 1])).astype(np.int64) % NB
            roi = masks[yr, xr]
            idx = ((i * n_r + roi) * 2 + k) * NCAN + can
            H += np.bincount(idx, minlength=len(H))
    return H.reshape(n_p, n_r, 2, NCAN).astype(np.int32)


def main():
    root = C.racine()
    D = C.dirs(root)
    P = C.charger_passes(D)
    masks = np.load(root / "masks" / "cellpose_masks.npy")
    n_roi = int(masks.max())
    t0s = time.time()
    p_rois = D["mesures"] / "rois_passes.npz"
    if p_rois.exists() and np.load(p_rois)["n_roi"] == n_roi:
        H = np.load(p_rois)["H"]
    else:
        H = ranger(root, D, P, masks)
        np.savez_compressed(p_rois, H=H, n_roi=n_roi)
    tot = H.sum((1, 3))
    ok_c = np.array_equal(tot[:, 0], P["H_cfp"][:, 0].sum(1)) and np.array_equal(tot[:, 1], P["H_yfp"][:, 0].sum(1))
    print(f"  photons rangés par ROI ({time.time() - t0s:.0f} s) ; contrôle contre l'étape 1 : "
          f"{'identique' if ok_c else 'DIFFÉRENT'} ; part des photons dans les ROIs : "
          f"{100 * H[:, 1:].sum() / H.sum():.1f} %")

    cal = pd.read_csv(CAL / "results/calibration.csv").set_index("lecture")
    pR = cal.loc["R", ["X0", "Xinf", "K_mM"]].values.astype(float)
    pT = cal.loc["tau_mono", ["X0", "Xinf", "K_mM"]].values.astype(float)
    DS = pd.read_csv(D["results"] / "declins_session.csv")
    s_cfp = float(DS[(DS.canal == "cfp") & (DS.modele == "1 exp.")].s_ps.iloc[0]) / 1000
    B = pd.read_csv(D["results"] / "bruit.csv")
    sd_tau_champ = float(B[B.lecture == "tau_cfp"].sd_passe.iloc[0])
    n_champ = float(np.median(P["H_cfp"][:, 0].sum(1)))
    M = flim.Modele(C.irf(root, "cfp"))
    t, stim = P["t_s"], P["stim"].astype(bool)
    dt = float(np.median(np.diff(t)))
    blocs = blocs_stim(stim)
    anc = t < t[blocs[0][0]]
    fen = [(t[a], t[b - 1] + dt) for a, b in blocs if t[b - 1] + dt + C.POST_S <= FIN_S]
    tau_fen = lambda h: M.ajuster_mono(h.astype(float), s=s_cfp, tau0=1.5)[0]

    lignes_cl, lignes_aa, series = [], [], {}
    cy_cx = center_of_mass(np.ones_like(masks), masks, range(1, n_roi + 1))
    for r in range(1, n_roi + 1):
        nC, nY = H[:, r, 0].sum(1).astype(float), H[:, r, 1].sum(1).astype(float)
        R = nY / np.maximum(nC, 1)
        dC = H[:, r, 0].astype(float)
        N_anc = dC[anc].sum()
        tau_a = tau_fen(dC[anc].sum(0))
        se_tau = sd_tau_champ * np.sqrt(n_champ / max(N_anc / anc.sum(), 1)) / np.sqrt(anc.sum())
        c_tau = inverser(tau_a, *pT)
        se_c = abs((inverser(tau_a + 1e-4, *pT) - inverser(tau_a - 1e-4, *pT)) / 2e-4) * se_tau
        R_a = nY[anc].sum() / nC[anc].sum()
        eta = R_a / iso(c_tau, *pR)
        Rf = lambda m: nY[m].sum() / nC[m].sum()          # rapport des photons sommés sur la fenêtre m
        cf = lambda m: inverser(Rf(m) / eta, *pR)          # sa [Cl⁻] (pas la moyenne des lectures par passe)
        series[r] = dict(R=R, nC=nC, nY=nY, eta=eta)
        lisible = bool(c_tau < C_MAX_LISIBLE)
        lignes_cl.append(dict(roi=r, y=cy_cx[r - 1][0], x=cy_cx[r - 1][1], aire_px=int((masks == r).sum()),
                              photons_cfp_par_passe=nC.mean(), photons_yfp_par_passe=nY.mean(), R_ancre=R_a,
                              tau_ancre_ns=tau_a, se_tau_ns=se_tau, c_ancre_mM=c_tau, se_c_mM=se_c, eta=eta,
                              lisible_mM=lisible))
        # avant / après, par stim
        for k, (a, b) in enumerate(fen, start=1):
            mb, ma = (t >= a - C.PRE_S) & (t < a), (t >= b) & (t < b + C.POST_S)
            lignes_aa.append(dict(roi=r, stim=k, lisible_mM=lisible, c_avant=cf(mb) - c_tau, c_apres=cf(ma) - c_tau,
                                  R_avant=100 * (Rf(mb) / R_a - 1), R_apres=100 * (Rf(ma) / R_a - 1),
                                  tau_avant_ps=1000 * (tau_fen(dC[mb].sum(0)) - tau_a),
                                  tau_apres_ps=1000 * (tau_fen(dC[ma].sum(0)) - tau_a)))
    CL = pd.DataFrame(lignes_cl)
    AA = pd.DataFrame(lignes_aa)
    for col in ("c", "R", "tau"):
        a, b = (f"{col}_avant", f"{col}_apres") if col != "tau" else ("tau_avant_ps", "tau_apres_ps")
        AA[f"d_{col}"] = AA[b] - AA[a]
    CL.to_csv(D["results"] / "rois_cl.csv", index=False)
    pd.set_option("display.width", 220)
    print(CL.round(3).to_string(index=False))

    # cinétique par ROI : toutes les stims ensemble ; écarts aux stims, dérive lente retirée
    lignes_k, traces, redr = [], {}, []
    for r in range(1, n_roi + 1):
        s = series[r]
        av = t < FIN_S
        x = np.log(s["nY"][av] / np.maximum(s["nC"][av], 1))
        sg = np.sqrt(1 / np.maximum(s["nY"][av], 1) + 1 / np.maximum(s["nC"][av], 1))
        m1 = ajuster(t[av], x, sg, stim[av].astype(float), "M1", dt)
        y = 100 * (x - (m1["pred"] - m1["reponse"]))
        for k, (a, b) in enumerate(fen, start=1):          # 10 s après − 10 s avant sur ln R redressé (étape 7)
            m0, ma = (t[av] >= a - C.PRE_S) & (t[av] < a), (t[av] >= b) & (t[av] < b + C.POST_S)
            mb = (t >= a - C.PRE_S) & (t < a)
            niveau = inverser(s["nY"][mb].sum() / s["nC"][mb].sum() / s["eta"], *pR)
            dln = y[ma].mean() - y[m0].mean()
            redr.append(dict(roi=r, stim=k, d_lnR_redresse_pct=dln, d_c_redresse=dln / 100 / dlnR_dc(niveau, pR)))
        segs = []
        for k, (a, b) in enumerate(blocs):
            t0, T = t[a], t[b - 1] - t[a] + dt
            fin_f = t[blocs[k + 1][0]] - C.PRE_S if k + 1 < len(blocs) else FIN_S
            m = (t[av] >= t0 - 20) & (t[av] < min(t0 + 110, fin_f))
            yy = y[m] - y[(t[av] >= t0 - C.PRE_S) & (t[av] < t0)].mean()
            segs.append((t[av][m] - t0, yy, 100 * sg[m], T))
        f = ajuster_conjoint(segs)
        ic = lambda tau, e: (tau * np.exp(-1.96 * e / tau), tau * np.exp(1.96 * e / tau)) if np.isfinite(e) else (np.nan, np.nan)
        ia, ir = ic(f["tau_act_s"], f["se_tau_act"]), ic(f["tau_ret_s"], f["se_tau_ret"])
        mes = lambda lo_hi: bool(f["A"] < 0 and 0.5 < lo_hi[0] and lo_hi[1] < 290)
        lignes_k.append(dict(roi=r, amplitude_pct=f["amplitude_pct"], latence_s=f["latence_s"], se_latence=f["se_latence"],
                             tau_act_s=f["tau_act_s"], tau_act_bas=ia[0], tau_act_haut=ia[1], act_mesurable=mes(ia),
                             tau_ret_s=f["tau_ret_s"], tau_ret_bas=ir[0], tau_ret_haut=ir[1], ret_mesurable=mes(ir),
                             chi2_ddl=f["chi2_ddl"]))
        traces[r] = (segs, f)
    K = pd.DataFrame(lignes_k)
    K.to_csv(D["results"] / "rois_cinetique.csv", index=False)
    print(K.round(2).to_string(index=False))
    AA = AA.merge(pd.DataFrame(redr), on=["roi", "stim"])
    AA.loc[~AA.lisible_mM, "d_c_redresse"] = np.nan
    AA.to_csv(D["results"] / "rois_avant_apres.csv", index=False)
    print(AA.groupby("stim")[["d_R", "d_tau"]].agg(["median", "mean"]).round(3).to_string())
    print(AA[AA.lisible_mM].groupby("stim")[["d_c", "d_c_redresse"]].agg(["median", "mean", "size"]).round(3)
          .to_string())

    figure_series(D, root, masks, CL, series, t, blocs, fen)
    figure_avant_apres(D, AA)
    figure_reunies(D, AA)
    figure_cinetique(D, K, traces)
    figure_ancrage_unique(D, masks, CL, series, AA, t, fen, pR, pT, dt)
    print(f"  ({time.time() - t0s:.0f} s)")


def figure_series(D, root, masks, CL, series, t, blocs, fen):
    z = np.load(D["mesures"] / "images.npz")
    img = np.sqrt(np.maximum(z["recale"], 0))
    n = len(CL)
    fig = plt.figure(figsize=(17, 8.6 + 0.62 * n))
    gs = fig.add_gridspec(2, 1, height_ratios=[8.6, 0.62 * n], hspace=0.16)
    ax = fig.add_subplot(gs[0])
    ax.imshow(img, cmap="gray", vmin=0, vmax=np.percentile(img, 99.7))
    cm = plt.get_cmap("plasma")
    ok = CL.lisible_mM.values
    lo, hi = np.nanpercentile(CL.c_ancre_mM[ok], [5, 95])
    for _, r in CL.iterrows():
        m = masks == r.roi
        col = cm((r.c_ancre_mM - lo) / max(hi - lo, 1e-9)) if r.lisible_mM else (0.6, 0.6, 0.6, 1)
        ax.contourf(m, levels=[0.5, 1.5], colors=[col], alpha=0.45)
        ax.contour(m, levels=[0.5], colors=[col], linewidths=1.2)
        lab = f"{r.c_ancre_mM:.0f} ± {1.96 * r.se_c_mM:.0f}" if r.lisible_mM else "n. l."
        ax.text(r.x, r.y, f"{int(r.roi)}\n{lab}", color="white", fontsize=9, ha="center", va="center", fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("ROIs (Cellpose sur l'image reconstruite) — numéro et [Cl⁻] de l'ancre (mM, IC 95 % ; τ CFP des 60 s "
                 "avant la 1re stim, isotherme du 2026-11-09 ; n. l. : τ trop près du plateau de l'isotherme)",
                 loc="left", fontsize=10.5, color=G.INK)
    ax = fig.add_subplot(gs[1])
    G.style(ax)
    G.bandes(ax, [(a, b) for a, b in fen])
    off = 0.0
    for i, r in enumerate(CL.roi):
        R = series[r]["R"]
        tt, rr = moyennes(t, R, pas=5.0)
        r0 = np.nanmean(R[t < t[blocs[0][0]]])
        ax.plot(tt, 100 * (rr / r0 - 1) + off, color=G.PALETTE[i % 6], lw=1.1)
        ax.text(-10, off, f"ROI {int(r)}", ha="right", va="center", fontsize=8.5, color=G.PALETTE[i % 6])
        off -= 4.0
    ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
    ax.set_xlim(-70, t.max() + 5)
    ax.set_yticks([])
    ax.plot([5, 5], [2.0, 4.0], color=G.INK, lw=2)
    ax.text(9, 3.0, "2 %", fontsize=8.5, va="center", color=G.INK)
    ax.set_xlabel("temps (s)", color=G.INK)
    G.titre(ax, "Rapport YFP/CFP de chaque ROI, ΔR/R₀ (moyennes de 5 s, décalées de 4 %)",
            "R₀ : les 60 s avant la 1re stim ; bandes jaunes : 1064 nm ; tirets : changement optique du champ")
    G.sauver(fig, D["results"] / "rois_series.png", dpi=120)


def figure_avant_apres(D, AA):
    rangs = [("c", "Δ[Cl⁻] (mM)", "c_avant", "c_apres", "mM"), ("R", "ΔR / R$_0$ (%)", "R_avant", "R_apres", "%"),
             ("tau", "Δτ CFP (ps)", "tau_avant_ps", "tau_apres_ps", "ps")]
    stims = sorted(AA.stim.unique())
    fig, axes = plt.subplots(3, len(stims), figsize=(2.7 * len(stims) + 1.2, 12), sharey="row", squeeze=False)
    for r, (cle, lab, a, b, unite) in enumerate(rangs):
        for c, k in enumerate(stims):
            ax = axes[r][c]
            G.style(ax)
            g = AA[(AA.stim == k) & (AA.lisible_mM | (cle != "c"))]
            w = g[[a, b]].to_numpy()
            for row in w:
                ax.plot([0, 1], row, color=G.TRACE, lw=0.8, alpha=0.8)
            for j in (0, 1):
                ax.scatter(np.full(len(w), j), w[:, j], s=14, color=G.PALETTE[c % 6], alpha=0.85,
                           edgecolors=G.SURFACE, linewidths=0.5, zorder=3)
            med = np.median(w, 0)
            ax.plot([0, 1], med, color=G.INK, lw=2.2, marker="o", ms=6, zorder=4)
            e = w[:, 1] - w[:, 0]
            p = stats.wilcoxon(e).pvalue if len(e) >= 6 else np.nan
            haut = int((e > 0).sum())
            ax.set_title(f"stim {k}\nécart méd. {np.median(e):+.2f} {unite} ; {haut}/{len(e)} ↑\np = {p:.1g}",
                         color=G.INK, fontsize=9.5)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"{C.PRE_S:.0f} s\navant", f"{C.POST_S:.0f} s\naprès"], color=G.INK, fontsize=9)
            ax.set_xlim(-0.35, 1.35)
            ax.axhline(0, color=G.INK_2, lw=0.8, ls="--")
        axes[r][0].set_ylabel(lab, color=G.INK, fontsize=11)
    fig.suptitle("Avant / après chaque stimulation, ROIs de l'image reconstruite — un trait par ROI, médiane en noir ; "
                 "p : Wilcoxon sur l'écart\nΔ[Cl⁻] : rapport de chaque ROI ancré sur son τ CFP avant la 1re stim, "
                 "par rapport à l'ancre (ROIs lisibles en mM) ; R₀, τ₀ : la fenêtre de l'ancre ; valeurs brutes, "
                 "dérive comprise", color=G.INK, fontsize=11,
                 x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    G.sauver(fig, D["results"] / "avant_apres_rois.png", dpi=160)


def figure_reunies(D, AA):
    """Un seul graphique : Δ[Cl⁻] avant / après, toutes les paires ROI × stim lisibles en mM. Comme pour le champ
    (avant_apres_reunies.png, étape 7), chaque stim est rapportée à une moyenne de ses 10 s avant : ici celle
    des ROIs. Le 0 est donc la moyenne des ROIs avant la stim ; chaque ROI garde son écart à cette moyenne
    avant comme après (pas ramené à 0 un par un). Les écarts après − avant, la moyenne et le test ne changent
    pas."""
    fig, ax = plt.subplots(figsize=(7.4, 6.8))
    G.style(ax)
    d = []
    exclus = sorted(AA[~AA.lisible_mM].roi.unique())
    for c, k in enumerate(sorted(AA.stim.unique())):
        g = AA[(AA.stim == k) & AA.lisible_mM]
        ref = g.c_avant.mean()
        a, b = (g.c_avant - ref).to_numpy(), (g.c_apres - ref).to_numpy()
        e = b - a
        d.append(e)
        dx = (c - 1.5) * 0.05
        for va, vb in zip(a, b):
            ax.plot([dx, 1 + dx], [va, vb], color=G.TRACE, lw=0.6, alpha=0.8, zorder=2)
        for j, v in ((0, a), (1, b)):
            ax.scatter(np.full(len(v), j + dx), v, s=16, color=G.PALETTE[c % 6], alpha=0.85, edgecolors=G.SURFACE,
                       linewidths=0.5, zorder=3, label=f"stim {k} : médiane {np.median(e):+.2f} mM" if j else None)
    e = np.concatenate(d)
    m, se = e.mean(), e.std(ddof=1) / np.sqrt(len(e))
    p = stats.wilcoxon(e).pvalue
    ax.errorbar([0, 1], [0, m], yerr=[0, 1.96 * se], color=G.INK, lw=2.6, marker="o", ms=8, capsize=5, zorder=5,
                label=f"moyenne : {m:+.2f} mM [{m - 1.96 * se:+.2f} ; {m + 1.96 * se:+.2f}]")
    ax.axhline(0, color=G.INK_2, lw=0.8, ls="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{C.PRE_S:.0f} s\navant", f"{C.POST_S:.0f} s\naprès"], color=G.INK, fontsize=10.5)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylabel("Δ[Cl⁻] (mM), écart à la moyenne des ROIs avant la stim", color=G.INK, fontsize=10.5)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.set_title(f"Avant / après, {len(e)} paires ROI × stim réunies", color=G.INK, fontsize=12, loc="left", pad=50)
    ax.text(0, 1.015, f"{int((e > 0).sum())}/{len(e)} ↑ ; p = {p:.1g} (Wilcoxon sur les écarts) ; un trait par ROI et "
            f"par stim\n0 : moyenne des ROIs dans les 10 s avant chaque stim ; ROI(s) "
            f"{', '.join(str(int(x)) for x in exclus) or 'aucun'} exclu(s) (non lisible(s) en mM)\n"
            "rapport de chaque ROI ancré sur son τ CFP avant la 1re stim ; dérive comprise", transform=ax.transAxes,
            color=G.INK_2, fontsize=9.2, va="bottom")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "avant_apres_rois_reunies.png", dpi=170)


def figure_cinetique(D, K, traces):
    n = len(K)
    cols = 4
    rows = int(np.ceil(n / cols)) + 1
    fig = plt.figure(figsize=(18, 3.1 * rows))
    gs = fig.add_gridspec(rows, cols, hspace=0.55, wspace=0.25)
    for i, r in enumerate(K.roi):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        G.style(ax)
        segs, f = traces[r]
        g = np.arange(-20, 110, 2.0)
        M = []
        for tp, yy, _, T in segs:
            M.append([yy[(tp >= u) & (tp < u + 2)].mean() if ((tp >= u) & (tp < u + 2)).any() else np.nan for u in g])
        ax.plot(g + 1, np.nanmean(M, 0), color=G.INK, lw=1.0)
        gg = np.linspace(-20, 110, 600)
        T = float(np.mean([s[3] for s in segs]))
        ax.plot(gg, modele([0, f["A"], f["latence_s"], f["tau_act_s"], f["tau_ret_s"]], gg, T), color=G.PALETTE[1], lw=1.5)
        ax.axvspan(0, T, color=G.STIM, alpha=0.45, lw=0, zorder=0)
        ax.axhline(0, color=G.TRACE, lw=0.8)
        k = K[K.roi == r].iloc[0]
        ta = f"{k.tau_act_s:.1f} s" if k.act_mesurable else "n.m."
        tr = f"{k.tau_ret_s:.0f} s" if k.ret_mesurable else "n.m."
        ax.set_title(f"ROI {int(r)} : {k.amplitude_pct:+.2f} % ; δ {k.latence_s:.1f} s ; τ_act {ta} ; τ_ret {tr}",
                     fontsize=9.2, color=G.INK, loc="left")
        if i % cols == 0:
            ax.set_ylabel("Δ ln R (%)", color=G.INK)
    ax = fig.add_subplot(gs[rows - 1, :])
    G.style(ax)
    x = np.arange(n)
    for j, (col, lo, hi, mes, c, nom) in enumerate((("tau_ret_s", "tau_ret_bas", "tau_ret_haut", "ret_mesurable", G.PALETTE[1], "τ de retour"),
                                                    ("tau_act_s", "tau_act_bas", "tau_act_haut", "act_mesurable", G.PALETTE[0], "τ d'activation"))):
        ok = K[mes].values.astype(bool)
        xx = x + (j - 0.5) * 0.25
        ax.errorbar(xx[ok], K[col][ok], yerr=[K[col][ok] - K[lo][ok], K[hi][ok] - K[col][ok]], fmt="o", color=c, ms=5,
                    capsize=2, lw=1, label=nom)
        if (~ok).any():
            ax.plot(xx[~ok], np.minimum(K[col][~ok], 78), "o", mfc="none", color=c, ms=6, label=f"{nom} non mesurable")
    ax.errorbar(x + 0.25, K.latence_s, yerr=1.96 * K.se_latence, fmt="s", color=G.PALETTE[2], ms=4, capsize=2, lw=1,
                label="latence δ")
    ax.set_xticks(x)
    ax.set_xticklabels([f"ROI {int(r)}" for r in K.roi], fontsize=9)
    ax.set_ylim(0, 80)
    ax.set_ylabel("s", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5, ncol=5)
    G.titre(ax, "Constantes de temps par ROI (toutes les stims ensemble)", "IC à 95 % ; cercles vides : non mesurable",
            fontsize=10.5)
    fig.suptitle("Réponse moyenne aux stims par ROI (ln R redressé, moyennes de 2 s) et modèle marche-retour avec latence",
                 color=G.INK, fontsize=12, x=0.01, ha="left")
    G.sauver(fig, D["results"] / "cinetique_rois.png", dpi=130)


def signe(v, d=1):
    """Écart signé, sans « -0.0 »."""
    return f"{0.0:.{d}f}" if round(v, d) == 0 else f"{v:+.{d}f}"


def figure_ancrage_unique(D, masks, CL, series, AA, t, fen, pR, pT, dt):
    """Comme hybride_ancrage_unique.png (étape 7), par ROI : [Cl⁻] lue au rapport de chaque soma, son échelle
    ancrée sur son τ CFP pendant les 60 s avant la 1re stim, puis fixée. Fenêtres de 2 s (trait fin) et
    glissantes de 10 s (trait épais), chacune inversée sur ses photons sommés. En haut : moyenne des ROIs
    lisibles en mM, IC du niveau tiré des ancres. Au-dessus des stims : 10 s après − 10 s avant, dérive lente
    retirée. Les panneaux des ROIs couvrent tous la même étendue en mM (seul le centre change)."""
    lis = CL[CL.lisible_mM]
    t_anc = fen[0][0]
    k10 = 2 * int(round(5.0 / dt)) + 1                     # passes d'une fenêtre glissante de ~10 s
    tr = {}
    for r in lis.itertuples():
        s = series[r.roi]
        t2, y2 = moyennes(t, s["nY"], 2.0)
        _, x2 = moyennes(t, s["nC"], 2.0)
        sy, sx = np.convolve(s["nY"], np.ones(k10), "same"), np.convolve(s["nC"], np.ones(k10), "same")
        tr[r.roi] = (inverser(y2 / x2 / s["eta"], *pR), inverser(sy / sx / s["eta"], *pR))
    bruit = lambda c2: np.nanstd(np.diff(c2[t2 < FIN_S]))  # écart des fenêtres de 2 s successives = bruit par seconde
    rep = AA[AA.lisible_mM].pivot(index="roi", columns="stim", values="d_c_redresse")
    cfp = G.COULEUR_CANAL["cfp"]
    nl = int(np.ceil(len(CL) / 3))
    fig = plt.figure(figsize=(18, 5.6 + 2.6 * nl))
    gs = fig.add_gridspec(1 + nl, 3, height_ratios=[1.9] + [1] * nl, hspace=0.5, wspace=0.13)

    # --- moyenne des ROIs lisibles
    ax = fig.add_subplot(gs[0, :2])
    G.style(ax)
    G.bandes(ax, fen)
    ax.axvspan(0, t_anc, ymin=0, ymax=0.03, color=cfp, alpha=0.6, lw=0)
    c2m = np.mean([tr[r][0] for r in lis.roi], axis=0)
    c10m = np.mean([tr[r][1] for r in lis.roi], axis=0)
    c_anc = lis.c_ancre_mM.mean()
    ic = 1.96 * np.sqrt((lis.se_c_mM ** 2).sum()) / len(lis)
    ax.plot(t2, c2m, color=G.TRACE, lw=0.7, label="fenêtres de 2 s")
    ax.fill_between(t, c10m - ic, c10m + ic, color=G.INK, alpha=0.10, lw=0,
                    label=f"IC 95 % du niveau (ancres : ± {ic:.1f} mM)")
    ax.plot(t, c10m, color=G.INK, lw=1.6, label="glissantes de 10 s")
    ax.errorbar([t_anc / 2], [c_anc], yerr=[ic], fmt="o", color=cfp, ms=7, capsize=4, lw=1.3,
                label=f"ancre : moyenne des τ CFP des ROIs, 0–{t_anc:.0f} s")
    for k, (t0, t1) in enumerate(fen, start=1):
        v = rep[k]
        ax.text((t0 + t1) / 2, 0.97, f"{v.mean():+.2f} mM\n{int((v > 0).sum())}/{len(v)} ↑",
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8.5, color=G.INK_2)
    ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
    ax.text(FIN_S + 3, 0.2, "changement optique du champ\n(CFP −5 %, τ stable)", transform=ax.get_xaxis_transform(),
            fontsize=8, color=G.INK_2, va="bottom")
    lo = min(np.nanpercentile(c2m, 0.3), np.nanmin(c10m - ic))
    hi = max(np.nanpercentile(c2m, 99.7), np.nanmax(c10m + ic))
    ax.set_ylim(lo - 0.22 * (hi - lo), hi + 0.2 * (hi - lo))
    ax.set_xlim(-8, t.max() + 5)
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.set_ylabel("[Cl⁻] (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.3, ncol=4, loc="lower left", bbox_to_anchor=(0, 0.02))
    E = max(np.nanpercentile(tr[r][0], 99.5) - np.nanpercentile(tr[r][0], 0.5) for r in lis.roi)
    ax.set_title("[Cl⁻] des somas lue au rapport YFP/CFP, échelle de chaque ROI ancrée sur son τ CFP avant la 1re "
                 "stim puis fixée", color=G.INK, fontsize=12, loc="left", pad=38)
    ax.text(0, 1.02, f"moyenne des {len(lis)} ROIs lisibles en mM : ancre {c_anc:.1f} ± {ic:.1f} mM, bruit "
            f"{bruit(c2m):.2f} mM par seconde ; au-dessus des stims : 10 s après − 10 s avant, dérive lente retirée\n"
            f"chaque fenêtre inversée sur ses photons sommés ; isothermes du 2026-11-09 ; panneaux des ROIs : même "
            f"étendue ({1.23 * E:.0f} mM) pour tous", transform=ax.transAxes, color=G.INK_2, fontsize=9.2, va="bottom")

    # --- carte des ROIs
    axm = fig.add_subplot(gs[0, 2])
    img = np.sqrt(np.maximum(np.load(D["mesures"] / "images.npz")["recale"], 0))
    axm.imshow(img, cmap="gray", vmin=0, vmax=np.percentile(img, 99.7))
    for r in CL.itertuples():
        axm.contour(masks == r.roi, levels=[0.5], colors=[G.PALETTE[1] if r.lisible_mM else "#9a9a9a"], linewidths=1.0)
        axm.text(r.x, r.y, str(int(r.roi)), color="white", fontsize=9, ha="center", va="center", fontweight="bold",
                 path_effects=[pe.withStroke(linewidth=2.2, foreground=G.INK)])
    axm.set_xticks([])
    axm.set_yticks([])
    axm.set_title("ROIs (Cellpose sur l'image reconstruite)\ngris : non lisible en mM", loc="left", fontsize=10,
                  color=G.INK)

    # --- un panneau par ROI
    for i, r in enumerate(CL.itertuples()):
        ax = fig.add_subplot(gs[1 + i // 3, i % 3])
        if not r.lisible_mM:
            ax.axis("off")
            ax.text(0.04, 0.5, f"ROI {int(r.roi)} : [Cl⁻] non lisible\n\nτ CFP de l'ancre {r.tau_ancre_ns:.2f} ns, trop "
                    f"près du plateau\nde l'isotherme ({pT[1]:.2f} ns) : {r.c_ancre_mM:.0f} ± {1.96 * r.se_c_mM:.0f} mM"
                    "\n\nson rapport : rois_series.png", transform=ax.transAxes, va="center", fontsize=9.5,
                    color=G.INK_2)
            continue
        c2, c10 = tr[r.roi]
        G.style(ax)
        G.bandes(ax, fen)
        ax.axvspan(0, t_anc, ymin=0, ymax=0.04, color=cfp, alpha=0.6, lw=0)
        ax.plot(t2, c2, color=G.TRACE, lw=0.6)
        ax.plot(t, c10, color=G.INK, lw=1.2)
        ax.errorbar([t_anc / 2], [r.c_ancre_mM], yerr=[1.96 * r.se_c_mM], fmt="o", color=cfp, ms=5, capsize=3, lw=1.1)
        ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
        for k, (t0, t1) in enumerate(fen, start=1):
            ax.text((t0 + t1) / 2, 0.97, signe(rep.loc[r.roi, k]), transform=ax.get_xaxis_transform(), ha="center",
                    va="top", fontsize=7.8, color=G.INK_2)
        mid = (np.nanpercentile(c2, 0.5) + np.nanpercentile(c2, 99.5)) / 2
        ax.set_ylim(mid - 0.55 * E, mid + 0.68 * E)
        ax.set_xlim(-8, t.max() + 5)
        ax.set_title(f"ROI {int(r.roi)} : ancre {r.c_ancre_mM:.1f} ± {1.96 * r.se_c_mM:.1f} mM ; bruit "
                     f"{bruit(c2):.1f} mM par seconde", loc="left", fontsize=9.5, color=G.INK)
        if i % 3 == 0:
            ax.set_ylabel("[Cl⁻] (mM)", color=G.INK)
        if i // 3 == nl - 1:
            ax.set_xlabel("temps (s)", color=G.INK)
        else:
            ax.tick_params(labelbottom=False)
    fig.subplots_adjust(left=0.045, right=0.99, top=0.93, bottom=0.045)
    G.sauver(fig, D["results"] / "hybride_ancrage_unique_rois.png", dpi=110)


if __name__ == "__main__":
    main()
