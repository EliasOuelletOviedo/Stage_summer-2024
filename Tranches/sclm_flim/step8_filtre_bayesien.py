"""Étape 8 — le filtre bayésien hybride (ratio + FLIM) de Simulateurs/bayesian_filter, sur les somas.

Usage :  python3 step8_filtre_bayesien.py <dossier>          (~1 à 5 min selon la session)

Le protocole est celui de `Simulateurs/bayesian_filter` (README de ce dossier), sans changer le filtre :
modèle « 1-2-3 » (η dérive, θ = (η0, β) partagé par les somas du champ, vitesse de c inconnue), 16 portes
temporelles, grille de c de 0,5 à 300 mM, filtre causal (temps réel) puis lissage hors ligne. Les constantes du
senseur viennent de la titration SCLM à 800 nm (`calibration_SCLM.json` : isothermes, λ, G, g2, κ).

Ce qu'il faut ajouter pour passer de la solution à la tranche, étape par étape :

1. **Instrument.** Les déclins de la QC-104 ont 256 canaux sur une période (48,8 ps), ceux de la titration
   (SPC-150N) 64 ps : l'IRF paramétrique est réajustée sur `irf.csv` de la session (même forme que pour la
   titration : gaussienne ⊗ traînes exponentielles). Fenêtre : canaux FENETRE (pleins, hors bords du TAC).
2. **Transfert de la calibration**, sur la somme de toute la session (hors stims), ajustement de Poisson
   conjoint de tous les somas :
   - partagés par le champ : le décalage d'horloge de chaque canal (le README du filtre demande de le mesurer à
     chaque session), η, et — variante « tranche », celle du filtre — les quatre durées de vie τ1, τ2, τD, τA ;
   - propres à chaque soma : c, A, le fond de chaque canal, et la YFP directe k0 (eNpHR est fusionnée à la YFP :
     de la YFP sans partenaire FRET, excitée directement ; c'est le terme k0 du modèle, nul en solution).
   La variante « solution » garde les durées de vie de la titration : elle ne sert que de contrôle (qc/).
3. **Filtre** sur les tranches de l'étape 5 (~3 s), une Grille par soma (sa YFP directe), son fond par tranche
   (celui de la session, réparti également). Les tranches de stim (fuite du 1064 nm dans ch2) et celles où le
   soma sort du champ ne sont pas lues : le filtre y prédit seulement. Le lissage de chaque soma utilise la loi
   de θ du champ (exact : sachant θ, les somas sont indépendants).
4. **Réponses**, comme à l'étape 5 : moyenne de c dans les 10 s après la fin de la stim − 10 s avant son début,
   et les mêmes témoins (fausses stims au milieu des intervalles).

La dérive de η est bornée par DERIVE_MAX_PCT_100S, en temps physique : la grille du protocole (±5 % par 1 000
trames) ne laisserait à η que ±0,9 % sur une session de 172 tranches de 3 s.

Sorties : qc/filtre_transfert.png, qc/filtre_transfert.json, qc/filtre_transfert_somas.csv ;
results/filtre_par_trame.csv, filtre_reponses.csv, filtre_resume.csv, filtre_bayesien.png, filtre_reponses.png.
"""
import copy
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

import config
from config import MIN_VALID_FRAC, PRE_WINDOW_S, RESPONSE_WINDOW_S
from figures import INK, INK_2, PALETTE, SURFACE, TRACE, bandes_stim, plt, sauver, style
from step5_analyse import COL_TEMOIN, regrouper, temoins, tranches

FILTRE = Path(__file__).resolve().parents[2] / "Simulateurs" / "bayesian_filter"
sys.path.insert(0, str(FILTRE))
from src.filtre import FiltreTempsReel, Grille, grille_c, lisseur_derive, portes_log  # noqa: E402
from src.modele import Senseur, charger_calibration, charger_params  # noqa: E402

PARAMS = "params_SCLM.toml"             # senseur SCLM régulier, 800 nm (comme les sessions en tranche)
FENETRE = (17, 254)                     # canaux pleins : 16 et 255 sont coupés par le TAC (QC-104, 0,82 ns à T)
DERIVE_MAX_PCT_100S = 5.0               # |β| maximal : 5 % de η par 100 s (±26 % sur 520 s)
ETA0 = (0.5, 2.0)                       # grille initiale de η0 (pas de 0,95 % comme le protocole) ; elle suit ensuite
TAUS = ("t1", "t2", "tD", "tA")
DEPART_TAUS = (0.1, 0.7, 3.05, 3.2)     # 3 exponentielles libres sur le bleu des somas du 2026-10-08
N_PART, N_ROI = 7, 5                    # paramètres partagés (4 τ, 2 décalages, ln η) ; par soma (ln c, ln A, k0, 2 fonds)


# ------------------------------------------------------------------ instrument
def instrument(root, cal_sol):
    """Calibration du filtre pour cette session : constantes du senseur de la titration, IRF et période de la
    QC-104 (IRF paramétrique ajustée sur irf.csv par la fonction de la titration)."""
    spec = importlib.util.spec_from_file_location(
        "sclm_hybride", FILTRE / charger_params(FILTRE / PARAMS)["donnees"]["dossier"] / "sclm_hybride.py")
    sh = importlib.util.module_from_spec(spec)
    sys.dont_write_bytecode = True
    spec.loader.exec_module(sh)
    I, dt = config.irf(root)
    n = config.CANAUX_TEMPS
    T = n * dt
    win = np.zeros(n, bool)
    win[FENETRE[0]:FENETRE[1] + 1] = True
    d = dict(T=T, NB=n, DT=dt, t=(np.arange(n) + 0.5) * dt, win={1: win, 2: win}, irf_raw=I)
    irfp = {k: sh.fit_irf(d, k, ntail=2 if k == 1 else 3) for k in (1, 2)}
    cal = copy.deepcopy(cal_sol)
    cal.update(periode_ns=T, n_cases=n, fenetre=list(FENETRE), decalage_moyen_ps=dict(bleu=0.0, jaune=0.0),
               irf={str(k): dict(mu=v["mu"], sig=v["sig"], poids=list(v["weights"]), taus=list(v["taus"]),
                                 dev_par_ddl=v["dev_per_dof"]) for k, v in irfp.items()})
    cal["constantes"]["k0"] = 0.0
    return cal, sh.dev_resid


# ------------------------------------------------------------------ transfert de la calibration
def senseur(cal, v, taus_libres, k0=0.0):
    c2 = copy.deepcopy(cal)
    if taus_libres:
        for i, nm in enumerate(TAUS):
            c2["constantes"][nm] = float(np.exp(v[i]))
    c2["constantes"]["k0"] = float(k0)
    return Senseur(c2, decalage_ps=(v[4], v[5]))


def transfert(cal, H, dres, taus_libres):
    """Ajustement de Poisson conjoint des somas sur la session (H : (M, 2, n_fenêtre)) ; voir l'en-tête."""
    M, _, n = H.shape
    k = cal["constantes"]

    def modele(v):
        S = senseur(cal, v, taus_libres)
        eta = np.exp(v[6])
        P = v[N_PART:].reshape(M, N_ROI)
        A = np.exp(P[:, 1])[:, None]
        mB, mY = S.mu(np.exp(P[:, 0]), A=A, eta=eta)
        mY = mY + eta * A * P[:, 2, None] * S.bYA[None, :]          # = le terme k0 du modèle, par soma
        return mB + P[:, 3, None], mY + P[:, 4, None]

    def residus(v):
        mB, mY = modele(v)
        return np.concatenate([np.r_[dres(H[m, 0], mB[m]), dres(H[m, 1], mY[m])] for m in range(M)])

    taus0 = DEPART_TAUS if taus_libres else [k[t] for t in TAUS]
    x0 = [*np.log(taus0), -185.0, -172.0, 0.0]
    lo = [*np.log([0.02, 0.2, 1.5, 1.5]), -500, -500, np.log(0.2)]
    hi = [*np.log([0.4, 1.5, 6.0, 6.0]), 200, 200, np.log(5.0)]
    if not taus_libres:                                              # τ fixés : bornes resserrées autour
        lo[:4] = [x - 1e-6 for x in x0[:4]]
        hi[:4] = [x + 1e-6 for x in x0[:4]]
    mB0, _ = Senseur(cal).mu(10.0)
    for m in range(M):
        x0 += [np.log(10.0), np.log(H[m, 0].sum() / mB0.sum()), 0.05, 0.1 * H[m, 0].min() + 1e-3,
               0.05 * H[m, 1].min() + 1e-3]
        lo += [np.log(0.1), -np.inf, 0, 0, 0]
        hi += [np.log(300.0), np.inf, 5, np.inf, np.inf]
    lo, hi = np.array(lo), np.array(hi)
    J = lil_matrix((M * 2 * n, len(x0)), dtype=int)
    J[:, :N_PART] = 1
    for m in range(M):
        J[m * 2 * n:(m + 1) * 2 * n, N_PART + m * N_ROI:N_PART + (m + 1) * N_ROI] = 1
    r = least_squares(residus, np.clip(x0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi), jac_sparsity=J,
                      x_scale="jac", max_nfev=3000)
    v = r.x
    mB, mY = modele(v)
    P = v[N_PART:].reshape(M, N_ROI)
    dev_roi = [float((np.r_[dres(H[m, 0], mB[m]), dres(H[m, 1], mY[m])] ** 2).sum() / (2 * n - N_ROI))
               for m in range(M)]
    return dict(v=v, taus=np.exp(v[:4]), decalage_ps=v[4:6], eta=float(np.exp(v[6])),
                c=np.exp(P[:, 0]), A=np.exp(P[:, 1]), k0=P[:, 2], fond=P[:, 3:5],
                dev_ddl=float((r.fun ** 2).sum() / (M * 2 * n - len(v))), dev_ddl_roi=np.array(dev_roi),
                mB=mB, mY=mY, taus_libres=taus_libres, message=r.message)


# ------------------------------------------------------------------ filtre
def parametres_filtre(dt_trame):
    f = dict(charger_params(FILTRE / PARAMS)["filtre"])
    f["derive_max"] = DERIVE_MAX_PCT_100S / 100 * 1000 * dt_trame / 100       # ln η par 1 000 trames
    f["eta0_min"], f["eta0_max"] = ETA0
    f["n_eta0"] = int(round(np.log(ETA0[1] / ETA0[0]) / np.log(1.5 / 0.7) * (f["n_eta0"] - 1))) + 1
    return f


def filtrer(cal, tr, H, lues, n_lues, f):
    """H : (M, T, 2, n_fenêtre) ; lues : (M, T) tranches lues. Filtre causal puis lissage de chaque soma."""
    M, T = H.shape[:2]
    c = grille_c(f)
    G = portes_log(H.shape[-1], f["portes"])
    grilles = [Grille(senseur(cal, tr["v"], tr["taus_libres"], tr["k0"][m]), c, [1.0], G) for m in range(M)]
    fond = tr["fond"] / n_lues[:, None]                         # photons par case et par tranche, par soma
    fl = FiltreTempsReel(grilles, M, f)
    zero = np.zeros((M, H.shape[-1]))
    rt = dict(c=np.full((T, M), np.nan), bas=np.full((T, M), np.nan), haut=np.full((T, M), np.nan),
              eta=np.empty(T), eta_sd=np.empty(T), poids=np.empty((T, len(f["s_c_candidats"]))))
    t0 = time.time()
    for k in range(T):
        B = np.where(lues[:, k, None], H[:, k, 0], 0)
        Y = np.where(lues[:, k, None], H[:, k, 1], 0)
        if lues[:, k].all():
            r = fl.mettre_a_jour(B, Y, fond)
        elif not lues[:, k].any():
            r = fl.mettre_a_jour(zero, zero)                    # aucune donnée : prédiction seule
        else:                                                   # somas sortis du champ : sans donnée ni fond
            r = fl.mettre_a_jour(B, Y, np.where(lues[:, k, None], fond, 0))
        for cle in ("c", "bas", "haut"):
            rt[cle][k] = r[cle]
        rt["eta"][k], rt["eta_sd"][k], rt["poids"][k] = r["eta"], r["eta_sd"], r["poids"]
    print(f"  filtre : {T} tranches × {M} somas en {time.time() - t0:.0f} s "
          f"({fl.n_regrilles} recentrages de θ) ; η final {rt['eta'][-1]:.3f} ± {100 * rt['eta_sd'][-1]:.2f} %")
    t0 = time.time()
    lis = dict(c=np.empty((T, M)), bas=np.empty((T, M)), haut=np.empty((T, M)))
    for m in range(M):
        g = grilles[m]
        Bg, Yg = g.regrouper(H[m, :, 0]), g.regrouper(H[m, :, 1])
        formes = np.zeros((T, len(c)))
        NB, NY = np.zeros(T), np.zeros(T)
        for k in np.flatnonzero(lues[m]):
            formes[k] = g.formes(Bg[k], Yg[k], fond[m])
            NB[k], NY[k] = g._sous_fond(Bg[k], Yg[k], fond[m])
        pc, le, sd, _ = lisseur_derive(fl, formes, NB, NY, grille=g)
        cdf = np.cumsum(pc, 1)
        lis["c"][:, m] = pc @ c
        lis["bas"][:, m] = c[np.argmax(cdf >= 0.025, 1)]
        lis["haut"][:, m] = c[np.argmax(cdf >= 0.975, 1)]
    lis["eta"], lis["eta_sd"] = np.exp(le), sd
    print(f"  lissage : {time.time() - t0:.0f} s")
    return rt, lis, fl


# ------------------------------------------------------------------ réponses
def reponses(tc, ok, C, fen, tem, rois, lecture):
    """Une ligne par soma × stim (ou témoin) : c moyen 10 s avant, 10 s après, écart (comme l'étape 5)."""
    rows = []
    for typ, wins in (("stim", fen), ("témoin", tem)):
        for i, (a, b) in enumerate(wins):
            av = (tc >= a - PRE_WINDOW_S) & (tc < a)
            ap = (tc > b) & (tc <= b + RESPONSE_WINDOW_S)
            for m, roi in enumerate(rois):
                x_av = C[av & ok[m], m].mean() if (av & ok[m]).any() else np.nan
                x_ap = C[ap & ok[m], m].mean() if (ap & ok[m]).any() else np.nan
                rows.append(dict(lecture=lecture, type=typ, numero=i + 1, t0_s=round(a, 2), t1_s=round(b, 2),
                                 roi=roi, c_avant_mM=x_av, c_apres_mM=x_ap, dc_mM=x_ap - x_av))
    return pd.DataFrame(rows)


def resumer(rep):
    rows = []
    for (lec, typ, num), g in rep.groupby(["lecture", "type", "numero"]):
        x = g["dc_mM"].dropna()
        p = stats.wilcoxon(x).pvalue if len(x) >= 6 and np.any(x != 0) else np.nan
        rows.append(dict(lecture=lec, type=typ, numero=num, n=len(x), dc_mM_median=x.median(),
                         n_hausse=int((x > 0).sum()), p_wilcoxon=p))
    return pd.DataFrame(rows)


def stims_contre_temoins(rep):
    """Par soma : moyenne des stims contre moyenne des témoins (Wilcoxon apparié), pour chaque lecture."""
    out = {}
    for lec, g in rep.groupby("lecture"):
        m = g.groupby(["roi", "type"])["dc_mM"].mean().unstack()
        if "témoin" not in m or m["témoin"].isna().all():
            continue
        m = m.dropna()
        out[lec] = dict(stim=float(m["stim"].median()), temoin=float(m["témoin"].median()), n=len(m),
                        p=float(stats.wilcoxon(m["stim"], m["témoin"]).pvalue) if len(m) >= 6 else np.nan)
    return out


# ------------------------------------------------------------------ figures
def figure_transfert(H, trs, n0, dt, out):
    """Somme des somas sur la session : données, modèle « solution » et modèle « tranche », résidus."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.6), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    t = (np.arange(H.shape[-1]) + n0 + 0.5) * dt
    cols = {False: PALETTE[1], True: PALETTE[0]}
    for j, (k, nom) in enumerate(((0, "bleu (ch1, CFP)"), (1, "jaune (ch2, YFP)"))):
        y = H[:, k].sum(0)
        ax, ar = axes[0, j], axes[1, j]
        style(ax), style(ar)
        ax.semilogy(t, y, color=INK, lw=0, marker=".", ms=2.5, label="somas, toute la session")
        for tr in trs:
            mu = (tr["mY"] if k else tr["mB"]).sum(0)
            lab = ("durées de vie en tranche" if tr["taus_libres"] else "durées de vie de la solution")
            ax.semilogy(t, mu, color=cols[tr["taus_libres"]], lw=1.3,
                        label=f"{lab} (déviance/ddl {tr['dev_ddl']:.0f})")
            ar.plot(t, 100 * (y / mu - 1), color=cols[tr["taus_libres"]], lw=0.9)
        ar.axhline(0, color=INK_2, lw=0.6)
        ar.set_ylim(-25, 25)
        ax.set_title(nom, loc="left", color=INK, fontsize=11)
        ar.set_xlabel("temps dans la période (ns)", color=INK)
        if j == 0:
            ax.set_ylabel("photons par canal", color=INK)
            ar.set_ylabel("données / modèle − 1 (%)", color=INK)
            ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    fig.suptitle("Transfert de la calibration en solution (SCLM, 800 nm) aux somas de la session",
                 x=0.01, ha="left", color=INK, fontsize=12.5)
    fig.tight_layout()
    sauver(fig, out)


def figure_traces(tc, fen, rt, lis, ok, out, titre):
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True, gridspec_kw=dict(height_ratios=[2, 2, 1.3]))
    for ax in axes:
        style(ax)
        bandes_stim(ax, fen)
    for ax, D, nom in ((axes[0], lis, "lissé (hors ligne)"), (axes[1], rt, "temps réel (causal)")):
        C = np.where(ok.T, D["c"], np.nan)
        for m in range(C.shape[1]):
            ax.plot(tc, C[:, m], color=TRACE, lw=0.6, alpha=0.7)
        med = np.nanmedian(C, 1)
        q1, q3 = np.nanpercentile(C, (25, 75), 1)
        ax.fill_between(tc, q1, q3, color=PALETTE[0], alpha=0.18, lw=0)
        ax.plot(tc, med, color=PALETTE[0], lw=1.8)
        ax.set_ylabel(f"[Cl⁻] {nom} (mM)", color=INK)
        hi = np.nanpercentile(C, 99)
        ax.set_ylim(0, hi * 1.08)
    axes[2].fill_between(tc, rt["eta"] * np.exp(-2 * rt["eta_sd"]), rt["eta"] * np.exp(2 * rt["eta_sd"]),
                         color=PALETTE[2], alpha=0.2, lw=0)
    axes[2].plot(tc, rt["eta"], color=PALETTE[2], lw=1.4, label="temps réel ± 2 sd")
    axes[2].plot(tc, lis["eta"], color=INK, lw=1.0, ls="--", label="lissé")
    axes[2].set_ylabel("η (jaune / bleu,\nrelatif à la solution)", color=INK)
    axes[2].legend(frameon=False, fontsize=8.5, loc="upper right")
    axes[2].set_xlabel("temps (s, horloge des cartes SPC)", color=INK)
    axes[0].set_title(titre, loc="left", color=INK, fontsize=12, pad=20)
    axes[0].text(0, 1.02, "gris : chaque soma ; bleu : médiane et écart interquartile ; bandes jaunes : stims",
                 transform=axes[0].transAxes, color=INK_2, fontsize=9, va="bottom")
    fig.tight_layout()
    sauver(fig, out)


def figure_reponses(rep, ratio, out):
    """Écart de c (lissé) par soma : chaque stim et chaque témoin ; puis contre ΔR/R0 de l'étape 5."""
    r = rep[rep["lecture"] == "lissé"]
    cols = [(t, n) for t, n in r[["type", "numero"]].drop_duplicates().itertuples(index=False)]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw=dict(width_ratios=[1.6, 1]))
    ax = axes[0]
    style(ax)
    rng = np.random.default_rng(0)
    for i, (t, n) in enumerate(cols):
        x = r[(r["type"] == t) & (r["numero"] == n)]["dc_mM"].dropna().values
        col = PALETTE[0] if t == "stim" else COL_TEMOIN
        ax.scatter(i + rng.uniform(-0.15, 0.15, len(x)), x, s=12, color=col, alpha=0.6, lw=0)
        ax.plot([i - 0.28, i + 0.28], [np.median(x)] * 2, color=INK, lw=2)
    ax.axhline(0, color=INK_2, lw=0.7)
    ax.set_xticks(range(len(cols)), [f"{'stim' if t == 'stim' else 'témoin'} {n}" for t, n in cols],
                  rotation=30, ha="right")
    ax.set_ylabel("Δ[Cl⁻] lissé, 10 s après − 10 s avant (mM)", color=INK)
    ax.set_title("Réponse de chaque soma (trait : médiane)", loc="left", color=INK, fontsize=11)
    ax = axes[1]
    style(ax)
    if ratio is not None and len(ratio):
        ax.scatter(ratio["dR_R0_pct"], ratio["dc_mM"], s=16, color=PALETTE[0], alpha=0.7, lw=0)
        rho, p = stats.spearmanr(ratio["dR_R0_pct"], ratio["dc_mM"], nan_policy="omit")
        ax.set_title(f"Par soma, moyenne des stims : ρ = {rho:.2f} (p = {p:.1g})", loc="left", color=INK,
                     fontsize=11)
    ax.axhline(0, color=INK_2, lw=0.7)
    ax.axvline(0, color=INK_2, lw=0.7)
    ax.set_xlabel("ΔR/R0 de l'étape 5 (%)", color=INK)
    ax.set_ylabel("Δ[Cl⁻] du filtre (mM)", color=INK)
    fig.tight_layout()
    sauver(fig, out)


# ------------------------------------------------------------------ principal
def main():
    warnings.filterwarnings("ignore", "All-NaN slice", RuntimeWarning)       # tranches de stim, sans lecture
    root = config.racine()
    D = config.dirs(root)
    info = json.loads((D["projections"] / "session.json").read_text())
    sel = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
    d = np.load(D["mesures"] / "declins.npz")
    fen = [(s["t0"], s["t1"]) for s in info["stims"]]
    kb, T, b0, b1, bst = tranches(info)
    tc = (b0 + b1) / 2
    dt_trame = float(np.median(b1 - b0))
    print(f"{root.name} — {len(sel)} somas, {T} tranches de {dt_trame:.1f} s")
    f0, f1 = FENETRE[0], FENETRE[1] + 1
    H = np.stack([regrouper(d["cellpose"][l], kb, T) for l in sel])[..., f0:f1].astype(float)   # (M, T, 2, n)
    frac = np.stack([regrouper(d["frac_cellpose"][:, l], kb, T, np.min) for l in sel])          # (M, T)
    lues = (frac >= MIN_VALID_FRAC) & ~bst[None, :]
    n_lues = lues.sum(1).astype(float)
    Hs = (H * lues[..., None, None]).sum(1)                                                     # (M, 2, n)

    cal_sol = charger_calibration(charger_params(FILTRE / PARAMS))
    cal, dres = instrument(root, cal_sol)
    print(f"  IRF QC-104 réajustée : déviance/ddl {cal['irf']['1']['dev_par_ddl']:.1f} (ch1), "
          f"{cal['irf']['2']['dev_par_ddl']:.1f} (ch2)")
    trs = []
    for libres in (False, True):
        t0 = time.time()
        tr = transfert(cal, Hs, dres, libres)
        trs.append(tr)
        print(f"  transfert, durées de vie {'en tranche' if libres else 'de la solution'} : déviance/ddl "
              f"{tr['dev_ddl']:.1f} ; τ {np.round(tr['taus'], 3)} ns ; décalages "
              f"{np.round(tr['decalage_ps'], 0)} ps ; η {tr['eta']:.3f} ; c médian {np.median(tr['c']):.1f} mM "
              f"({time.time() - t0:.0f} s)")
    tr = trs[1]
    figure_transfert(Hs, trs, f0, cal["periode_ns"] / cal["n_cases"], D["qc"] / "filtre_transfert.png")
    json.dump(dict(fenetre=list(FENETRE), irf=cal["irf"], periode_ns=cal["periode_ns"],
                   variantes={("tranche" if t["taus_libres"] else "solution"): dict(
                       taus_ns=dict(zip(TAUS, map(float, t["taus"]))),
                       decalage_ps=dict(bleu=float(t["decalage_ps"][0]), jaune=float(t["decalage_ps"][1])),
                       eta=t["eta"], deviance_par_ddl=t["dev_ddl"], c_median_mM=float(np.median(t["c"])))
                       for t in trs}),
              open(D["qc"] / "filtre_transfert.json", "w"), indent=2, ensure_ascii=False)
    pd.DataFrame({"roi": sel, **{f"{'tranche' if t['taus_libres'] else 'solution'}_{k}": v
                                  for t in trs for k, v in (
                                      ("c_mM", t["c"]), ("k0", t["k0"]),
                                      ("fond_bleu_pct", 100 * t["fond"][:, 0] * Hs.shape[-1] / Hs[:, 0].sum(1)),
                                      ("fond_jaune_pct", 100 * t["fond"][:, 1] * Hs.shape[-1] / Hs[:, 1].sum(1)),
                                      ("deviance_par_ddl", t["dev_ddl_roi"]))},
                  "photons_bleus_par_tranche": Hs[:, 0].sum(1) / n_lues,
                  "photons_jaunes_par_tranche": Hs[:, 1].sum(1) / n_lues}).to_csv(
        D["qc"] / "filtre_transfert_somas.csv", index=False)

    f = parametres_filtre(dt_trame)
    rt, lis, fl = filtrer(cal, tr, H, lues, n_lues, f)

    rows = []
    for m, roi in enumerate(sel):
        rows.append(pd.DataFrame(dict(roi=roi, tranche=np.arange(1, T + 1), t_s=np.round(tc, 3),
                                      lue=lues[m], c_temps_reel_mM=rt["c"][:, m], bas_temps_reel_mM=rt["bas"][:, m],
                                      haut_temps_reel_mM=rt["haut"][:, m], c_lisse_mM=lis["c"][:, m],
                                      bas_lisse_mM=lis["bas"][:, m], haut_lisse_mM=lis["haut"][:, m],
                                      eta_temps_reel=rt["eta"], eta_sd_temps_reel=rt["eta_sd"],
                                      eta_lisse=lis["eta"])))
    pd.concat(rows).to_csv(D["results"] / "filtre_par_trame.csv", index=False)

    titre = f"Filtre bayésien hybride (ratio + FLIM) — {len(sel)} somas, η partagé par le champ"
    figure_traces(tc, fen, rt, lis, lues, D["results"] / "filtre_bayesien.png", titre)

    ok = ~bst[None, :] & (frac >= MIN_VALID_FRAC)
    tem = temoins(fen, info["blocs"][-1]["t1"])
    rep = pd.concat([reponses(tc, ok, lis["c"], fen, tem, sel, "lissé"),
                     reponses(tc, ok, rt["c"], fen, tem, sel, "temps réel")], ignore_index=True)
    rep.to_csv(D["results"] / "filtre_reponses.csv", index=False)
    rs = resumer(rep)
    rs.to_csv(D["results"] / "filtre_resume.csv", index=False)
    sct = stims_contre_temoins(rep)

    ratio = None
    aa = D["results"] / "avant_apres_stim.csv"
    if aa.exists():
        a = pd.read_csv(aa)
        a = a[(a["famille"] == "cellpose") & (a["type"] == "stim")]
        a = a.groupby(a["roi"].astype(str))["dR_R0_pct_ecart"].mean()
        c_ = rep[(rep["lecture"] == "lissé") & (rep["type"] == "stim")]
        c_ = c_.groupby(c_["roi"].astype(str))["dc_mM"].mean()
        ratio = pd.DataFrame({"dR_R0_pct": a, "dc_mM": c_}).dropna()
    figure_reponses(rep, ratio, D["results"] / "filtre_reponses.png")

    print("\n  Réponses (Δc, mM, médiane des somas ; hausses / n ; p de Wilcoxon) :")
    for r in rs.itertuples():
        print(f"    {r.lecture:10s} {r.type:6s} {r.numero} : {r.dc_mM_median:+6.2f}  ({r.n_hausse}/{r.n}, "
              f"p = {r.p_wilcoxon:.1g})")
    for lec, s in sct.items():
        print(f"  {lec} : par soma, stims {s['stim']:+.2f} contre témoins {s['temoin']:+.2f} mM "
              f"(n = {s['n']}, p = {s['p']:.2g})")
    if ratio is not None and len(ratio) > 2:
        rho, p = stats.spearmanr(ratio["dR_R0_pct"], ratio["dc_mM"])
        print(f"  Δc (lissé) contre ΔR/R0 de l'étape 5, par soma : ρ = {rho:.2f} (p = {p:.2g})")
    print(f"  η : {rt['eta'][0]:.3f} (1re tranche) → {rt['eta'][-1]:.3f} ± {100 * rt['eta_sd'][-1]:.2f} % ; "
          f"poids des vitesses de c {np.round(rt['poids'][-1], 3)} (s_c = {f['s_c_candidats']})")
    print(f"  c lissé, médiane des somas sur la session : {np.nanmedian(np.where(lues.T, lis['c'], np.nan)):.1f} mM")
    print(f"Résultats -> {D['results']}")


if __name__ == "__main__":
    main()
