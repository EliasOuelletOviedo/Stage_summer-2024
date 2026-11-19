"""Étape 4 — réponses aux stims : taille, témoins, cinétique, ordre de grandeur en mM.

Usage :  python3 step4_reponses.py <dossier>          (quelques secondes ; après l'étape 3)

Lectures (étape 3, photons d'acquisition) : ln R, ln photons CFP, ln photons YFP, τ CFP
mono, tm CFP. Les passes après FIN_S (le rapport y monte de ~5 % sans stim) sont écartées.

1. Modèles cinétiques ajustés sur toute la série (base quadratique dans le temps, poids du
   bruit de photons), u = 1 pendant les passes de stim :
     M1  pompe + relaxation :      dx/dt = g·u − x/τ
     M2  avec retard :             u → passe-bas (τa) → dx/dt = g·y − x/τ   (M1 si τa → 0) ;
                                   τa = v·τ avec v ≤ 1 (les deux constantes jouent le même
                                   rôle dans la forme : on appelle τa la plus courte)
     M3  avec part permanente :    M1 sur α·g, plus (1 − α)·g·u intégré sans retour
   Comparaison par AIC (M2 et M3 seulement pour le rapport et le YFP : sur les lectures plus
   bruitées, ils servent à imiter la tendance lente) ; puis M1 avec un gain par stim (τ
   commun ; rapport et photons seulement) : la réponse se répète-t-elle ? La base quadratique de M1, ajustée avec la
   réponse, sert ensuite à redresser chaque lecture. (Une pente linéaire estimée sur les
   moments calmes, −1,2 %/100 s pour le rapport, est trop raide pour les intervalles entre
   stims et y fabrique un faux dépassement de +0,5 % 40 à 100 s après les stims ; la base
   quadratique n'en laisse rien.)
2. Sans modèle, pour chaque stim, sur la lecture redressée : Δ = moyenne des 10 s après la
   fin − moyenne des 10 s avant le début (convention des tranches), et Δ en fin de stim
   (5 dernières s). Témoins : la même mesure autour de fausses stims de 10 s placées là où
   rien ne se passe — deux dans la base (début à 12 et 35 s), une à la fin de chaque
   intervalle (95 s après la fin de la stim, avant la suivante). z = Δ / erreur type de
   bruit de photons (étape 3).
3. Réponse moyenne alignée sur le début des stims (tranches de 2 s), stim par stim.
4. Ordre de grandeur en mM, avec la calibration en tranche du 2026-11-09 (autre tranche,
   imagerie au polygone) : sensibilité relative du rapport d ln R / dc et pentes de τ et tm
   à une [Cl⁻] basale supposée (0, 4,4 ou 10 mM). Seulement indicatif (voir le README).

Sorties : results/reponses_stims.csv, results/reponses_modeles.csv, results/reponses_gains.csv,
results/reponses_mM.csv, results/reponses.png, results/reponses_temoins.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.signal import lfilter

import config as C
import figures as G
from figures import plt
from step3_series import FIN_S, blocs_stim

LECTURES = {                      # colonne, transformation, facteur d'affichage, unité, nom
    "lnR": ("R", np.log, 100, "%", "rapport YFP/CFP"),
    "lnC": ("n_cfp", np.log, 100, "%", "photons CFP"),
    "lnY": ("n_yfp", np.log, 100, "%", "photons YFP"),
    "tauC": ("tau_cfp", None, 1000, "ps", "τ CFP mono"),
    "tmC": ("tm_cfp", None, 1000, "ps", "tm CFP"),
}
MODELES = {"lnR": ("M1", "M2", "M3"), "lnY": ("M1", "M2", "M3"), "lnC": ("M1", "M3"),
           "tauC": ("M1",), "tmC": ("M1",)}   # M2/M3 seulement là où le rapport signal/bruit les porte
TEMOINS_BASE = (12.0, 35.0)       # débuts des fausses stims dans la base (s)
TEMOIN_APRES = 95.0               # fausse stim à cette distance de la fin de chaque stim (s)
CALIB_0911 = "../SCLM_calibration_2026-11-09/results/calibration.csv"
BASAUX = (0.0, 4.4, 10.0)


def serie(T, cle):
    col, f, k, unite, nom = LECTURES[cle]
    x = T[col].values.astype(float)
    return f(x) if f else x


def sigma_photons(T, B, cle):
    if cle == "lnR":
        return np.sqrt(1 / T.n_yfp + 1 / T.n_cfp).values
    if cle == "lnC":
        return 1 / np.sqrt(T.n_cfp.values)
    if cle == "lnY":
        return 1 / np.sqrt(T.n_yfp.values)
    col = LECTURES[cle][0]
    return np.full(len(T), float(B[B.lecture == col].sd_photons.iloc[0]))


def moyenne(t, x, a, b):
    m = (t >= a) & (t < b)
    return x[m].mean(), m.sum()


def mesures(t, x, s, debut, fin):
    """Δ après (10 s après − 10 s avant) et Δ fin de stim (5 dernières s − 10 s avant), erreurs types."""
    b0, nb = moyenne(t, x, debut - C.PRE_S, debut)
    a0, na = moyenne(t, x, fin, fin + C.POST_S)
    f0, nf = moyenne(t, x, fin - 5, fin)
    sp = np.nanmedian(s)
    return dict(delta_apres=a0 - b0, se_apres=sp * np.sqrt(1 / nb + 1 / na),
                delta_fin_stim=f0 - b0, se_fin_stim=sp * np.sqrt(1 / nb + 1 / nf))


# ------------------------------------------------------------------ modèles cinétiques
def passe_bas(u, tau, dt):
    a = np.exp(-dt / tau)
    return lfilter([1 - a], [1, -a], u)


def integrateur(u, tau, dt):
    a = np.exp(-dt / tau)
    return lfilter([(1 - a) * tau], [1, -a], u)


def reponse(p, mod, u, dt):
    if mod == "M1":
        g, tau = p
        return integrateur(g * u, tau, dt)
    if mod == "M2":
        g, tau, v = p
        return integrateur(g * passe_bas(u, v * tau, dt), tau, dt)
    if mod == "M3":
        g, tau, alpha = p
        return integrateur(alpha * g * u, tau, dt) + np.cumsum((1 - alpha) * g * u) * dt
    raise ValueError(mod)


def ajuster(t, x, s, u, mod, dt, gains=False, blocs=None):
    """Base quadratique + réponse ; rend paramètres, erreurs types, χ²/ddl, AIC, prédiction."""
    tt = (t - t.mean()) / 100
    if gains:
        us = [np.where((np.arange(len(u)) >= a) & (np.arange(len(u)) < b), 1.0, 0.0) for a, b in blocs]
    p0 = {"M1": [-0.001, 20.0], "M2": [-0.001, 20.0, 0.3], "M3": [-0.001, 20.0, 0.7]}[mod]
    lb = {"M1": [-1, 0.5], "M2": [-1, 0.5, 0.005], "M3": [-1, 0.5, 0.0]}[mod]
    ub = {"M1": [1, 3000], "M2": [1, 3000, 1.0], "M3": [1, 3000, 1.0]}[mod]
    if gains:                        # un gain par stim, τ communs
        p0 = [p0[0]] * len(blocs) + p0[1:]
        lb = [lb[0]] * len(blocs) + lb[1:]
        ub = [ub[0]] * len(blocs) + ub[1:]

    def pred(q):
        base = q[0] + q[1] * tt + q[2] * tt ** 2
        if gains:
            nb = len(blocs)
            return base + sum(reponse([q[3 + k]] + list(q[3 + nb:]), mod, us[k], dt) for k in range(nb))
        return base + reponse(q[3:], mod, u, dt)

    def res(q):
        return (x - pred(q)) / s

    q0 = [np.mean(x), 0.0, 0.0] + p0
    i_tau = 3 + (len(blocs) if gains else 1)          # position de τ (relaxation) dans q
    best = None
    for f in (1.0, 0.3, 3.0):                           # quelques départs de τ
        qq = list(q0)
        qq[i_tau] = q0[i_tau] * f
        r = least_squares(res, qq, bounds=([-np.inf] * 3 + lb, [np.inf] * 3 + ub), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    r = best
    n, k = len(x), len(r.x)
    chi2 = float((r.fun ** 2).sum())
    try:
        cov = np.linalg.pinv(r.jac.T @ r.jac) * chi2 / (n - k)
        se = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        cov = np.full((k, k), np.nan)
        se = np.full(k, np.nan)
    return dict(q=r.x, se=se, cov=cov, chi2_ddl=chi2 / (n - k), aic=n * np.log(chi2 / n) + 2 * k, pred=pred(r.x),
                reponse=pred(r.x) - (r.x[0] + r.x[1] * tt + r.x[2] * tt ** 2))


def main():
    root = C.racine()
    D = C.dirs(root)
    T = pd.read_csv(D["results"] / "series_passes.csv")
    B = pd.read_csv(D["results"] / "bruit.csv")
    ST = pd.read_csv(D["results"] / "stims.csv")
    garde = T.t_s.values < FIN_S
    T = T[garde].reset_index(drop=True)
    t = T.t_s.values
    stim = T.stim.values.astype(bool)
    blocs = blocs_stim(stim)
    dt = float(np.median(np.diff(t)))
    u = stim.astype(float)
    fen = [(r.t_debut_s, r.t_fin_s) for r in ST.itertuples() if r.t_fin_s < FIN_S]
    temoins = [(a, a + 10.0) for a in TEMOINS_BASE] + [(b + TEMOIN_APRES, b + TEMOIN_APRES + 10.0)
                                                       for a, b in fen if b + TEMOIN_APRES + 10 + C.POST_S < min(
                                                           [x for x, _ in fen if x > b] + [FIN_S])]
    pd.set_option("display.width", 230)

    # --- modèles cinétiques (base quadratique + réponse)
    series_x, sig, mods, preds = {}, {}, [], {}
    for cle in LECTURES:
        x, s = serie(T, cle), sigma_photons(T, B, cle)
        series_x[cle], sig[cle] = x, s
        k = LECTURES[cle][2]
        for mod in MODELES[cle]:
            r = ajuster(t, x, s, u, mod, dt)
            preds[(cle, mod)] = r
            q, se, cov = r["q"][3:], r["se"][3:], r["cov"][3:, 3:]
            resp = r["reponse"]
            i_min = int(np.argmax(np.abs(resp[blocs[0][0]:blocs[0][0] + 200]))) + blocs[0][0]
            taua = se_taua = np.nan
            if mod == "M2":                                   # τa = v·τ, erreur type par la covariance
                taua = q[2] * q[1]
                gr = np.array([0.0, q[2], q[1]])
                se_taua = float(np.sqrt(max(gr @ cov @ gr, 0)))
            mods.append(dict(lecture=cle, modele=mod, chi2_ddl=r["chi2_ddl"], aic=r["aic"],
                             g_par_s=k * q[0], se_g=k * se[0], tau_s=q[1], se_tau=se[1],
                             tau_a_s=taua, se_tau_a=se_taua,
                             alpha=q[2] if mod == "M3" else np.nan, se_alpha=se[2] if mod == "M3" else np.nan,
                             regime_permanent=k * q[0] * q[1],
                             pic=k * resp[i_min], t_pic_apres_fin_s=t[i_min] - t[blocs[0][1] - 1]))
    MO = pd.DataFrame(mods)
    MO["delta_aic"] = MO.aic - MO.groupby("lecture").aic.transform("min")
    MO.to_csv(D["results"] / "reponses_modeles.csv", index=False)
    print(MO.round(4).to_string(index=False))

    # --- un gain par stim (M1, τ commun) ; pas pour τ et tm, trop bruités pour cinq paramètres de réponse
    gl = []
    for cle in ("lnR", "lnC", "lnY"):
        r = ajuster(t, series_x[cle], sig[cle], u, "M1", dt, gains=True, blocs=blocs)
        k = LECTURES[cle][2]
        for i in range(len(blocs)):
            gl.append(dict(lecture=cle, stim=i + 1, g_par_s=k * r["q"][3 + i], se_g=k * r["se"][3 + i],
                           tau_s=r["q"][3 + len(blocs)], se_tau=r["se"][3 + len(blocs)]))
    GL = pd.DataFrame(gl)
    GL.to_csv(D["results"] / "reponses_gains.csv", index=False)
    print(GL.round(4).to_string(index=False))

    # --- sans modèle, sur les lectures redressées par la base de M1
    redresse = {cle: series_x[cle] - (preds[(cle, "M1")]["pred"] - preds[(cle, "M1")]["reponse"]) for cle in LECTURES}
    lignes = []
    for cle, (col, f, k, unite, nom) in LECTURES.items():
        x, s = redresse[cle], sig[cle]
        for genre, liste in (("stim", fen), ("témoin", temoins)):
            for i, (a, b) in enumerate(liste, start=1):
                m = mesures(t, x, s, a, b)
                lignes.append(dict(lecture=cle, nom=nom, unite=unite, genre=genre, numero=i, debut_s=a,
                                   delta_apres=k * m["delta_apres"], se_apres=k * m["se_apres"],
                                   z_apres=m["delta_apres"] / m["se_apres"],
                                   delta_fin_stim=k * m["delta_fin_stim"], z_fin_stim=m["delta_fin_stim"] / m["se_fin_stim"]))
    R = pd.DataFrame(lignes)
    R.to_csv(D["results"] / "reponses_stims.csv", index=False)
    print(R.round(3).to_string(index=False))
    synth = R.groupby(["lecture", "genre"]).agg(delta=("delta_apres", "mean"), sd=("delta_apres", "std"),
                                                 z_moyen=("z_apres", "mean"), n=("z_apres", "size")).reset_index()
    print(synth.round(3).to_string(index=False))

    # --- ordre de grandeur en mM (calibration du 2026-11-09)
    cal = pd.read_csv(Path(__file__).parent / CALIB_0911).set_index("lecture")
    conv = []
    m1 = MO[MO.modele == "M1"].set_index("lecture")
    sd_passe = {"lnR": float(B[B.lecture == "R"].sd_passe.iloc[0] / T.R.mean()),
                "tauC": float(B[B.lecture == "tau_cfp"].sd_passe.iloc[0]),
                "tmC": float(B[B.lecture == "tm_cfp"].sd_passe.iloc[0])}
    for cb in BASAUX:
        for cle, lec in (("lnR", "R"), ("tauC", "tau_mono"), ("tmC", "tau_m")):
            X0, Xi, K = cal.loc[lec, ["X0", "Xinf", "K_mM"]]
            pente = (Xi - X0) * K / (K + cb) ** 2                       # dX/dc
            if cle == "lnR":
                pente = pente / (Xi + (X0 - Xi) * K / (K + cb))         # d ln R / dc
            k = LECTURES[cle][2]
            d_obs = R[(R.lecture == cle) & (R.genre == "stim")].delta_apres.mean() / k
            conv.append(dict(basal_suppose_mM=cb, lecture=cle, pente_par_mM=k * pente,
                             delta_apres_mM=d_obs / pente, pic_M1_mM=m1.loc[cle, "pic"] / k / pente,
                             g_M1_mM_par_s=m1.loc[cle, "g_par_s"] / k / pente,
                             sd_passe_mM=sd_passe[cle] / abs(pente), sd_1s_mM=sd_passe[cle] / abs(pente) / np.sqrt(1 / dt)))
    CV = pd.DataFrame(conv)
    CV.to_csv(D["results"] / "reponses_mM.csv", index=False)
    print(CV.round(3).to_string(index=False))

    figure(D, t, redresse, sig, blocs, fen, preds, R)
    figure_temoins(D, R)


def aligner(t, x, debuts, avant=30, apres=110, pas=1.0):
    g = np.arange(-avant, apres, pas)
    out = []
    for d in debuts:
        b = x[(t >= d - C.PRE_S) & (t < d)].mean()
        out.append([np.nanmean(x[(t >= d + a) & (t < d + a + pas)]) - b if ((t >= d + a) & (t < d + a + pas)).any()
                    else np.nan for a in g])
    return g + pas / 2, np.array(out)


def figure(D, t, series_x, sig, blocs, fen, preds, R):
    fig, axes = plt.subplots(1, 5, figsize=(26, 5))
    debuts = [a for a, _ in fen]
    for ax, cle in zip(axes, LECTURES):
        G.style(ax)
        col, f, k, unite, nom = LECTURES[cle]
        g, A = aligner(t, series_x[cle], debuts, pas=2.0)
        for i, a in enumerate(A):
            ax.plot(g, k * a, color=G.PALETTE[i % 6], lw=0.7, alpha=0.45, label=f"stim {i + 1}")
        m = np.nanmean(A, 0)
        se = np.nanstd(A, 0, ddof=1) / np.sqrt(A.shape[0])
        ax.fill_between(g, k * (m - se), k * (m + se), color=G.INK, alpha=0.12, lw=0)
        ax.plot(g, k * m, color=G.INK, lw=1.6, label="moyenne ± ET")
        for mod, ls in (("M1", "--"), ("M2", "-")):
            if (cle, mod) not in preds:
                continue
            r = preds[(cle, mod)]
            gg, AA = aligner(t, r["reponse"], debuts, pas=2.0)
            ax.plot(gg, k * np.nanmean(AA, 0), color=G.PALETTE[1] if mod == "M2" else G.PALETTE[3], lw=1.3, ls=ls,
                    label=f"{mod} (χ²/ddl {r['chi2_ddl']:.2f})")
        ax.axvspan(0, fen[0][1] - fen[0][0], color=G.STIM, alpha=0.45, lw=0, zorder=0)
        ax.axhline(0, color=G.TRACE, lw=0.8)
        ax.set_xlabel("temps depuis le début de la stim (s)", color=G.INK)
        ax.set_ylabel(f"Δ ({unite})", color=G.INK)
        ax.legend(frameon=False, fontsize=7.5, ncol=2)
        r = R[(R.lecture == cle) & (R.genre == "stim")]
        G.titre(ax, nom, f"10 s après − 10 s avant : {r.delta_apres.mean():+.2f} {unite} (z moyen {r.z_apres.mean():+.1f})",
                fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "reponses.png")


def figure_temoins(D, R):
    fig, ax = plt.subplots(figsize=(11, 4.6))
    G.style(ax)
    for j, cle in enumerate(LECTURES):
        for genre, c, dx in (("stim", G.PALETTE[0], -0.12), ("témoin", G.TRACE, 0.12)):
            z = R[(R.lecture == cle) & (R.genre == genre)].z_apres.values
            ax.plot(np.full(len(z), j + dx), z, "o", color=c, ms=6, label=("stims" if genre == "stim" else "témoins")
                    if j == 0 else None)
    ax.axhline(0, color=G.TRACE, lw=0.8)
    for y in (-2, 2):
        ax.axhline(y, color=G.TRACE, lw=0.8, ls=":")
    ax.set_xticks(range(len(LECTURES)))
    ax.set_xticklabels([LECTURES[c][4] for c in LECTURES])
    ax.set_ylabel("z = Δ / erreur type de bruit de photons", color=G.INK)
    ax.legend(frameon=False, fontsize=9)
    G.titre(ax, "10 s après − 10 s avant, en unités de bruit de photons",
            "témoins : fausses stims dans la base et à la fin de chaque intervalle ; base de M1 retirée", fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "reponses_temoins.png")


if __name__ == "__main__":
    main()
