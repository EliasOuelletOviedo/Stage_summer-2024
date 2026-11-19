"""Étape 3 — lectures passe par passe, dérive, bruit, et ce que FLIMApp suivait en direct.

Usage :  python3 step3_series.py <dossier>          (~1 min ; après l'étape 1)

Chaque passe (210 ms de scan, dont la moitié d'acquisition : 105 ms) donne, sur les seuls
photons d'acquisition (code 1, ceux que FLIMApp analyse) :
  - le rapport R = YFP / CFP ;
  - le temps moyen d'arrivée tm de chaque canal, dans la fenêtre [17, 255) (sans modèle) ;
  - la durée de vie mono-exponentielle τ de chaque canal (flim.py : IRF de FLIMApp décalée,
    fond libre, vraisemblance de Poisson). Le décalage de l'IRF est instrumental : ajusté
    une fois sur le déclin de toute la session, puis fixé passe par passe ;
  - la durée de vie de FLIMApp en direct (frames.csv), pour comparaison.

Déclins de toute la session : 1, 2 et 3 exponentielles (décalage libre) — ce que contient
le déclin (composante prompte, décalage de l'IRF).

Dérive : ln x = a_k + β·t sur les moments calmes (avant la 1re stim ; de la fin de chaque
stim + EXCLURE_APRES_STIM_S jusqu'à la stim suivante), un niveau libre par intervalle (les
marches des stims restent hors de la pente), une pente commune. La fin de session, où le
rapport monte de ~5 %, est décrite à part (FIN_S).

Bruit : résidus passe par passe autour d'une moyenne glissante centrée de 41 passes
(~10 s, stims comprises) comparés au bruit de photons attendu — rapport : R²(1/N_Y + 1/N_C) ;
tm : variance des temps d'arrivée / N ; τ : Monte-Carlo (déclins de Poisson tirés du
déclin moyen de la session, au nombre de photons d'une passe, réajustés). Écart d'Allan
de ln R selon la durée moyennée.

Sorties : results/series_passes.csv, results/declins_session.csv, results/derive.csv,
results/bruit.csv, results/series.png, results/bruit.png, results/declins.png
"""
import numpy as np
import pandas as pd

import config as C
import figures as G
import flim
import spc as S
from figures import plt

FIN_S = 495.0           # après : le rapport monte de ~5 % sans stim (CFP −6 %) — décrit à part
MC = 400                # tirages Monte-Carlo pour le bruit de photons de τ
GRAINE = 20261111


def blocs_stim(stim):
    d = np.diff(np.r_[0, stim.astype(int), 0])
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def calmes(t, stim, blocs):
    """Masque des passes calmes et indice d'intervalle (0 : avant la 1re stim)."""
    ok = np.zeros(len(t), bool)
    inter = np.full(len(t), -1)
    ok[:blocs[0][0]] = True
    inter[:blocs[0][0]] = 0
    for k, (a, b) in enumerate(blocs, start=1):
        fin = blocs[k][0] if k < len(blocs) else len(t)
        m = (t >= t[b - 1] + C.EXCLURE_APRES_STIM_S) & (np.arange(len(t)) < fin) & (t < FIN_S)
        ok |= m
        inter[m] = k
    return ok & (t < FIN_S), inter


def derive(t, y, ok, inter):
    """ln y = a_k + β·t : β et son erreur type (moindres carrés, un niveau par intervalle)."""
    ks = np.unique(inter[ok])
    X = np.c_[np.stack([(inter[ok] == k).astype(float) for k in ks], 1), t[ok]]
    ly = np.log(y[ok])
    coef, *_ = np.linalg.lstsq(X, ly, rcond=None)
    r = ly - X @ coef
    cov = np.linalg.inv(X.T @ X) * r.var(ddof=X.shape[1])
    return coef[-1], np.sqrt(cov[-1, -1])


def glissante(x, n=41):
    k = np.ones(n) / n
    pad = np.r_[np.full(n // 2, x[:n // 2].mean()), x, np.full(n // 2, x[-(n // 2):].mean())]
    return np.convolve(pad, k, mode="valid")


def allan(x, ms):
    """Écart d'Allan (non chevauchant) de la série x pour des moyennes de m points."""
    out = []
    for m in ms:
        n = len(x) // m
        if n < 3:
            out.append(np.nan)
            continue
        y = x[:n * m].reshape(n, m).mean(1)
        out.append(np.sqrt(0.5 * np.mean(np.diff(y) ** 2)))
    return np.array(out)


def main():
    root = C.racine()
    D = C.dirs(root)
    P = C.charger_passes(D)
    t, stim = P["t_s"], P["stim"]
    blocs = blocs_stim(stim)
    rng = np.random.default_rng(GRAINE)
    M = {nom: flim.Modele(C.irf(root, nom)) for nom in ("cfp", "yfp")}
    Hc = {nom: P[f"H_{nom}"][:, S.ACQ].astype(float) for nom in ("cfp", "yfp")}

    # --- déclins de toute la session
    lignes_dec, session = [], {}
    for nom in ("cfp", "yfp"):
        y = Hc[nom].sum(0)
        tau, s, b, dv = M[nom].ajuster_mono(y, tau0=1.5 if nom == "cfp" else 2.9)
        session[nom] = dict(tau=tau, s=s, b=b, y=y)
        lignes_dec.append(dict(canal=nom, modele="1 exp.", taus_ns=f"{tau:.3f}", amplitudes="1", s_ps=1000 * s, fond=b,
                               dev_ddl=dv, part_prompte=np.nan))
        for n in (2, 3):
            r = M[nom].ajuster_multi(y, n, depart=dict(taus=np.geomspace(0.3, 3.5, n), amps=np.full(n, 1 / n), s=s, b=b))
            i = r["amps"] * r["taus"]
            prompte = float(i[0] / i.sum()) if r["taus"][0] < 0.05 else 0.0
            session[nom][f"multi{n}"] = r
            lignes_dec.append(dict(canal=nom, modele=f"{n} exp.", taus_ns=" / ".join(f"{x:.3f}" for x in r["taus"]),
                                   amplitudes=" / ".join(f"{x:.3f}" for x in r["amps"]), s_ps=1000 * r["s"], fond=r["b"],
                                   dev_ddl=r["dev_ddl"], part_prompte=prompte))
    DS = pd.DataFrame(lignes_dec)
    DS.to_csv(D["results"] / "declins_session.csv", index=False)
    pd.set_option("display.width", 220)
    print(DS.round(4).to_string(index=False))

    # --- lectures passe par passe
    T = pd.DataFrame(dict(passe=np.arange(1, len(t) + 1), t_s=t, stim=stim))
    for nom in ("cfp", "yfp"):
        H = Hc[nom]
        T[f"n_{nom}"] = H.sum(1)
        T[f"tm_{nom}"] = [M[nom].temps_moyen(h) for h in H]
        s, tau0 = session[nom]["s"], session[nom]["tau"]
        fits = [M[nom].ajuster_mono(h, s=s, tau0=tau0) for h in H]
        T[f"tau_{nom}"] = [f[0] for f in fits]
        T[f"fond_{nom}"] = [f[2] for f in fits]
        T[f"flimapp_tau_{nom}"] = P[f"flimapp_tau_{nom}"]
        print(f"  {nom.upper()} : τ passe par passe ajustés")
    T["R"] = T.n_yfp / T.n_cfp
    T.to_csv(D["results"] / "series_passes.csv", index=False)

    # --- dérive sur les moments calmes
    ok, inter = calmes(t, stim, blocs)
    der = []
    for col in ("R", "n_cfp", "n_yfp", "tm_cfp", "tm_yfp", "tau_cfp", "tau_yfp", "flimapp_tau_cfp"):
        b, se = derive(t, T[col].values, ok, inter)
        der.append(dict(lecture=col, pente_pct_par_100s=1e4 * b, se=1e4 * se, moyenne=T[col][ok].mean()))
    fin = t >= FIN_S
    avant = (t >= FIN_S - 30) & (t < FIN_S)
    for col in ("R", "n_cfp", "n_yfp", "tm_cfp", "tau_cfp"):
        a0, a1 = T[col][avant].mean(), T[col][t >= t[-1] - 30].mean()
        der.append(dict(lecture=f"{col} (fin : {FIN_S:.0f} s → fin)", pente_pct_par_100s=np.nan, se=np.nan,
                        moyenne=a1, variation_pct=100 * (a1 / a0 - 1)))
    DR = pd.DataFrame(der)
    DR.to_csv(D["results"] / "derive.csv", index=False)
    print(DR.round(4).to_string(index=False))

    # --- bruit
    m_cfp = session["cfp"]["y"] / session["cfp"]["y"].sum()
    lo, hi = C.FENETRE
    tc = C.temps_canaux()
    bruit = []
    for col, attendu in (
            ("R", T.R.values * np.sqrt(1 / T.n_yfp + 1 / T.n_cfp)),
            ("tm_cfp", np.sqrt(np.cov(tc[lo:hi], aweights=Hc["cfp"].sum(0)[lo:hi]) / Hc["cfp"][:, lo:hi].sum(1))),
            ("tm_yfp", np.sqrt(np.cov(tc[lo:hi], aweights=Hc["yfp"].sum(0)[lo:hi]) / Hc["yfp"][:, lo:hi].sum(1)))):
        x = T[col].values
        r = (x - glissante(x))[t < FIN_S]
        bruit.append(dict(lecture=col, sd_passe=r.std() * np.sqrt(41 / 40), sd_photons=np.median(attendu),
                          exces=r.std() * np.sqrt(41 / 40) / np.median(attendu)))
    # τ : Monte-Carlo au nombre de photons d'une passe
    for nom in ("cfp", "yfp"):
        n = int(T[f"n_{nom}"].median())
        p = session[nom]["y"] / session[nom]["y"].sum()
        sims = [M[nom].ajuster_mono(rng.multinomial(n, p), s=session[nom]["s"], tau0=session[nom]["tau"])[0] for _ in range(MC)]
        for col, nm in ((f"tau_{nom}", "hors ligne"), (f"flimapp_tau_{nom}", "FLIMApp")):
            x = T[col].values
            r = (x - glissante(x))[t < FIN_S]
            bruit.append(dict(lecture=col, sd_passe=r.std() * np.sqrt(41 / 40), sd_photons=np.std(sims),
                              exces=r.std() * np.sqrt(41 / 40) / np.std(sims)))
    B = pd.DataFrame(bruit)
    B.to_csv(D["results"] / "bruit.csv", index=False)
    print(B.round(5).to_string(index=False))
    # FLIMApp contre hors ligne
    for nom in ("cfp", "yfp"):
        d = T[f"flimapp_tau_{nom}"] - T[f"tau_{nom}"]
        print(f"  τ {nom.upper()} FLIMApp − hors ligne : {1000 * d.mean():+.1f} ps (écart-type {1000 * d.std():.1f} ps), "
              f"corrélation {np.corrcoef(T[f'flimapp_tau_{nom}'], T[f'tau_{nom}'])[0, 1]:.3f}")

    figures_series(D, T, blocs, ok)
    figure_bruit(D, T, B)
    figure_declins(D, M, session)


def moyennes(t, x, pas=5.0):
    e = np.arange(0, t.max() + pas, pas)
    k = np.digitize(t, e) - 1
    tt = np.array([t[k == i].mean() for i in np.unique(k)])
    xx = np.array([np.nanmean(x[k == i]) for i in np.unique(k)])
    return tt, xx


def figures_series(D, T, blocs, ok):
    t = T.t_s.values
    fen = [(t[a], t[b - 1]) for a, b in blocs]
    fig, axes = plt.subplots(4, 1, figsize=(15, 13), sharex=True)
    for ax in axes:
        G.style(ax)
        G.bandes(ax, fen)
        ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
    ax = axes[0]
    ax.plot(t, T.R, color=G.TRACE, lw=0.5, alpha=0.7)
    ax.plot(*moyennes(t, T.R.values), color=G.INK, lw=1.4)
    ax.set_ylabel("YFP / CFP", color=G.INK)
    G.titre(ax, "Rapport YFP / CFP, passe par passe (gris) et moyennes de 5 s",
            f"photons d'acquisition seulement ; bandes jaunes : 1064 nm à 100 % (PWM 50 %) ; tirets : {FIN_S:.0f} s")
    ax = axes[1]
    for nom in ("cfp", "yfp"):
        x = T[f"n_{nom}"].values
        tt, xx = moyennes(t, x)
        ax.plot(tt, 100 * (xx / np.mean(x[t < 60]) - 1), color=G.COULEUR_CANAL[nom], lw=1.4, label=C.NOMS[nom])
    ax.set_ylabel("photons, écart aux 60 premières s (%)", color=G.INK)
    ax.legend(frameon=False, fontsize=9)
    G.titre(ax, "Photons d'acquisition par canal (moyennes de 5 s)")
    ax = axes[2]
    for col, c, lab in (("tau_cfp", G.COULEUR_CANAL["cfp"], "τ CFP mono (hors ligne)"),
                        ("flimapp_tau_cfp", G.PALETTE[3], "τ CFP de FLIMApp en direct")):
        ax.plot(*moyennes(t, T[col].values), color=c, lw=1.4, label=lab)
    ax2 = ax.twinx()
    ax2.plot(*moyennes(t, T.tm_cfp.values), color=G.PALETTE[2], lw=1.1, alpha=0.8, label="tm CFP (axe de droite)")
    ax2.tick_params(colors=G.INK_2, length=0)
    ax2.set_ylabel("tm CFP (ns)", color=G.INK)
    ax.set_ylabel("τ CFP (ns)", color=G.INK)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8.5, ncol=3)
    G.titre(ax, "Durée de vie et temps moyen du CFP (moyennes de 5 s)")
    ax = axes[3]
    for col, c, lab in (("tau_yfp", G.COULEUR_CANAL["yfp"], "τ YFP mono (hors ligne)"),
                        ("flimapp_tau_yfp", G.PALETTE[1], "τ « ch1 » de FLIMApp en direct (= YFP)")):
        ax.plot(*moyennes(t, T[col].values), color=c, lw=1.4, label=lab)
    ax.set_ylabel("τ YFP (ns)", color=G.INK)
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Durée de vie du YFP (moyennes de 5 s)", "ch1 = IN1 = YFP dans cette session : c'est elle que FLIMApp suivait sur ch1")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "series.png")


def figure_bruit(D, T, B):
    t = T.t_s.values
    fig, axes = plt.subplots(1, 3, figsize=(19, 5))
    for ax in axes:
        G.style(ax)
    ax = axes[0]
    ms = 2 ** np.arange(0, 8)
    dt = np.median(np.diff(t))
    for col, c, lab, rel in (("R", G.INK, "ln R", True), ("tm_cfp", G.PALETTE[2], "tm CFP (ns)", False),
                             ("tau_cfp", G.COULEUR_CANAL["cfp"], "τ CFP (ns)", False)):
        x = T[col].values[t < FIN_S]
        x = np.log(x) if rel else x
        a = allan(x - glissante(x, 401), ms)
        s1 = B[B.lecture == col].sd_photons.iloc[0] / (T.R.mean() if rel else 1)
        ax.loglog(ms * dt, a, "o-", color=c, ms=4, label=lab)
        ax.loglog(ms * dt, s1 / np.sqrt(ms), ":", color=c, lw=1)
    ax.set_xlabel("durée moyennée (s)", color=G.INK)
    ax.set_ylabel("écart d'Allan", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Écart d'Allan (points) et bruit de photons (pointillés)", "passes avant 495 s, tendance lente retirée (401 passes)", fontsize=10.5)
    ax = axes[1]
    z = (T.R - glissante(T.R.values)) / (T.R * np.sqrt(1 / T.n_yfp + 1 / T.n_cfp))
    z = z[t < FIN_S]
    ax.hist(z, bins=60, density=True, color=G.TRACE)
    g = np.linspace(-4, 4, 200)
    ax.plot(g, np.exp(-g ** 2 / 2) / np.sqrt(2 * np.pi), color=G.INK, lw=1.2)
    ax.set_xlabel("résidu du rapport / bruit de photons attendu", color=G.INK)
    G.titre(ax, "Rapport : résidus passe par passe", f"écart-type {z.std():.3f} (1 = bruit de photons seul)", fontsize=10.5)
    ax = axes[2]
    for nom in ("cfp", "yfp"):
        ax.plot(T[f"tau_{nom}"], T[f"flimapp_tau_{nom}"], ".", ms=2, color=G.COULEUR_CANAL[nom], alpha=0.5, label=C.NOMS[nom])
    lim = [min(T.tau_cfp.min(), T.flimapp_tau_cfp.min()), max(T.tau_yfp.max(), T.flimapp_tau_yfp.max())]
    ax.plot(lim, lim, color=G.TRACE, lw=0.8)
    ax.set_xlabel("τ hors ligne, passe par passe (ns)", color=G.INK)
    ax.set_ylabel("τ de FLIMApp en direct (ns)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "FLIMApp en direct contre l'ajustement hors ligne", "même passe, mêmes photons", fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "bruit.png")


def figure_declins(D, M, session):
    fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.2]))
    tc = C.temps_canaux()
    for j, nom in enumerate(("cfp", "yfp")):
        m = M[nom]
        y = session[nom]["y"]
        yw = m.fen(y)
        x = tc[m.lo:m.hi]
        ax = axes[0, j]
        G.style(ax)
        ax.semilogy(x, yw, ".", ms=3, color=G.COULEUR_CANAL[nom], label="session entière (acquisition)")
        r3 = session[nom]["multi3"]
        mu3 = m.attendu(yw, r3["taus"], r3["amps"], r3["s"], r3["b"])
        mu1 = m.attendu(yw, [session[nom]["tau"]], [1.0], session[nom]["s"], session[nom]["b"])
        ax.semilogy(x, mu1, color=G.PALETTE[1], lw=1, label=f"1 exp. ({session[nom]['tau']:.2f} ns)")
        ax.semilogy(x, mu3, color=G.INK, lw=1, label="3 exp.")
        g0 = m.irf_fenetre(0.0)
        g1 = m.irf_fenetre(session[nom]["s"])
        ax.semilogy(x, g0 / g0.max() * yw.max(), color=G.TRACE, lw=1, ls="--", label="IRF de FLIMApp")
        ax.semilogy(x, g1 / g1.max() * yw.max(), color=G.INK_2, lw=1, ls=":",
                    label=f"IRF décalée de {1000 * session[nom]['s']:+.0f} ps")
        ax.set_ylim(yw.max() * 3e-3, yw.max() * 1.6)
        ax.legend(frameon=False, fontsize=8)
        G.titre(ax, f"{C.NOMS[nom]} : déclin de toute la session", f"{yw.sum() / 1e6:.1f} M photons dans la fenêtre", fontsize=10.5)
        ax = axes[1, j]
        G.style(ax)
        ax.plot(x, (yw - mu1) / np.sqrt(mu1), color=G.PALETTE[1], lw=0.8)
        ax.plot(x, (yw - mu3) / np.sqrt(mu3), color=G.INK, lw=0.8)
        ax.set_xlabel("temps (ns)", color=G.INK)
        ax.set_ylabel("résidus (σ)", color=G.INK)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "declins.png")


if __name__ == "__main__":
    main()
