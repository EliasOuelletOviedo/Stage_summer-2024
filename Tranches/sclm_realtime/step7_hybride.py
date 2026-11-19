"""Étape 7 — lecture hybride : le niveau absolu du temps de vie, la précision du rapport.

Usage :  python3 step7_hybride.py <dossier>          (quelques secondes ; après les étapes 3 et 4)

Le rapport YFP/CFP est ~10 fois plus précis que τ CFP, mais son échelle η (trajet optique,
YFP sans FRET…) n'est pas celle de la calibration du 2026-11-09 : lu avec elle, il donne des
concentrations négatives. τ CFP, lu avec la même calibration, donne des valeurs plausibles.
La lecture hybride prend le niveau de τ et les variations du rapport :

  1. ancres : dans chaque fenêtre calme (avant la 1re stim ; de la fin de chaque stim +
     EXCLURE_APRES_STIM_S jusqu'à la stim suivante ou la fin), τ CFP moyen → c_τ par
     l'isotherme de τ mono du 2026-11-09, et ln R moyen ; l'échelle du rapport y vaut
     η = R / R_iso(c_τ), R_iso l'isotherme du rapport du 2026-11-09 ;
  2. (A) ancrage continu : ln η interpolé linéairement entre les ancres (constant avant la
     première et après la dernière) — le niveau lent suit τ ;
     (B) ancrage unique : η de la première ancre (avant toute stim) pour toute la session —
     le niveau lent suit le rapport ;
  3. chaque passe : c = R_iso⁻¹(R / η(t)), concentrations négatives gardées.

Diagnostic de la dérive lente : Δτ contre Δ ln R entre ancres, comparés à deux directions —
un changement de Cl⁻ (pentes des deux isothermes du 2026-11-09 au niveau de la 1re ancre) et
une perte d'accepteur (YFP blanchie : ces capteurs n'ont plus de FRET ; modèle simple à une
efficacité E = 1 − τ_DA/τ_D, τ_D = 2,72 ns mesuré en solution : Δτ / Δ ln R = −τ_D·E).

Réponses aux stims en mM, dérive lente retirée : ln R redressé par la base quadratique de M1
(étape 4), converti avec la pente locale de l'isotherme du rapport au niveau de la lecture.

Sorties : results/hybride.csv, results/hybride_ancres.csv, results/hybride_reponses.csv,
results/hybride.png, results/hybride_ancrage_unique.png (le rapport seul, ancré une fois),
results/avant_apres_stim.png et .csv (avant / après chaque stim : [Cl⁻] ancrée, ΔR/R₀, Δτ CFP),
results/avant_apres_reunies.png (Δ[Cl⁻] avant / après, les stims réunies)
"""
import numpy as np
import pandas as pd

import config as C
import figures as G
from figures import plt
from step3_series import FIN_S, blocs_stim, glissante, moyennes
from step4_reponses import TEMOIN_APRES, TEMOINS_BASE, ajuster
from step5_absolu import CAL, inverser

TAU_D_NS = 2.72          # donneur sans FRET, mesuré en solution (README du 2026-10-09)


def iso(c, X0, Xi, K):
    return Xi + (X0 - Xi) * K / (K + c)


def dlnR_dc(c, pR):
    return (np.log(iso(c + 1e-3, *pR)) - np.log(iso(c - 1e-3, *pR))) / 2e-3


def main():
    root = C.racine()
    D = C.dirs(root)
    T = pd.read_csv(D["results"] / "series_passes.csv")
    B = pd.read_csv(D["results"] / "bruit.csv")
    cal = pd.read_csv(CAL / "results/calibration.csv").set_index("lecture")
    pR = cal.loc["R", ["X0", "Xinf", "K_mM"]].values.astype(float)
    pT = cal.loc["tau_mono", ["X0", "Xinf", "K_mM"]].values.astype(float)
    t, R, tau = T.t_s.values, T.R.values, T.tau_cfp.values
    stim = T.stim.values.astype(bool)
    blocs = blocs_stim(stim)
    sd_tau = float(B[B.lecture == "tau_cfp"].sd_passe.iloc[0])
    dt = float(np.median(np.diff(t)))

    # --- ancres : fenêtres calmes
    fen = [(0.0, t[blocs[0][0]])]
    for k, (a, b) in enumerate(blocs):
        fin = t[blocs[k + 1][0]] if k + 1 < len(blocs) else t[-1] + 1e-6
        fen.append((t[b - 1] + C.EXCLURE_APRES_STIM_S, fin))
    anc = []
    for a, b in fen:
        m = (t >= a) & (t < b)
        if m.sum() < 20:
            continue
        tb, lr = tau[m].mean(), np.log(R[m]).mean()
        c_tau = inverser(tb, *pT)
        se_tb = sd_tau / np.sqrt(m.sum())
        se_c = abs((inverser(tb + 1e-4, *pT) - inverser(tb - 1e-4, *pT)) / 2e-4) * se_tb
        anc.append(dict(debut_s=a, fin_s=b, t_s=t[m].mean(), passes=int(m.sum()), tau_ns=tb, se_tau_ns=se_tb,
                        c_tau_mM=c_tau, se_c_mM=se_c, R=np.exp(lr), se_lnR=float(np.log(R[m]).std() / np.sqrt(m.sum())),
                        eta=np.exp(lr - np.log(iso(c_tau, *pR))), se_ln_eta=abs(dlnR_dc(c_tau, pR)) * se_c,
                        apres_fin=bool(a >= FIN_S - 1)))
    A = pd.DataFrame(anc)
    A.to_csv(D["results"] / "hybride_ancres.csv", index=False)

    # --- lectures
    ln_eta_A = np.interp(t, A.t_s, np.log(A.eta))
    c_A = inverser(R / np.exp(ln_eta_A), *pR)                       # ancrage continu
    c_B = inverser(R / A.eta.iloc[0], *pR)                          # ancrage unique (avant toute stim)
    c_tau = inverser(tau, *pT)
    se_niveau = np.interp(t, A.t_s, A.se_c_mM)
    pd.DataFrame(dict(passe=T.passe, t_s=t, stim=stim, R=R, tau_cfp=tau, eta_continu=np.exp(ln_eta_A),
                      c_hybride_mM=c_A, c_hybride_ancre_unique_mM=c_B, c_tau_mM=c_tau, se_niveau_mM=se_niveau)) \
        .to_csv(D["results"] / "hybride.csv", index=False)

    # --- diagnostic de la dérive : Δτ contre Δ ln R
    c1 = A.c_tau_mM.iloc[0]
    dtau_dc = (iso(c1 + 1e-3, *pT) - iso(c1 - 1e-3, *pT)) / 2e-3            # ns / mM
    pente_cl = 1000 * dtau_dc / (100 * dlnR_dc(c1, pR))                     # ps par % de ln R
    E = 1 - A.tau_ns.iloc[0] / TAU_D_NS
    pente_acc = -1000 * TAU_D_NS * E / 100                                  # ps par % de ln R
    A["dlnR_pct"] = 100 * (np.log(A.R) - np.log(A.R.iloc[0]))
    A["dtau_ps"] = 1000 * (A.tau_ns - A.tau_ns.iloc[0])

    # --- réponses aux stims, dérive lente retirée (base de M1 sur ln R)
    avant = t < FIN_S
    x = np.log(R[avant])
    s = np.sqrt(1 / T.n_yfp.values[avant] + 1 / T.n_cfp.values[avant])
    m1 = ajuster(t[avant], x, s, stim[avant].astype(float), "M1", dt)
    y = x - (m1["pred"] - m1["reponse"])                                  # ln R redressé
    rep, alig = [], []
    for k, (a, b) in enumerate(blocs, start=1):
        t0, t1 = t[a], t[b - 1] + dt
        if t1 + C.POST_S > FIN_S:
            continue
        m0 = (t[avant] >= t0 - C.PRE_S) & (t[avant] < t0)
        ma = (t[avant] >= t1) & (t[avant] < t1 + C.POST_S)
        dln = y[ma].mean() - y[m0].mean()
        lig = dict(stim=k, delta_lnR_pct=100 * dln)
        for nom, c in (("A", c_A), ("B", c_B)):
            niveau = c[(t >= t0 - C.PRE_S) & (t < t0)].mean()
            lig[f"niveau_{nom}_mM"] = niveau
            lig[f"delta_{nom}_mM"] = dln / dlnR_dc(niveau, pR)
        rep.append(lig)
        g = np.arange(-20, 80, 2.0)
        niveau = c_A[(t >= t0 - C.PRE_S) & (t < t0)].mean()
        base = y[m0].mean()
        tt = t[avant]
        alig.append([(y[(tt >= t0 + u) & (tt < t0 + u + 2)].mean() - base) / dlnR_dc(niveau, pR)
                     if ((tt >= t0 + u) & (tt < t0 + u + 2)).any() else np.nan for u in g])
    RP = pd.DataFrame(rep)
    RP.to_csv(D["results"] / "hybride_reponses.csv", index=False)

    # --- précision
    prec = {}
    for nom, c in (("A", c_A), ("B", c_B), ("tau", c_tau)):
        r = (c - glissante(c))[avant]
        prec[nom] = r.std() * np.sqrt(41 / 40)
    pd.set_option("display.width", 230)
    print(A.round(4).to_string(index=False))
    print(RP.round(3).to_string(index=False))
    print(f"  pentes attendues : Cl⁻ {pente_cl:.1f} ps par % de ln R ; perte d'accepteur {pente_acc:.1f} (E = {E:.2f}) ; "
          f"observées entre ancres : {np.round(A.dtau_ps[1:] / A.dlnR_pct[1:], 1).tolist()}")
    print(f"  précision par passe : hybride {prec['A']:.2f} mM (ancrage continu), {prec['B']:.2f} mM (ancrage unique), "
          f"τ seul {prec['tau']:.1f} mM ; par seconde : {prec['A'] * np.sqrt(dt):.2f}, {prec['B'] * np.sqrt(dt):.2f}, "
          f"{prec['tau'] * np.sqrt(dt):.1f} mM")
    print(f"  niveau : A {c_A[t < 60].mean():.2f} → {c_A[(t > 420) & (t < 450)].mean():.2f} mM ; "
          f"B {c_B[t < 60].mean():.2f} → {c_B[(t > 420) & (t < 450)].mean():.2f} mM")
    figure(D, t, c_A, c_B, tau, pT, A, blocs, prec, RP, se_niveau, alig, pente_cl, pente_acc, E, dt)
    figure_unique(D, t, c_B, A, blocs, prec, RP, dt)
    avant_apres(D, t, c_B, R, tau, A, blocs, dt)
    avant_apres_reunies(D, t, c_B, blocs, dt)


def avant_apres_reunies(D, t, c_B, blocs, dt):
    """Un seul graphique : Δ[Cl⁻] (rapport ancré une fois sur τ) 10 s avant / 10 s après, les stims réunies.
    Chaque stim est rapportée à la moyenne de ses 10 s avant ; valeurs brutes, dérive comprise (comme
    avant_apres_stim.png). Points : moyennes de 1 s, la k-ième seconde avant reliée à la k-ième seconde après par
    un trait mince ; noir : moyenne des stims, IC à 95 % de Student sur les écarts des stims."""
    from scipy import stats
    fen = [(t[a], t[b - 1] + dt) for a, b in blocs if t[b - 1] + dt + C.POST_S <= FIN_S]
    fig, ax = plt.subplots(figsize=(7.2, 6.8))
    G.style(ax)
    ecarts = []
    for i, (t0, t1) in enumerate(fen):
        mb, ma = (t >= t0 - C.PRE_S) & (t < t0), (t >= t1) & (t < t1 + C.POST_S)
        ref = c_B[mb].mean()
        e = c_B[ma].mean() - ref
        ecarts.append(e)
        col = G.PALETTE[i % len(G.PALETTE)]
        secs = []
        for m in (mb, ma):
            k = np.floor(t[m] - t[m].min()).astype(int)
            secs.append(np.array([c_B[m][k == u].mean() for u in np.unique(k)]) - ref)
        n = min(len(secs[0]), len(secs[1]))
        dx = (i - 1.5) * 0.06 + np.linspace(-0.025, 0.025, n)
        for kk in range(n):                                   # k-ième seconde avant ↔ k-ième seconde après
            ax.plot([dx[kk], 1 + dx[kk]], [secs[0][kk], secs[1][kk]], color=G.TRACE, lw=0.6, alpha=0.7, zorder=2)
        for j in (0, 1):
            ax.scatter(j + dx, secs[j][:n], s=16, color=col, alpha=0.8, edgecolors=G.SURFACE, linewidths=0.5, zorder=3,
                       label=f"stim {i + 1} : {e:+.2f} mM" if j == 0 else None)
    ecarts = np.array(ecarts)
    m, se = ecarts.mean(), ecarts.std(ddof=1) / np.sqrt(len(ecarts))
    ic = stats.t.ppf(0.975, len(ecarts) - 1) * se
    p = stats.ttest_1samp(ecarts, 0).pvalue
    ax.errorbar([0, 1], [0, m], yerr=[0, ic], color=G.INK, lw=2.6, marker="o", ms=8, capsize=5, zorder=5,
                label=f"moyenne : {m:+.2f} mM [{m - ic:+.2f} ; {m + ic:+.2f}]")
    ax.axhline(0, color=G.INK_2, lw=0.8, ls="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{C.PRE_S:.0f} s\navant", f"{C.POST_S:.0f} s\naprès"], color=G.INK, fontsize=10.5)
    ax.set_xlim(-0.45, 1.45)
    ax.set_ylabel("Δ[Cl⁻] (mM), par rapport aux 10 s avant chaque stim", color=G.INK, fontsize=10.5)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title(f"Avant / après, les {len(ecarts)} stims réunies", color=G.INK, fontsize=12, loc="left", pad=36)
    ax.text(0, 1.015, f"{int((ecarts > 0).sum())}/{len(ecarts)} ↑ ; p = {p:.1g} (Student sur les {len(ecarts)} écarts)\n"
            "moyennes de 1 s, la k-ième s avant reliée à la k-ième s après ; rapport ancré sur τ avant la 1re stim",
            transform=ax.transAxes, color=G.INK_2, fontsize=9.2, va="bottom")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "avant_apres_reunies.png", dpi=170)
    print(f"  stims réunies : écart {m:+.3f} mM, IC 95 % [{m - ic:+.3f} ; {m + ic:+.3f}], p = {p:.2g}")


def avant_apres(D, t, c_B, R, tau, A, blocs, dt):
    """Avant / après chaque stim, comme avant_apres_stim.png d'octobre, sur la série unique (champ entier) :
    [Cl⁻] du rapport ancré une fois sur τ, ΔR/R₀ et Δτ CFP (référence : la fenêtre de l'ancre, 60 s avant la
    1re stim). Valeurs brutes, dérive comprise. Points : moyennes de 1 s ; noir : moyenne et IC à 95 % des
    passes de la fenêtre ; z = écart / erreur type des passes. Témoins (fausses stims de l'étape 4) : dans le
    CSV seulement."""
    a0 = A.iloc[0]
    ref = (t >= a0.debut_s) & (t < a0.fin_s)
    series = [("c", "[Cl⁻] (mM)", c_B, "mM"), ("R", "ΔR / R$_0$ (%)", 100 * (R / R[ref].mean() - 1), "%"),
              ("tau", "Δτ CFP (ps)", 1000 * (tau - tau[ref].mean()), "ps")]
    fen = [(t[a], t[b - 1] + dt) for a, b in blocs if t[b - 1] + dt + C.POST_S <= FIN_S]
    debuts = [x for x, _ in fen]
    tem = [(x, x + 10.0) for x in TEMOINS_BASE] + [
        (b + TEMOIN_APRES, b + TEMOIN_APRES + 10.0) for _, b in fen
        if b + TEMOIN_APRES + 10.0 + C.POST_S < min([d for d in debuts if d > b] + [FIN_S])]

    def fenetres(x, t0, t1):
        mb, ma = (t >= t0 - C.PRE_S) & (t < t0), (t >= t1) & (t < t1 + C.POST_S)
        return x[mb], x[ma], t[mb], t[ma]

    def par_seconde(x, tt):
        k = np.floor(tt - tt.min()).astype(int)
        return np.array([x[k == i].mean() for i in np.unique(k)])

    lignes = []
    groupes = [(f"stim {i + 1}", [w]) for i, w in enumerate(fen)]          # les témoins : dans le CSV seulement
    fig, axes = plt.subplots(3, len(groupes), figsize=(2.6 * len(groupes) + 1.3, 11.5), sharey="row", squeeze=False)
    for r, (cle, lab, x, unite) in enumerate(series):
        for c, (nom, ws) in enumerate(groupes):
            ax = axes[r][c]
            G.style(ax)
            temoin = nom.startswith("témoin")
            couleur = G.TRACE if temoin else G.PALETTE[c % len(G.PALETTE)]
            moy = []
            for (t0, t1) in ws:
                xb, xa, tb, ta = fenetres(x, t0, t1)
                e = xa.mean() - xb.mean()
                z = e / np.sqrt(xb.var(ddof=1) / len(xb) + xa.var(ddof=1) / len(xa))
                lignes.append(dict(lecture=cle, groupe="témoin" if temoin else nom, debut_s=t0, avant=xb.mean(),
                                   apres=xa.mean(), ecart=e, z=z))
                moy.append((xb.mean(), xa.mean(), xb.std(ddof=1) / np.sqrt(len(xb)), xa.std(ddof=1) / np.sqrt(len(xa))))
                if temoin:
                    ax.plot([0, 1], [xb.mean(), xa.mean()], color=G.TRACE, lw=0.9)
                    ax.scatter([0, 1], [xb.mean(), xa.mean()], s=16, color=G.INK_2, zorder=3)
                else:
                    for j, (v, tt) in enumerate(((xb, tb), (xa, ta))):
                        sec = par_seconde(v, tt)
                        ax.scatter(j + np.linspace(-0.12, 0.12, len(sec)), sec, s=14, color=couleur, alpha=0.75,
                                   edgecolors=G.SURFACE, linewidths=0.5, zorder=3)
            m = np.array(moy)
            if temoin:
                mb, ma = m[:, 0].mean(), m[:, 1].mean()
                ax.plot([0, 1], [mb, ma], color=G.INK, lw=2.2, marker="o", ms=6, zorder=4)
                e = m[:, 1] - m[:, 0]
                titre = f"{nom}\nécart moy. {e.mean():+.2f} {unite}\nde {e.min():+.2f} à {e.max():+.2f} ({len(e)})"
            else:
                (mb, ma, sb, sa), = m
                ax.errorbar([0, 1], [mb, ma], yerr=[1.96 * sb, 1.96 * sa], color=G.INK, lw=2.2, marker="o", ms=6,
                            capsize=4, zorder=4)
                l = lignes[-1]
                titre = f"{nom}\nécart {l['ecart']:+.2f} {unite}\nz = {l['z']:+.1f}"
            ax.set_title(titre, color=G.INK, fontsize=9.5)
            ax.set_xticks([0, 1])
            ax.set_xticklabels([f"{C.PRE_S:.0f} s\navant", f"{C.POST_S:.0f} s\naprès"], color=G.INK, fontsize=9)
            ax.set_xlim(-0.4, 1.4)
            if cle != "c":
                ax.axhline(0, color=G.INK_2, lw=0.8, ls="--")
        axes[r][0].set_ylabel(lab, color=G.INK, fontsize=11)
    fig.suptitle("Avant / après chaque stimulation — champ entier (pas de ROI dans cette session)\n"
                 "[Cl⁻] : rapport ancré sur τ CFP avant la 1re stim, puis fixé ; R₀, τ₀ : 60 s avant la 1re stim\n"
                 "points : moyennes de 1 s ; noir : moyenne et IC à 95 % des passes ; z : écart / erreur type des "
                 "passes ; valeurs brutes, dérive comprise", color=G.INK, fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    G.sauver(fig, D["results"] / "avant_apres_stim.png", dpi=170)
    for cle, _, x, _ in series:                                            # témoins, pour le CSV
        for (t0, t1) in tem:
            xb, xa, _, _ = fenetres(x, t0, t1)
            e = xa.mean() - xb.mean()
            lignes.append(dict(lecture=cle, groupe="témoin", debut_s=t0, avant=xb.mean(), apres=xa.mean(), ecart=e,
                               z=e / np.sqrt(xb.var(ddof=1) / len(xb) + xa.var(ddof=1) / len(xa))))
    pd.DataFrame(lignes).to_csv(D["results"] / "avant_apres_stim.csv", index=False)
    L = pd.DataFrame(lignes)
    pd.set_option("display.width", 200)
    print(L.round(3).to_string(index=False))


def figure_unique(D, t, c_B, A, blocs, prec, RP, dt):
    """Le rapport seul, son échelle ancrée sur τ pendant les 60 s avant la 1re stim, puis fixée."""
    a0 = A.iloc[0]
    fen = [(t[a], t[b - 1]) for a, b in blocs]
    fig, ax = plt.subplots(figsize=(15, 6))
    G.style(ax)
    G.bandes(ax, fen)
    ax.axvspan(a0.debut_s, a0.fin_s, ymin=0, ymax=0.03, color=G.COULEUR_CANAL["cfp"], alpha=0.6, lw=0)
    tt, cc = moyennes(t, c_B, pas=2.0)
    ic = 1.96 * a0.se_c_mM
    ax.plot(t, c_B, color=G.TRACE, lw=0.4, alpha=0.6, label="passe par passe (240 ms)")
    ax.fill_between(tt, cc - ic, cc + ic, color=G.INK, alpha=0.10, lw=0, label=f"IC 95 % du niveau (ancre : ± {ic:.2f} mM)")
    ax.plot(tt, cc, color=G.INK, lw=1.6, label="moyennes de 2 s")
    ax.errorbar([a0.t_s], [a0.c_tau_mM], yerr=[ic], fmt="o", color=G.COULEUR_CANAL["cfp"], ms=7, capsize=4, lw=1.3,
                label=f"ancre : τ CFP moyen, 0–{a0.fin_s:.0f} s")
    lo, hi = np.percentile(c_B, [0.3, 99.7])
    ax.set_ylim(lo - 0.6, hi + 0.9)
    for (t0, t1), r in zip(fen, RP.itertuples()):
        ax.text((t0 + t1) / 2, 0.97, f"{r.delta_B_mM:+.2f} mM", transform=ax.get_xaxis_transform(), ha="center",
                va="top", fontsize=8.5, color=G.INK_2)
    ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
    ax.text(FIN_S + 3, 0.05, "changement optique du champ\n(CFP −5 %, τ stable)", transform=ax.get_xaxis_transform(),
            fontsize=8, color=G.INK_2, va="bottom")
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.set_ylabel("[Cl⁻] (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5, ncol=4, loc="lower left", bbox_to_anchor=(0, 0.04))
    G.titre(ax, "[Cl⁻] lue au rapport YFP/CFP, échelle ancrée sur τ CFP avant la 1re stim puis fixée",
            f"ancre (0–{a0.fin_s:.0f} s) : {a0.c_tau_mM:.1f} ± {ic:.1f} mM, η = {a0.eta:.2f} ; précision "
            f"{prec['B']:.2f} mM par passe, {prec['B'] * np.sqrt(dt):.2f} mM par seconde ; au-dessus des stims : "
            "10 s après − 10 s avant, dérive lente retirée ; isothermes du 2026-11-09")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "hybride_ancrage_unique.png")


def figure(D, t, c_A, c_B, tau, pT, A, blocs, prec, RP, se_niveau, alig, pente_cl, pente_acc, E, dt):
    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.55, 1], hspace=0.42, wspace=0.3)
    fen = [(t[a], t[b - 1]) for a, b in blocs]

    ax = fig.add_subplot(gs[0, :])
    G.style(ax)
    G.bandes(ax, fen)
    for _, r in A.iterrows():
        ax.axvspan(r.debut_s, r.fin_s, ymin=0, ymax=0.02, color=G.COULEUR_CANAL["cfp"], alpha=0.5, lw=0)
    tt, cA = moyennes(t, c_A, pas=2.0)
    _, ss = moyennes(t, se_niveau, pas=2.0)
    _, cB = moyennes(t, c_B, pas=2.0)
    ax.plot(t, c_A, color=G.TRACE, lw=0.4, alpha=0.55, label="hybride, passe par passe (240 ms)")
    ax.fill_between(tt, cA - 1.96 * ss, cA + 1.96 * ss, color=G.INK, alpha=0.10, lw=0, label="IC 95 % du niveau (ancres)")
    ax.plot(tt, cA, color=G.INK, lw=1.6, label="hybride, ancrage continu sur τ (moyennes de 2 s)")
    ax.plot(tt, cB, color=G.PALETTE[1], lw=1.3, ls="--", label="hybride, ancrage unique avant la 1re stim")
    t10, tau10 = moyennes(t, tau, pas=10.0)
    ax.plot(t10, inverser(tau10, *pT), color=G.COULEUR_CANAL["cfp"], lw=1.0, alpha=0.75,
            label="τ CFP seul (τ moyen de 10 s)")
    ax.errorbar(A.t_s, A.c_tau_mM, yerr=1.96 * A.se_c_mM, fmt="o", color=G.COULEUR_CANAL["cfp"], ms=6, capsize=3,
                lw=1.2, label="ancres : τ moyen d'une fenêtre calme (IC 95 %)")
    ax.axvline(FIN_S, color=G.TRACE, lw=0.8, ls="--")
    ax.text(FIN_S + 3, 0.05, "changement optique du champ\n(CFP −5 %, τ stable)", transform=ax.get_xaxis_transform(),
            fontsize=8, color=G.INK_2, va="bottom")
    ax.set_ylim(0, 20)
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.set_ylabel("[Cl⁻] (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.3, ncol=3, loc="upper left")
    G.titre(ax, "[Cl⁻] dans le temps : niveau du temps de vie, variations du rapport",
            f"précision par passe (240 ms) : hybride {prec['A']:.2f} mM contre {prec['tau']:.0f} mM pour τ seul ; "
            "isothermes de la calibration en tranche du 2026-11-09 ; bandes jaunes : 1064 nm")

    ax = fig.add_subplot(gs[1, 0])
    G.style(ax)
    xl = np.array([min(A.dlnR_pct.min(), -5) - 0.5, 0.5])
    ax.plot(xl, pente_cl * xl, color=G.PALETTE[2], lw=1.3, label=f"Cl⁻ (calibration) : {pente_cl:.1f} ps par %")
    ax.plot(xl, pente_acc * xl, color=G.PALETTE[4], lw=1.3, label=f"perte d'accepteur (E = {E:.2f}) : {pente_acc:.0f} ps par %")
    ok = ~A.apres_fin.values
    ax.errorbar(A.dlnR_pct[ok], A.dtau_ps[ok], xerr=1.96 * 100 * A.se_lnR[ok], yerr=1.96 * 1000 * A.se_tau_ns[ok],
                fmt="o", color=G.INK, ms=5, capsize=2, label="ancres (dérive lente)")
    if (~ok).any():
        ax.plot(A.dlnR_pct[~ok], A.dtau_ps[~ok], "o", mfc="none", color=G.INK, ms=6, label="ancre après 495 s")
    r = RP.delta_lnR_pct.mean()
    st = pd.read_csv(D["results"] / "reponses_stims.csv")
    dtau = st[(st.lecture == "tauC") & (st.genre == "stim")].delta_apres
    ax.errorbar([r], [dtau.mean()], yerr=[1.96 * dtau.std(ddof=1) / np.sqrt(len(dtau))], fmt="s", color=G.PALETTE[1],
                ms=6, capsize=2, label="réponse aux stims (moyenne)")
    ax.axhline(0, color=G.TRACE, lw=0.8)
    ax.axvline(0, color=G.TRACE, lw=0.8)
    ax.set_xlabel("Δ ln R (%) depuis la 1re ancre", color=G.INK)
    ax.set_ylabel("Δ τ CFP (ps)", color=G.INK)
    ax.legend(frameon=False, fontsize=7.8)
    G.titre(ax, "La dérive lente n'a pas la signature du Cl⁻", "τ monte trop pour la baisse du rapport", fontsize=10.5)

    ax = fig.add_subplot(gs[1, 1])
    G.style(ax)
    G.bandes(ax, fen)
    ax.plot(t, np.exp(np.interp(t, A.t_s, np.log(A.eta))), color=G.INK, lw=1.2, label="η, ancrage continu")
    ax.errorbar(A.t_s, A.eta, yerr=1.96 * A.se_ln_eta * A.eta, fmt="o", color=G.COULEUR_CANAL["cfp"], ms=5,
                capsize=3, label="ancres (IC 95 %)")
    ax.axhline(A.eta.iloc[0], color=G.PALETTE[1], lw=1.1, ls="--", label="η, ancrage unique")
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.set_ylabel("η (rapport mesuré / attendu)", color=G.INK)
    ax.legend(frameon=False, fontsize=8)
    G.titre(ax, "Échelle du rapport", f"ancrage continu : de {A.eta.iloc[0]:.2f} à {A.eta.iloc[-1]:.2f}", fontsize=10.5)

    ax = fig.add_subplot(gs[1, 2])
    G.style(ax)
    g = np.arange(-20, 80, 2.0) + 1
    for i, a in enumerate(alig):
        ax.plot(g, a, color=G.PALETTE[i % 6], lw=0.8, alpha=0.6, label=f"stim {RP.stim.iloc[i]}")
    ax.plot(g, np.nanmean(alig, 0), color=G.INK, lw=1.8, label="moyenne")
    ax.axvspan(0, fen[0][1] - fen[0][0] + dt, color=G.STIM, alpha=0.45, lw=0, zorder=0)
    ax.axhline(0, color=G.TRACE, lw=0.8)
    ax.set_xlabel("temps depuis le début de la stim (s)", color=G.INK)
    ax.set_ylabel("Δ[Cl⁻] (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    G.titre(ax, "Réponse aux stims, dérive lente retirée",
            f"10 s après − 10 s avant : {RP.delta_A_mM.mean():+.2f} mM (ancrage continu), "
            f"{RP.delta_B_mM.mean():+.2f} (unique)", fontsize=10.5)
    G.sauver(fig, D["results"] / "hybride.png")


if __name__ == "__main__":
    main()
