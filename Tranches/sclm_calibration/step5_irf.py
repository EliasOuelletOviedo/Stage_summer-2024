"""Étape 5 — déclin CFP et IRF : décalages, contrôles, trois populations à τ communes.

Usage :  python3 step5_irf.py <dossier>

IRF : celle mesurée par SPCM pour le jeu (data/FLIM/irf*.sdt, `config.charger_irf`), utilisée
telle quelle (fond plat retiré). Données : déclin CFP réuni des somas gardés de chaque FOV
(10⁶ à 10⁷ photons).

1. Par FOV : 3 exponentielles, décalage libre, puis avec un élargissement gaussien libre de
   l'IRF. Contrôles (rapportés, pas de choix à faire) : le décalage doit être proche de 0
   (l'IRF est dans le repère des données : mêmes décalages du TDC) et stable d'un FOV à
   l'autre ; l'élargissement voulu par les données doit rester petit devant l'IRF.
2. Trois durées de vie **communes à tous les FOV**, amplitudes, fond et décalage par FOV
   (le modèle « 1-2-3 » de Simulateurs/bayesian_filter : τ₁ très courte, τ₂ ~0,5 ns, τ_D le
   donneur sans FRET). Les τ servent à l'étape 6 (amplitudes par soma), les amplitudes de
   chaque FOV à l'étape 7.
3. Somas (≥ 20 000 photons) : déviance/ddl avec ces τ communes, et en mono-exponentielle.

Sorties : qc/irf.json, qc/irf_fovs.csv, qc/irf.png, qc/irf_somas.png
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import config as C
import figures as G
from figures import plt
from flim import Modele
from mesures import selection

MIN_PHOTONS_SOMA = 20_000


def declins_fov(D, a):
    """(déclin réuni des somas gardés (2, 256), déclins de chaque soma (n, 2, 256), labels)."""
    m = np.load(D["mesures"] / f"{a['nom']}.npz")
    n = int(m["n_rois"])
    dec = m["declins"][m["garde"]].sum(0)              # (n + 2, 2, 256)
    k = selection(D, a["nom"], n)
    return dec[k].sum(0).astype(float), dec[k].astype(float), k


def global_taus(M, Y, taus0, s0):
    """τ communes à tous les déclins Y, amplitudes / fond / décalage par déclin."""
    def interne(taus):
        tot, res = 0.0, []
        for y, s in zip(Y, s0):
            N = y[M.lo:M.hi].sum()

            def f(p):
                z = np.r_[p[2:], 0.0]
                a = np.exp(z - z.max())
                return M.deviance(y, M.attendu(y, taus, a / a.sum(), p[0], p[1])) / N
            r = minimize(f, [s, 0.02, 0.0, 0.0], method="L-BFGS-B",
                         bounds=[(s - 3, s + 3), (0, 0.5), (-12, 12), (-12, 12)])
            z = np.r_[r.x[2:], 0.0]
            a = np.exp(z - z.max())
            res.append(dict(s=float(r.x[0]), b=float(r.x[1]), amps=a / a.sum(), deviance=r.fun * N))
            tot += r.fun * N
        return tot, res
    r = minimize(lambda lt: interne(np.exp(lt))[0] / 1e6, np.log(taus0), method="Nelder-Mead",
                 options=dict(xatol=1e-3, fatol=1e-6, maxiter=200))
    taus = np.exp(r.x)
    tot, res = interne(taus)
    o = np.argsort(taus)[::-1]
    for d in res:
        d["amps"] = d["amps"][o]
    return taus[o], res, tot


def main():
    root = C.racine()
    D = C.dirs(root)
    A = [a for a in C.acquisitions(root) if a["type"] == "fov" and (D["mesures"] / f"{a['nom']}.npz").exists()]
    if not A:
        sys.exit("Aucune mesure — lance d'abord step4_declins.py")
    donnees = {a["nom"]: declins_fov(D, a) for a in A}
    M = Modele(C.charger_irf("cfp", root))
    lignes = []
    for a in A:
        y = donnees[a["nom"]][0][0]
        r = M.ajuster_multi(y, 3)
        rb = M.ajuster_multi(y, 3, s=r["s"], sigma_b_libre=True)
        lignes.append(dict(acquisition=a["nom"], condition=a["condition"], photons=int(y[M.lo:M.hi].sum()),
                           dev_ddl=r["dev_ddl"], s_canaux=r["s"], s_ps=r["s"] * C.DT_NS * 1000,
                           taus=np.round(r["taus"], 4).tolist(), amps=np.round(r["amps"], 4).tolist(), b=r["b"],
                           sigma_b_ps=rb["sigma_b"] * 1000, dev_ddl_elargie=rb["dev_ddl"]))
    T = pd.DataFrame(lignes)
    T.to_csv(D["qc"] / "irf_fovs.csv", index=False)
    print(T[["acquisition", "photons", "dev_ddl", "s_ps", "sigma_b_ps", "taus"]].round(3).to_string(index=False))

    Y = [donnees[a["nom"]][0][0] for a in A]
    taus0 = np.median(np.stack(T.taus.values), axis=0)
    taus, res, _ = global_taus(M, Y, taus0, T.s_canaux.values)
    print(f"  τ communes (CFP) : {np.round(taus, 4)} ns")
    soma_rows = []
    for a, r in zip(A, res):
        _, dec, k = donnees[a["nom"]]
        for lab, d in zip(k, dec):
            y = d[0]
            N = y[M.lo:M.hi].sum()
            if N < MIN_PHOTONS_SOMA:
                continue
            _, _, dd = M.amplitudes(y, taus, r["s"])
            soma_rows.append(dict(acquisition=a["nom"], soma=int(lab), photons=int(N), dev_ddl_3pop=dd,
                                  dev_ddl_mono=M.ajuster_mono(y)[3]))
    S = pd.DataFrame(soma_rows)
    controles = dict(s_median_ps=float(np.median([r["s"] for r in res]) * C.DT_NS * 1000),
                     s_ecart_type_ps=float(np.std([r["s"] for r in res]) * C.DT_NS * 1000),
                     sigma_b_median_ps=float(T.sigma_b_ps.median()),
                     dev_ddl_somas_3pop=float(S.dev_ddl_3pop.median()) if len(S) else np.nan,
                     dev_ddl_somas_mono=float(S.dev_ddl_mono.median()) if len(S) else np.nan,
                     n_somas=len(S))
    print(f"  contrôles : {json.dumps({k: round(v, 3) for k, v in controles.items()})}")
    out = dict(irf=str(C.irf_sdt(root).name), voie_cfp="cfp", ligne_cfp=C.IRF_LIGNE["cfp"],
               controles=controles, taus_cfp=taus.tolist(),
               fovs={a["nom"]: dict(s=r["s"], b=r["b"], amps=r["amps"].tolist(),
                                    dev_ddl=r["deviance"] / (M.n_fen - 6)) for a, r in zip(A, res)},
               s_median=float(np.median([r["s"] for r in res])))
    (D["qc"] / "irf.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print("  -> qc/irf.json, qc/irf_fovs.csv")
    figure(D, A, M, donnees, taus, res, T, S, controles)


def figure(D, A, M, donnees, taus, res, T, S, controles):
    t = C.temps_canaux()
    conds = [c for c in C.CONDITIONS_FIXES if any(a["condition"] == c for a in A)]
    fig, axes = plt.subplots(2, len(conds), figsize=(4.4 * len(conds), 6.6), sharex=True,
                             gridspec_kw=dict(height_ratios=[3, 1.2]), squeeze=False)
    for j, cond in enumerate(conds):
        i = [k for k, a in enumerate(A) if a["condition"] == cond][0]
        y = donnees[A[i]["nom"]][0][0]
        r = res[i]
        m = M.attendu(y, taus, r["amps"], r["s"], r["b"])
        ax, axr = axes[0, j], axes[1, j]
        for x in (ax, axr):
            G.style(x)
        sl = slice(M.lo, M.hi)
        ax.semilogy(t[sl], y[sl], ".", ms=2.5, color=G.COULEUR_CANAL["cfp"], label="données")
        ax.semilogy(t[sl], m[sl], color=G.INK, lw=1, label="3 populations, τ communes")
        g = M.irf_decalee(r["s"])
        ax.semilogy(t, g / g.max() * y.max(), color=G.PALETTE[1], lw=0.9, label="IRF (décalée)")
        ax.set_ylim(max(y[sl].min(), 1) * 0.5, y.max() * 2)
        rz = (y - m) / np.sqrt(np.maximum(m, 1))
        axr.plot(t[sl], rz[sl], color=G.INK_2, lw=0.7)
        axr.axhline(0, color=G.TRACE, lw=0.8)
        G.titre(ax, f"{A[i]['nom']} — {y[sl].sum() / 1e6:.1f} M photons",
                f"déviance/ddl {r['deviance'] / (M.n_fen - 6):.1f} ; décalage {r['s'] * C.DT_NS * 1000:+.0f} ps",
                fontsize=10.5)
        axr.set_xlabel("temps (ns)", color=G.INK)
        if j == 0:
            ax.set_ylabel("photons (CFP, somas réunis)", color=G.INK)
            axr.set_ylabel("résidus (σ)", color=G.INK)
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"Déclins CFP réunis et IRF du jeu — τ communes {np.round(taus, 3)} ns", color=G.INK, fontsize=12)
    fig.tight_layout()
    G.sauver(fig, D["qc"] / "irf.png")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax in axes:
        G.style(ax)
    axes[0].plot(T.photons / 1e6, T.s_ps, "o", color=G.COULEUR_CANAL["cfp"])
    axes[0].axhline(0, color=G.TRACE)
    axes[0].set_xlabel("photons du déclin réuni (M)", color=G.INK)
    axes[0].set_ylabel("décalage ajusté de l'IRF (ps)", color=G.INK)
    G.titre(axes[0], "Décalage de l'IRF par FOV", f"médiane {controles['s_median_ps']:+.0f} ps, "
            f"écart-type {controles['s_ecart_type_ps']:.0f} ps (attendu ≈ 0 : même repère)")
    axes[1].plot(T.photons / 1e6, T.sigma_b_ps, "o", color=G.COULEUR_CANAL["cfp"])
    axes[1].set_xlabel("photons (M)", color=G.INK)
    axes[1].set_ylabel("élargissement voulu σ_b (ps)", color=G.INK)
    G.titre(axes[1], "Élargissement libre de l'IRF", f"médiane {controles['sigma_b_median_ps']:.0f} ps")
    if len(S):
        axes[2].hist(S.dev_ddl_3pop, bins=np.linspace(0.6, 3, 40), color=G.COULEUR_CANAL["cfp"], alpha=0.8,
                     label="3 populations, τ communes")
        axes[2].hist(S.dev_ddl_mono, bins=np.linspace(0.6, 3, 40), histtype="step", color=G.INK, lw=1.2,
                     label="mono-exponentielle")
        axes[2].legend(frameon=False, fontsize=8)
    axes[2].set_xlabel("déviance/ddl d'un soma (CFP)", color=G.INK)
    G.titre(axes[2], f"Somas (≥ {MIN_PHOTONS_SOMA // 1000} k photons), n = {len(S)}",
            f"médianes : 3 populations {controles['dev_ddl_somas_3pop']:.2f}, mono {controles['dev_ddl_somas_mono']:.2f}")
    fig.tight_layout()
    G.sauver(fig, D["qc"] / "irf_somas.png")


if __name__ == "__main__":
    main()
