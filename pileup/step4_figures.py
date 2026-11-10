"""Étape 4 — figures (resultats/figures/).

1_taux.png             taux en fonction de la puissance (deux photons : pente 2), pertes de la 150N
2_temps_morts.png      intervalles entre photons successifs : temps mort de chaque carte
3_declins.png          déclins au plus faible et au plus fort taux local, et leur rapport (mesure, simulation)
4_biais_tranches.png   τ des tranches de taux local en fonction de μ_enr, et la simulation
5_biais_detecteur.png  biais simulé et taux enregistré en fonction des photons au détecteur par impulsion
6_balayage_150N.png    150N : biais et temps mort effectif simulés selon D, et les mesures
7_tau_puissance.png    τ du déclin global de chaque acquisition, par puissance

Durée : ~30 s.
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from config import (F_LASER, FENETRE, FIGURES, MESURES, PERIODE_NS, RESULTATS, SERIES, SOLUTIONS,  # noqa: E402
                    cle, preparer, temps_canaux)
from simulation import gabarit, simuler  # noqa: E402
from step1_photons import BORDS_INTERVALLES_150N, BORDS_INTERVALLES_QC  # noqa: E402
from step2_durees import EXCLUES  # noqa: E402

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
COULEUR = {"SPC-QC-104": "#2a78d6", "SPC-150N": "#eb6834"}
SOMBRE = {"SPC-QC-104": "#104281", "SPC-150N": "#8a3412"}
CLAIR = {"SPC-QC-104": "#9ec5f4", "SPC-150N": "#f6b99c"}
MARQUEUR = {"SPC-QC-104": "o", "SPC-150N": "s"}
NOM = {"SPC-QC-104": "SPC-QC-104", "SPC-150N": "SPC-150N"}
plt.rcParams.update({"font.size": 10, "axes.titlesize": 11.5, "axes.titlelocation": "left",
                     "axes.labelcolor": INK, "text.color": INK, "legend.frameon": False})


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_2, length=0)


def rampe(carte, n):
    """Couleurs de la puissance la plus faible (claire) à la plus forte (foncée)."""
    cm = LinearSegmentedColormap.from_list("r", [CLAIR[carte], COULEUR[carte], SOMBRE[carte]])
    return [cm(x) for x in np.linspace(0, 1, n)] if n > 1 else [COULEUR[carte]]


def sauver(fig, nom):
    fig.savefig(FIGURES / nom, dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {nom}")


def lire(nom):
    return list(csv.DictReader(open(RESULTATS / nom)))


def f(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return np.nan


# ----------------------------------------------------------------------------- 1
def fig_taux(acq):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), facecolor=SURFACE)
    for ax, sol in zip(axes, ("low", "high")):
        style(ax)
        for carte in ("SPC-QC-104", "SPC-150N"):
            rr = sorted((r for r in acq if r["serie"] == "spc" and r["carte"] == carte
                         and r["solution"] == sol and r["canal"] == "1"), key=lambda r: int(r["puissance"]))
            rr = [r for r in rr if (carte, sol, int(r["puissance"]), 1) not in EXCLUES]
            p = np.array([int(r["puissance"]) for r in rr])
            enr = np.array([f(r, "mu_enr_image") for r in rr]) * F_LASER / 1e6
            ax.plot(p, enr, marker=MARQUEUR[carte], color=COULEUR[carte], lw=2, ms=7,
                    label=f"{NOM[carte]} : enregistré")
            if carte == "SPC-QC-104":
                ancre = (p[0], enr[0])
            if carte == "SPC-150N":
                # taux au détecteur pendant le balayage : CFD de StopInfo × (balayage / moyenne)
                cfd = np.array([f(r, "taux_cfd_MHz") / f(r, "taux_stopinfo_MHz") for r in rr]) * enr
                ax.plot(p, cfd, marker="s", mfc=SURFACE, color=COULEUR[carte], lw=1.5, ls="--", ms=7,
                        label="SPC-150N : au détecteur (CFD)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        xr = np.geomspace(ancre[0], max(int(r["puissance"]) for r in acq if r["solution"] == sol) * 0.9, 20)
        ax.plot(xr, ancre[1] * (xr / ancre[0]) ** 2, color=INK_2, lw=0.8, ls=":",
                label="pente 2 (absorption à deux photons)")
        ax.set_xlabel("puissance (consigne, nom du fichier)")
        ax.set_ylabel("taux pendant le balayage, canal 1 (MHz)")
        ax.set_title(f"MQAE, {SOLUTIONS[sol]}")
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xticks([5, 10, 15, 20, 25] if sol == "low" else [10, 20, 30, 40, 60])
        ax.set_xticklabels([str(v) for v in ax.get_xticks()])
        ax.legend(fontsize=8.5, loc="upper left")
    fig.suptitle("Taux de photons : la 150N plafonne vers 5 MHz, la QC-104 suit jusqu'à 14 MHz",
                 x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout()
    sauver(fig, "1_taux.png")


# ----------------------------------------------------------------------------- 2
def fig_temps_morts(tm):
    d_qc = f(next(r for r in tm if r["carte"] == "SPC-QC-104" and r["canal"] == "1"), "temps_mort_ns")
    d150 = {s: np.mean([f(r, "temps_mort_ns") for r in tm if r["carte"] == "SPC-150N" and r["solution"] == s])
            for s in ("low", "high")}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), facecolor=SURFACE)
    for ax in axes:
        style(ax)
    h = np.load(MESURES / f"{cle('SPC-QC-104', 'low', 20)}_ch1.npz")["hist_int"].astype(float)
    c = 0.5 * (BORDS_INTERVALLES_QC[1:] + BORDS_INTERVALLES_QC[:-1])
    axes[0].fill_between(c, h / h[(c > 8) & (c < 12)].mean(), step="mid", color=COULEUR["SPC-QC-104"], lw=0)
    axes[0].axvline(d_qc, color=INK, lw=1)
    axes[0].axvline(PERIODE_NS, color=INK_2, lw=0.8, ls=":")
    axes[0].text(d_qc - 0.4, 1.5, f"temps mort\n{d_qc:.1f} ns".replace(".", ","), fontsize=9, ha="right")
    axes[0].text(PERIODE_NS + 0.4, 1.5, "période laser\n12,49 ns", fontsize=9, color=INK_2)
    axes[0].set_xlim(0, 40)
    axes[0].set_ylim(0, 1.7)
    axes[0].set_xlabel("intervalle entre deux photons successifs de l'entrée 1 (ns)")
    axes[0].set_ylabel("nombre (relatif)")
    axes[0].set_title("SPC-QC-104 (faible Cl⁻, P20) : intervalles exacts")
    h = np.load(MESURES / f"{cle('SPC-150N', 'low', 20)}_ch1.npz")["hist_int"].astype(float)
    c = 0.5 * (BORDS_INTERVALLES_150N[1:] + BORDS_INTERVALLES_150N[:-1])
    axes[1].bar(c, h / h.max(), width=23, color=COULEUR["SPC-150N"], lw=0)
    axes[1].axvline(d150["low"], color=INK, lw=1)
    axes[1].text(260, 0.8, f"temps mort effectif\n{d150['low']:.0f} ns (faible Cl⁻)\n"
                 f"{d150['high']:.0f} ns (fort Cl⁻)", fontsize=9)
    axes[1].set_xlim(0, 600)
    axes[1].set_xlabel("intervalle entre deux photons successifs du module 1 (ns, pas de 25 ns)")
    axes[1].set_title("SPC-150N (faible Cl⁻, P20) : macrotemps de 25 ns")
    fig.suptitle("Temps mort : 7,5 ns pour la QC-104, plus court que la période ; ~9 périodes pour la 150N",
                 x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout()
    sauver(fig, "2_temps_morts.png")


# ----------------------------------------------------------------------------- 3
def fig_declins(sim, choix150):
    t = temps_canaux()
    fig, axes = plt.subplots(2, 4, figsize=(17, 6.6), facecolor=SURFACE,
                             gridspec_kw={"height_ratios": [2, 1.2]}, sharex=True)
    for j, (carte, sol, _, ps) in enumerate(SERIES):
        lo, hi = FENETRE[carte]
        bas = np.load(MESURES / f"{cle(carte, sol, ps[0])}_ch1.npz")
        p_max = max(p for p in ps if np.load(MESURES / f"{cle(carte, sol, p)}_ch1.npz")["dec_global"].sum() > 1e6)
        haut = np.load(MESURES / f"{cle(carte, sol, p_max)}_ch1.npz")
        d0 = bas["dec_global"]
        d1, mu1 = haut["dec_tranches"][-1], haut["mu_tranches"][-1]
        n0, n1 = d0 / d0[lo:hi].sum(), d1 / d1[lo:hi].sum()
        m = slice(lo, hi)
        ax = axes[0, j]
        style(ax)
        ax.semilogy(t[m], n0[m], color=CLAIR[carte], lw=2, label=f"P{ps[0]}, μ_enr = {bas['mu_tranches'].mean():.3f}")
        ax.semilogy(t[m], n1[m], color=SOMBRE[carte], lw=1.2, label=f"P{p_max}, tranche la plus vive, μ_enr = {mu1:.3f}")
        ax.set_title(f"{NOM[carte]}, {SOLUTIONS[sol]}")
        ax.legend(fontsize=8, loc="upper right")
        if j == 0:
            ax.set_ylabel("photons (normalisé)")
        # simulation au même μ_enr : gabarit = déclin de la plus faible puissance
        mode = "photon" if carte == "SPC-QC-104" else choix150[0]
        rr = [r for r in sim if r["carte"] == carte and r["solution"] == sol and r["mode"] == mode]
        me = np.array([f(r, "mu_enr") for r in rr])
        md = np.array([f(r, "mu_det") for r in rr])
        Dn = f(rr[0], "temps_mort_ns")
        mu_det = float(np.interp(mu1, me, md))
        hs, _, _ = simuler(gabarit(d0, FENETRE[carte]), mu_det, Dn, mode, photons=6_000_000, graine=3)
        hs0, _, _ = simuler(gabarit(d0, FENETRE[carte]), 0.002, Dn, mode, photons=6_000_000, graine=4)
        rs = (hs / hs[m].sum()) / (hs0 / hs0[m].sum())
        r = n1 / np.where(n0 > 0, n0, np.nan)
        k = np.ones(5) / 5
        lisse = lambda v: np.convolve(np.nan_to_num(v, nan=1.0), k, "same")
        ax = axes[1, j]
        style(ax)
        ax.axhline(1, color=INK_2, lw=0.8)
        ax.plot(t[m], lisse(r)[m], color=COULEUR[carte], lw=1.4, label="mesure (moyenne sur 5 canaux)")
        ax.plot(t[m], lisse(rs)[m], color=INK, lw=1.2, ls="--", label="simulation")
        ax.set_ylim(0.8, 1.2)
        ax.set_xlabel("temps dans la période (ns)")
        if j == 0:
            ax.set_ylabel("rapport fort / faible taux")
        ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("Forme des déclins : au fort taux, la QC-104 perd des photons tardifs ; la 150N presque pas",
                 x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout()
    sauver(fig, "3_declins.png")


# ----------------------------------------------------------------------------- 4
def fig_tranches(dv, pentes, sim, choix150):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURFACE)
    for ax, sol in zip(axes, ("low", "high")):
        style(ax)
        ax.axhline(0, color=INK_2, lw=0.8)
        for carte, s, _, ps in SERIES:
            if s != sol:
                continue
            tau0 = f(next(r for r in pentes if r["carte"] == carte and r["solution"] == sol), "tau0_ns")
            cols = rampe(carte, len(ps))
            for p, c in zip(ps, cols):
                rr = [r for r in dv if r["serie"] == "spc" and r["carte"] == carte and r["solution"] == sol
                      and r["canal"] == "1" and r["puissance"] == str(p) and r["tranche"] != "global"]
                if not rr:
                    continue
                mu = np.array([f(r, "mu_enr") for r in rr])
                b = 100 * (np.array([f(r, "tau_ns") for r in rr]) / tau0 - 1)
                e = 100 * np.array([f(r, "err_tau_ns") for r in rr]) / tau0
                ax.errorbar(mu, b, yerr=e, color=c, marker=MARQUEUR[carte], ms=4.5, lw=1, elinewidth=0.7,
                            capsize=0, mec=SURFACE, mew=0.5)
            for mode, ls in (("photon", "-"), ("synchro", "--")):
                rr = [r for r in sim if r["carte"] == carte and r["solution"] == sol and r["mode"] == mode]
                if not rr:
                    continue
                me = [0] + [f(r, "mu_enr") for r in rr]
                bs = [0] + [100 * f(r, "biais_tau") for r in rr]
                lab = (f"{NOM[carte]} simulé, temps mort {f(rr[0], 'temps_mort_ns'):.1f} ns"
                       + ("" if carte == "SPC-QC-104" else f" (« {mode} »)"))
                ax.plot(me, bs, color=SOMBRE[carte], lw=1.6, ls=ls, label=lab)
        ax.set_xlim(0, 0.23)
        ax.set_ylim(-17 if sol == "low" else -8, 4)
        ax.set_xlabel("photons enregistrés par impulsion laser, taux local μ_enr\n(× 80 MHz = taux en MHz)")
        ax.set_ylabel("biais de τ (%)")
        ax.set_title(f"MQAE, {SOLUTIONS[sol]}")
        ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("τ par tranche de taux local (8 tranches par image ; une ligne par acquisition, foncé = puissance forte)",
                 x=0.01, ha="left", fontsize=12.5)
    fig.text(0.01, -0.02, "Points : τ ajusté (exponentielle incomplète ⊛ IRF, MLE) relatif à τ0, ordonnée à "
             "l'origine de la droite τ(μ_enr) de la carte. Ronds : QC-104, carrés : 150N.", fontsize=8.5, color=INK_2)
    fig.tight_layout()
    sauver(fig, "4_biais_tranches.png")


# ----------------------------------------------------------------------------- 5
def fig_detecteur(sim, choix150):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), facecolor=SURFACE)
    for ax, sol in zip(axes[:2], ("low", "high")):
        style(ax)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.axhspan(-1, 1, color=GRID, alpha=0.6, lw=0)
        for carte, mode in (("SPC-QC-104", "photon"), ("SPC-150N", choix150[0])):
            rr = [r for r in sim if r["carte"] == carte and r["solution"] == sol and r["mode"] == mode]
            md = [0] + [f(r, "mu_det") for r in rr]
            bs = [0] + [100 * f(r, "biais_tau") for r in rr]
            ax.plot(md, bs, color=COULEUR[carte], marker=MARQUEUR[carte], ms=5, lw=2, label=NOM[carte])
        ax.set_xlabel("photons au détecteur par impulsion, μ_det")
        ax.set_ylabel("biais de τ simulé (%)")
        ax.set_title(f"{SOLUTIONS[sol]} (τ0 ≈ {f(rr[0], 'tau0_ns'):.2f} ns)")
        ax.legend(fontsize=8.5, loc="lower left")
        ax.text(0.005, 1.3, "±1 %", fontsize=8, color=INK_2)
    ax = axes[2]
    style(ax)
    for carte, mode in (("SPC-QC-104", "photon"), ("SPC-150N", choix150[0])):
        rr = [r for r in sim if r["carte"] == carte and r["solution"] == "low" and r["mode"] == mode]
        md = np.array([0] + [f(r, "mu_det") for r in rr])
        me = np.array([0] + [f(r, "mu_enr") for r in rr])
        ax.plot(md * F_LASER / 1e6, me * F_LASER / 1e6, color=COULEUR[carte], marker=MARQUEUR[carte], ms=5,
                lw=2, label=NOM[carte])
    x = np.array([0, 25])
    ax.plot(x, x, color=INK_2, lw=0.8, ls=":", label="sans perte")
    ax.set_xlabel("taux au détecteur (MHz)")
    ax.set_ylabel("taux enregistré (MHz)")
    ax.set_title("Photons enregistrés (faible Cl⁻)")
    ax.legend(fontsize=8.5, loc="upper left")
    fig.suptitle("Simulation : biais de τ et taux enregistré, à même flux de photons au détecteur",
                 x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout()
    sauver(fig, "5_biais_detecteur.png")


# ----------------------------------------------------------------------------- 6
def fig_balayage(bal):
    fig, axes = plt.subplots(2, 2, figsize=(14, 6.8), facecolor=SURFACE, sharex=True)
    for j, sol in enumerate(("low", "high")):
        for i in range(2):
            style(axes[i, j])
        for mode, ls in (("photon", "-"), ("synchro", "--")):
            rr = [r for r in bal if r["mode"] == mode and r["solution"] == sol]
            D = np.array([f(r, "D_ns") for r in rr])
            axes[0, j].plot(D, [100 * f(r, "biais_simule") for r in rr], color=COULEUR["SPC-150N"], ls=ls, lw=1.6,
                            label=f"simulé, variante « {mode} »")
            axes[1, j].plot(D, [f(r, "D_effectif_ns") for r in rr], color=COULEUR["SPC-150N"], ls=ls, lw=1.6)
        r0 = [r for r in bal if r["solution"] == sol]
        bm = 100 * np.mean([f(r, "biais_mesure") for r in r0])
        be = 100 * np.mean([f(r, "biais_mesure_entre") for r in r0])
        axes[0, j].axhspan(bm - 0.7, bm + 0.7, color=GRID, lw=0)
        axes[0, j].axhline(bm, color=INK, lw=1, label="mesuré (pente intra-acquisition)")
        axes[0, j].axhline(be, color=INK_2, lw=1, ls=":", label="mesuré (pente entre acquisitions)")
        de = f(r0[0], "D_effectif_mesure_ns")
        axes[1, j].axhspan(de - 1.5, de + 1.5, color=GRID, lw=0)
        axes[1, j].axhline(de, color=INK, lw=1)
        for k in range(8, 11):
            for i in range(2):
                axes[i, j].axvline(k * PERIODE_NS, color=INK_2, lw=0.6, ls=":")
        axes[0, j].set_title(f"{SOLUTIONS[sol]}, μ_det = 0,15")
        axes[1, j].set_xlabel("temps mort D (ns) ; pointillés : 8, 9, 10 périodes laser")
        if j == 1:
            axes[0, j].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    axes[0, 0].set_ylabel("biais de τ (%)")
    axes[1, 0].set_ylabel("temps mort effectif\n1/r_enr − 1/r_det (ns)")
    axes[1, 0].text(0.01, 0.97, "trait noir et bande grise : mesure ± incertitude", transform=axes[1, 0].transAxes,
                    fontsize=8, color=INK_2, va="top")
    fig.suptitle("SPC-150N : le biais dépend de la phase où la carte se libère (D modulo la période)",
                 x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout()
    sauver(fig, "6_balayage_150N.png")


# ----------------------------------------------------------------------------- 7
def fig_puissance(dv):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), facecolor=SURFACE)
    for ax, sol in zip(axes, ("low", "high")):
        style(ax)
        for carte in ("SPC-QC-104", "SPC-150N"):
            for k, ls, mfc in (("1", "-", None), ("2", ":", SURFACE)):
                rr = sorted((r for r in dv if r["serie"] == "spc" and r["carte"] == carte and r["solution"] == sol
                             and r["canal"] == k and r["tranche"] == "global"), key=lambda r: int(r["puissance"]))
                ax.errorbar([int(r["puissance"]) for r in rr], [f(r, "tau_ns") for r in rr],
                            yerr=[f(r, "err_tau_ns") for r in rr], color=COULEUR[carte], marker=MARQUEUR[carte],
                            ls=ls, lw=1.6, ms=6, mfc=mfc or COULEUR[carte], capsize=0,
                            label=f"{NOM[carte]}, canal {k}" + (" (peu de photons)" if k == "2" else ""))
        ax.set_xlabel("puissance (consigne)")
        ax.set_ylabel("τ du déclin global (ns)")
        ax.set_title(f"MQAE, {SOLUTIONS[sol]}")
        ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("τ de chaque acquisition : ce que donnerait l'analyse de l'image entière",
                 x=0.01, ha="left", fontsize=12.5)
    fig.tight_layout()
    sauver(fig, "7_tau_puissance.png")


def main():
    preparer()
    acq, dv, tm = lire("acquisitions.csv"), lire("durees_de_vie.csv"), lire("temps_morts.csv")
    pentes, sim, bal = lire("pentes.csv"), lire("simulation.csv"), lire("balayage_150N.csv")
    meilleur = min(bal, key=lambda r: f(r, "chi2"))
    choix150 = (meilleur["mode"], f(meilleur, "D_ns"))
    fig_taux(acq)
    fig_temps_morts(tm)
    fig_declins(sim, choix150)
    fig_tranches(dv, pentes, sim, choix150)
    fig_detecteur(sim, choix150)
    fig_balayage(bal)
    fig_puissance(dv)


if __name__ == "__main__":
    main()
