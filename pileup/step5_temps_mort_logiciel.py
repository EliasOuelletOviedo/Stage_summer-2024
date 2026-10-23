"""Étape 5 — QC-104 : temps mort allongé en logiciel pour supprimer l'empilement.

Le flux .spc donne le temps exact de chaque photon : le macrotemps (2,048 ns) dit à
quelle impulsion il appartient, le microtemps donne sa phase. On peut donc imposer après
coup un temps mort D' ≥ 7,5 ns (le temps mort natif) : tout photon arrivé moins de D'
après le dernier photon gardé est jeté. Avec D' juste au-dessus d'une période laser, la
carte se libère dans l'impulsion suivante juste après la phase du photon qui l'a bloquée :
les photons précoces qu'elle perd compensent les photons tardifs que l'empilement fait
perdre (le mécanisme qui protège la 150N, voir étape 3).

- simulation : pour chaque solution et chaque ppp (photons par pixel et par impulsion au
  détecteur, μ_det), un même flux de photons poissoniens subit tous les D' de la grille ;
  biais de τ = τ(D', ppp) / τ(ppp → 0) − 1 ;
- mesure : chaque acquisition de la QC-104 (canal 1) subit les D' de MESURES_D ; τ des 8
  tranches de taux local (mêmes tranches que l'étape 1, définies avant le temps mort) ;
  biais d'une tranche = τ / τ(0) − 1, τ(0) étant l'ordonnée à l'origine de τ(ppp) au
  temps mort natif, toutes tranches réunies (une valeur par solution, comme la
  simulation) ; moyenné par classe de ppp.
  Le ppp d'une tranche est son taux enregistré natif ramené au détecteur par la
  simulation (pertes de 1 à 8 %).

Sorties : resultats/temps_mort_logiciel_simulation.csv, temps_mort_logiciel_mesure.csv,
figures/8_biais_temps_mort.png. Durée : ~15 min.
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

from ajustement import Modele  # noqa: E402
from config import (FENETRE, FIGURES, LISSAGE_PX, MESURES, N_CANAUX, N_TRANCHES, PERIODE_NS,  # noqa: E402
                    PIXEL_NS, RESULTATS, SERIES, SOLUTIONS, cle, dossier_donnees, preparer)
from lecture import irf, lire_spc, pixels  # noqa: E402
from simulation import flux_poisson, gabarit, garder  # noqa: E402

PPP = (0.02, 0.05, 0.10, 0.15, 0.20)
GRILLE_D = np.round(np.r_[7.54, np.arange(8.0, 40.01, 0.25)], 2)
MESURES_D = (7.54, 9.0, 10.0, 11.0, 12.0, 12.5, 13.0, 13.5, 14.0, 15.0, 16.0, 17.5, 19.0,
             21.0, 23.0, 24.5, 25.5, 26.5, 28.0, 31.0, 34.0, 37.5, 38.5)
PHOTONS = 3_000_000
D_NATIF = 7.54
LARGEUR_CLASSE = 0.025        # ppp ± 0,025 : tranches mesurées rattachées à chaque courbe

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
RAMPE = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]     # ppp croissant


def ecrire(nom, lignes):
    with open(RESULTATS / nom, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lignes[0]))
        w.writeheader()
        for l in lignes:
            w.writerow({c: (f"{v:.6g}" if isinstance(v, float) else v) for c, v in l.items()})


def temps_exacts(f):
    """Temps absolu exact (ns) de chaque photon d'une entrée de la QC-104 : intervalle
    entre photons successifs = k · P + Δmicrotemps, k entier tiré du macrotemps."""
    dT, du = np.diff(f.t_ns), np.diff(f.micro_ns)
    k = np.round((dT - du) / PERIODE_NS)
    return f.micro_ns[0] + np.r_[0.0, np.cumsum(k * PERIODE_NS + du)]


def simulation(modele):
    lignes = []
    for carte, sol, _, ps in SERIES:
        if carte != "SPC-QC-104":
            continue
        M = modele
        g = np.load(MESURES / f"{cle(carte, sol, ps[0])}_ch1.npz")["dec_global"]
        gb = gabarit(g, FENETRE[carte])
        tau_g, s0, _, _, _ = M.ajuster(g)
        T, c, _ = flux_poisson(gb, 0.002, PHOTONS, graine=1)
        tau0 = M.ajuster(np.bincount(c[garder(T, D_NATIF)], minlength=N_CANAUX).astype(float),
                         s=s0, tau0=tau_g)[0]
        for ppp in PPP:
            T, c, n_imp = flux_poisson(gb, ppp, PHOTONS, graine=2)
            Tl = T.tolist()
            for Dp in GRILLE_D:
                g_ = garder(T, Dp, Tl)
                h = np.bincount(c[g_], minlength=N_CANAUX).astype(float)
                tau = M.ajuster(h, s=s0, tau0=tau0)[0]
                lignes.append({"solution": sol, "ppp": ppp, "D_ns": float(Dp), "tau_ns": tau,
                               "biais": tau / tau0 - 1, "ppp_enregistre": g_.sum() / n_imp,
                               "photons_gardes": g_.sum() / len(T)})
            print(f"  simulation {sol} ppp {ppp:.2f} : biais {100 * lignes[-len(GRILLE_D)]['biais']:+.1f} % "
                  f"à {D_NATIF} ns, min |biais| {100 * min(abs(l['biais']) for l in lignes[-len(GRILLE_D):]):.2f} %",
                  flush=True)
    return lignes


def mesure(D, modele, sim):
    """τ par tranche de taux local pour chaque D' de MESURES_D, acquisitions de la QC-104."""
    tranches = []
    for carte, sol, sous, ps in SERIES:
        if carte != "SPC-QC-104":
            continue
        # ppp enregistré natif -> ppp au détecteur (simulation à D' = 7,54 ns)
        rr = sorted((l for l in sim if l["solution"] == sol and l["D_ns"] == D_NATIF), key=lambda l: l["ppp"])
        enr = np.r_[0, [l["ppp_enregistre"] for l in rr]]
        det = np.r_[0, [l["ppp"] for l in rr]]
        for p in ps:
            f = lire_spc(carte, D / sous, p)[1]
            if len(f.t_ns) < 1_000_000:          # détecteur coupé (P25, P70) ou trop peu de photons
                continue
            T = temps_exacts(f)
            y, x, expo = pixels(f)
            dans = y >= 0
            C = np.zeros((len(expo), 1024))
            np.add.at(C, (y[dans], x[dans]), 1)
            E = np.repeat(expo[:, None], 1024, axis=1)
            mu = (ndi.gaussian_filter(C, LISSAGE_PX, mode="nearest")
                  / np.maximum(ndi.gaussian_filter(E, LISSAGE_PX, mode="nearest"), 1e-9)
                  / PIXEL_NS * PERIODE_NS)[y[dans], x[dans]]
            bords = np.quantile(mu, np.linspace(0, 1, N_TRANCHES + 1))
            bords[-1] += 1e-12
            tr = np.clip(np.searchsorted(bords, mu, side="right") - 1, 0, N_TRANCHES - 1)
            canal = f.canal[dans]
            g = np.bincount(canal, minlength=N_CANAUX).astype(float)
            tau_g, s0, _, _, _ = modele.ajuster(g)
            Tl = T.tolist()
            for Dp in MESURES_D:
                k = garder(T, Dp, Tl)[dans] if Dp > D_NATIF else np.ones(dans.sum(), bool)
                for j in range(N_TRANCHES):
                    m = tr == j
                    h = np.bincount(canal[k & m], minlength=N_CANAUX).astype(float)
                    tau, _, _, _, _ = modele.ajuster(h, s=s0, tau0=tau_g)
                    mu_j = float(mu[m].mean())
                    tranches.append({"solution": sol, "puissance": p, "D_ns": Dp, "tranche": j,
                                     "ppp_enregistre_natif": mu_j, "ppp": float(np.interp(mu_j, enr, det)),
                                     "tau_ns": tau, "photons_gardes": float(k[m].mean())})
            print(f"  mesure {sol} P{p} : {len(MESURES_D)} temps morts", flush=True)
    # biais : τ / τ(0) − 1, avec τ(0) une seule valeur par solution, comme la simulation :
    # ordonnée à l'origine de τ en fonction du ppp au temps mort natif, toutes tranches réunies
    for sol in ("low", "high"):
        nat = [t for t in tranches if t["solution"] == sol and t["D_ns"] == D_NATIF]
        tau0 = np.polyfit([t["ppp"] for t in nat], [t["tau_ns"] for t in nat], 1)[1]
        for t in tranches:
            if t["solution"] == sol:
                t["tau0_ns"] = tau0
                t["biais"] = t["tau_ns"] / tau0 - 1
    return tranches


def figure(sim, mes):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), facecolor=SURFACE)
    for ax, sol in zip(axes, ("low", "high")):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_2, length=0)
        ax.axhspan(-1, 1, color=GRID, alpha=0.7, lw=0)
        ax.axhline(0, color=INK_2, lw=0.8)
        for n in (1, 2, 3):
            ax.axvline(n * PERIODE_NS, color=INK_2, lw=0.7, ls=":")
        ax.axvline(D_NATIF, color=INK, lw=0.9, ls="--")
        for ppp, col in zip(PPP, RAMPE):
            rr = sorted((l for l in sim if l["solution"] == sol and l["ppp"] == ppp), key=lambda l: l["D_ns"])
            ax.plot([l["D_ns"] for l in rr], [100 * l["biais"] for l in rr], color=col, lw=2,
                    label=f"{ppp:.2f} ppp")
            pts = []
            for Dp in MESURES_D:
                b = [100 * t["biais"] for t in mes if t["solution"] == sol and t["D_ns"] == Dp
                     and abs(t["ppp"] - ppp) <= LARGEUR_CLASSE]
                if len(b) >= 3:
                    pts.append((Dp, np.mean(b), np.std(b) / np.sqrt(len(b))))
            if pts:
                pts = np.array(pts)
                ax.errorbar(pts[:, 0], pts[:, 1], yerr=pts[:, 2], fmt="o", color=col, ms=5.5, mec=SURFACE,
                            mew=0.8, elinewidth=0.8, capsize=0, zorder=5)
        ax.set_xlim(7, 40)
        ax.set_xlabel("temps mort imposé D' (ns) ; pointillés : 1, 2, 3 périodes laser (12,49 ns)")
        ax.set_ylabel("biais de τ (%)")
        tau0 = next(t["tau0_ns"] for t in mes if t["solution"] == sol)
        ax.set_title(f"MQAE, {SOLUTIONS[sol]} (τ ≈ {tau0:.2f} ns)".replace(".", ","), loc="left", fontsize=11.5)
        ax.text(D_NATIF + 0.3, ax.get_ylim()[0] * 0.97, "natif\n7,5 ns", fontsize=8.5, color=INK, va="bottom")
    axes[1].legend(title="photons par pixel\net par impulsion\n(au détecteur)", fontsize=9, title_fontsize=9,
                   frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle("SPC-QC-104 : biais de τ selon le temps mort imposé en logiciel", x=0.01, ha="left", fontsize=13)
    fig.text(0.01, -0.03, "Lignes : simulation (photons poissoniens, déclin mesuré au plus bas taux, ajustement MLE "
             "exponentielle incomplète ⊛ IRF). Points : données de la QC-104 (tranches de taux local à ± 0,025 ppp "
             "de la courbe ; moyenne ± erreur-type).\nBande grise : ± 1 %.", fontsize=8.5, color=INK_2)
    fig.tight_layout()
    fig.savefig(FIGURES / "8_biais_temps_mort.png", dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("  -> 8_biais_temps_mort.png")


def main():
    D = dossier_donnees()
    preparer()
    modele = Modele(irf(D, "SPC-QC-104", 1), "SPC-QC-104")
    sim = simulation(modele)
    ecrire("temps_mort_logiciel_simulation.csv", sim)
    mes = mesure(D, modele, sim)
    ecrire("temps_mort_logiciel_mesure.csv", mes)
    figure(sim, mes)


if __name__ == "__main__":
    main()
