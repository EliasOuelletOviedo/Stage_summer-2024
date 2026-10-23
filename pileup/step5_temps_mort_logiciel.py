"""Étape 5 — QC-104 : temps mort allongé en logiciel pour supprimer l'empilement.

Mesures seulement, sur la MQAE (canal 1 de la QC-104, faible et fort Cl⁻).

Le flux .spc donne le temps exact de chaque photon : le macrotemps (2,048 ns) dit à
quelle impulsion il appartient, le microtemps donne sa phase. On peut donc imposer après
coup un temps mort D' ≥ 7,5 ns (le temps mort natif) : tout photon arrivé moins de D'
après le dernier photon gardé est jeté. Avec D' juste au-dessus d'une période laser, la
carte se libère dans l'impulsion suivante juste après la phase du photon qui l'a bloquée :
les photons précoces qu'elle perd compensent les photons tardifs que l'empilement fait
perdre.

Pour chaque acquisition et chaque D' de GRILLE_D : τ des 8 tranches de taux local (mêmes
tranches que l'étape 1, définies avant le temps mort ; ajustement MLE, exponentielle
incomplète ⊛ IRF, décalage de l'IRF fixé à celui du déclin global natif).

- ppp d'une tranche : photons enregistrés par pixel et par impulsion au temps mort natif
  (taux local mesuré dans l'image) ;
- biais = τ / τ(0) − 1, τ(0) étant l'ordonnée à l'origine de τ(ppp) au temps mort natif,
  toutes tranches réunies (une valeur par solution) ;
- courbes : moyenne des tranches de chaque classe de ppp (± LARGEUR_CLASSE).

Sorties : resultats/temps_mort_logiciel.csv, figures/8_biais_temps_mort.png.
Durée : ~10 min.
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

from ajustement import Modele  # noqa: E402
from config import (FIGURES, LISSAGE_PX, N_CANAUX, N_TRANCHES, PERIODE_NS, PIXEL_NS,  # noqa: E402
                    RESULTATS, SERIES, SOLUTIONS, dossier_donnees, preparer)
from lecture import irf, lire_spc, pixels  # noqa: E402
from simulation import garder  # noqa: E402

D_NATIF = 7.54
GRILLE_D = np.unique(np.round(np.r_[D_NATIF, np.arange(8, 40.01, 1.0), np.arange(11.5, 15.01, 0.5),
                                    np.arange(24, 27.51, 0.5), np.arange(36.5, 40.01, 0.5)], 2))
PPP = (0.02, 0.05, 0.10, 0.15, 0.20)
LARGEUR_CLASSE = 0.025
MIN_PHOTONS = 1_000_000          # détecteur coupé (P25, P70) ou acquisition trop pauvre : écartée

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
RAMPE = dict(zip(PPP, ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]))     # ppp croissant


def temps_exacts(f):
    """Temps absolu exact (ns) de chaque photon d'une entrée de la QC-104 : intervalle
    entre photons successifs = k · P + Δmicrotemps, k entier tiré du macrotemps."""
    dT, du = np.diff(f.t_ns), np.diff(f.micro_ns)
    k = np.round((dT - du) / PERIODE_NS)
    return f.micro_ns[0] + np.r_[0.0, np.cumsum(k * PERIODE_NS + du)]


def mesurer(D, modele):
    tranches = []
    for carte, sol, sous, ps in SERIES:
        if carte != "SPC-QC-104":
            continue
        for p in ps:
            f = lire_spc(carte, D / sous, p)[1]
            if len(f.t_ns) < MIN_PHOTONS:
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
            tau_g, s0, _, _, _ = modele.ajuster(np.bincount(canal, minlength=N_CANAUX).astype(float))
            Tl = T.tolist()
            for Dp in GRILLE_D:
                k = garder(T, Dp, Tl)[dans] if Dp > D_NATIF else np.ones(dans.sum(), bool)
                for j in range(N_TRANCHES):
                    m = tr == j
                    h = np.bincount(canal[k & m], minlength=N_CANAUX).astype(float)
                    tau, _, _, _, err = modele.ajuster(h, s=s0, tau0=tau_g)
                    tranches.append({"solution": sol, "puissance": p, "D_ns": float(Dp), "tranche": j,
                                     "ppp": float(mu[m].mean()), "tau_ns": tau, "err_tau_ns": err,
                                     "photons_gardes": float(k[m].mean())})
            print(f"  {sol} P{p} : {len(GRILLE_D)} temps morts", flush=True)
    for sol in ("low", "high"):
        nat = [t for t in tranches if t["solution"] == sol and t["D_ns"] == D_NATIF]
        tau0 = np.polyfit([t["ppp"] for t in nat], [t["tau_ns"] for t in nat], 1)[1]
        for t in tranches:
            if t["solution"] == sol:
                t["tau0_ns"] = tau0
                t["biais"] = t["tau_ns"] / tau0 - 1
    return tranches


def classes(mes, sol):
    """{ppp: tableau (D', biais moyen, erreur-type, part des photons gardés)} par classe."""
    out = {}
    for ppp in PPP:
        pts = []
        for Dp in GRILLE_D:
            rr = [t for t in mes if t["solution"] == sol and t["D_ns"] == Dp
                  and abs(t["ppp"] - ppp) <= LARGEUR_CLASSE]
            if len(rr) >= 3:
                b = np.array([t["biais"] for t in rr])
                pts.append((Dp, b.mean(), b.std() / np.sqrt(len(b)),
                            np.mean([t["photons_gardes"] for t in rr])))
        if pts:
            out[ppp] = np.array(pts)
    return out


def figure(mes):
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
        for ppp, pts in classes(mes, sol).items():
            col = RAMPE[ppp]
            ax.fill_between(pts[:, 0], 100 * (pts[:, 1] - pts[:, 2]), 100 * (pts[:, 1] + pts[:, 2]),
                            color=col, alpha=0.25, lw=0)
            ax.plot(pts[:, 0], 100 * pts[:, 1], color=col, lw=2, marker="o", ms=3.5, mec=SURFACE, mew=0.5,
                    label=f"{ppp:.2f} ppp")
        ax.set_xlim(7, 40)
        ax.set_xlabel("temps mort imposé D' (ns) ; pointillés : 1, 2, 3 périodes laser (12,49 ns)")
        ax.set_ylabel("biais de τ (%)")
        tau0 = next(t["tau0_ns"] for t in mes if t["solution"] == sol)
        ax.set_title(f"MQAE, {SOLUTIONS[sol]} (τ(0) = {tau0:.2f} ns)".replace(".", ","), loc="left", fontsize=11.5)
        ax.text(D_NATIF + 0.3, ax.get_ylim()[0], "natif\n7,5 ns", fontsize=8.5, color=INK, va="bottom")
    axes[1].legend(title="photons enregistrés\npar pixel et par impulsion\n(temps mort natif)", fontsize=9,
                   title_fontsize=9, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle("SPC-QC-104, données de MQAE : biais de τ selon le temps mort imposé en logiciel",
                 x=0.01, ha="left", fontsize=13)
    fig.text(0.01, -0.03, "Tranches de taux local des acquisitions (8 par image) à ± 0,025 ppp de la valeur ; "
             "moyenne ± erreur-type (bande). τ : MLE, exponentielle incomplète ⊛ IRF mesurée. "
             "τ(0) : τ extrapolé à ppp → 0 au temps mort natif.\nBande grise : ± 1 %.", fontsize=8.5, color=INK_2)
    fig.tight_layout()
    fig.savefig(FIGURES / "8_biais_temps_mort.png", dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("  -> 8_biais_temps_mort.png")


def main():
    D = dossier_donnees()
    preparer()
    mes = mesurer(D, Modele(irf(D, "SPC-QC-104", 1), "SPC-QC-104"))
    with open(RESULTATS / "temps_mort_logiciel.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mes[0]))
        w.writeheader()
        for l in mes:
            w.writerow({c: (f"{v:.6g}" if isinstance(v, float) else v) for c, v in l.items()})
    figure(mes)
    for sol in ("low", "high"):
        for ppp, pts in classes(mes, sol).items():
            i0 = np.argmin(np.abs(pts[:, 0] - D_NATIF))
            i = np.argmin(np.abs(pts[:, 1]))
            print(f"  {sol} {ppp:.2f} ppp : biais natif {100 * pts[i0, 1]:+.1f} %, minimal {100 * pts[i, 1]:+.2f} % "
                  f"à D' = {pts[i, 0]:.2f} ns ({100 * pts[i, 3]:.0f} % des photons gardés)")


if __name__ == "__main__":
    main()
