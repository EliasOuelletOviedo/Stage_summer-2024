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
- courbes : moyenne des tranches de chaque classe de ppp (± LARGEUR_CLASSE) ;
- pointillés : biais de la SPC-150N à son temps mort natif (~115 ns), même définition, à
  partir de ses tranches (durees_de_vie.csv, étape 2). Son temps mort plafonne son taux
  vers 0,08 ppp enregistré : pas de ligne au-delà.

`--figure` : refait la figure depuis resultats/temps_mort_logiciel.csv, sans mesurer.

Figure 9 : biais de τ et taux enregistré selon les photons au détecteur par impulsion,
QC-104 au temps mort natif et à D' = 12,5 ns (une période laser), 150N au sien.

Sorties : resultats/temps_mort_logiciel.csv, figures/8_biais_temps_mort.png,
figures/9_biais_detecteur_corrige.png.
Durée : ~10 min.
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

from ajustement import Modele  # noqa: E402
from config import (F_LASER, FIGURES, LISSAGE_PX, N_CANAUX, N_TRANCHES, PERIODE_NS, PIXEL_NS,  # noqa: E402
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


def biais_150n():
    """{solution: {ppp: biais moyen}} de la SPC-150N, canal 1 : τ / τ(0) − 1 sur ses
    tranches de taux local, τ(0) ordonnée à l'origine de τ(ppp) ; classes à ± LARGEUR_CLASSE."""
    rows = [r for r in csv.DictReader(open(RESULTATS / "durees_de_vie.csv"))
            if r["carte"] == "SPC-150N" and r["canal"] == "1" and r["serie"] == "spc" and r["tranche"] != "global"]
    out = {}
    for sol in ("low", "high"):
        rr = [r for r in rows if r["solution"] == sol]
        mu = np.array([float(r["mu_enr"]) for r in rr])
        tau = np.array([float(r["tau_ns"]) for r in rr])
        b = tau / np.polyfit(mu, tau, 1)[1] - 1
        out[sol] = {ppp: float(b[np.abs(mu - ppp) <= LARGEUR_CLASSE].mean()) for ppp in PPP
                    if (np.abs(mu - ppp) <= LARGEUR_CLASSE).sum() >= 3}
    return out


def figure(mes):
    b150 = biais_150n()
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
            if ppp in b150[sol]:
                ax.axhline(100 * b150[sol][ppp], color=col, lw=1.6, ls=(0, (4, 3)))
        presents = ", ".join(f"{q:.2f}".replace(".", ",") for q in b150[sol])
        ax.plot([], [], color=INK_2, lw=1.6, ls=(0, (4, 3)),
                label=f"SPC-150N, temps mort ~115 ns\n(lignes à {presents} ppp seulement :\nson taux plafonne vers 0,08 ppp)")
        ax.set_xlim(7, 40)
        ax.set_xlabel("temps mort imposé D' (ns) ; pointillés : 1, 2, 3 périodes laser (12,49 ns)")
        ax.set_ylabel("biais de τ (%)")
        tau0 = next(t["tau0_ns"] for t in mes if t["solution"] == sol)
        ax.set_title(f"MQAE, {SOLUTIONS[sol]} (τ(0) = {tau0:.2f} ns)".replace(".", ","), loc="left", fontsize=11.5)
        ax.text(D_NATIF + 0.3, ax.get_ylim()[0], "natif\n7,5 ns", fontsize=8.5, color=INK, va="bottom")
    axes[1].legend(title="photons enregistrés\npar pixel et par impulsion\n(temps mort natif)", fontsize=9,
                   title_fontsize=9, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle("Données de MQAE : biais de τ de la QC-104 selon le temps mort imposé en logiciel, "
                 "et de la 150N (pointillés)",
                 x=0.01, ha="left", fontsize=13)
    fig.text(0.01, -0.03, "Tranches de taux local des acquisitions (8 par image) à ± 0,025 ppp de la valeur ; "
             "moyenne ± erreur-type (bande). τ : MLE, exponentielle incomplète ⊛ IRF mesurée. "
             "τ(0) : τ extrapolé à ppp → 0 au temps mort natif.\nBande grise : ± 1 %.", fontsize=8.5, color=INK_2)
    fig.tight_layout()
    fig.savefig(FIGURES / "8_biais_temps_mort.png", dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("  -> 8_biais_temps_mort.png")


D_CORRECTION = 12.5              # D' retenu : une période laser (12,49 ns)
BORDS_DET = np.arange(0, 0.2501, 0.025)


def temps_morts_natifs():
    """Temps mort natif de chaque carte (canal 1) : QC-104 par les intervalles, 150N
    effectif (1/r_TAC − 1/r_CFD) moyen par solution ; temps_morts.csv, étape 2."""
    rows = list(csv.DictReader(open(RESULTATS / "temps_morts.csv")))
    qc = float(next(r["temps_mort_ns"] for r in rows if r["carte"] == "SPC-QC-104" and r["canal"] == "1"))
    n150 = {s: np.mean([float(r["temps_mort_ns"]) for r in rows
                        if r["carte"] == "SPC-150N" and r["canal"] == "1" and r["solution"] == s])
            for s in ("low", "high")}
    return qc, n150


def vers_detecteur(mu_enr, d_ns):
    """ppp au détecteur depuis le ppp enregistré, formule non paralysable :
    r_det = r_enr / (1 − r_enr · d)."""
    return mu_enr / (1 - mu_enr * d_ns / PERIODE_NS)


def par_classe(x, y, bords=BORDS_DET):
    """(x moyen, y moyen, erreur-type) par classe de x, au moins 2 tranches."""
    out = []
    for a, b in zip(bords[:-1], bords[1:]):
        m = (x >= a) & (x < b)
        if m.sum() >= 2:
            out.append((x[m].mean(), y[m].mean(), y[m].std() / np.sqrt(m.sum())))
    return np.array(out)


def figure_detecteur(mes):
    """Figure 9 : biais de τ et taux enregistré selon les photons au détecteur par impulsion,
    QC-104 au temps mort natif et à D' = 12,5 ns, 150N au sien. Données de MQAE seulement."""
    d_qc, d_150 = temps_morts_natifs()
    rows150 = [r for r in csv.DictReader(open(RESULTATS / "durees_de_vie.csv"))
               if r["carte"] == "SPC-150N" and r["canal"] == "1" and r["serie"] == "spc" and r["tranche"] != "global"]
    series = {}
    for sol in ("low", "high"):
        nat = [t for t in mes if t["solution"] == sol and t["D_ns"] == D_NATIF]
        cor = {(t["puissance"], t["tranche"]): t for t in mes if t["solution"] == sol and t["D_ns"] == D_CORRECTION}
        mu = np.array([t["ppp"] for t in nat])
        det = vers_detecteur(mu, d_qc)
        series[("natif", sol)] = (det, np.array([t["biais"] for t in nat]), mu)
        c = [cor[(t["puissance"], t["tranche"])] for t in nat]
        series[("12,5", sol)] = (det, np.array([t["biais"] for t in c]), mu * np.array([t["photons_gardes"] for t in c]))
        rr = [r for r in rows150 if r["solution"] == sol]
        m5 = np.array([float(r["mu_enr"]) for r in rr])
        tau = np.array([float(r["tau_ns"]) for r in rr])
        series[("150N", sol)] = (vers_detecteur(m5, d_150[sol]), tau / np.polyfit(m5, tau, 1)[1] - 1, m5)
    styles = {"natif": ("SPC-QC-104, temps mort natif (7,5 ns)", "#86b6ef", "o", "--"),
              "12,5": ("SPC-QC-104, temps mort 12,5 ns (= période laser)", "#104281", "o", "-"),
              "150N": (f"SPC-150N, temps mort natif (~{np.mean(list(d_150.values())):.0f} ns)", "#eb6834", "s", "-")}

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8), facecolor=SURFACE)
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_2, length=0)
    for ax, sol in zip(axes[:2], ("low", "high")):
        ax.axhspan(-1, 1, color=GRID, alpha=0.7, lw=0)
        ax.axhline(0, color=INK_2, lw=0.8)
        for cle_s, (lab, col, mk, ls) in styles.items():
            det, b, _ = series[(cle_s, sol)]
            pts = par_classe(det, 100 * b)
            ax.fill_between(pts[:, 0], pts[:, 1] - pts[:, 2], pts[:, 1] + pts[:, 2], color=col, alpha=0.25, lw=0)
            ax.plot(pts[:, 0], pts[:, 1], color=col, marker=mk, ms=5, lw=2, ls=ls, mec=SURFACE, mew=0.6, label=lab)
        tau0 = next(t["tau0_ns"] for t in mes if t["solution"] == sol)
        ax.set_title(f"{SOLUTIONS[sol]} (τ(0) ≈ {tau0:.2f} ns)".replace(".", ","), loc="left", fontsize=11.5)
        ax.set_xlabel("photons au détecteur par pixel et par impulsion")
        ax.set_ylabel("biais de τ mesuré (%)")
        ax.set_xlim(0, 0.25)
        ax.legend(fontsize=8, loc="lower left", frameon=False)
    ax = axes[2]
    for cle_s, (lab, col, mk, ls) in styles.items():
        det, _, enr = series[(cle_s, "low")]
        pts = par_classe(det, enr)
        ax.plot(pts[:, 0] * F_LASER / 1e6, pts[:, 1] * F_LASER / 1e6, color=col, marker=mk, ms=5, lw=2, ls=ls,
                mec=SURFACE, mew=0.6, label=lab.split(",")[0] + ("" if cle_s == "150N" else
                                                                  (", natif" if cle_s == "natif" else ", 12,5 ns")))
    ax.plot([0, 24], [0, 24], color=INK_2, lw=0.8, ls=":", label="sans perte")
    ax.set_xlabel("taux au détecteur (MHz)")
    ax.set_ylabel("taux enregistré (MHz)")
    ax.set_title("Photons enregistrés (faible Cl⁻)", loc="left", fontsize=11.5)
    ax.legend(fontsize=8, loc="upper left", frameon=False)
    fig.suptitle("Données de MQAE : biais de τ et taux enregistré, à même flux de photons au détecteur",
                 x=0.01, ha="left", fontsize=13)
    fig.text(0.01, -0.04, "Tranches de taux local (8 par image), regroupées par classes de 0,025 ppp (au moins 2 tranches) ; moyenne ± "
             "erreur-type. Biais = τ / τ(0) − 1, τ(0) extrapolé à ppp → 0 au temps mort natif de chaque carte. "
             "τ : MLE, exponentielle incomplète ⊛ IRF mesurée.\nppp au détecteur = ppp enregistré / (1 − ppp "
             "enregistré · temps mort / période) avec les temps morts mesurés (QC-104 : 7,5 ns ; 150N : 116 ns "
             "faible Cl⁻, 120 ns fort Cl⁻). Bande grise : ± 1 %.", fontsize=8.5, color=INK_2)
    fig.tight_layout()
    fig.savefig(FIGURES / "9_biais_detecteur_corrige.png", dpi=170, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print("  -> 9_biais_detecteur_corrige.png")


def main():
    import sys
    if "--figure" in sys.argv:
        mes = [{k: (v if k == "solution" else float(v)) for k, v in r.items()}
               for r in csv.DictReader(open(RESULTATS / "temps_mort_logiciel.csv"))]
        figure(mes)
        figure_detecteur(mes)
        return
    D = dossier_donnees()
    preparer()
    mes = mesurer(D, Modele(irf(D, "SPC-QC-104", 1), "SPC-QC-104"))
    with open(RESULTATS / "temps_mort_logiciel.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mes[0]))
        w.writeheader()
        for l in mes:
            w.writerow({c: (f"{v:.6g}" if isinstance(v, float) else v) for c, v in l.items()})
    figure(mes)
    figure_detecteur(mes)
    for sol in ("low", "high"):
        for ppp, pts in classes(mes, sol).items():
            i0 = np.argmin(np.abs(pts[:, 0] - D_NATIF))
            i = np.argmin(np.abs(pts[:, 1]))
            print(f"  {sol} {ppp:.2f} ppp : biais natif {100 * pts[i0, 1]:+.1f} %, minimal {100 * pts[i, 1]:+.2f} % "
                  f"à D' = {pts[i, 0]:.2f} ns ({100 * pts[i, 3]:.0f} % des photons gardés)")


if __name__ == "__main__":
    main()
