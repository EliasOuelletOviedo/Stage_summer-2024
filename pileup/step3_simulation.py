"""Étape 3 — biais mesuré et simulation de l'empilement (canal 1).

1. Biais mesuré : pour chaque carte et solution, τ de toutes les tranches de taux de
   toutes les acquisitions, ajusté par τ = τ0 · (1 + β · μ_enr) (moindres carrés pondérés
   par l'erreur de τ). β = biais relatif par photon enregistré par impulsion. Comme
   contrôle, β « intra-acquisition » : même pente, une ordonnée par acquisition (insensible
   aux écarts d'une acquisition à l'autre).
2. Temps mort de la 150N. Le biais simulé oscille avec D, de période P : il s'annule
   quand la carte se libère à la phase du photon qui l'a bloquée (D ≈ multiple de P).
   D est donc choisi sur deux mesures à la fois : les temps morts effectifs
   1/r_TAC − 1/r_CFD (116 ns faible Cl⁻, 120 ns fort Cl⁻) et le biais mesuré au même
   μ_enr, avec la pente intra-acquisition : pour la 150N, la pente entre acquisitions
   est presque nulle alors que celle de l'intérieur des images ne l'est pas, et seule
   la seconde compare des photons enregistrés dans les mêmes conditions. Balayage de D
   pour les deux variantes (« photon », « synchro »).
3. Courbes : biais de τ en fonction de μ_det pour la QC-104 (temps mort des intervalles)
   et pour la 150N (D retenu de chaque variante).

Sorties : resultats/pentes.csv, resultats/balayage_150N.csv, resultats/simulation.csv.
Durée : ~15 min.
"""
import csv

import numpy as np

from ajustement import Modele
from config import FENETRE, MESURES, PERIODE_NS, RESULTATS, SERIES, cle, dossier_donnees, preparer
from lecture import irf
from simulation import gabarit, simuler, temps_mort_effectif

MU_DET = (0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25, 0.30)
MU_BALAYAGE = 0.15            # μ_det des acquisitions de la 150N au plus fort taux (~12 MHz CFD)
GRILLE_D = np.arange(98, 128.01, 0.5)
PHOTONS_BALAYAGE = 1_500_000
PHOTONS_COURBES = 3_000_000
SIGMA_D = 1.5                 # ns : dispersion des temps morts effectifs mesurés
SIGMA_BIAIS = 0.007           # bruit de la simulation et de la mesure sur le biais


def ecrire(nom, lignes):
    with open(RESULTATS / nom, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lignes[0]))
        w.writeheader()
        for l in lignes:
            w.writerow({c: (f"{v:.6g}" if isinstance(v, float) else v) for c, v in l.items()})


def pentes():
    """{(carte, sol): dict(tau0, beta, err_beta, beta_intra, mu_max)} sur les tranches du canal 1."""
    rows = [r for r in csv.DictReader(open(RESULTATS / "durees_de_vie.csv"))
            if r["serie"] == "spc" and r["canal"] == "1" and r["tranche"] != "global"]
    out = {}
    for carte, sol, _, _ in SERIES:
        rr = [r for r in rows if r["carte"] == carte and r["solution"] == sol
              and np.isfinite(float(r["err_tau_ns"]))]
        mu = np.array([float(r["mu_enr"]) for r in rr])
        tau = np.array([float(r["tau_ns"]) for r in rr])
        w = 1 / np.array([float(r["err_tau_ns"]) for r in rr]) ** 2
        A = np.c_[np.ones_like(mu), mu]
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        a, b = cov @ (A.T @ (w * tau))
        # intra-acquisition : une ordonnée par puissance, pente commune
        p = np.array([int(r["puissance"]) for r in rr])
        ps = np.unique(p)
        B = np.c_[(p[:, None] == ps[None, :]).astype(float), mu]
        cB = np.linalg.inv(B.T @ (B * w[:, None]))
        cb = cB @ (B.T @ (w * tau))
        tau0 = float(a)
        out[(carte, sol)] = {"carte": carte, "solution": sol, "tau0_ns": tau0,
                             "beta": float(b) / tau0, "err_beta": float(np.sqrt(cov[1, 1])) / tau0,
                             "beta_intra": float(cb[-1]) / tau0,
                             "err_beta_intra": float(np.sqrt(cB[-1, -1])) / tau0,
                             "mu_enr_max": float(mu.max()), "tranches": len(rr)}
        o = out[(carte, sol)]
        print(f"  {carte:10s} {sol:4s} : τ0 = {tau0:.3f} ns, biais {100 * o['beta'] * 0.1:+.2f} % "
              f"par 0,1 photon/impulsion (intra-acquisition {100 * o['beta_intra'] * 0.1:+.2f} %)")
    return out


def main():
    D = dossier_donnees()
    preparer()
    P = pentes()
    ecrire("pentes.csv", list(P.values()))

    tm = list(csv.DictReader(open(RESULTATS / "temps_morts.csv")))
    d_qc = float(next(r["temps_mort_ns"] for r in tm if r["carte"] == "SPC-QC-104" and r["canal"] == "1"))
    deff = {s: np.mean([float(r["temps_mort_ns"]) for r in tm
                        if r["carte"] == "SPC-150N" and r["solution"] == s]) for s in ("low", "high")}

    gab, ref = {}, {}
    for carte, sol, _, ps in SERIES:
        g = np.load(MESURES / f"{cle(carte, sol, ps[0])}_ch1.npz")["dec_global"]
        gab[(carte, sol)] = gabarit(g, FENETRE[carte])
        M = Modele(irf(D, carte, 1), carte)
        tau0, s0, _, _, _ = M.ajuster(g)
        ref[(carte, sol)] = (M, tau0, s0)

    def biais(carte, sol, h):
        M, tau0, s0 = ref[(carte, sol)]
        return M.ajuster(h, s=s0, tau0=tau0)[0] / tau0 - 1

    # --- balayage du temps mort de la 150N -------------------------------------------------
    bal = []
    for mode in ("photon", "synchro"):
        for Dn in GRILLE_D:
            chi2 = 0.0
            for sol in ("low", "high"):
                h, me, md = simuler(gab[("SPC-150N", sol)], MU_BALAYAGE, Dn, mode,
                                    photons=PHOTONS_BALAYAGE, graine=7)
                bs = biais("SPC-150N", sol, h)
                bm = P[("SPC-150N", sol)]["beta_intra"] * me
                de = temps_mort_effectif(me, md)
                chi2 += ((de - deff[sol]) / SIGMA_D) ** 2 + ((bs - bm) / SIGMA_BIAIS) ** 2
                bal.append({"mode": mode, "D_ns": Dn, "solution": sol, "mu_det": md, "mu_enr": me,
                            "D_effectif_ns": de, "D_effectif_mesure_ns": deff[sol],
                            "biais_simule": bs, "biais_mesure": bm,
                            "biais_mesure_entre": P[("SPC-150N", sol)]["beta"] * me})
            for l in bal[-2:]:
                l["chi2"] = chi2
        print(f"  balayage « {mode} » fait")
    ecrire("balayage_150N.csv", bal)
    choix = {}
    for mode in ("photon", "synchro"):
        b = min((l for l in bal if l["mode"] == mode), key=lambda l: l["chi2"])
        choix[mode] = b["D_ns"]
        print(f"  150N « {mode} » : D = {b['D_ns']:.1f} ns (χ² = {b['chi2']:.1f} pour 4 mesures)")

    # --- courbes ------------------------------------------------------------------------------
    lignes = []
    cas = [("SPC-QC-104", "photon", d_qc)] + [("SPC-150N", m, choix[m]) for m in ("photon", "synchro")]
    for carte, mode, Dn in cas:
        for sol in ("low", "high"):
            _, tau0, _ = ref[(carte, sol)]
            for mu in MU_DET:
                h, me, md = simuler(gab[(carte, sol)], mu, Dn, mode, photons=PHOTONS_COURBES, graine=11)
                b = biais(carte, sol, h)
                lignes.append({"carte": carte, "mode": mode, "temps_mort_ns": Dn, "solution": sol,
                               "mu_det": md, "mu_enr": me, "taux_enr_MHz": me / PERIODE_NS * 1e3,
                               "tau0_ns": tau0, "tau_ns": tau0 * (1 + b), "biais_tau": b})
            f = lignes[-1]
            print(f"  {carte:10s} {mode:7s} {sol:4s} : à μ_det = {f['mu_det']:.2f} → μ_enr = {f['mu_enr']:.3f} "
                  f"({f['taux_enr_MHz']:.1f} MHz), biais {100 * f['biais_tau']:+.1f} %")
    ecrire("simulation.csv", lignes)
    print("  -> pentes.csv, balayage_150N.csv, simulation.csv")


if __name__ == "__main__":
    main()
