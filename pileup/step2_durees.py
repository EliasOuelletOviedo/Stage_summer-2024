"""Étape 2 — durées de vie et temps morts.

- τ : exponentielle incomplète reconvoluée par l'IRF mesurée, ajustée au maximum de
  vraisemblance de Poisson (ajustement.py), pour le déclin global de chaque acquisition
  (décalage de l'IRF libre) et pour ses tranches de taux local (décalage fixé) ;
- temps mort de la QC-104 : demi-montée de l'histogramme des intervalles entre photons
  successifs d'une même entrée (toutes les acquisitions réunies) ;
- temps mort effectif de la 150N : 1/r_TAC − 1/r_CFD (StopInfo), acquisition par acquisition.

Sorties : resultats/durees_de_vie.csv, resultats/temps_morts.csv.
Durée : ~1 min.
"""
import csv

import numpy as np

from ajustement import Modele
from config import (CANAUX, MESURES, MIN_PHOTONS, RESULTATS, SERIE_150N_SDT, SERIES, cle,
                    dossier_donnees, preparer)
from lecture import irf
from step1_photons import BORDS_INTERVALLES_QC

# Acquisitions sans fluorescence, écartées : le détecteur s'est coupé (protection contre
# la surcharge) : le taux s'effondre et le déclin devient plat (bruit de fond seul).
EXCLUES = {("SPC-QC-104", "low", 25, 1): "détecteur ch1 coupé (233 000 photons, taux 0,09 MHz)",
           ("SPC-QC-104", "high", 70, 1): "détecteurs coupés : déclin plat",
           ("SPC-QC-104", "high", 70, 2): "détecteurs coupés : déclin plat"}


def lire_acquisitions():
    rows = list(csv.DictReader(open(RESULTATS / "acquisitions.csv")))
    return {(r["serie"], r["carte"], r["solution"], int(r["puissance"]), int(r["canal"])): r
            for r in rows}


def demi_montee(h, bords):
    """Temps (ns) où l'histogramme des intervalles atteint la moitié de son plateau
    (médiane de 8 à 12 ns), par interpolation linéaire."""
    c = 0.5 * (bords[1:] + bords[:-1])
    plateau = np.median(h[(c > 8) & (c < 12)])
    i = int(np.argmax(h >= plateau / 2))
    return float(np.interp(plateau / 2, [h[i - 1], h[i]], [c[i - 1], c[i]]))


def main():
    D = dossier_donnees()
    preparer()
    acq = lire_acquisitions()
    modeles = {(c, k): Modele(irf(D, c, k), c) for c in ("SPC-150N", "SPC-QC-104") for k in CANAUX}
    lignes = []

    def ajouter(serie, carte, sol, p, k, z):
        M = modeles[(carte, k)]
        g = z["dec_global"]
        exclue = EXCLUES.get((carte, sol, p, k), "")
        if exclue or g.sum() < MIN_PHOTONS:
            print(f"  {carte} {sol} {p} ch{k} écartée : {exclue or 'trop peu de photons'}")
            return
        tau, s, b, dev, err = M.ajuster(g)
        a = acq.get((serie, carte, sol, p, k), {})
        base = {"serie": serie, "carte": carte, "solution": sol, "puissance": p, "canal": k}
        lignes.append({**base, "tranche": "global", "mu_enr": a.get("mu_enr_image", ""),
                       "photons": int(g.sum()), "tau_ns": tau, "err_tau_ns": err,
                       "fond": b, "decalage": s,
                       "deviance_reduite": dev})
        if "dec_tranches" in z:
            for j, (mu, d) in enumerate(zip(z["mu_tranches"], z["dec_tranches"])):
                t, _, bt, dv, e = M.ajuster(d, s=s, tau0=tau)
                lignes.append({**base, "tranche": j, "mu_enr": mu, "photons": int(d.sum()),
                               "tau_ns": t, "err_tau_ns": e, "fond": bt, "decalage": s, "deviance_reduite": dv})
        print(f"  {carte:10s} {sol:4s} P{p:<3} ch{k} : τ = {tau:.3f} ± {err:.3f} ns, "
              f"fond {b:.3f}, déviance réduite {dev:.2f}")

    for carte, sol, _, puissances in SERIES:
        for p in puissances:
            for k in CANAUX:
                ajouter("spc", carte, sol, p, k, np.load(MESURES / f"{cle(carte, sol, p)}_ch{k}.npz"))
    carte, sol, _, puissances = SERIE_150N_SDT
    for p in puissances:
        for k in CANAUX:
            f = MESURES / f"150N-sdt_{sol}_{p}_ch{k}.npz"
            if f.exists():
                ajouter("sdt", carte, sol, p, k, np.load(f))

    with open(RESULTATS / "durees_de_vie.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lignes[0]))
        w.writeheader()
        for l in lignes:
            w.writerow({c: (f"{v:.6g}" if isinstance(v, float) else v) for c, v in l.items()})
    print(f"  -> durees_de_vie.csv ({len(lignes)} lignes)")

    # --- temps morts ---------------------------------------------------------------------
    tm = []
    for k in CANAUX:
        h = sum(np.load(MESURES / f"{cle('SPC-QC-104', sol, p)}_ch{k}.npz")["hist_int"]
                for _, sol, _, ps in SERIES[:2] for p in ps
                if ("SPC-QC-104", sol, p, k) not in EXCLUES)
        d = demi_montee(h, BORDS_INTERVALLES_QC)
        tm.append({"carte": "SPC-QC-104", "canal": k, "solution": "toutes", "puissance": "toutes",
                   "methode": "demi-montée des intervalles entre photons", "temps_mort_ns": d})
        print(f"  QC-104 ch{k} : temps mort {d:.2f} ns (intervalles)")
    for (serie, carte, sol, p, k), a in sorted(acq.items(), key=str):
        if carte != "SPC-150N" or not a.get("taux_cfd_MHz"):
            continue
        c, t = float(a["taux_cfd_MHz"]), float(a["taux_stopinfo_MHz"])
        if c < 1.0:            # sous 1 MHz, les pertes (< 10 %) sont trop petites pour D
            continue
        tm.append({"carte": carte, "canal": k, "solution": sol, "puissance": p,
                   "methode": f"1/r_TAC − 1/r_CFD ({serie})", "temps_mort_ns": (1 / t - 1 / c) * 1e3})
    with open(RESULTATS / "temps_morts.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(tm[0]))
        w.writeheader()
        for l in tm:
            w.writerow({c: (f"{v:.4g}" if isinstance(v, float) else v) for c, v in l.items()})
    for sol in ("low", "high"):
        v = [l["temps_mort_ns"] for l in tm if l["carte"] == "SPC-150N" and l["solution"] == sol]
        print(f"  150N {sol} : temps mort effectif {np.mean(v):.1f} ± {np.std(v):.1f} ns ({len(v)} mesures)")
    print("  -> temps_morts.csv")


if __name__ == "__main__":
    main()
