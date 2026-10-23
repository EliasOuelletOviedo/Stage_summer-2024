"""Étape 1 — des photons aux déclins, rangés par taux local.

Pour chaque acquisition (carte, solution, puissance) et chaque canal :

- taux : photons détectés / durée du flux, et taux de StopInfo du .sdt (CFD = détecteur,
  TAC/ADC = conversions, pour la 150N) ;
- carte d'intensité, lissée (LISSAGE_PX), divisée par l'exposition de chaque ligne :
  taux local enregistré, en photons par impulsion laser (μ_enr) ;
- déclins (256 canaux sur une période) : global, et par tranche de taux local
  (N_TRANCHES tranches de même nombre de photons). La solution est uniforme : seul le
  taux change d'une tranche à l'autre (vignettage, bord sombre en bas des images
  fort Cl⁻ de la QC-104) ;
- intervalles entre photons successifs d'une même entrée (temps mort). QC-104 : le
  macrotemps (2,048 ns) dit à quelle impulsion appartient chaque photon, le microtemps
  donne l'intervalle exact. 150N : macrotemps de 25 ns seulement.

Sortie : resultats/mesures/<clé>.npz, resultats/acquisitions.csv.
Durée : ~1 min.
"""
import csv

import numpy as np
from scipy import ndimage as ndi

from config import (CARTES, LISSAGE_PX, MESURES, MIN_PHOTONS, N_CANAUX, N_TRANCHES,
                    PERIODE_NS, PIXEL_NS, RESULTATS, SERIE_150N_SDT, SERIES, cle,
                    dossier_donnees, preparer)
from lecture import info_mesure, lire_sdt, lire_spc, pixels

BORDS_INTERVALLES_QC = np.arange(0, 60.01, 0.1)       # ns
BORDS_INTERVALLES_150N = np.arange(0, 1000.1, 25.0)   # ns, un tic de macrotemps


def intervalles(f, carte):
    if carte == "SPC-QC-104":
        dT = np.diff(f.t_ns)
        du = np.diff(f.micro_ns)
        k = np.round((dT - du) / PERIODE_NS)        # impulsions entre les deux photons
        vrai = k * PERIODE_NS + du
        return np.histogram(vrai, BORDS_INTERVALLES_QC)[0], float(np.mean(k == 0))
    return np.histogram(np.diff(f.t_ns), BORDS_INTERVALLES_150N)[0], 0.0


def tranches(mu, canal):
    """Déclins par tranche de taux local : N_TRANCHES tranches de même nombre de photons."""
    ok = canal >= 0
    bords = np.quantile(mu[ok], np.linspace(0, 1, N_TRANCHES + 1))
    bords[-1] += 1e-12
    i = np.clip(np.searchsorted(bords, mu, side="right") - 1, 0, N_TRANCHES - 1)
    dec = np.zeros((N_TRANCHES, N_CANAUX))
    np.add.at(dec, (i[ok], canal[ok]), 1)
    mu_moy = np.array([mu[ok & (i == j)].mean() for j in range(N_TRANCHES)])
    return mu_moy, dec


def traiter(carte, sol, dossier, p, sdt):
    flux = lire_spc(carte, dossier, p)
    info = {k: info_mesure(sdt, 0 if carte == "SPC-QC-104" else k - 1) for k in (1, 2)} \
        if sdt is not None else {}
    lignes = []
    for k, f in flux.items():
        n = len(f.t_ns)
        val = f.canal >= 0
        duree = (f.t_ns[-1] - f.t_ns[0]) * 1e-9 if n > 1 else np.nan
        dec_global = np.bincount(f.canal[val], minlength=N_CANAUX).astype(float)
        hist_int, meme_impulsion = intervalles(f, carte)
        ligne = {"carte": carte, "solution": sol, "puissance": p, "canal": k,
                 "photons": n, "photons_valides": int(val.sum()), "duree_s": duree,
                 "taux_enr_MHz": n / duree / 1e6 if n > 1 else np.nan,
                 "meme_impulsion": meme_impulsion}
        if info:
            i = info[k]
            # QC-104 : StopInfo range le taux de l'entrée 1 dans « cfd », celui de l'entrée 2
            # dans « tac ». 150N : un module par canal, cfd = détecteur, tac = conversions.
            if carte == "SPC-QC-104":
                ligne["taux_stopinfo_MHz"] = (i["cfd"] if k == 1 else i["tac"]) / 1e6
            else:
                ligne["taux_stopinfo_MHz"] = i["tac"] / 1e6
                ligne["taux_cfd_MHz"] = i["cfd"] / 1e6
        out = {"dec_global": dec_global, "hist_int": hist_int}
        if n >= MIN_PHOTONS:
            y, x, expo = pixels(f)
            dans = y >= 0
            C = np.zeros((len(expo), 1024))
            np.add.at(C, (y[dans], x[dans]), 1)
            E = np.repeat(expo[:, None], 1024, axis=1)
            Cs = ndi.gaussian_filter(C, LISSAGE_PX, mode="nearest")
            Es = ndi.gaussian_filter(E, LISSAGE_PX, mode="nearest")
            mu_carte = Cs / np.maximum(Es, 1e-9) / PIXEL_NS * PERIODE_NS
            mu = mu_carte[y[dans], x[dans]]
            mu_moy, dec = tranches(mu, f.canal[dans])
            ligne.update({"n_trames": len(f.trames_ns), "photons_image": int(dans.sum()),
                          "mu_enr_image": float(mu.mean()),
                          "mu_enr_min": float(mu_moy[0]), "mu_enr_max": float(mu_moy[-1])})
            out.update({"mu_tranches": mu_moy, "dec_tranches": dec,
                        "carte_mu": mu_carte[::8, ::8]})
        np.savez_compressed(MESURES / f"{cle(carte, sol, p)}_ch{k}.npz", **out)
        lignes.append(ligne)
        print(f"  {cle(carte, sol, p)} ch{k} : {n:>9} photons, {ligne['taux_enr_MHz']:.2f} MHz"
              + (f", μ_enr image {ligne['mu_enr_image']:.4f}" if "mu_enr_image" in ligne else ""))
    return lignes


def serie_sdt(D):
    """1re série de la 150N : déclins globaux et taux du .sdt seulement."""
    carte, sol, sous, puissances = SERIE_150N_SDT
    lignes = []
    for p in puissances:
        s = lire_sdt(D / sous / f"{p}.sdt")
        for b, d in enumerate(s.data):
            i = info_mesure(s, b)
            m = s.measure_info[b]
            m = m.flatten()[0] if hasattr(m, "flatten") else m
            k = int(m["FCSInfo"]["mod"]) + 1            # module 0 -> canal 1
            dec = d.sum(axis=(0, 1)).astype(float)
            np.savez_compressed(MESURES / f"150N-sdt_{sol}_{p}_ch{k}.npz", dec_global=dec)
            lignes.append({"carte": carte, "solution": sol, "puissance": p, "canal": k,
                           "serie": "sdt", "photons": int(dec.sum()), "duree_s": i["duree_s"],
                           "taux_stopinfo_MHz": i["tac"] / 1e6, "taux_cfd_MHz": i["cfd"] / 1e6})
            print(f"  150N-sdt {sol} {p} ch{k} : {int(dec.sum())} photons, "
                  f"CFD {i['cfd'] / 1e6:.2f} MHz, TAC {i['tac'] / 1e6:.2f} MHz")
    return lignes


def main():
    D = dossier_donnees()
    preparer()
    lignes = []
    for carte, sol, sous, puissances in SERIES:
        assert carte in CARTES
        for p in puissances:
            sdt_path = D / sous / f"{p}.sdt"
            sdt = lire_sdt(sdt_path) if sdt_path.exists() else None
            for l in traiter(carte, sol, D / sous, p, sdt):
                l["serie"] = "spc"
                lignes.append(l)
    lignes += serie_sdt(D)
    cles = []
    for l in lignes:
        cles += [c for c in l if c not in cles]
    with open(RESULTATS / "acquisitions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cles)
        w.writeheader()
        for l in lignes:
            w.writerow({c: (f"{v:.6g}" if isinstance(v, float) else v) for c, v in l.items()})
    print(f"  -> {RESULTATS / 'acquisitions.csv'}")


if __name__ == "__main__":
    main()
