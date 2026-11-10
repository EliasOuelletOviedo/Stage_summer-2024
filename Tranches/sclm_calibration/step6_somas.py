"""Étape 6 — mesures par soma dans les FOV des conditions fixes.

Usage :  python3 step6_somas.py <dossier>

Pour chaque ROI de chaque FOV (cycles gardés réunis) :
  - photons CFP et YFP, rapport YFP/CFP **brut** (photons enregistrés, sans fond retiré) ;
  - le même rapport corrigé du temps mort de la QC-104 (`mesures.comptes_cycles`) et le taux
    enregistré par impulsion μ de chaque canal ;
  - durée de vie CFP : mono-exponentielle (décalage libre), temps moyen τm, amplitudes des
    trois populations à τ communes et leur moyenne pondérée τ_3pop (étape 5) ;
  - stabilité : pente du rapport sur les 6 cycles (%/min).
Mêmes mesures pour le neuropile et le champ entier de chaque FOV.

Sorties : results/somas.csv (une ligne par ROI × FOV ; colonne `garde` = tri de l'étape 3),
results/somas_cycles.csv (rapport cycle par cycle).
"""
import json
import sys

import numpy as np
import pandas as pd

import config as C
from flim import Modele
from mesures import comptes_cycles, lectures_cfp, selection


def main():
    root = C.racine()
    D = C.dirs(root)
    p_irf = D["qc"] / "irf.json"
    if not p_irf.exists():
        sys.exit("qc/irf.json absent : lance d'abord step5_irf.py")
    irf = json.loads(p_irf.read_text())
    M = Modele(C.charger_irf(irf["voie_cfp"], root))
    taus = np.array(irf["taus_cfp"])
    A = [a for a in C.acquisitions(root) if a["type"] == "fov" and (D["mesures"] / f"{a['nom']}.npz").exists()]
    lignes, cyc_lignes = [], []
    for a in A:
        z = np.load(D["projections"] / f"{a['nom']}.npz")
        m = np.load(D["mesures"] / f"{a['nom']}.npz")
        n = int(m["n_rois"])
        g = m["garde"]
        garde = set(selection(D, a["nom"], n).tolist())
        s_fov = irf["fovs"][a["nom"]]["s"]
        cc = comptes_cycles(z, m["etiquettes"], n + 1, m["decalages"])[g]     # (k, n+2, 2, 3)
        dec = m["declins"][g]                                               # (k, n+2, 2, 256)
        aire = np.bincount(m["etiquettes"].ravel(), minlength=n + 2)
        t_min = (m["debut"][g] - m["debut"][g][0]) / 60
        # champ entier : déclins de l'étape 1, comptes sur toute l'image
        dec_champ = z["declin"][g]
        cc_champ = comptes_cycles(z, np.ones(m["etiquettes"].shape, np.int32), 1, m["decalages"])[g][:, 1]
        roi_ids = list(range(1, n + 2)) + ["champ"]
        for roi in roi_ids:
            if roi == "champ":
                d_k, c_k, typ, ar = dec_champ, cc_champ, "champ", int(np.prod(m["etiquettes"].shape))
            else:
                d_k, c_k = dec[:, roi], cc[:, roi]
                typ, ar = ("neuropile", int(aire[roi])) if roi == n + 1 else ("soma", int(aire[roi]))
            tot = d_k.sum(0).astype(float)                                   # (2, 256)
            n_c, n_y = tot[0].sum(), tot[1].sum()
            cs = c_k.sum(0)                                                  # (2, 3)
            r_k = d_k[:, 1].sum(1) / np.maximum(d_k[:, 0].sum(1), 1)
            pente = np.polyfit(t_min, r_k / np.mean(r_k), 1)[0] * 100 if len(r_k) > 2 and n_c > 0 else np.nan
            row = dict(acquisition=a["nom"], condition=a["condition"], cl_mM=a["cl_mM"], type=typ,
                       roi=roi if roi != "champ" else 0, garde=(typ == "soma" and roi in garde),
                       aire_px=ar, cycles=int(g.sum()), photons_cfp=n_c, photons_yfp=n_y,
                       R=n_y / max(n_c, 1), R_corrige=cs[1, 1] / max(cs[0, 1], 1e-9),
                       mu_cfp=cs[0, 2] / max(cs[0, 0], 1e-9), mu_yfp=cs[1, 2] / max(cs[1, 0], 1e-9),
                       pente_R_pct_min=pente)
            row.update(lectures_cfp(M, tot[0], s_fov, taus))
            lignes.append(row)
            if typ == "soma":
                for k, rk in enumerate(r_k):
                    cyc_lignes.append(dict(acquisition=a["nom"], roi=roi, cycle=k, t_min=t_min[k],
                                           photons_cfp=int(d_k[k, 0].sum()), photons_yfp=int(d_k[k, 1].sum()), R=rk))
        S = pd.DataFrame([l for l in lignes if l["acquisition"] == a["nom"] and l["type"] == "soma" and l["garde"]])
        print(f"  {a['nom']:>12} : {len(S)} somas gardés / {n} ; R médian {S.R.median():.3f}, "
              f"τ mono {S.tau_mono.median():.3f} ns, τm {S.tau_m.median():.3f} ns, "
              f"μ YFP médian {S.mu_yfp.median():.4f}")
    pd.DataFrame(lignes).to_csv(D["results"] / "somas.csv", index=False)
    pd.DataFrame(cyc_lignes).to_csv(D["results"] / "somas_cycles.csv", index=False)
    print("  -> results/somas.csv, results/somas_cycles.csv")


if __name__ == "__main__":
    main()
