"""Mesures communes aux étapes 6 et 8 : sélection, comptes corrigés du temps mort, lectures
de la durée de vie CFP.

**Temps mort.** La QC-104 a un temps mort de 7,5 ns par entrée, plus court que la période
(12,46 ns) mais plus long que le gros d'un déclin : au plus un photon par impulsion et par
entrée est enregistré (pileup/, 2026-10-23). Pour un pixel balayé E = passages × 4,16 fois
(impulsions), n photons enregistrés valent −E·ln(1 − n/E) photons arrivés (Poisson). La
correction ne dépend que du taux local : elle se calcule pixel par pixel, cycle par cycle,
puis se somme sur le ROI. μ = n/E est le taux enregistré par impulsion ; sa moyenne pondérée
par les photons dit à quel point un ROI est exposé à l'empilement.
"""
import json

import numpy as np

import config as C
import recalage


def selection(D, nom, n):
    """Labels des ROIs gardés (tous si aucun tri)."""
    p = D["selection"] / f"{nom}_selection.json"
    if not p.exists():
        return np.arange(1, n + 1)
    return np.array(json.loads(p.read_text())["kept"], int)


def comptes_cycles(z, lab, n_lab, decalages):
    """Par cycle et par label : Σn, Σn_corrigé, Σn²/E pour chaque canal.

    z : projections/<acq>.npz (images (cycles, 2, 512, 1024), passages (cycles, 512)) ;
    lab : étiquettes (0..n_lab) dans le repère de la somme ; decalages : (cycles, 2) entiers,
    ceux de l'étape 4. Rend (cycles, n_lab + 1, 2, 3)."""
    ims, pas = z["images"], z["passages"]
    out = np.zeros((len(ims), n_lab + 1, 2, 3))
    for k in range(len(ims)):
        lk = recalage.deplacer(lab, decalages[k]) if tuple(decalages[k]) != (0, 0) else lab
        lk = lk.ravel().astype(np.int64)
        E = (pas[k].astype(float) * C.IMPULSIONS_PAR_PIXEL)[:, None] * np.ones((1, C.PIXELS_LIGNE))
        E = np.maximum(E.ravel(), 1e-9)
        for c in range(2):
            n = ims[k, c].ravel().astype(float)
            p = np.minimum(n / E, 0.999)
            corr = -E * np.log1p(-p)
            for j, v in enumerate((n, corr, n * n / E)):
                out[k, :, c, j] = np.bincount(lk, weights=v, minlength=n_lab + 1)
    return out


def lectures_cfp(M, y, s_fov, taus):
    """Lectures de la durée de vie sur un déclin CFP (256 canaux).

    tau_mono : mono-exponentielle, décalage libre (le modèle de FLIMApp) ;
    tau_m    : temps moyen − celui de l'IRF au décalage du FOV, sans modèle (fond compris) ;
    a1..a3   : amplitudes des trois populations à τ communes (étape 5), fond libre ;
    tau_3pop : leur durée de vie moyenne pondérée par l'intensité.
    """
    y = np.asarray(y, float)
    N = float(y[M.lo:M.hi].sum())
    out = dict(photons_fenetre=N)
    tau, s, b, dd = M.ajuster_mono(y)
    out.update(tau_mono=tau, s_mono=s, b_mono=b, dev_ddl_mono=dd)
    out["tau_m"] = M.temps_moyen(y, s_fov) if N >= C.MIN_PHOTONS_FIT else np.nan
    a, b3, dd3 = M.amplitudes(y, taus, s_fov)
    out.update({f"a{i + 1}": float(v) for i, v in enumerate(a)})
    out.update(b_3pop=b3, dev_ddl_3pop=dd3)
    if np.all(np.isfinite(a)):
        t = np.asarray(taus)
        out["tau_3pop"] = float((a * t ** 2).sum() / (a * t).sum())
    else:
        out["tau_3pop"] = np.nan
    return out
