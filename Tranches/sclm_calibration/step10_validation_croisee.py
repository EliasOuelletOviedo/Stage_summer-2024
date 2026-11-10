"""Étape 10 — validation croisée des lectures de [Cl⁻] : calibrer sans un FOV, lire ce FOV.

Usage :  python3 step10_validation_croisee.py <dossier>              (~10 min) tirage aléatoire des ROIs
         python3 step10_validation_croisee.py <dossier> --par-fov    un FOV de chaque condition par pli

**Tirage aléatoire (par défaut, sorties validation_aleatoire_*)** : dans chaque [Cl⁻] imposée, les ROIs sont
mélangés (graine fixe) et répartis en 5 plis ; le pli k calibre sur les 4/5 des ROIs et lit le 1/5 restant. Les mêmes
champs servent à calibrer et à lire : on mesure la justesse sur de nouvelles cellules des mêmes champs, pas sur un
nouveau champ. Critère : médiane des ROIs lus de chaque pli contre la [Cl⁻] imposée (biais, erreur quadratique sur
les 5 plis). **Par FOV (--par-fov, sorties validation_croisee_*)** : ci-dessous.

Aux étapes 7 et 9, les somas relus à 0, 20 et 131 mM sont ceux qui ont servi à calibrer : leur relecture dit la
dispersion, pas l'exactitude (le τ mono retombe même exactement sur les concentrations imposées, par construction).
Ici, 5 plis : le pli k retire le k-ième FOV de chaque [Cl⁻] imposée, recalibre sur les 4 autres FOV de chaque
condition, puis lit les 3 FOV retirés. Chaque FOV est lu une fois, par une calibration qui ne l'a pas vu. Le basal,
jamais utilisé pour calibrer, est lu par chaque pli (sensibilité de sa lecture à la calibration).

Lectures (mêmes définitions qu'aux étapes 7 et 9, sélection de l'étape 3) :
  - isotherme sur les médianes des somas (étape 7) : rapport brut, rapport corrigé du temps mort, τ mono, τm,
    τ des 3 populations ;
  - modèle FLIM à deux canaux recalibré à chaque pli (étape 9), 2 populations en FRET et 1 seule : temps de vie
    unique, ratio (η = 1), FLIM du système, hybride (η partagé par le FOV). L'IRF paramétrique est ajustée une fois
    (elle ne dépend pas des somas).

Critères, fixés avant de regarder les résultats : erreur de la médiane de chaque FOV lu par rapport à la [Cl⁻]
imposée (biais moyen et erreur quadratique moyenne sur les 5 FOV de la condition) ; dispersion entre somas
(1,48·MAD) ; part de lectures hors courbe (+∞, comptées dans les médianes). Les [Cl⁻] négatives sont gardées.

Sorties : results/validation_croisee_somas.csv, results/validation_croisee_fovs.csv,
results/validation_croisee_resume.csv, results/validation_croisee_basal.csv, results/validation_croisee.png
"""
import re
import sys
import time

import numpy as np
import pandas as pd

import config as C
import deux_canaux as DC
import figures as G
from figures import plt
from step7_calibration import ajuster, inverser
from step9_flim_deux_canaux import C_POSITIFS, CONDITIONS, FENETRE, LN_ETA, TAU_GRILLE, somas

LECT_ISO = ("R", "R_corrige", "tau_mono", "tau_m", "tau_3pop")
NOMS = {"R": "rapport brut (isotherme)", "R_corrige": "rapport corrigé (isotherme)", "tau_mono": "τ mono (isotherme)",
        "tau_m": "τm (isotherme)", "tau_3pop": "τ 3 populations (isotherme)"}
for nt in (2, 1):
    for c, n in (("tau_unique", "temps de vie unique"), ("ratio", "ratio, η = 1"), ("flim", "FLIM du système"),
                 ("hybride", "hybride, η par FOV")):
        NOMS[f"{c}_{nt}p"] = f"{n} ({nt} pop.)"


def num_fov(nom):
    return int(re.search(r"fov(\d+)$", nom).group(1))


def modele(root, B, Y, conc, irfp, nt):
    """Calibration conjointe et lecteur (étape 9) sur les déclins réunis B, Y de chaque [Cl⁻]."""
    d = DC.donnees(root, B, Y, conc, FENETRE)
    d["irfp"] = irfp
    R = DC.calibrer(d, irfp, nt)
    iso = DC.isothermes(R, conc)
    S = DC.senseur(d, R, iso)
    cg = np.r_[np.linspace(DC.c_min_modele(iso), 0, 200, endpoint=False), C_POSITIFS]
    L = DC.LecteurSoma(S, cg, TAU_GRILLE)
    w = d["win"][1]
    tau_cal = np.array([L.tau_unique(b[w]) for b in d["B"]])
    return L, w, tau_cal, R["dev"] / R["dof"]


def c_depuis_tau(t, tau_cal, conc):
    """Courbe τ(c) des déclins de calibration (étape 9) : interpolation monotone, droites au-delà, pas de borne."""
    from scipy.interpolate import PchipInterpolator
    if not np.all(np.diff(tau_cal) > 0):
        return np.nan
    if t < tau_cal[0]:
        return conc[0] + (t - tau_cal[0]) * (conc[1] - conc[0]) / (tau_cal[1] - tau_cal[0])
    if t > tau_cal[-1]:
        return conc[-1] + (t - tau_cal[-1]) * (conc[-1] - conc[-2]) / (tau_cal[-1] - tau_cal[-2])
    return float(PchipInterpolator(tau_cal, conc)(t))


def lire_fov(L, w, tau_cal, conc, labs, Bs, Ys):
    """Lectures de l'étape 9 pour les somas d'un FOV : {label: {lecture: c}}."""
    haut = L.c[-1] * (1 - 1e-9)
    out, profs = {}, []
    for lab, b, y in zip(labs, Bs, Ys):
        b, y = b[w], y[w]
        if b.sum() < C.MIN_PHOTONS_FIT or y.sum() < C.MIN_PHOTONS_FIT:
            continue
        pr = L.lire(b, y)
        quant = -0.5 * (pr["ln_q"] - L.lnR1) ** 2 / pr["v"]
        cf, cr = DC.argmax_affine(pr["formes"], L.c), DC.argmax_affine(quant, L.c)
        out[int(lab)] = dict(tau_unique=c_depuis_tau(L.tau_unique(b), tau_cal, conc),
                             ratio=np.inf if cr >= haut else cr, flim=np.inf if cf >= haut else cf)
        profs.append((int(lab), pr))
    if profs:
        prof_fov = lambda g: sum(np.array([L.hybride(pr, x).max() for x in g]) for _, pr in profs)
        le0 = LN_ETA[int(np.argmax(prof_fov(LN_ETA)))]
        fine = le0 + np.linspace(-0.05, 0.05, 101)
        le = DC.argmax_affine(prof_fov(fine), fine)
        for lab, pr in profs:
            ch = DC.argmax_affine(L.hybride(pr, le), L.c)
            out[lab]["hybride"] = np.inf if ch >= haut else ch
    return out


def mad(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v))) if len(v) else np.nan


PAR_FOV, PREFIXE = False, "validation_aleatoire"
N_PLIS, GRAINE = 5, 20261110


def main():
    global PAR_FOV, PREFIXE
    if "--par-fov" in sys.argv[1:]:
        PAR_FOV, PREFIXE = True, "validation_croisee"
    root = C.racine()
    D = C.dirs(root)
    t0 = time.time()
    A = [a for a in C.acquisitions(root) if a["type"] == "fov" and (D["mesures"] / f"{a['nom']}.npz").exists()]
    cond_de = {a["nom"]: a["condition"] for a in A}
    so = somas(D, A)
    T6 = pd.read_csv(D["results"] / "somas.csv")
    T6 = T6[(T6.type == "soma") & T6.garde]
    conc = [C.CL_MM[c] for c in CONDITIONS]
    imposes = [n for n in so if cond_de[n] in CONDITIONS]
    if PAR_FOV:            # pli k : le k-ième FOV de chaque condition
        plis = {k: {(n, int(l)) for n in imposes if num_fov(n) == k for l in so[n][0]}
                for k in sorted({num_fov(n) for n in imposes})}
    else:                  # tirage aléatoire des ROIs, stratifié par condition
        rng = np.random.default_rng(GRAINE)
        plis = {k: set() for k in range(1, N_PLIS + 1)}
        for c in CONDITIONS:
            rois = [(n, int(l)) for n in imposes if cond_de[n] == c for l in so[n][0]]
            for i, j in enumerate(rng.permutation(len(rois))):
                plis[i % N_PLIS + 1].add(rois[j])
    tous = {(n, int(l)) for n in imposes for l in so[n][0]}
    # IRF paramétrique : une fois, sur les déclins réunis de toutes les conditions imposées
    B0 = [sum(so[n][1].sum(0) for n in so if cond_de[n] == c) for c in CONDITIONS]
    Y0 = [sum(so[n][2].sum(0) for n in so if cond_de[n] == c) for c in CONDITIONS]
    irfp = DC.ajuster_irf(DC.donnees(root, B0, Y0, conc, FENETRE))
    basal = [n for n in so if cond_de[n] == "basal"]
    cle6 = list(zip(T6.acquisition, T6.roi.astype(int)))

    rows, cal_rows = [], []
    for k, test in plis.items():
        train = tous - test
        # --- isothermes sur les médianes (étape 7)
        dans_train = np.array([c in train for c in cle6])
        lu = np.array([c in test or cond_de.get(c[0]) == "basal" for c in cle6])
        tr6, s6 = T6[dans_train], T6[lu]
        for col in LECT_ISO:
            f = ajuster(conc, [tr6[tr6.condition == c][col].median() for c in CONDITIONS])
            v = inverser(s6[col].values, *f) if f else np.full(len(s6), np.nan)
            cal_rows.append(dict(pli=k, lecture=col, K_mM=f[2] if f else np.nan))
            for (_, r), c in zip(s6.iterrows(), v):
                rows.append(dict(pli=k, acquisition=r.acquisition, condition=r.condition, roi=int(r.roi),
                                 lecture=col, c_mM=float(c)))
        # --- modèle à deux canaux, recalibré sur les ROIs d'entraînement
        def reunir(i):
            return [sum(so[n][i][j] for n in imposes if cond_de[n] == c for j, l in enumerate(so[n][0])
                        if (n, int(l)) in train) for c in CONDITIONS]
        B, Y = reunir(1), reunir(2)
        for nt in (2, 1):
            L, w, tau_cal, dev = modele(root, B, Y, conc, irfp, nt)
            cal_rows.append(dict(pli=k, lecture=f"modele_{nt}p", deviance_ddl=dev))
            for n in imposes + basal:
                labs, Bs, Ys = so[n]
                garde = np.array([cond_de[n] == "basal" or (n, int(l)) in test for l in labs], bool)
                if not garde.any():
                    continue
                for lab, lec in lire_fov(L, w, tau_cal, conc, labs[garde], Bs[garde], Ys[garde]).items():
                    for nom, c in lec.items():
                        rows.append(dict(pli=k, acquisition=n, condition=cond_de[n], roi=lab, lecture=f"{nom}_{nt}p",
                                         c_mM=float(c)))
        print(f"  pli {k} : {len(test)} ROIs lus hors calibration ({time.time() - t0:.0f} s)")
    V = pd.DataFrame(rows)
    V.to_csv(D["results"] / f"{PREFIXE}_somas.csv", index=False)
    pd.DataFrame(cal_rows).to_csv(D["results"] / f"{PREFIXE}_calibrations.csv", index=False)

    # --- groupes lus hors calibration (un FOV, ou les ROIs d'un pli) : médiane contre la [Cl⁻] imposée
    imp = V[V.condition.isin(CONDITIONS)].copy()
    imp["groupe"] = imp.acquisition if PAR_FOV else "pli" + imp.pli.astype(str)
    F = (imp.groupby(["lecture", "condition", "groupe"])
         .agg(mediane=("c_mM", "median"), dispersion=("c_mM", mad),
              hors_courbe_pct=("c_mM", lambda v: 100 * float(np.isinf(v).mean())), somas=("c_mM", "size")).reset_index())
    F["cl_mM"] = F.condition.map(C.CL_MM)
    F["erreur"] = F.mediane - F.cl_mM
    F.to_csv(D["results"] / f"{PREFIXE}_groupes.csv", index=False)
    res = []
    for (lec, cond), g in F.groupby(["lecture", "condition"]):
        e = g.erreur.replace([np.inf, -np.inf], np.nan)
        res.append(dict(lecture=lec, nom=NOMS.get(lec, lec), condition=cond, cl_mM=C.CL_MM[cond],
                        biais=float(e.mean()) if e.notna().all() else np.inf,
                        rmse=float(np.sqrt((e ** 2).mean())) if e.notna().all() else np.inf,
                        dispersion_somas=float(g.dispersion.median()), hors_courbe_pct=float(g.hors_courbe_pct.mean())))
    RS = pd.DataFrame(res)
    RS.to_csv(D["results"] / f"{PREFIXE}_resume.csv", index=False)
    # --- basal : lecture par chaque pli
    Bv = V[V.condition == "basal"].groupby(["lecture", "pli"]).c_mM.median().reset_index()
    BS = Bv.groupby("lecture").c_mM.agg(["median", "min", "max"]).reset_index()
    BS.to_csv(D["results"] / f"{PREFIXE}_basal.csv", index=False)
    for q in ("rmse", "biais", "dispersion_somas", "hors_courbe_pct"):
        print(f"\n{q}")
        print(RS.pivot(index="lecture", columns="condition", values=q).round(1).to_string())
    print("\nbasal (médiane des somas, par pli) :")
    print(BS.round(2).to_string(index=False))
    figure(D, F, RS, BS)
    print(f"  -> results/{PREFIXE}_* ({time.time() - t0:.0f} s)")


def figure(D, F, RS, BS):
    lects = [l for l in NOMS if l in set(F.lecture)]
    fig, axes = plt.subplots(1, 4, figsize=(21, 6.4), gridspec_kw=dict(width_ratios=[1, 1, 1, 0.9]))
    for ax, cond in zip(axes, CONDITIONS):
        G.style(ax)
        for j, l in enumerate(lects):
            g = F[(F.lecture == l) & (F.condition == cond)]
            v = g.mediane.replace(np.inf, np.nan)
            ax.plot(np.full(len(v), j) + np.linspace(-0.2, 0.2, len(v)), v, "o", ms=4,
                    color=C.COULEUR_CONDITION[cond], alpha=0.8)
            n_inf = int(np.isinf(g.mediane).sum())
            if n_inf:
                ax.text(j, 1.0, f"{n_inf}∞", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7,
                        color=G.PALETTE[4])
        ax.axhline(C.CL_MM[cond], color=G.INK, ls="--", lw=1)
        ax.set_xticks(range(len(lects)))
        ax.set_xticklabels([NOMS[l] for l in lects], rotation=70, ha="right", fontsize=7.5)
        lo, hi = {"0 mM": (-15, 25), "20 mM": (-10, 90), "131 mM": (0, 400)}[cond]
        ax.set_ylim(lo, hi)
        quoi = "chaque FOV" if PAR_FOV else "chaque pli (ROIs tirés au hasard)"
        G.titre(ax, f"{cond} — médiane de {quoi}, lue hors calibration",
                "tireté : [Cl⁻] imposée ; « n∞ » : médiane hors courbe", fontsize=10)
    axes[0].set_ylabel("[Cl⁻] lue (mM)", color=G.INK)
    ax = axes[3]
    G.style(ax)
    b = BS.set_index("lecture").loc[lects]
    y = np.arange(len(lects))
    ax.errorbar(b["median"], y, xerr=[b["median"] - b["min"], b["max"] - b["median"]], fmt="o",
                color=C.COULEUR_CONDITION["basal"], capsize=3)
    ax.set_yticks(y)
    ax.set_yticklabels([NOMS[l] for l in lects], fontsize=7.5)
    ax.invert_yaxis()
    ax.axvline(0, color=G.TRACE)
    ax.set_xlabel("[Cl⁻]i basale (mM) : médiane des 5 plis, barres min–max", color=G.INK)
    G.titre(ax, "Basal lu par chaque pli", "sensibilité à la calibration", fontsize=10)
    fig.tight_layout()
    G.sauver(fig, D["results"] / f"{PREFIXE}.png")


if __name__ == "__main__":
    main()
