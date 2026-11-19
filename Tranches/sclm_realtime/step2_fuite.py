"""Étape 2 — ce que le laser à 1064 nm ajoute aux photons, et ce qui en passe dans la mesure.

Usage :  python3 step2_fuite.py <dossier>          (quelques secondes ; après l'étape 1)

Pendant une stim, le PWM alterne les échantillons de 100 µs : phase opto (1064 nm à sa
commande, code 0, photons jetés par FLIMApp) puis phase d'acquisition (1064 nm éteint,
code 1, photons analysés). Hors stim, le même code alterne, laser opto éteint. Le laser à
850 nm, lui, ne s'arrête jamais (ni pendant les pauses de 30 ms entre les passes).

Trois mesures, pour chaque canal :
  - excès dans les phases opto : photons opto / photons des pauses (aucun laser opto
    pendant les pauses : la référence suit la fluorescence au même moment, y compris la
    réponse au Cl⁻ et le blanchiment), stim contre les 40 passes avant et après ;
  - ce qui passe dans les phases d'acquisition : profil des photons dans la période du PWM
    (tranches de 1 µs depuis le M0 de la passe), chaque profil rapporté à son plateau
    d'acquisition (µs 30 à 99 de la période) ; l'excès des premières µs d'acquisition en
    stim est la fuite du 1064 nm dans la mesure ;
  - forme des photons ajoutés : déclin de l'excès opto (stim − référence), décomposé en un
    fond plat (excitation non synchrone du 850 nm) et la forme du déclin à 850 nm.

Sorties : results/fuite.csv, results/fuite_profil.csv, results/fuite.png
"""
import numpy as np
import pandas as pd

import config as C
import figures as G
import spc as S
from figures import plt

VOISINS = 40            # passes de référence avant et après chaque stim
PLATEAU = (30, 100)     # µs de la période du PWM où l'acquisition est sans fuite


def blocs_stim(stim):
    d = np.diff(np.r_[0, stim.astype(int), 0])
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def phases(prof, n_phase):
    """Quelles tranches de la période sont l'acquisition (code ≠ 0) : d'après le profil hors stim."""
    acq = prof[0, 0] > prof[0, 1]
    return acq


def main():
    root = C.racine()
    D = C.dirs(root)
    P = C.charger_passes(D)
    stim, t = P["stim"], P["t_s"]
    n_phase = int(P["n_phase"])
    us = (np.arange(n_phase) + 0.5) * 2 * P["tics_echantillon"] * P["tic_s"] * 1e6 / n_phase
    blocs = blocs_stim(stim)
    lignes, profils, excedents = [], {}, {}
    for nom in ("yfp", "cfp"):
        H = P[f"H_{nom}"].astype(float)
        nA, nO, nP = H[:, S.ACQ].sum(1), H[:, S.OPTO].sum(1), H[:, S.PAUSE].sum(1)
        nP = np.where(nP > 0, nP, np.nan)                   # la dernière passe n'a pas de pause
        q_o, q_a = nO / nP, nA / nP
        ex_dec = np.zeros(H.shape[2])
        ref_dec = np.zeros(H.shape[2])
        for k, (a, b) in enumerate(blocs, start=1):
            ref = np.r_[max(a - VOISINS, 0):a, b:min(b + VOISINS, len(t))]
            ref = ref[~stim[ref] & np.isfinite(nP[ref])]
            r_o, r_a = q_o[ref].mean(), q_a[ref].mean()
            e_o = q_o[a:b] / r_o - 1
            e_a = q_a[a:b] / r_a - 1
            se_o = np.sqrt(e_o.var(ddof=1) / len(e_o) + (q_o[ref] / r_o).var(ddof=1) / len(ref))
            se_a = np.sqrt(e_a.var(ddof=1) / len(e_a) + (q_a[ref] / r_a).var(ddof=1) / len(ref))
            pente = np.polyfit(np.arange(b - a) * np.median(np.diff(t)), e_o, 1)[0]
            lignes.append(dict(canal=nom, stim=k, passes=b - a, exces_opto_pct=100 * e_o.mean(), se_opto_pct=100 * se_o,
                               pente_exces_opto_pct_par_s=100 * pente,
                               exces_acq_pct=100 * e_a.mean(), se_acq_pct=100 * se_a,
                               photons_ajoutes_par_passe=float((e_o.mean()) * r_o * nP[a:b].mean())))
            # déclin des photons ajoutés : opto en stim − opto de référence ramené aux pauses
            ex_dec += H[a:b, S.OPTO].sum(0) - H[ref, S.OPTO].sum(0) * np.nansum(nP[a:b]) / np.nansum(nP[ref])
            ref_dec += H[ref, S.OPTO].sum(0) * np.nansum(nP[a:b]) / np.nansum(nP[ref])
        excedents[nom] = (ex_dec, ref_dec)
        # profil dans la période du PWM, par passe, chaque profil sur son plateau d'acquisition
        prof = P[f"prof_{nom}"].astype(float)
        acq = phases(prof, n_phase)
        tot = prof.sum(1)                                     # [stim, phase] : codes réunis
        n_s, n_h = stim.sum(), (~stim).sum()
        plat = slice(*PLATEAU)
        norm_s, norm_h = tot[1, plat].mean(), tot[0, plat].mean()
        exces = tot[1] / norm_s - tot[0] / norm_h                # en part du débit d'acquisition, par µs
        profils[nom] = dict(us=us, stim=tot[1] / n_s, hors=tot[0] / n_h, exces=exces, acq=acq,
                            code0_stim=prof[1, 1] / n_s, code0_hors=prof[0, 1] / n_h)
        # fuite dans l'acquisition : excès sur les tranches d'acquisition hors plateau
        hors_plateau = acq.copy()
        hors_plateau[plat] = False
        fuite = exces[hors_plateau].sum() / acq.sum()
        # retard d'allumage : première tranche opto où l'excès dépasse la moitié du plateau opto
        opto = ~acq
        e_opto = exces[opto]
        mi = 0.5 * np.median(e_opto[len(e_opto) // 2:])
        i0 = np.flatnonzero(opto)[0]
        retard = (np.flatnonzero((exces > mi) & opto)[0] - i0) if mi > 0.005 else np.nan
        k_acq0 = np.flatnonzero(acq)[0]
        lignes.append(dict(canal=nom, stim=0, passes=int(n_s),
                           exces_opto_pct=100 * e_opto.mean(), exces_acq_pct=100 * fuite,
                           retard_allumage_us=float(retard) * (us[1] - us[0]) if np.isfinite(retard) else np.nan,
                           fuite_premiere_us_pct=100 * exces[k_acq0] / acq.sum()))
    L = pd.DataFrame(lignes)
    # forme des photons ajoutés : fond plat + forme à 850 nm (moindres carrés pondérés, fenêtre utile)
    lo, hi = C.FENETRE
    formes = []
    for nom in ("yfp", "cfp"):
        ex, ref = excedents[nom]
        X = np.c_[np.ones(hi - lo) / (hi - lo), ref[lo:hi] / ref[lo:hi].sum()]
        w = 1 / np.maximum(ref[lo:hi] + np.abs(ex[lo:hi]), 1)
        A = (X * w[:, None]).T @ X
        coef = np.linalg.solve(A, (X * w[:, None]).T @ ex[lo:hi])
        cov = np.linalg.inv(A)
        formes.append(dict(canal=nom, exces_total=float(ex[lo:hi].sum()), plat=float(coef[0]), forme_850=float(coef[1]),
                           part_plate=float(coef[0] / coef.sum()), se_part_plate=float(np.sqrt(cov[0, 0]) / abs(coef.sum()))))
    F = pd.DataFrame(formes)
    L.to_csv(D["results"] / "fuite.csv", index=False)
    pd.DataFrame({f"{nom}_{k}": v for nom, d in profils.items() for k, v in d.items() if k != "us"} | {"us": us}) \
        .to_csv(D["results"] / "fuite_profil.csv", index=False)
    pd.set_option("display.width", 220)
    print(L.round(3).to_string(index=False))
    print(F.round(4).to_string(index=False))
    figure(D, P, profils, excedents, L, F, blocs)


def figure(D, P, profils, excedents, L, F, blocs):
    fig, axes = plt.subplots(1, 4, figsize=(22, 4.8))
    for ax in axes:
        G.style(ax)
    t, stim = P["t_s"], P["stim"]
    ax = axes[0]
    a, b = blocs[0]
    sl = slice(max(a - 60, 0), b + 60)
    for nom in ("yfp", "cfp"):
        H = P[f"H_{nom}"].astype(float)
        nP = H[:, S.PAUSE].sum(1)
        nP = np.where(nP > 0, nP, np.nan)
        q = H[:, S.OPTO].sum(1) / nP
        ref = np.nanmedian(q[~stim])
        ax.plot(t[sl], 100 * (q[sl] / ref - 1), color=G.COULEUR_CANAL[nom], lw=1, label=f"{C.NOMS[nom]}, phases opto")
        qa = H[:, S.ACQ].sum(1) / nP
        ax.plot(t[sl], 100 * (qa[sl] / np.nanmedian(qa[~stim]) - 1), color=G.COULEUR_CANAL[nom], lw=0.8, alpha=0.45,
                ls=":", label=f"{C.NOMS[nom]}, phases d'acquisition")
    G.bandes(ax, [(t[a], t[b - 1])])
    ax.set_xlabel("temps (s)", color=G.INK)
    ax.set_ylabel("photons / photons des pauses, écart (%)", color=G.INK)
    ax.legend(frameon=False, fontsize=8)
    G.titre(ax, "Stim 1 : photons ajoutés, passe par passe", "référence : les pauses (850 nm seul), même passe", fontsize=10.5)

    ax = axes[1]
    for nom in ("yfp", "cfp"):
        p = profils[nom]
        ax.plot(p["us"], 100 * p["exces"], color=G.COULEUR_CANAL[nom], lw=1.2, label=C.NOMS[nom])
    o = np.flatnonzero(~profils["yfp"]["acq"])
    ax.axvspan(p["us"][o[0]] - 0.5, p["us"][o[-1]] + 0.5, color=G.STIM, alpha=0.3, lw=0)
    ax.text(p["us"][o[0]] + 2, 0.92, "phase opto (code 0)", fontsize=8.5, color=G.INK_2, transform=ax.get_xaxis_transform())
    ax.set_xlabel("temps dans la période du PWM (µs depuis le début de la passe, modulo 200)", color=G.INK)
    ax.set_ylabel("excès stim − hors stim (% du débit d'acquisition)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Dans la période du PWM", "chaque profil rapporté à son plateau d'acquisition (µs 30–99)", fontsize=10.5)

    ax = axes[2]
    tc = C.temps_canaux()
    lo, hi = C.FENETRE
    ex, ref = excedents["yfp"]
    ax.plot(tc[lo:hi], ex[lo:hi] / 4, color=G.COULEUR_CANAL["yfp"], lw=1, label="photons YFP ajoutés (phases opto)")
    ax.plot(tc[lo:hi], ref[lo:hi] / 4 * (ex[lo:hi].sum() / ref[lo:hi].sum()), color=G.INK_2, lw=1, ls="--",
            label="déclin à 850 nm, même nombre de photons")
    ax.set_yscale("log")
    ax.set_xlabel("temps (ns)", color=G.INK)
    ax.set_ylabel("photons (4 stims, par canal)", color=G.INK)
    ax.legend(frameon=False, fontsize=8)
    fy = F[F.canal == "yfp"].iloc[0]
    G.titre(ax, "Forme des photons ajoutés", f"part plate (non synchrone du 850 nm) : {100 * fy.part_plate:.0f} %", fontsize=10.5)

    ax = axes[3]
    x = np.arange(4)
    for j, nom in enumerate(("yfp", "cfp")):
        l = L[(L.canal == nom) & (L.stim > 0)]
        ax.bar(x + (j - 0.5) * 0.38, l.exces_opto_pct, 0.36, yerr=1.96 * l.se_opto_pct, color=G.COULEUR_CANAL[nom],
               label=f"{C.NOMS[nom]}, phases opto", error_kw=dict(ecolor=G.INK_2, lw=0.8))
        ax.plot(x + (j - 0.5) * 0.38, l.exces_acq_pct, "o", color=G.INK, ms=4,
                label="phases d'acquisition" if j == 0 else None)
    ax.axhline(0, color=G.TRACE, lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"stim {k}" for k in range(1, 5)])
    ax.set_ylabel("photons ajoutés (%)", color=G.INK)
    ax.legend(frameon=False, fontsize=8)
    G.titre(ax, "Par stim", "barres : IC à 95 % ; points : phases d'acquisition", fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "fuite.png")


if __name__ == "__main__":
    main()
