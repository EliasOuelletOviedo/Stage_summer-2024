#!/usr/bin/env python3
"""
calibration_trois_lectures.py — SCLM à 850 nm (titration du 2026-10-16) : trois façons de mesurer [Cl⁻].

    f   modèle de sclm_hybride.py, deux canaux, même amplitude A (fit_joint)
    α   d_D = α e^{-t/τ_DA} + (1 - α) e^{-t/τ_D}, canal bleu seul (τ_DA, τ_D communs aux 8 concentrations)
    τ   durée de vie du donneur : une exponentielle sur le canal bleu (τ libre par déclin)
Chaque lecture : calibrée sur les sommes (10 acquisitions), isotherme X(c) = X∞ + (X₀ - X∞)/(1 + c/K), puis lue sur
les sommes amincies à N photons bleus (fond et décalage connus, amplitude libre ; valeurs négatives gardées).
Erreur = √(écart-type robuste² + biais médian²).

Usage : python calibration_trois_lectures.py [--N 1e4] [--tirages 300]
Sorties : calibration_trois_lectures.png, calibration_trois_lectures_850nm.csv
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sclm_hybride import load, fit_irf, fit_joint, Fine, CONC, dev_resid
from comparaison_canaux import bases, lire
from calibration_f_1_vs_2_canaux import fit_bleu

LAMBDA = 850
COUL = {"f": "#008300", "alpha": "#2a78d6", "tau": "#eb6834"}
NOM = {"f": "modèle f (2 canaux, même amplitude)", "alpha": "α (canal bleu seul)", "tau": "τ du donneur (mono-exponentielle)"}


def iso(c, X0, Xinf, K):
    return Xinf + (X0 - Xinf) / (1.0 + np.asarray(c, float) / K)


def ajuster_iso(x):
    pp, _ = curve_fit(iso, CONC, x, p0=[x[0], x[-1], 50.0], bounds=([-10, -10, 0.1], [10, 10, 1e4]), maxfev=20000)
    return pp


def inv_iso(X, p):
    X0, Xinf, K = p
    return K * (X0 - X) / (X - Xinf)


def mono(E, y, bg, w):
    """τ d'une exponentielle ⊗ IRF (déjà chargée dans E), amplitude libre, fond connu."""
    def r(th):
        m = np.exp(th[1]) * E.conv_exp(th[0])[w] + bg
        return dev_resid(y, np.maximum(m, 1e-9))
    A0 = max(y.sum() - bg.sum(), 1.0) / E.conv_exp(1.8)[w].sum()
    return least_squares(r, [1.8, np.log(A0)], bounds=([0.05, -np.inf], [10.0, np.inf]), x_scale='jac').x[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=float, default=1e4)
    ap.add_argument('--tirages', type=int, default=300)
    a = ap.parse_args()
    ici = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.default_rng(20261124)
    d = load(os.path.join(ici, f'irf_bruts_purs_{LAMBDA}nm.csv'))
    irfp = {1: fit_irf(d, 1, 2), 2: fit_irf(d, 2, 3)}
    w = d['win'][1]
    NC = len(CONC)
    R2, RB = fit_joint(d, irfp), fit_bleu(d, irfp)
    p2, pB = R2['p'], RB['p']
    E = Fine(d['T'], d['NB'])

    def irf_bleu(ic):
        i = irfp[1]
        E.set_irf(E.irf(i['mu'] + pB[f's1_{ic}'], i['sig'], i['weights'], i['taus']))

    # calibrations : valeur de chaque somme, puis isotherme
    val = {"f": np.array([p2[f'f1_{c}'] for c in range(NC)]), "alpha": np.array([pB[f'f1_{c}'] for c in range(NC)])}
    taus = []
    for ic in range(NC):
        irf_bleu(ic)
        taus.append(mono(E, d['B'][ic][w], np.full(w.sum(), pB[f'bgB{ic}']), w))
    val["tau"] = np.array(taus)
    P = {k: ajuster_iso(v) for k, v in val.items()}
    for k in NOM:
        print(f"{NOM[k]:38s} : {np.round(val[k], 3)}  K = {P[k][2]:.1f} mM")

    lignes = []
    for ic, c in enumerate(CONC):
        yB, yY = d['B'][ic][w], d['Y'][ic][w]
        q = a.N / yB.sum()
        tir = [(rng.binomial(yB.astype(np.int64), q).astype(float), rng.binomial(yY.astype(np.int64), q).astype(float))
               for _ in range(a.tirages)]
        b2 = {k: v[w] for k, v in bases(d, irfp, p2, p2[f's1_{ic}'], p2[f's2_{ic}']).items()}
        pBb = dict(pB, tA=3.0, lam=0.0, G=0.0, k1=0.0, k0=0.0, **{f's2_{ic}': 0.0})
        bB = {k: v[w] for k, v in bases(d, irfp, pBb, pB[f's1_{ic}'], 0.0).items()}
        g2B, g2Y, gB = (np.full(w.sum(), x) for x in (p2[f'bgB{ic}'], p2[f'bgY{ic}'], pB[f'bgB{ic}']))
        irf_bleu(ic)
        lus = {"f": [], "alpha": [], "tau": []}
        for tB, tY in tir:
            lus["f"].append(inv_iso(lire(tB, tY, "deux", b2, q * g2B, q * g2Y), P["f"]))
            lus["alpha"].append(inv_iso(lire(tB, tY, "bleu", bB, q * gB, q * gB), P["alpha"]))
            lus["tau"].append(inv_iso(mono(E, tB, q * gB, w), P["tau"]))
        for k, v in lus.items():
            v = np.asarray(v, float)
            v = v[np.isfinite(v)]
            sd = (np.percentile(v, 75) - np.percentile(v, 25)) / 1.349
            bm = float(np.median(v) - c)
            lignes.append(dict(cl_mM=c, lecture=NOM[k], N_bleu=a.N, ecart_type_robuste=sd, biais_median=bm, erreur=float(np.hypot(sd, bm))))
        t = pd.DataFrame(lignes[-3:]).set_index('lecture')
        print(f"  {c:5.1f} mM : " + " · ".join(f"{k} σ {t.loc[NOM[k], 'ecart_type_robuste']:.2f}, biais {t.loc[NOM[k], 'biais_median']:+.2f}"
                                              for k in NOM))
    res = pd.DataFrame(lignes)
    res.to_csv(os.path.join(ici, f'calibration_trois_lectures_{LAMBDA}nm.csv'), index=False)

    # figure
    virg = lambda x: f'{x:g}'.replace('.', ',')
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    cc = np.linspace(0, 131, 400)

    def axe_c(ax):
        ax.set_xscale('symlog', linthresh=2.5, linscale=0.6)
        ax.set_xticks(CONC); ax.set_xticklabels([virg(c) for c in CONC])
        ax.set_xlabel('[Cl⁻] imposée (mM)'); ax.grid(True, alpha=0.3)
    for ax, k, lab, tit in ((axes[0, 0], "f", "f", "(a) Modèle f, deux canaux"), (axes[0, 1], "alpha", "α", "(b) α, canal bleu seul"),
                            (axes[1, 0], "tau", "τ du donneur (ns)", "(c) Durée de vie du donneur, une exponentielle")):
        ax.plot(CONC, val[k], 'o', color=COUL[k], ms=7)
        ax.plot(cc, iso(cc, *P[k]), color=COUL[k], lw=1.6, label=f"K = {virg(round(P[k][2], 1))} mM")
        axe_c(ax); ax.set_ylabel(lab); ax.legend(frameon=False, fontsize=10)
        ax.set_title(tit, loc='left', fontsize=12)
    ax = axes[1, 1]
    for k in NOM:
        s = res[res.lecture == NOM[k]]
        ax.plot(s.cl_mM, s.erreur, 'o-', color=COUL[k], ms=6, lw=1.8, label=NOM[k])
    axe_c(ax); ax.set_yscale('log')
    yt = [v for v in (0.2, 0.5, 1, 2, 5, 10, 20, 50, 100) if res.erreur.min() / 1.5 <= v <= res.erreur.max() * 1.5]
    ax.set_yticks(yt); ax.set_yticklabels([virg(v) for v in yt]); ax.minorticks_off()
    ax.set_ylabel('erreur sur la [Cl⁻] lue (mM)')
    ax.set_title('(d) Erreur à 10⁴ photons bleus (bruit et biais)', loc='left', fontsize=12)
    ax.legend(frameon=False, fontsize=9)
    fig.suptitle('SCLM à 850 nm : modèle f, α et durée de vie du donneur', x=0.01, ha='left', fontsize=14)
    fig.text(0.01, 0.005, "(a–c) valeur de chaque somme de 10 acquisitions et isotherme ; (d) sommes amincies, "
             f"{a.tirages} lectures : √(écart-type robuste² + biais médian²)", fontsize=9, color='#52514e')
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(os.path.join(ici, 'calibration_trois_lectures.png'), dpi=150)
    print('Figure écrite : calibration_trois_lectures.png')


if __name__ == '__main__':
    main()
