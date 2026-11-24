#!/usr/bin/env python3
"""
comparaison_canaux.py — mesure de [Cl⁻] avec le canal bleu seul ou avec les deux canaux (modèle à une population FRET).

Pour chaque longueur d'onde (800, 850 nm) : ajustement conjoint de la titration (sclm_hybride.fit_joint, une
population FRET), puis lecture de [Cl⁻] histogramme par histogramme avec le même modèle, constantes figées :

    bleu  : A [ f e^{-t/τ1} + (1 - f) e^{-t/τD} ] ⊗ IRF + fond                         (affine en f)
    jaune : η A [ λ·donneur(t) + G f k (montée) + (k1 f + k0) e^{-t/τA} ] ⊗ IRF + fond   (affine en f)

Lecteurs (maximum de vraisemblance de Poisson, f libre de signe et de borne, fond et décalages connus) :
    bleu seul           f, A              — forme du déclin du donneur
    jaune seul          f, A              — forme du canal jaune
    2 canaux, η libre   f, A, η           — les deux formes, sans le rapport d'intensité
    2 canaux, η connu   f, A (η = 1)      — les deux formes et le rapport jaune / bleu
f est converti en [Cl⁻] par l'isotherme f(c) = fmin + (fmax - fmin)/(1 + c/K) inversée (valeurs négatives gardées).

Deux mesures de la précision, pour un même nombre de photons bleus N (le jaune vient avec, au rapport de l'échantillon) :
    - borne de Cramér-Rao du modèle (information de Fisher, fond compris) ;
    - écart-type robuste des lectures sur les histogrammes réels amincis (tirage binomial des sommes de 10 acquisitions).
Et le biais systématique de chaque lecteur : lecture de la somme complète − [Cl⁻] imposée.

Usage : python comparaison_canaux.py [--N 1e4] [--tirages 300]
Sorties : comparaison_canaux_<λ>nm.csv, comparaison_canaux.png
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sclm_hybride import load, fit_irf, fit_joint, fit_isotherms, Fine, CONC, dev_resid

LECTEURS = {"bleu": "bleu seul (forme du donneur)", "jaune": "jaune seul (forme)",
            "deux_eta": "2 canaux, η libre (formes seules)", "deux": "2 canaux, η connu (formes + rapport J/B)"}
COUL = {"bleu": "#2a78d6", "jaune": "#c99a06", "deux_eta": "#4a3aa7", "deux": "#008300"}


def bases(d, irfp, p, s1, s2):
    """Vecteurs u, v (μ = A (u + f v)) des deux canaux pour des décalages donnés."""
    E = Fine(d['T'], d['NB'])
    E.set_irf(E.irf(irfp[1]['mu'] + s1, irfp[1]['sig'], irfp[1]['weights'], irfp[1]['taus']))
    b1, bD = E.conv_exp(p['t1']), E.conv_exp(p['tD'])
    E.set_irf(E.irf(irfp[2]['mu'] + s2, irfp[2]['sig'], irfp[2]['weights'], irfp[2]['taus']))
    y1, yD, yA = E.conv_exp(p['t1']), E.conv_exp(p['tD']), E.conv_exp(p['tA'])
    rise = p['t1'] * p['tA'] / (p['tA'] - p['t1']) * (yA - y1)
    kT = 1 / p['t1'] - 1 / p['tD']
    return dict(uB=bD, vB=b1 - bD, uY=p['lam'] * yD + p['k0'] * yA,
                vY=p['lam'] * (y1 - yD) + p['G'] * kT * rise + p['k1'] * yA)


def modele(th, lect, b, bgB, bgY):
    f, A = th[0], np.exp(th[1])
    eta = np.exp(th[2]) if lect == "deux_eta" else 1.0
    out = []
    if lect != "jaune":
        out.append(A * (b['uB'] + f * b['vB']) + bgB)
    if lect != "bleu":
        out.append(A * eta * (b['uY'] + f * b['vY']) + bgY)
    return np.maximum(np.concatenate(out), 1e-9)


def lire(yB, yY, lect, b, bgB, bgY):
    y = np.concatenate([v for v, k in ((yB, "B"), (yY, "Y")) if (k == "B" and lect != "jaune") or (k == "Y" and lect != "bleu")])
    ref = b['uB'] + 0.5 * b['vB'] if lect != "jaune" else b['uY'] + 0.5 * b['vY']
    obs = yB if lect != "jaune" else yY
    A0 = max(obs.sum() - (bgB if lect != "jaune" else bgY).sum(), 1.0) / ref.sum()
    th0 = [0.5, np.log(A0)] + ([0.0] if lect == "deux_eta" else [])
    r = least_squares(lambda th: dev_resid(y, modele(th, lect, b, bgB, bgY)), th0, x_scale='jac', max_nfev=400)
    return r.x[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=float, default=1e4, help='photons bleus par acquisition (fenêtre)')
    ap.add_argument('--tirages', type=int, default=300)
    ap.add_argument('--sans-calcul', action='store_true', help='refait seulement la figure depuis les CSV')
    a = ap.parse_args()
    ici = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.default_rng(20261123)
    resultats = {}
    if a.sans_calcul:
        resultats = {lam: pd.read_csv(os.path.join(ici, f'comparaison_canaux_{lam}nm.csv')) for lam in (800, 850)}
    for lam in ((800, 850) if not a.sans_calcul else ()):
        d = load(os.path.join(ici, f'irf_bruts_purs_{lam}nm.csv'))
        irfp = {1: fit_irf(d, 1, 2), 2: fit_irf(d, 2, 3)}
        R = fit_joint(d, irfp)
        p = R['p']
        iso = fit_isotherms(R)['f1']
        fmax, fmin, K = iso['fmax'], iso['fmin'], iso['K']
        f_de = lambda c: fmin + (fmax - fmin) / (1 + c / K)
        c_de = lambda f: K * (fmax - f) / (f - fmin)
        pente = lambda c: (fmax - fmin) / K / (1 + c / K) ** 2          # |df/dc|
        w = d['win'][1]
        print(f"\n{lam} nm : τ1 = {p['t1']:.3f}, τD = {p['tD']:.3f}, τA = {p['tA']:.3f} ns, λ = {p['lam']:.3f}, G = {p['G']:.3f}, "
              f"k1 = {p['k1']:.3f}, k0 = {p['k0']:.3f} ; déviance/ddl {R['dev'] / R['dof']:.1f} ; f : {fmax:.3f} → {fmin:.3f}, K = {K:.1f} mM")
        lignes = []
        for ic, c in enumerate(CONC):
            b = {k: v[w] for k, v in bases(d, irfp, p, p[f's1_{ic}'], p[f's2_{ic}']).items()}
            yB, yY = d['B'][ic][w], d['Y'][ic][w]
            bgB, bgY = np.full(w.sum(), p[f'bgB{ic}']), np.full(w.sum(), p[f'bgY{ic}'])
            # biais systématique : lecture de la somme complète
            plein = {k: c_de(lire(yB, yY, k, b, bgB, bgY)) - c for k in LECTEURS}
            # amincissement à N photons bleus
            q = a.N / yB.sum()
            lus = {k: [] for k in LECTEURS}
            for _ in range(a.tirages):
                tB, tY = rng.binomial(yB.astype(np.int64), q).astype(float), rng.binomial(yY.astype(np.int64), q).astype(float)
                for k in LECTEURS:
                    lus[k].append(c_de(lire(tB, tY, k, b, q * bgB, q * bgY)))
            # Cramér-Rao (fond compris), mêmes N
            A, f0 = p[f'A{ic}'] * q, R['p'][f'f1_{ic}']
            cr = {}
            for k in LECTEURS:
                blocs = []
                if k != "jaune":
                    m = A * (b['uB'] + f0 * b['vB']) + q * bgB
                    blocs.append((m, np.c_[A * b['vB'], b['uB'] + f0 * b['vB'], np.zeros_like(m)][:, :3 if k == "deux_eta" else 2]))
                if k != "bleu":
                    m = A * (b['uY'] + f0 * b['vY']) + q * bgY
                    blocs.append((m, np.c_[A * b['vY'], b['uY'] + f0 * b['vY'], A * (b['uY'] + f0 * b['vY'])][:, :3 if k == "deux_eta" else 2]))
                I = sum((D.T / m) @ D for m, D in blocs)
                cr[k] = np.sqrt(np.linalg.inv(I)[0, 0]) / pente(max(c, 0))
            for k in LECTEURS:
                v = np.asarray(lus[k], float)
                v = v[np.isfinite(v)]
                lignes.append(dict(longueur_onde=lam, cl_mM=c, lecteur=k, N_bleu=a.N, photons_jaunes=float(q * yY.sum()),
                                   ecart_type_robuste=(np.percentile(v, 75) - np.percentile(v, 25)) / 1.349,
                                   biais_median=float(np.median(v) - c), cramer_rao=cr[k], biais_somme_complete=plein[k]))
            t = pd.DataFrame(lignes[-4:]).set_index('lecteur')
            print(f"  {c:5.1f} mM (J/B {yY.sum() / yB.sum():.1f}) : " + " · ".join(
                f"{k} σ {t.loc[k, 'ecart_type_robuste']:.2f} (CR {t.loc[k, 'cramer_rao']:.2f}, biais {t.loc[k, 'biais_somme_complete']:+.2f})"
                for k in LECTEURS))
        df = pd.DataFrame(lignes)
        df.to_csv(os.path.join(ici, f'comparaison_canaux_{lam}nm.csv'), index=False)
        resultats[lam] = df

    # figure : écart-type de la [Cl⁻] lue, réel aminci (points) et borne de Cramér-Rao (traits)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, lam in zip(axes, (800, 850)):
        df = resultats[lam]
        for k, nom in LECTEURS.items():
            s = df[df.lecteur == k]
            # bleu seul et 2 canaux à η libre coïncident : le bleu en trait large et anneaux, l'autre par-dessus
            large = k == "bleu"
            ax.plot(s.cl_mM, s.cramer_rao, '-', color=COUL[k], lw=5 if large else 1.8, alpha=0.35 if large else 0.85,
                    zorder=1 if large else 2)
            ax.plot(s.cl_mM, s.ecart_type_robuste, 'o', ms=12 if large else 6, label=nom, zorder=1 if large else 3,
                    markerfacecolor='none' if large else COUL[k], markeredgecolor=COUL[k], markeredgewidth=2 if large else 0)
        ax.set_xscale('symlog', linthresh=2.5, linscale=0.6)
        ax.set_xticks(CONC)
        ax.set_xticklabels([f'{c:g}'.replace('.', ',') for c in CONC])
        ax.set_yscale('log')
        ax.set_xlabel('[Cl⁻] (mM)')
        ax.set_title(f'SCLM, {lam} nm · {a.N:.0e} photons bleus par acquisition'.replace('1e+04', '10⁴'), loc='left', fontsize=12)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel('écart-type de la [Cl⁻] lue (mM)')
    axes[0].legend(fontsize=9, frameon=False)
    fig.text(0.01, 0.005, 'points : histogrammes réels amincis (écart-type robuste) ; traits : borne de Cramér-Rao du modèle à une population FRET ; '
             'bleu seul (anneaux) et 2 canaux à η libre (points pleins) se superposent',
             fontsize=9, color='#52514e')
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(ici, 'comparaison_canaux.png'), dpi=150)
    print('\nFigure écrite : comparaison_canaux.png')


if __name__ == '__main__':
    main()
