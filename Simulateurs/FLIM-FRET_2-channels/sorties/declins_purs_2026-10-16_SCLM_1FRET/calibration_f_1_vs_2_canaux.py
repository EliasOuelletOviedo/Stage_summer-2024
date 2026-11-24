#!/usr/bin/env python3
"""
calibration_f_1_vs_2_canaux.py — calibration de SCLM avec le canal bleu seul contre avec les deux canaux, modèle
de sclm_hybride.py (une population FRET, variable f), titration du 2026-10-16, 800 et 850 nm.

1 canal  : ajustement des 8 déclins bleus seuls, τ1 et τD communs ; f, amplitude, fond et décalage par concentration :
               bleu = A [ f e^{-t/τ1} + (1 - f) e^{-t/τD} ] ⊗ IRF + fond
2 canaux : l'ajustement conjoint de sclm_hybride.fit_joint, même amplitude A dans les deux canaux :
               bleu  = A [ f e^{-t/τ1} + (1 - f) e^{-t/τD} ]
               jaune = A [ λ·donneur(t) + G f k (e^{-t/τ1} ⊛ e^{-t/τA}) + (k1 f + k0) e^{-t/τA} ],  k = 1/τ1 - 1/τD
Pour chacune : isotherme f(c) = fmin + (fmax - fmin)/(1 + c/K) sur les 8 fractions.

Mesure de [Cl⁻] : chaque somme (10 acquisitions) est amincie à N photons bleus (le jaune suit) ; f est lu avec les
constantes de sa calibration (maximum de vraisemblance de Poisson, amplitude libre, fond et décalage connus) puis
converti par son isotherme (valeurs négatives gardées). Écart-type robuste des lectures, borne de Cramér-Rao, écart
de la somme complète à l'isotherme.

(Recréé le 2026-11-24 : calibration_1_vs_2_canaux.py a depuis été réécrit pour le modèle α, β.)

Usage : python calibration_f_1_vs_2_canaux.py [--N 1e4] [--tirages 300]
Sorties : calibration_f_1_vs_2_canaux.png, calibration_f_1_vs_2_canaux_<λ>nm.csv, calibration_f_1_vs_2_canaux_parametres.csv
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, curve_fit
from scipy.optimize._numdiff import approx_derivative
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sclm_hybride import load, fit_irf, fit_joint, Fine, CONC, Params, dev_resid, isotherm
from comparaison_canaux import bases, lire

COUL = {"1": "#2a78d6", "2": "#008300"}
NOM = {"1": "1 canal (bleu)", "2": "2 canaux (bleu + jaune)"}
BASE = 'calibration_f_1_vs_2_canaux'


def fit_bleu(d, irfp):
    """Calibration sur le canal bleu seul : τ1, τD communs ; décalage, amplitude, f et fond par concentration."""
    E = Fine(d['T'], d['NB'])
    NC = len(CONC)
    P = Params()
    P.add('t1', 0.5, 0.02, 3.0); P.add('tD', 3.75, 2.5, 6.0)            # mêmes bornes que l'ajustement conjoint
    for c in range(NC):
        P.add(f's1_{c}', 0.0, -0.3, 0.3); P.add(f'A{c}', d['B'][c].max() / d['DT'], 0, np.inf)
        P.add(f'f1_{c}', 0.4, 0, 1); P.add(f'bgB{c}', 300.0, 0, np.inf)
    w = d['win'][1]

    def resid(v):
        p = P.unpack(v)
        out = []
        for c in range(NC):
            E.set_irf(E.irf(irfp[1]['mu'] + p[f's1_{c}'], irfp[1]['sig'], irfp[1]['weights'], irfp[1]['taus']))
            m = p[f'A{c}'] * (p[f'f1_{c}'] * E.conv_exp(p['t1']) + (1 - p[f'f1_{c}']) * E.conv_exp(p['tD'])) + p[f'bgB{c}']
            out.append(dev_resid(d['B'][c][w], m[w]))
        return np.concatenate(out)
    r = least_squares(resid, P.start, bounds=P.bounds, x_scale='jac', max_nfev=20000, xtol=1e-14, ftol=1e-14, gtol=1e-12)
    dev = float(np.sum(r.fun ** 2)); dof = w.sum() * NC - len(r.x)
    J = approx_derivative(resid, r.x, method='3-point', bounds=P.bounds)
    cov = np.linalg.pinv(J.T @ J) * dev / dof
    return dict(P=P, x=r.x, p=P.unpack(r.x), cov=cov, dev=dev, dof=int(dof))


def isotherme_de(R):
    idx = {n: i for i, n in enumerate(R['P'].names)}
    f = np.array([R['p'][f'f1_{c}'] for c in range(len(CONC))])
    s = np.array([np.sqrt(R['cov'][idx[f'f1_{c}'], idx[f'f1_{c}']]) for c in range(len(CONC))])
    pp, cc = curve_fit(isotherm, CONC, f, p0=[f[0], 0.0, 60], sigma=s, bounds=([0, -0.5, 1], [1, 1, 1e4]))
    return dict(f=f, se=s, fmax=pp[0], fmin=pp[1], K=pp[2], se_K=float(np.sqrt(cc[2, 2])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=float, default=1e4, help='photons bleus par acquisition')
    ap.add_argument('--tirages', type=int, default=300)
    a = ap.parse_args()
    ici = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.default_rng(20261124)
    res, params, isos = {}, [], {}
    for lam in (800, 850):
        d = load(os.path.join(ici, f'irf_bruts_purs_{lam}nm.csv'))
        irfp = {1: fit_irf(d, 1, 2), 2: fit_irf(d, 2, 3)}
        w = d['win'][1]
        cal = {"1": fit_bleu(d, irfp), "2": fit_joint(d, irfp)}
        isos[lam] = {k: isotherme_de(R) for k, R in cal.items()}
        for k, R in cal.items():
            p, iso = R['p'], isos[lam][k]
            params.append(dict(longueur_onde=lam, calibration=NOM[k], tau1_ns=p['t1'], tauD_ns=p['tD'],
                               deviance_par_ddl=R['dev'] / R['dof'], f_0mM=iso['f'][0], f_131mM=iso['f'][-1],
                               fmax=iso['fmax'], fmin=iso['fmin'], K_mM=iso['K'], K_erreur_type=iso['se_K']))
            print(f"{lam} nm, {NOM[k]:24s} : τ1 = {p['t1']:.3f}, τD = {p['tD']:.3f} ns, déviance/ddl {R['dev'] / R['dof']:.1f} ; "
                  f"f {iso['f'][0]:.3f} → {iso['f'][-1]:.3f}, K = {iso['K']:.1f} ± {iso['se_K']:.1f} mM")
        lignes = []
        for ic, c in enumerate(CONC):
            yB, yY = d['B'][ic][w], d['Y'][ic][w]
            q = a.N / yB.sum()
            tirages = [(rng.binomial(yB.astype(np.int64), q).astype(float), rng.binomial(yY.astype(np.int64), q).astype(float))
                       for _ in range(a.tirages)]
            for k, R in cal.items():
                p, iso = dict(R['p']), isos[lam][k]
                if k == "1":                               # le jaune n'est pas lu : constantes jaunes sans effet
                    p.update(tA=3.0, lam=0.0, G=0.0, k1=0.0, k0=0.0, **{f's2_{ic}': 0.0, f'bgY{ic}': 0.0})
                b = {kk: v[w] for kk, v in bases(d, irfp, p, p[f's1_{ic}'], p[f's2_{ic}']).items()}
                bgB, bgY = np.full(w.sum(), p[f'bgB{ic}']), np.full(w.sum(), p[f'bgY{ic}'])
                lecteur = "bleu" if k == "1" else "deux"
                c_de = lambda f: iso['K'] * (iso['fmax'] - f) / (f - iso['fmin'])
                pente = (iso['fmax'] - iso['fmin']) / iso['K'] / (1 + c / iso['K']) ** 2
                ecart = c_de(lire(yB, yY, lecteur, b, bgB, bgY)) - c
                lus = np.array([c_de(lire(tB, tY, lecteur, b, q * bgB, q * bgY)) for tB, tY in tirages])
                lus = lus[np.isfinite(lus)]
                A, f0 = p[f'A{ic}'] * q, p[f'f1_{ic}']
                blocs = [(A * (b['uB'] + f0 * b['vB']) + q * bgB, np.c_[A * b['vB'], b['uB'] + f0 * b['vB']])]
                if k == "2":
                    blocs.append((A * (b['uY'] + f0 * b['vY']) + q * bgY, np.c_[A * b['vY'], b['uY'] + f0 * b['vY']]))
                I = sum((D.T / m) @ D for m, D in blocs)
                lignes.append(dict(longueur_onde=lam, cl_mM=c, calibration=NOM[k], N_bleu=a.N,
                                   ecart_type_robuste=(np.percentile(lus, 75) - np.percentile(lus, 25)) / 1.349,
                                   cramer_rao=np.sqrt(np.linalg.inv(I)[0, 0]) / pente, ecart_isotherme=ecart))
            t = pd.DataFrame(lignes[-2:]).set_index('calibration')
            print(f"   {c:5.1f} mM : " + " · ".join(f"{n} σ {t.loc[n, 'ecart_type_robuste']:.2f} (CR {t.loc[n, 'cramer_rao']:.2f}, "
                                                    f"écart {t.loc[n, 'ecart_isotherme']:+.2f})" for n in t.index))
        res[lam] = pd.DataFrame(lignes)
        res[lam].to_csv(os.path.join(ici, f'{BASE}_{lam}nm.csv'), index=False)
    pd.DataFrame(params).to_csv(os.path.join(ici, f'{BASE}_parametres.csv'), index=False)

    virg = lambda x: f'{x:g}'.replace('.', ',')
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    cc = np.linspace(0, 131, 400)
    for i, lam in enumerate((800, 850)):
        ax = axes[i, 0]
        for k in ("1", "2"):
            iso = isos[lam][k]
            ax.errorbar(CONC, iso['f'], yerr=iso['se'], fmt='o', color=COUL[k], ms=6, capsize=3,
                        label=f"{NOM[k]} : K = {virg(round(iso['K'], 1))} mM")
            ax.plot(cc, isotherm(cc, iso['fmax'], iso['fmin'], iso['K']), color=COUL[k], lw=1.6)
        ax.set_xscale('symlog', linthresh=2.5, linscale=0.6)
        ax.set_xticks(CONC); ax.set_xticklabels([virg(c) for c in CONC])
        ax.set_xlabel('[Cl⁻] imposée (mM)'); ax.set_ylabel('fraction FRET f')
        ax.set_title(f'Calibration, {lam} nm', loc='left', fontsize=12)
        ax.legend(frameon=False, fontsize=9); ax.grid(True, alpha=0.3)
        ax = axes[i, 1]
        for k in ("1", "2"):
            s = res[lam][res[lam].calibration == NOM[k]]
            ax.plot(s.cl_mM, s.cramer_rao, '-', color=COUL[k], lw=1.6, alpha=0.8)
            ax.plot(s.cl_mM, s.ecart_type_robuste, 'o', color=COUL[k], ms=6, label=NOM[k])
        ax.set_xscale('symlog', linthresh=2.5, linscale=0.6); ax.set_yscale('log')
        ax.set_xticks(CONC); ax.set_xticklabels([virg(c) for c in CONC])
        yt = [0.2, 0.5, 1, 2, 5, 10, 20, 50]
        ax.set_yticks(yt); ax.set_yticklabels([virg(v) for v in yt]); ax.minorticks_off()
        ax.set_xlabel('[Cl⁻] imposée (mM)'); ax.set_ylabel('écart-type de la [Cl⁻] lue (mM)')
        ax.set_title(f'Précision à 10⁴ photons bleus, {lam} nm', loc='left', fontsize=12)
        ax.legend(frameon=False, fontsize=9); ax.grid(True, alpha=0.3)
    fig.text(0.01, 0.005, "gauche : fraction FRET de chaque concentration (± erreur-type) et isotherme ajustée ; droite : points = "
             "sommes réelles amincies (écart-type robuste de 300 lectures), traits = borne de Cramér-Rao", fontsize=9, color='#52514e')
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    fig.savefig(os.path.join(ici, f'{BASE}.png'), dpi=150)
    print(f'Figure écrite : {BASE}.png')


if __name__ == '__main__':
    main()
