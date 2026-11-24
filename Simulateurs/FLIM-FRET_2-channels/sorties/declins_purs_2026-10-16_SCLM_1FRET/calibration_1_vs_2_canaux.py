#!/usr/bin/env python3
"""
calibration_1_vs_2_canaux.py — calibration de SCLM à 850 nm (titration du 2026-10-16) : canal bleu seul, canal jaune
seul, ou les deux canaux.

Modèle de chaque canal (déclin ⊛ IRF mesurée du canal, décalée ; amplitude, décalage et fond par acquisition) :
    bleu  (donneur)    d_D(t) = α e^{-t/τ_DA} + (1 - α) e^{-t/τ_D}
    jaune (accepteur)  d_A(t) = β [e^{-t/τ_DA} ⊛ e^{-t/τ_A}] + (1 - β) e^{-t/τ_A}       (⊛ : convolution, t en ns)
    canal jaune mesuré = A_A d_A + λ · (photons du donneur) d_D     (fuite du donneur dans le jaune, λ commun)

Calibrations, chacune avec ses propres durées de vie :
    bleu seul   τ_DA, τ_D communs ; α par concentration
    jaune seul  τ_DA, τ_D, τ_A, λ communs ; β et α par concentration (fuite rapportée à A_A, faute de bleu)
    2 canaux    τ_DA commun aux deux canaux, τ_D, τ_A, λ ; α, β et amplitudes A_D, A_A par concentration ;
                la fuite vaut λ A_D d_D (le α du bleu)
Isotherme de chaque variable X (α, β, et ρ = A_A / A_D pour les deux canaux) : X(c) = X∞ + (X₀ - X∞)/(1 + c/K).

Mesure de [Cl⁻] : chaque somme (10 acquisitions) est amincie à N photons bleus (le jaune suit), puis lue par maximum
de vraisemblance de Poisson avec les constantes de sa calibration, fond et décalages connus, valeurs négatives gardées :
    bleu seul               α (amplitude libre) → [Cl⁻] par l'isotherme de α
    jaune seul              β (amplitude libre) → [Cl⁻] par l'isotherme de β
    2 canaux : α et β       [Cl⁻] directement, α(c) et β(c), amplitudes des deux canaux libres
    2 canaux : α, β et ρ    [Cl⁻] directement, et l'amplitude du jaune suit celle du bleu : A_A = ρ(c) A_D
Écart-type robuste des lectures, borne de Cramér-Rao, écart de la somme complète à l'isotherme.

Usage : python calibration_1_vs_2_canaux.py [--N 1e4] [--tirages 300]
Sorties : calibration_1_vs_2_canaux.png, calibration_1_vs_2_canaux_850nm.csv, calibration_1_vs_2_canaux_parametres.csv
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

from sclm_hybride import load, fit_irf, Fine, CONC, Params, dev_resid

LAMBDA = 850
NC = len(CONC)
COUL = {"bleu": "#2a78d6", "jaune": "#c99a06", "ab": "#4a3aa7", "abr": "#008300"}
NOM = {"bleu": "1 canal : bleu seul (α)", "jaune": "1 canal : jaune seul (β)",
       "ab": "2 canaux : α et β", "abr": "2 canaux : α, β et intensité du jaune"}


def iso(c, X0, Xinf, K):
    return Xinf + (X0 - Xinf) / (1.0 + np.asarray(c, float) / K)


def inv_iso(X, p):
    X0, Xinf, K = p
    return K * (X0 - X) / (X - Xinf)


def diso(c, p):
    X0, Xinf, K = p
    return -(X0 - Xinf) / K / (1.0 + c / K) ** 2


class Formes:
    """Déclins ⊛ IRF des deux canaux pour des durées de vie et des décalages donnés (vecteurs sur toute la grille)."""
    def __init__(self, d, irfp):
        self.d, self.irfp, self.E = d, irfp, Fine(d['T'], d['NB'])

    def bleu(self, tDA, tD, s):
        i = self.irfp[1]
        self.E.set_irf(self.E.irf(i['mu'] + s, i['sig'], i['weights'], i['taus']))
        return self.E.conv_exp(tDA), self.E.conv_exp(tD)                      # FRET, sans FRET

    def jaune(self, tDA, tA, tD, s):
        i = self.irfp[2]
        self.E.set_irf(self.E.irf(i['mu'] + s, i['sig'], i['weights'], i['taus']))
        yA, yDA, yD = self.E.conv_exp(tA), self.E.conv_exp(tDA), self.E.conv_exp(tD)
        return tDA * tA / (tA - tDA) * (yA - yDA), yA, yDA, yD               # FRET (convolution), directe, fuite (DA, D)


def calibrer(d, irfp, canaux):
    """Ajustement d'un ou des deux canaux sur les 8 concentrations ; renvoie paramètres, covariance, déviance par canal."""
    F, w = Formes(d, irfp), d['win'][1]
    P = Params()
    P.add('tDA', 0.5, 0.02, 2.0)
    if 'B' in canaux:
        P.add('tD', 2.8, 1.0, 6.0)
    if 'Y' in canaux:
        P.add('tA', 3.2, 1.0, 6.0); P.add('lam', 0.1, 0, 5)
        if 'B' not in canaux:
            P.add('tD', 2.8, 1.0, 6.0)
    for c in range(NC):
        if 'B' in canaux:
            P.add(f'sB{c}', 0.0, -0.3, 0.3); P.add(f'AB{c}', d['B'][c].max() / d['DT'], 0, np.inf)
            P.add(f'al{c}', 0.5, 0, 1); P.add(f'bgB{c}', 300.0, 0, np.inf)
        if 'Y' in canaux:
            P.add(f'sY{c}', 0.0, -0.3, 0.3); P.add(f'AY{c}', d['Y'][c].max() / d['DT'], 0, np.inf)
            P.add(f'be{c}', 0.3, 0, 1); P.add(f'bgY{c}', 300.0, 0, np.inf)
            if 'B' not in canaux:
                P.add(f'al{c}', 0.5, 0, 1)

    def resid(v):
        p = P.unpack(v)
        rB, rY = [], []
        for c in range(NC):
            if 'B' in canaux:
                b1, bD = F.bleu(p['tDA'], p['tD'], p[f'sB{c}'])
                m = p[f'AB{c}'] * (p[f'al{c}'] * b1 + (1 - p[f'al{c}']) * bD) + p[f'bgB{c}']
                rB.append(dev_resid(d['B'][c][w], m[w]))
            if 'Y' in canaux:
                ri, yA, yDA, yD = F.jaune(p['tDA'], p['tA'], p['tD'], p[f'sY{c}'])
                fuite = p['lam'] * p[f'AB{c}' if 'B' in canaux else f'AY{c}'] * (p[f'al{c}'] * yDA + (1 - p[f'al{c}']) * yD)
                m = p[f'AY{c}'] * (p[f'be{c}'] * ri + (1 - p[f'be{c}']) * yA) + fuite + p[f'bgY{c}']
                rY.append(dev_resid(d['Y'][c][w], m[w]))
        return np.concatenate(rB + rY)
    r = least_squares(resid, P.start, bounds=P.bounds, x_scale='jac', max_nfev=20000, xtol=1e-12, ftol=1e-12)
    dev, dof = float(np.sum(r.fun ** 2)), len(r.fun) - len(r.x)
    J = approx_derivative(resid, r.x, method='3-point', bounds=P.bounds)
    cov = np.linalg.pinv(J.T @ J) * dev / dof
    return dict(P=P, p=P.unpack(r.x), cov=cov, idx={n: i for i, n in enumerate(P.names)}, dev_ddl=dev / dof)


def isotherme(R, nom):
    """Isotherme de α (nom 'al'), β ('be') ou ρ = A_A/A_D ('rho'), pondérée par les erreurs-types de l'ajustement."""
    p, cov, idx = R['p'], R['cov'], R['idx']
    if nom == 'rho':
        x = np.array([p[f'AY{c}'] / p[f'AB{c}'] for c in range(NC)])
        se = []
        for c in range(NC):
            iY, iB = idx[f'AY{c}'], idx[f'AB{c}']
            g = np.array([1 / p[f'AB{c}'], -p[f'AY{c}'] / p[f'AB{c}'] ** 2])
            C = cov[np.ix_([iY, iB], [iY, iB])]
            se.append(float(np.sqrt(max(g @ C @ g, 1e-30))))
        se = np.array(se)
    else:
        x = np.array([p[f'{nom}{c}'] for c in range(NC)])
        se = np.array([np.sqrt(max(cov[idx[f'{nom}{c}'], idx[f'{nom}{c}']], 1e-30)) for c in range(NC)])
    se = np.maximum(se, 1e-3 * max(np.ptp(x), 1e-6))                          # valeurs collées à une borne
    try:
        pp, cc = curve_fit(iso, CONC, x, p0=[x[0], x[-1], 30.0], sigma=se, bounds=([-2, -2, 0.1], [5, 5, 1e4]), maxfev=20000)
        return dict(x=x, se=se, p=pp, se_K=float(np.sqrt(cc[2, 2])))
    except Exception:
        return dict(x=x, se=se, p=None, se_K=np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--N', type=float, default=1e4, help='photons bleus par acquisition')
    ap.add_argument('--tirages', type=int, default=300)
    a = ap.parse_args()
    ici = os.path.dirname(os.path.abspath(__file__))
    rng = np.random.default_rng(20261124)
    d = load(os.path.join(ici, f'irf_bruts_purs_{LAMBDA}nm.csv'))
    irfp = {1: fit_irf(d, 1, 2), 2: fit_irf(d, 2, 3)}
    w = d['win'][1]
    F = Formes(d, irfp)

    cal = {"bleu": calibrer(d, irfp, "B"), "jaune": calibrer(d, irfp, "Y"), "deux": calibrer(d, irfp, "BY")}
    I = {"bleu": {"al": isotherme(cal["bleu"], "al")}, "jaune": {"be": isotherme(cal["jaune"], "be")},
         "deux": {k: isotherme(cal["deux"], k) for k in ("al", "be", "rho")}}
    lignes_p = []
    for k, R in cal.items():
        p = R['p']
        txt = ", ".join(f"{n} = {p[n]:.3f} ns" for n in ("tDA", "tD", "tA") if n in p) + (f", λ = {p['lam']:.3f}" if 'lam' in p else "")
        print(f"{k:5s} : {txt} ; déviance/ddl {R['dev_ddl']:.1f}")
        for v, s in I[k].items():
            K = s['p'][2] if s['p'] is not None else np.nan
            print(f"       {v:3s} : {np.round(s['x'], 3)}  K = {K:.1f} mM")
            lignes_p.append(dict(calibration=k, variable={'al': 'alpha', 'be': 'beta', 'rho': 'A_A/A_D'}[v],
                                 tau_DA=p['tDA'], tau_D=p.get('tD', np.nan), tau_A=p.get('tA', np.nan), deviance_ddl=R['dev_ddl'],
                                 **{f'{c:g}mM': x for c, x in zip(CONC, s['x'])},
                                 X0=s['p'][0] if s['p'] is not None else np.nan, Xinf=s['p'][1] if s['p'] is not None else np.nan,
                                 K_mM=K, K_erreur_type=s['se_K']))
    pd.DataFrame(lignes_p).to_csv(os.path.join(ici, 'calibration_1_vs_2_canaux_parametres.csv'), index=False)

    # --- mesure de [Cl⁻] sur les sommes amincies
    pB, pY, p2 = cal["bleu"]['p'], cal["jaune"]['p'], cal["deux"]['p']
    ia, ib = I["bleu"]["al"]['p'], I["jaune"]["be"]['p']
    a2, b2, r2 = I["deux"]["al"]['p'], I["deux"]["be"]['p'], I["deux"]["rho"]['p']
    Kmin = min(x[2] for x in (a2, b2, r2))
    lignes = []
    for ic, c in enumerate(CONC):
        yB, yY = d['B'][ic][w], d['Y'][ic][w]
        q = a.N / yB.sum()
        tir = [(rng.binomial(yB.astype(np.int64), q).astype(float), rng.binomial(yY.astype(np.int64), q).astype(float))
               for _ in range(a.tirages)] + [(yB, yY)]                     # le dernier : somme complète
        # bases de chaque calibration pour cette acquisition
        b1, bD = (v[w] for v in F.bleu(pB['tDA'], pB['tD'], pB[f'sB{ic}']))
        XB, gB = np.c_[bD, b1 - bD], np.full(w.sum(), pB[f'bgB{ic}'])
        ri, yA, yDA, yD = (v[w] for v in F.jaune(pY['tDA'], pY['tA'], pY['tD'], pY[f'sY{ic}']))
        lY = pY['lam']                                     # μ = A_A [ (yA + λ yD) + β (ri − yA) + α λ (yDA − yD) ]
        XY, gY = np.c_[yA + lY * yD, ri - yA, lY * (yDA - yD)], np.full(w.sum(), pY[f'bgY{ic}'])
        b1_2, bD_2 = (v[w] for v in F.bleu(p2['tDA'], p2['tD'], p2[f'sB{ic}']))
        ri_2, yA_2, yDA_2, yD_2 = (v[w] for v in F.jaune(p2['tDA'], p2['tA'], p2['tD'], p2[f'sY{ic}']))
        l2 = p2['lam']
        gB2, gY2 = np.full(w.sum(), p2[f'bgB{ic}']), np.full(w.sum(), p2[f'bgY{ic}'])

        def mu2(cc, AB, AY, s):
            al, be = iso(cc, *a2), iso(cc, *b2)
            return np.r_[AB * (al * b1_2 + (1 - al) * bD_2) + s * gB2,
                         AY * (be * ri_2 + (1 - be) * yA_2) + l2 * AB * (al * yDA_2 + (1 - al) * yD_2) + s * gY2]

        def lire2(tB, tY, s, avec_rho):
            y = np.r_[tB, tY]
            AB0 = max(tB.sum() - s * gB2.sum(), 1.0) / (0.5 * (b1_2 + bD_2)).sum()
            AY0 = max(tY.sum() - s * gY2.sum(), 1.0) / yA_2.sum()
            if avec_rho:
                f = lambda th: dev_resid(y, np.maximum(mu2(th[0], np.exp(th[1]), iso(th[0], *r2) * np.exp(th[1]), s), 1e-9))
                th0, lo, hi = [20.0, np.log(AB0)], [-0.9 * Kmin, -np.inf], [3000.0, np.inf]
            else:
                f = lambda th: dev_resid(y, np.maximum(mu2(th[0], np.exp(th[1]), np.exp(th[2]), s), 1e-9))
                th0, lo, hi = [20.0, np.log(AB0), np.log(AY0)], [-0.9 * Kmin, -np.inf, -np.inf], [3000.0, np.inf, np.inf]
            return least_squares(f, th0, bounds=(lo, hi), x_scale='jac', max_nfev=400).x[0]

        lus = {k: [] for k in NOM}
        for j, (tB, tY) in enumerate(tir):
            s = 1.0 if j == len(tir) - 1 else q
            # un canal : μ = θ₁ u + θ₂ v + fond connu ; α (ou β) = θ₂ / θ₁
            thB = _lin_fond(tB, XB, s * gB)
            thY = _lin_fond(tY, XY, s * gY)
            lus["bleu"].append(inv_iso(thB[1] / thB[0], ia))
            lus["jaune"].append(inv_iso(thY[1] / thY[0], ib))
            lus["ab"].append(lire2(tB, tY, s, False))
            lus["abr"].append(lire2(tB, tY, s, True))
        # bornes de Cramér-Rao (information de Fisher, fond compris), au point de l'isotherme
        cr = {}
        AB_, AY_ = pB[f'AB{ic}'] * q, pY[f'AY{ic}'] * q
        al_, be_ = iso(c, *ia), iso(c, *ib)
        XYf = np.c_[XY[:, 0] + pY[f'al{ic}'] * XY[:, 2], XY[:, 1]]
        for k, (X, A_, x_, g_, dx) in {"bleu": (XB, AB_, al_, q * gB, diso(c, ia)), "jaune": (XYf, AY_, be_, q * gY, diso(c, ib))}.items():
            m = A_ * (X[:, 0] + x_ * X[:, 1]) + g_
            D = np.c_[A_ * X[:, 1], X[:, 0] + x_ * X[:, 1]]
            cr[k] = np.sqrt(np.linalg.inv((D.T / m) @ D)[0, 0]) / abs(dx)
        AB2, AY2 = p2[f'AB{ic}'] * q, p2[f'AY{ic}'] * q
        for k, avec_rho in (("ab", False), ("abr", True)):
            h = 1e-3 * max(c, 1.0)
            if avec_rho:
                g = lambda th: mu2(th[0], th[1], iso(th[0], *r2) * th[1], q)
                th0 = np.array([c, AB2])
            else:
                g = lambda th: mu2(th[0], th[1], th[2], q)
                th0 = np.array([c, AB2, AY2])
            m = g(th0)
            D = np.array([(g(th0 + e) - g(th0 - e)) / (2 * e[i]) for i, e in enumerate(np.diag([h] + [1e-4 * x for x in th0[1:]]))]).T
            cr[k] = np.sqrt(np.linalg.inv((D.T / m) @ D)[0, 0])
        for k in NOM:
            v = np.asarray(lus[k][:-1], float)
            v = v[np.isfinite(v)]
            sd = (np.percentile(v, 75) - np.percentile(v, 25)) / 1.349
            bm = float(np.median(v) - c)
            lignes.append(dict(cl_mM=c, lecture=NOM[k], N_bleu=a.N, ecart_type_robuste=sd, biais_median=bm,
                               erreur_totale=float(np.hypot(sd, bm)), cramer_rao=cr[k], ecart_isotherme=lus[k][-1] - c))
        t = pd.DataFrame(lignes[-4:]).set_index('lecture')
        print(f"  {c:5.1f} mM : " + " · ".join(f"{k} σ {t.loc[NOM[k], 'ecart_type_robuste']:.2f}, biais {t.loc[NOM[k], 'biais_median']:+.2f}, "
                                              f"erreur {t.loc[NOM[k], 'erreur_totale']:.2f}" for k in NOM))
    res = pd.DataFrame(lignes)
    res.to_csv(os.path.join(ici, f'calibration_1_vs_2_canaux_{LAMBDA}nm.csv'), index=False)
    figure(ici, cal, I, res, a)


def _lin_fond(y, X, fond, it=60):
    """Poisson, lien identité, fond connu : μ = X θ + fond (θ de signe libre)."""
    th = np.linalg.lstsq(X, np.maximum(y - fond, 0.5), rcond=None)[0]
    for _ in range(it):
        mu = np.maximum(X @ th + fond, 1e-9)
        wX = X / mu[:, None]
        step = np.linalg.solve(X.T @ wX, wX.T @ (y - mu))
        a = 1.0
        while np.any(X @ (th + a * step) + fond <= 0) and a > 1e-4:
            a /= 2
        th = th + a * step
        if np.max(np.abs(a * step)) < 1e-10 * max(1.0, np.max(np.abs(th))):
            break
    return th


def figure(ici, cal, I, res, a):
    virg = lambda x: f'{x:g}'.replace('.', ',')
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.8))
    cc = np.linspace(0, 131, 400)

    def courbe(ax, s, col, nom):
        ax.errorbar(CONC, s['x'], yerr=s['se'], fmt='o', color=col, ms=6, capsize=3,
                    label=nom + (f" : K = {virg(round(s['p'][2], 1))} mM" if s['p'] is not None else ""))
        if s['p'] is not None:
            ax.plot(cc, iso(cc, *s['p']), color=col, lw=1.6)

    def axe_c(ax):
        ax.set_xscale('symlog', linthresh=2.5, linscale=0.6)
        ax.set_xticks(CONC); ax.set_xticklabels([virg(c) for c in CONC])
        ax.set_xlabel('[Cl⁻] imposée (mM)'); ax.grid(True, alpha=0.3)

    ax = axes[0, 0]
    courbe(ax, I["bleu"]["al"], COUL["bleu"], "bleu seul")
    courbe(ax, I["deux"]["al"], COUL["abr"], "2 canaux")
    axe_c(ax); ax.set_ylabel('α'); ax.legend(frameon=False, fontsize=9)
    ax.set_title('(a) α, canal bleu', loc='left', fontsize=12)
    ax = axes[0, 1]
    courbe(ax, I["deux"]["be"], COUL["abr"], "2 canaux")
    courbe(ax, I["jaune"]["be"], COUL["jaune"], "jaune seul")          # par-dessus : ses points à β = 0 restent visibles
    axe_c(ax); ax.set_ylabel('β'); ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(-0.1, 1.1)                       # le jaune seul a des erreurs-types de ±16 : β n'y est pas déterminé
    ax.set_title('(b) β, canal jaune', loc='left', fontsize=12)
    ax = axes[1, 0]
    courbe(ax, I["deux"]["rho"], COUL["abr"], "2 canaux")
    axe_c(ax); ax.set_ylabel('amplitude du jaune / amplitude du bleu'); ax.legend(frameon=False, fontsize=9)
    ax.set_title('(c) Intensité du jaune relative au bleu', loc='left', fontsize=12)
    ax = axes[1, 1]
    for k, nom in NOM.items():
        s = res[res.lecture == nom]
        ax.plot(s.cl_mM, s.erreur_totale, 'o-', color=COUL[k], ms=6, lw=1.8, label=nom)
    axe_c(ax); ax.set_yscale('log')
    lo, hi = np.nanmin(res.erreur_totale.values), np.nanmax(res.erreur_totale.values)
    yt = [v for v in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500) if lo / 1.5 <= v <= hi * 1.5]
    ax.set_yticks(yt); ax.set_yticklabels([virg(v) for v in yt]); ax.minorticks_off()
    ax.set_ylabel('erreur sur la [Cl⁻] lue (mM)')
    ax.set_title(f'(d) Erreur à 10⁴ photons bleus (bruit et biais)', loc='left', fontsize=12)
    ax.legend(frameon=False, fontsize=9)
    p = {k: R['p'] for k, R in cal.items()}
    fig.suptitle('Calibration de SCLM à 850 nm : un canal ou deux (avec la fuite λ du donneur dans le jaune)', x=0.01, ha='left', fontsize=14)
    fig.text(0.01, 0.005,
             "d_D = α e^(−t/τ_DA) + (1 − α) e^(−t/τ_D) ;  d_A = β [e^(−t/τ_DA) ⊛ e^(−t/τ_A)] + (1 − β) e^(−t/τ_A) ;  "
             "canal jaune = A_A d_A + λ × (donneur) d_D.\n"
             f"τ (ns) — bleu seul : τ_DA {virg(round(p['bleu']['tDA'], 3))}, τ_D {virg(round(p['bleu']['tD'], 3))} ; "
             f"jaune seul : τ_DA {virg(round(p['jaune']['tDA'], 3))}, τ_D {virg(round(p['jaune']['tD'], 3))}, τ_A {virg(round(p['jaune']['tA'], 3))}, λ {virg(round(p['jaune']['lam'], 3))} ; "
             f"2 canaux : τ_DA {virg(round(p['deux']['tDA'], 3))}, τ_D {virg(round(p['deux']['tD'], 3))}, τ_A {virg(round(p['deux']['tA'], 3))}, λ {virg(round(p['deux']['lam'], 3))}.\n"
             "(a–c) variable de chaque concentration ± erreur-type, et isotherme ; (d) sommes réelles amincies, "
             f"{a.tirages} lectures : √(écart-type robuste² + biais médian²).", fontsize=8.6, color='#52514e')
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(os.path.join(ici, 'calibration_1_vs_2_canaux.png'), dpi=150)
    print('Figure écrite : calibration_1_vs_2_canaux.png')


if __name__ == '__main__':
    main()
