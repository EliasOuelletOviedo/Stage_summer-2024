#!/usr/bin/env python3
"""
histogrammes.py — Décomposition des histogrammes de calibration SuperClomeleon en composantes.

Données : titration KIRIN du 2026-10-16 (somme des 10 acquisitions de chaque concentration et IRF,
irf_bruts_KIRIN_<λ>nm.csv, produit par step18_bruts_KIRIN.jl), lue par défaut dans le dossier du script.
Le chargement (sclm_hybride.load) ne garde qu'une période laser des 256 canaux du TAC.

Grille 4 × 2 :
  colonne de gauche  (canal bleu)  : τ1 (FRET rapide), τ2 (FRET lent), τ_D (donneurs sans FRET), somme
  colonne de droite  (canal jaune) : excitation directe, émission sensibilisée (FRET), fuite du donneur, somme
Chaque panneau superpose les 8 concentrations. Dans les panneaux « somme », le trait est le modèle complet
(composantes + fond) et les points sont les photons mesurés.

Usage :
    python histogrammes.py                          # KIRIN 800 nm, échelle log
    python histogrammes.py --lin                    # échelle linéaire
    python histogrammes.py autre.csv --out decomposition.png

Nécessite sclm_hybride.py dans le même dossier (ajustement conjoint, ~1 min).
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # pour trouver sclm_hybride.py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

from sclm_hybride import load, fit_irf, fit_joint, Fine, CONC


def composantes(d, irfp, p):
    """Recalcule, pour chaque concentration, les composantes du modèle conjoint ajusté (comptes par canal)."""
    E = Fine(d['T'], d['NB'])
    taus = [p['t1'], p['t2']]
    Gs = [p['G'], p['G'] * p.get('G2r', 1.0)]
    kT = [1 / t - 1 / p['tD'] for t in taus]
    out = []
    for c in range(len(CONC)):
        f1, f2 = p[f'f1_{c}'], p[f'f2_{c}']
        fD = max(1 - f1 - f2, 0.0)
        A = p[f'A{c}']
        # canal bleu
        E.set_irf(E.irf(irfp[1]['mu'] + p[f's1_{c}'], irfp[1]['sig'], irfp[1]['weights'], irfp[1]['taus']))
        b1 = A * f1 * E.conv_exp(taus[0])
        b2 = A * f2 * E.conv_exp(taus[1])
        bD = A * fD * E.conv_exp(p['tD'])
        bgB = p[f'bgB{c}']
        # canal jaune
        E.set_irf(E.irf(irfp[2]['mu'] + p[f's2_{c}'], irfp[2]['sig'], irfp[2]['weights'], irfp[2]['taus']))
        y1, y2, yD, yA = (E.conv_exp(t) for t in (taus[0], taus[1], p['tD'], p['tA']))
        fuite = p['lam'] * A * (f1 * y1 + f2 * y2 + fD * yD)
        rises = [t * p['tA'] / (p['tA'] - t) * (yA - y) for t, y in zip(taus, (y1, y2))]
        fret = A * sum(G_ * f * k * r for G_, f, k, r in zip(Gs, (f1, f2), kT, rises))
        directe = A * (p['k1'] * (f1 + f2) + p['k0']) * yA
        bgY = p[f'bgY{c}']
        out.append(dict(b1=b1, b2=b2, bD=bD, bsum=b1 + b2 + bD + bgB,
                        directe=directe, fret=fret, fuite=fuite, ysum=directe + fret + fuite + bgY))
    return out


def figure_contributions(comp, w, p, out, comptes=False):
    """Part de chaque composante dans le signal de son canal (fond exclu), sommée sur la fenêtre, contre [Cl⁻].
    comptes=True : photons de chaque composante dans la fenêtre, sans normalisation."""
    canaux = (('Canal bleu (1)', (('b1', f"τ₁ = {p['t1']:.3f} ns (FRET rapide)"), ('b2', f"τ₂ = {p['t2']:.3f} ns (FRET lent)"),
                                           ('bD', f"τ_D = {p['tD']:.3f} ns (sans FRET)"))),
              ('Canal jaune (2)', (('directe', 'Excitation directe'), ('fret', 'Émission sensibilisée (FRET)'),
                                       ('fuite', f"Fuite du donneur (λ = {p['lam']:.2f})"))))
    lignes = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=not comptes)
    for ax, (titre, cles) in zip(axes, canaux):
        tot = np.array([sum(comp[c][k][w].sum() for k, _ in cles) for c in range(len(CONC))])
        for (k, nom), couleur in zip(cles, ('#2a78d6', '#e07b39', '#3a9e5c')):
            ph = np.array([comp[c][k][w].sum() for c in range(len(CONC))])
            ax.plot(CONC, ph if comptes else 100 * ph / tot, 'o-', color=couleur, lw=1.8, ms=5, label=nom)
            lignes += [dict(canal=titre.split()[1], composante=k, cl_mM=c, photons=ph[i], part_pct=100 * ph[i] / tot[i])
                       for i, c in enumerate(CONC)]
        ax.set_xscale('symlog', linthresh=1.0, linscale=0.5)
        ax.set_xticks(CONC); ax.set_xticklabels([f'{c:g}' for c in CONC])
        ax.set_xlabel('[Cl⁻] (mM), échelle log (0 à part)'); ax.set_title(titre, fontsize=12, loc='left', fontweight='bold')
        ax.grid(True, alpha=0.3); ax.legend(fontsize=9)
    axes[0].set_ylabel('photons dans la fenêtre (somme des acquisitions)' if comptes else 'part du signal du canal (%), fond exclu')
    fig.tight_layout(); fig.savefig(out, dpi=150)
    if not comptes:
        import pandas as pd
        pd.DataFrame(lignes).to_csv(os.path.splitext(out)[0] + '.csv', index=False)
    print(f"Figure écrite : {out}")


def main():
    ap = argparse.ArgumentParser()
    ici = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument('csv', nargs='?', default=os.path.join(ici, 'irf_bruts_KIRIN_800nm.csv'),
                    help='sommes et IRF au format de step18_bruts_KIRIN.jl (défaut : KIRIN 800 nm)')
    ap.add_argument('--out', default=os.path.join(ici, 'decomposition_histogrammes.png'))
    ap.add_argument('--out-contrib', default=os.path.join(ici, 'contributions_concentration.png'))
    ap.add_argument('--lin', action='store_true', help='axe vertical linéaire (log par défaut)')
    a = ap.parse_args()

    d = load(a.csv)
    irfp = {1: fit_irf(d, 1, 2), 2: fit_irf(d, 2, 3)}
    R = fit_joint(d, irfp, G2=True)
    p = R['p']
    comp = composantes(d, irfp, p)

    t = d['t']
    w = d['win'][1] & d['win'][2]
    couleurs = cm.viridis(np.linspace(0, 0.92, len(CONC)))
    labels = [f'{c:g} mM' for c in CONC]

    lignes = [
        (('b1', f"τ₁ = {p['t1']:.3f} ns (FRET rapide)"), ('directe', 'Excitation directe de l’accepteur')),
        (('b2', f"τ₂ = {p['t2']:.3f} ns (FRET lent)"), ('fret', 'Émission sensibilisée (FRET)')),
        (('bD', f"τ_D = {p['tD']:.3f} ns (donneurs sans FRET)"), ('fuite', f"Fuite du donneur (λ = {p['lam']:.2f})")),
        (('bsum', 'Somme (modèle + fond) et données'), ('ysum', 'Somme (modèle + fond) et données')),
    ]

    fig, axes = plt.subplots(4, 2, figsize=(13, 15), sharex=True)
    for i, ((kb, tb), (ky, ty)) in enumerate(lignes):
        for j, (k, titre, data) in enumerate(((kb, tb, d['B']), (ky, ty, d['Y']))):
            ax = axes[i, j]
            for c in range(len(CONC)):
                y = comp[c][k]
                ax.plot(t[w], y[w], color=couleurs[c], lw=1.5, label=labels[c])
                if k in ('bsum', 'ysum'):
                    ax.plot(t[w], data[c][w], 'o', ms=1.8, color=couleurs[c], alpha=0.5, mec='none')
            ax.set_title(titre, fontsize=11, loc='left')
            if not a.lin:
                ax.set_yscale('log')
                ymax = max(comp[c][k][w].max() for c in range(len(CONC)))
                ax.set_ylim(ymax * 1e-3, ymax * 1.6)
            ax.grid(True, alpha=0.3)
            if j == 0:
                ax.set_ylabel('photons / canal')
            if i == 3:
                ax.set_xlabel('temps (ns)')
    axes[0, 0].annotate('Canal bleu (1)', xy=(0.5, 1.22), xycoords='axes fraction',
                        ha='center', fontsize=13, fontweight='bold')
    axes[0, 1].annotate('Canal jaune (2)', xy=(0.5, 1.22), xycoords='axes fraction',
                        ha='center', fontsize=13, fontweight='bold')
    axes[0, 1].legend(title='[Cl⁻]', fontsize=8, ncol=2, loc='upper right')
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    figure_contributions(comp, w, p, a.out_contrib)
    figure_contributions(comp, w, p, os.path.splitext(a.out_contrib)[0] + '_comptes.png', comptes=True)
    print(f"Figure écrite : {a.out}  (déviance/ddl = {R['dev'] / R['dof']:.3f})")


if __name__ == '__main__':
    main()