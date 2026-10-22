"""figures.py — rendu après simulation (matplotlib, sans affichage)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLEU, ORANGE, GRIS, NOIR, VERT = '#2a6fdb', '#e8871e', '#9a9a9a', '#111111', '#2e9e5b'


def donnees(S, mes, cal, p, rng, chemin):
    """Les données d'origine : histogrammes mesurés contre le modèle calibré, et une trame simulée de chaque sorte."""
    from .simulation import trame_modele
    from .modele import Senseur
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    t = (np.arange(S.n_fen) + cal['fenetre'][0] + 0.5) * cal['periode_ns'] / cal['n_cases']
    dm = cal['decalage_moyen_ps']
    for j, (canal, H) in enumerate((('bleu', mes.B), ('jaune', mes.Y))):
        a, r = ax[j, 0], ax[j, 1]
        for c, coul in ((0.0, BLEU), (20.0, VERT), (131.0, ORANGE)):
            i = mes.indice(c); n = H[i]
            ti = cal['titration'][i]
            # chaque acquisition avec son propre décalage d'horloge ajusté (section 3)
            Si = Senseur(cal, decalage_ps=(ti['decalage_bleu_ps'] - dm['bleu'], ti['decalage_jaune_ps'] - dm['jaune']))
            muB, muY = Si.mu(c, 1.0); m = muB if j == 0 else muY
            fond = ti['fond_bleu_par_case'] if j == 0 else ti['fond_jaune_par_case']
            m = m * (n.sum() - fond * len(n)) / m.sum() + fond
            a.semilogy(t, n, '.', ms=2.5, color=coul, label=f'{c:g} mM mesuré')
            a.semilogy(t, m, '-', lw=0.8, color=coul)
            r.plot(t, (n - m) / np.sqrt(m), '-', lw=0.6, color=coul, label=f'{c:g} mM')
        a.set_title(f"canal {canal} : titration {p['donnees']['nom']} (points), modèle (traits)")
        a.set_xlabel('t (ns)'); a.set_ylabel('photons par case'); a.legend(fontsize=8)
        r.set_title(f'canal {canal} : résidus (n − μ)/√μ, modèle aux isothermes', fontsize=10)
        r.set_xlabel('t (ns)'); r.axhline(0, color=NOIR, lw=0.5); r.legend(fontsize=8)
        s = ax[j, 2]
        Bm, Ym = mes.trame(10.0, p['exemple']['photons'], rng)
        Bs, Ys = trame_modele(S, 10.0, p['exemple']['photons'], rng)
        s.step(t, Bm if j == 0 else Ym, color=GRIS, lw=0.8, label='tirée des mesures (10 mM)')
        s.step(t, Bs if j == 0 else Ys, color=BLEU if j == 0 else ORANGE, lw=0.8, label='tirée du modèle (10 mM)')
        s.set_yscale('symlog', linthresh=1); s.set_xlabel('t (ns)')
        s.set_title(f'une trame de {p["exemple"]["photons"]} photons, canal {canal}'); s.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(chemin, dpi=130); plt.close(fig)


def exemple(D, p, chemin):
    """Section 15 : panneau a, [Cl⁻] ; panneau b, η."""
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    k = D.trame
    a = ax[0]
    a.plot(k, np.clip(D.c_flim, 0, 80), '.', ms=3, color=GRIS, label='durée de vie seule, trame par trame')
    a.plot(k, D.c_ratio_eta1, '-', lw=1, color=ORANGE, label='ratio lu avec la calibration (η = 1)')
    a.fill_between(k, D.lisse_bas, D.lisse_haut, color=BLEU, alpha=0.2, lw=0, label='lissé, intervalle à 95 %')
    a.plot(k, D.c_filtre, '-', lw=0.8, color=VERT, label='filtre hybride, temps réel')
    a.plot(k, D.c_lisse, '-', lw=1.4, color=BLEU, label='filtre hybride, lissé')
    a.plot(k, D.c_vrai, '-', lw=1.2, color=NOIR, label='vraie valeur')
    a.set_ylim(0, 45); a.set_ylabel('[Cl⁻] (mM)'); a.legend(fontsize=8, ncol=2, loc='upper right')
    e = p['exemple']
    a.set_title(f"a · un soma, η = {e['eta']} inconnu du filtre, {e['photons']} photons par trame, A −30 % (photoblanchiment)")
    b = ax[1]
    b.fill_between(k, D.eta_filtre * np.exp(-2 * D.eta_filtre_sd), D.eta_filtre * np.exp(2 * D.eta_filtre_sd),
                   color=ORANGE, alpha=0.2, lw=0)
    b.plot(k, D.eta_filtre, color=ORANGE, label='temps réel (± 2 sd)')
    b.plot(k, D.eta_lisse, color=BLEU, label='après coup')
    b.axhline(e['eta'], color=NOIR, lw=0.8, ls='--', label='vrai')
    b.set_ylim(0.95, 1.6); b.set_ylabel('η'); b.set_xlabel('trame'); b.legend(fontsize=8); b.set_title('b · η estimé')
    fig.tight_layout(); fig.savefig(chemin, dpi=130); plt.close(fig)


def cramer_rao(lignes, nom, chemin):
    D = pd.DataFrame(lignes)
    fig, a = plt.subplots(figsize=(8, 5))
    for col, coul, ls in (('forme jaune seule', GRIS, ':'), ('formes 2 canaux (eta inconnu)', ORANGE, '-'),
                          ('ratio (eta connu)', BLEU, '-'), ('ratio + formes (eta connu)', VERT, '--')):
        a.loglog(D.c_mM, D[col], ls, color=coul, marker='o', ms=3, label=col)
    a.set_xlabel('[Cl⁻] (mM)'); a.set_ylabel('σ_c minimal (mM)')
    a.set_title(f'borne de Cramér-Rao, modèle calibré sur {nom}, 10 000 photons')
    a.grid(True, which='both', alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(chemin, dpi=130); plt.close(fig)


def temps_reel(D, chemin):
    fig, a = plt.subplots(figsize=(8, 5))
    for (nom, x), coul in zip(D.groupby('variante', sort=False), (ORANGE, BLEU, NOIR)):
        a.plot(x.photons, x.rmse_mM, 'o-', color=coul, label=nom)
    a.set_xscale('log'); a.set_xlabel('photons par trame et par ROI'); a.set_ylabel('erreur typique après la trame 100 (mM)')
    a.set_title('filtre temps réel, η = 1,25 inconnu, [Cl⁻] propre à chaque ROI'); a.grid(True, alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(chemin, dpi=130); plt.close(fig)


def amincissement(trace, D, nom, chemin):
    fig, ax = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw=dict(width_ratios=[2, 1]))
    a = ax[0]; k = trace['trame']
    a.plot(k, np.clip(trace['durée de vie seule'], 0, 120), '.', ms=3, color=GRIS, label='durée de vie seule')
    a.plot(k, trace['ratio (η = 1)'], '-', lw=1, color=ORANGE, label='ratio (η = 1)')
    a.plot(k, trace['filtre, fond ignoré'], '-', lw=1, color=VERT, label='filtre, fond ignoré')
    a.plot(k, trace['filtre, fond connu'], '-', lw=1.4, color=BLEU, label='filtre, fond connu')
    a.plot(k, trace['c_vrai'], '-', color=NOIR, lw=1.2, label='vraie valeur (titration)')
    a.set_ylim(0, 80); a.set_xlabel('trame'); a.set_ylabel('[Cl⁻] (mM)'); a.legend(fontsize=8)
    a.set_title(f'trames tirées des histogrammes mesurés ({nom}, 2026-10-16), η = 1,25 inconnu')
    b = ax[1]
    meth = [m for m in D.methode.unique() if m != 'durée de vie seule']
    pal = sorted(D.palier_mM.unique()); w = 0.8 / len(meth)
    for i, (m, coul) in enumerate(zip(meth, (ORANGE, BLEU, VERT))):
        x = D[D.methode == m].set_index('palier_mM').loc[pal]
        b.bar(np.arange(len(pal)) + i * w, x.biais_mM, w, yerr=x.sd_mM, color=coul, label=m, capsize=2)
    b.set_xticks(np.arange(len(pal)) + w); b.set_xticklabels([f'{c:g} mM' for c in pal]); b.axhline(0, color=NOIR, lw=0.5)
    b.set_ylabel('biais (mM), barres : sd par trame'); b.set_title('biais par palier (moyenne des répétitions)'); b.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(chemin, dpi=130); plt.close(fig)
