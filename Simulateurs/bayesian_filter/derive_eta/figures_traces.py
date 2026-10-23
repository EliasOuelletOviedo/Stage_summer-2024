"""Figures des 5 variantes du test de dérive, dans le style de fig_exemple_eta (src/figures.py, exemple), à partir des
traces de traces.py (sorties/traces/) → sorties/figures/.   python3 figures_traces.py [dossier des traces]"""
import sys, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ICI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(ICI))   # bayesian_filter
from src.figures import GRIS, ORANGE, BLEU, NOIR

DOS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ICI, 'sorties', 'traces')
SOR = os.path.join(os.path.dirname(DOS), 'figures') if len(sys.argv) == 1 else os.path.join(DOS, 'figures'); os.makedirs(SOR, exist_ok=True)
VARIANTES = [v for v in ('ref', '1', '2', '12', '123') if os.path.exists(os.path.join(DOS, f'{v}.npz'))]
com = np.load(os.path.join(DOS, 'commun.npz')); C, eta_vrai, c_ratio, c_flim = com['C'], com['eta_vrai'], com['c_ratio'], com['c_flim']
M, T = C.shape; k = np.arange(1, T + 1)
SCEN = f'10 ROI, 1 000 photons par ROI et par trame, η 1 → 1,1 (inconnu), A −30 %, {T} trames'
fr = lambda x, n=2: f'{x:.{n}f}'.replace('.', ',')


def bilan(d):
    err = d['c'] - C; cov = np.mean((d['bas'] <= C) & (C <= d['haut']))
    return f"10 ROI : rms {fr(np.sqrt(np.mean(err ** 2)))} mM, biais {fr(err.mean())} mM, couverture à 95 % : {100 * cov:.0f} %"


def poids(d):
    w, et = d['poids'], str(d['etiquette_poids'])
    if len(w) == 0:
        return ''
    val = [0.0, 1e-4, 3e-4, 1e-3, 3e-3] if et == 's_η' else [0.004, 0.013, 0.04]
    return (' — poids des modèles à la dernière trame, ' + et + ' par trame : '
            + ', '.join(('0 (η constant)' if s == 0 else f'{fr(100 * s, 2 if s < 1e-3 else 1)} %') + f' : {fr(x)}'
                        for s, x in zip(val, w) if x >= 0.005))


def detail(v, d, m):
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    rms = lambda x: f" — rms {fr(np.sqrt(np.mean((x - C[m]) ** 2)))} mM"
    a = ax[0]
    a.plot(k, np.clip(c_flim[m], 0, 80), '.', ms=1, color=GRIS, label='durée de vie seule, trame par trame' + rms(c_flim[m]))
    a.plot(k, c_ratio[m], '-', lw=0.5, color=ORANGE, label='ratio lu avec la calibration (η = 1)' + rms(c_ratio[m]))
    a.fill_between(k, d['bas'][m], d['haut'][m], color=BLEU, alpha=0.2, lw=0, label='temps réel, intervalle à 95 %')
    a.plot(k, d['c'][m], '-', lw=1.2, color=BLEU, label='filtre hybride, temps réel' + rms(d['c'][m]))
    a.plot(k, C[m], '-', lw=1.2, color=NOIR, label='vraie valeur')
    a.set_ylim(0, 45); a.set_ylabel('[Cl⁻] (mM)'); a.legend(fontsize=8, loc='upper right')
    a.set_title(f"{d['titre']} — ROI {m + 1}\n{SCEN}\n{bilan(d)}", fontsize=10)
    b = ax[1]; i = m if d['eta'].shape[0] > 1 else 0; e, sd = d['eta'][i], d['eta_sd'][i]
    b.fill_between(k, e * np.exp(-2 * sd), e * np.exp(2 * sd), color=BLEU, alpha=0.2, lw=0)
    b.plot(k, e, color=BLEU, label='temps réel (± 2 sd)' + ('' if d['eta'].shape[0] > 1 else ', partagé par les 10 ROI'))
    b.plot(k, eta_vrai, color=NOIR, lw=0.8, ls='--', label='vrai')
    b.set_ylim(0.9, 1.2); b.set_ylabel('η'); b.set_xlabel('trame'); b.legend(fontsize=8, loc='upper left')
    b.set_title('b · η estimé' + (f' (ROI {m + 1})' if d['eta'].shape[0] > 1 else '') + poids(d), fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(SOR, f'fig_derive_{v}_roi{m + 1}.png'), dpi=130); plt.close(fig)


def mosaique(v, d):
    fig = plt.figure(figsize=(15, 15)); gs = fig.add_gridspec(6, 2, height_ratios=[1, 1, 1, 1, 1, 1.3])
    for m in range(M):
        a = fig.add_subplot(gs[m // 2, m % 2])
        a.plot(k, c_ratio[m], '-', lw=0.3, color=ORANGE)
        a.fill_between(k, d['bas'][m], d['haut'][m], color=BLEU, alpha=0.2, lw=0)
        a.plot(k, d['c'][m], '-', lw=0.8, color=BLEU); a.plot(k, C[m], '-', lw=1.0, color=NOIR)
        a.set_ylim(0, 45); a.text(0.01, 0.95, f"ROI {m + 1} — rms {fr(np.sqrt(np.mean((d['c'][m] - C[m]) ** 2)))} mM",
                                  transform=a.transAxes, va='top', fontsize=8)
        if m % 2 == 0: a.set_ylabel('[Cl⁻] (mM)')
        a.tick_params(labelsize=7)
    b = fig.add_subplot(gs[5, :])
    for i in range(d['eta'].shape[0]):
        b.plot(k, d['eta'][i], color=BLEU, lw=0.6 if d['eta'].shape[0] > 1 else 1.2, alpha=0.6 if d['eta'].shape[0] > 1 else 1)
    if d['eta'].shape[0] == 1:
        e, sd = d['eta'][0], d['eta_sd'][0]; b.fill_between(k, e * np.exp(-2 * sd), e * np.exp(2 * sd), color=BLEU, alpha=0.2, lw=0)
    b.plot(k, eta_vrai, color=NOIR, lw=0.8, ls='--')
    b.set_ylim(0.9, 1.2); b.set_ylabel('η'); b.set_xlabel('trame')
    b.set_title('η estimé, temps réel' + (' (une courbe par ROI)' if d['eta'].shape[0] > 1 else ' (± 2 sd), partagé') + poids(d), fontsize=9)
    fig.suptitle(f"{d['titre']}\n{SCEN} — noir : vraie valeur, bleu : filtre temps réel (95 %), orange : ratio (η = 1)\n{bilan(d)}",
                 fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(SOR, f'fig_derive_{v}_10roi.png'), dpi=110); plt.close(fig)


D = {v: dict(np.load(os.path.join(DOS, f'{v}.npz'))) for v in VARIANTES}
for v, d in D.items():
    detail(v, d, 0); mosaique(v, d)
# comparaison : η et erreur sur c des 5 variantes
fig, ax = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
couleurs = dict(ref='#999999', **{'1': '#e69f00', '2': '#56b4e9', '12': '#009e73', '123': '#0b4f9c'})
for v, d in D.items():
    nom = str(d['titre']).split(' :')[0]
    e = d['eta'].mean(0) if d['eta'].shape[0] > 1 else d['eta'][0]
    ax[0].plot(k, e, color=couleurs[v], lw=1.2, label=nom + (' (moyenne des 10 η par ROI)' if d['eta'].shape[0] > 1 else ''))
    err = np.sqrt(np.mean((d['c'] - C) ** 2, 0)); n_l = max(1, min(50, T // 5)); lis = np.convolve(err, np.ones(n_l) / n_l, mode="same")
    ax[1].plot(k, lis, color=couleurs[v], lw=1.2, label=f"{nom} — rms {fr(np.sqrt(np.mean((d['c'] - C) ** 2)))} mM")
ax[0].plot(k, eta_vrai, color=NOIR, ls='--', lw=0.8, label='vrai'); ax[0].set_ylim(0.9, 1.2); ax[0].set_ylabel('η')
ax[0].legend(fontsize=8, loc='upper left'); ax[0].set_title(f'a · η estimé en temps réel — {SCEN}', fontsize=10)
ax[1].set_ylabel('erreur rms sur c, 10 ROI (mM)\n(moyenne glissante de 50 trames)'); ax[1].set_xlabel('trame'); ax[1].legend(fontsize=8)
ax[1].set_title('b · erreur sur [Cl⁻]', fontsize=10); ax[1].set_ylim(0, None)
fig.tight_layout(); fig.savefig(os.path.join(SOR, 'fig_derive_comparaison.png'), dpi=130); plt.close(fig)
print(sorted(os.listdir(SOR)))
