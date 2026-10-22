#!/usr/bin/env python3
"""
step0_calibration.py — le pont unique avec les données de titration du 2026-10-16.

Ajuste la calibration des sections 3 à 6 du document sur les 16 histogrammes sommés (8 [Cl⁻] × 2 canaux) d'un
jeu de données (SCLM_mCer3 par défaut, extrait par step17_bruts_mCer3.jl ; SCLM régulier avec params_SCLM.toml,
extrait par step16_declins_purs.jl), avec le module de calibration qui accompagne ces données (sclm_hybride.py :
IRF paramétrique, ajustement conjoint de Poisson, isothermes ; seul son lecteur de CSV change d'un jeu à l'autre).
Écrit la calibration (calibration.json par défaut), que lisent le modèle direct, le simulateur et le filtre (src/). Rien d'autre dans ce dossier ne lit les données
brutes, sauf le générateur « amincissement » de src/simulation.py, qui tire des trames des histogrammes mesurés.

Lancer depuis ce dossier :  python3 step0_calibration.py  [--params params_SCLM.toml]        (~1 min)
"""
import importlib.util, json, os, sys, time
import numpy as np

ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ICI)
from src.modele import charger_params, chemin_params


def module_calibration(p):
    sys.dont_write_bytecode = True                    # ne pas réécrire le __pycache__ versionné du dossier des données
    chemin = os.path.join(ICI, p['donnees']['dossier'], p['donnees']['module_calibration'])
    spec = importlib.util.spec_from_file_location('sclm_hybride', chemin)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def calibrer(p, csv, sortie):
    """Ajuste la calibration sur les histogrammes de `csv` (format du jeu de `p`) et l'écrit dans `sortie` (json)."""
    t0 = time.time()
    sh = module_calibration(p)
    d = sh.load(csv)
    irfp = {ch: sh.fit_irf(d, ch, ntail=2 if ch == 1 else 3) for ch in (1, 2)}
    R = sh.fit_joint(d, irfp, G2=True)
    iso = sh.fit_isotherms(R)
    q, names = R['p'], R['P'].names
    se = dict(zip(names, np.sqrt(np.maximum(np.diag(R['cov']), 0))))
    cles = ['t1', 't2', 'tD', 'tA', 'lam', 'G', 'G2r', 'k1', 'k0']
    NC, w = len(sh.CONC), d['win'][1]
    assert np.array_equal(d['win'][1], d['win'][2]), "fenêtres bleue et jaune différentes"
    # composition de chaque histogramme de la titration : photons du senseur et du fond dans la fenêtre
    mB, mY = R['model_fn'](R['x'])
    titration = []
    for c in range(NC):
        titration.append(dict(
            c_mM=float(sh.CONC[c]), A=float(q[f'A{c}']), f1=float(q[f'f1_{c}']), f2=float(q[f'f2_{c}']),
            decalage_bleu_ps=float(1e3 * q[f's1_{c}']), decalage_jaune_ps=float(1e3 * q[f's2_{c}']),
            fond_bleu_par_case=float(q[f'bgB{c}']), fond_jaune_par_case=float(q[f'bgY{c}']),
            photons_bleus=float(d['B'][c][w].sum()), photons_jaunes=float(d['Y'][c][w].sum()),
            part_fond_bleu=float(q[f'bgB{c}'] * w.sum() / mB[c][w].sum()),
            part_fond_jaune=float(q[f'bgY{c}'] * w.sum() / mY[c][w].sum())))
    out = dict(
        jeu=p['donnees']['nom'], source=os.path.relpath(csv, ICI),
        periode_ns=float(d['T']), n_cases=int(d['NB']), fenetre=[int(i) for i in np.where(w)[0][[0, -1]]],
        constantes={k: float(q[k]) for k in cles}, erreurs_types={k: float(se[k]) for k in cles},
        isothermes={k: dict(fmax=float(v['fmax']), fmin=float(v['fmin']), K=float(v['K']),
                            valeurs=[float(x) for x in v['values']]) for k, v in iso.items()},
        irf={str(ch): dict(mu=float(v['mu']), sig=float(v['sig']), poids=[float(x) for x in v['weights']],
                           taus=[float(x) for x in v['taus']], dev_par_ddl=float(v['dev_per_dof'])) for ch, v in irfp.items()},
        decalage_moyen_ps=dict(bleu=float(np.mean([t['decalage_bleu_ps'] for t in titration])),
                               jaune=float(np.mean([t['decalage_jaune_ps'] for t in titration]))),
        deviance=float(R['dev']), ddl=int(R['dof']), titration=titration,
        note="η = 1 par définition dans les conditions de la calibration ; λ, G et κ sont relatifs à ces conditions.")
    with open(sortie, 'w') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    k = out['constantes']
    print(f"déviance/ddl = {R['dev']:.0f}/{R['dof']} = {R['dev'] / R['dof']:.3f}")
    print("τ1 {t1:.4f}  τ2 {t2:.4f}  τD {tD:.4f}  τA {tA:.4f} ns | λ {lam:.3f}  G {G:.3f}  g2 {G2r:.3f}  κ {k1:.3f}".format(**k))
    for n, v in out['isothermes'].items():
        print(f"{n:2s} : max {v['fmax']:.3f}  min {v['fmin']:.3f}  K {v['K']:.1f} mM")
    print("part du fond dans la fenêtre (bleu / jaune, %) :",
          " ".join(f"{t['c_mM']:g}:{100 * t['part_fond_bleu']:.1f}/{100 * t['part_fond_jaune']:.1f}" for t in titration))
    print(f"{os.path.relpath(sortie, ICI)} écrit ({time.time() - t0:.0f} s)")
    return out


def main():
    p = charger_params(chemin_params(sys.argv))
    calibrer(p, os.path.join(ICI, p['donnees']['dossier'], p['donnees']['titration']),
             os.path.join(ICI, p['donnees']['calibration']))


if __name__ == '__main__':
    main()
