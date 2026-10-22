"""
exemple.py — une course de la section 15 : un soma, trames simulées le long de la trajectoire de l'exemple, lues par le
ratio (η = 1), la durée de vie seule, et le filtre hybride (temps réel et lissé) sur la grille (c, η) complète.
Sert à run.py (partie 3) et à exemples.py (la même course sous d'autres conditions).

η peut évoluer d'une trame à l'autre : ln η_k = ln η_{k−1} + N(0, s_η²). s_η est inconnu ; chaque candidat de
[filtre] s_eta_candidats (0 = η constant) donne un modèle, et les modèles sont moyennés selon leur vraisemblance
marginale — celle des trames 1..k pour le filtre temps réel (causal), celle de toute la course pour le lisseur.
"""
import time
import numpy as np
import pandas as pd

from .filtre import Grille, grille_eta, transition, filtre_lisseur, moyenne_de_modeles, resume_c, resume_eta
from .simulation import trame_modele, trajectoire_exemple

METHODES_EXEMPLE = (('ratio avec la calibration (η = 1)', 'c_ratio_eta1'), ('durée de vie seule, trame par trame', 'c_flim'),
                    ('filtre hybride, temps réel', 'c_filtre'), ('filtre hybride, lissé après coup', 'c_lisse'))


def grille_lisseur(S, g, f):
    """La grille de c de g, avec la grille de η plus fine du lisseur ([filtre] n_eta_lisseur)."""
    return Grille(S, g.c, grille_eta(dict(f, n_eta=f['n_eta_lisseur'])))


def etirer(e, facteur):
    """La même trajectoire sur `facteur` fois plus de trames : montée, descente et constantes de temps à l'échelle."""
    e = dict(e)
    e['trames'] = facteur * e['trames']
    e['trame_montee'] = facteur * (e['trame_montee'] - 1) + 1
    e['trame_descente'] = facteur * (e['trame_descente'] - 1) + 1
    e['tau_montee_trames'] *= facteur; e['tau_descente_trames'] *= facteur
    return e


def course_exemple(S, g, f, e, rng, eta=None, lisser=True):
    """
    Une course : e est le bloc [exemple] (éventuellement étiré), f le bloc [filtre] ; eta, s'il est donné, est le η vrai
    de chaque trame (sinon e['eta'] partout). lisser=False : filtre temps réel seul, pas de colonnes « lisse ».
    Renvoie (table par trame, histogrammes, durée du filtre et du lisseur).
    """
    T = e['trames']
    C = trajectoire_exemple(e)
    eta = np.full(T, e['eta']) if eta is None else np.asarray(eta, float)
    A = 1 + e['derive_A'] * np.arange(T) / (T - 1)            # photoblanchiment : le filtre ne le connaît pas
    H = [trame_modele(S, C[k], e['photons'] * A[k], rng, eta[k]) for k in range(T)]
    t0 = time.perf_counter()
    gl = grille_lisseur(S, g, f)
    formes = np.array([gl.formes(B, Y) for B, Y in H])                   # (T, n_c) : le morceau « quand »
    NB = np.array([B.sum() for B, _ in H]); NY = np.array([Y.sum() for _, Y in H])
    lB, lY = gl.logs_ratio()                                               # le morceau « combien » : NB·lB + NY·lY
    loglik = lambda k: formes[k][:, None] + NB[k] * lB + NY[k] * lY
    Kc = transition(np.log(gl.c), f['s'], f['s_saut'], f['w'])
    courses = [filtre_lisseur(loglik, T, Kc, transition(np.log(gl.eta), s), lisser) for s in f['s_eta_candidats']]
    pcf, pef, wf = moyenne_de_modeles(courses, causal=True)
    fm, flo, fhi = resume_c(pcf, gl.c); ef, ef_sd = resume_eta(pef, gl.eta)
    D = pd.DataFrame(dict(trame=np.arange(1, T + 1), c_vrai=C, eta_vrai=eta, N_B=NB, N_Y=NY, R=NY / NB,
                          c_ratio_eta1=[g.lire_ratio(B, Y) for B, Y in H], c_flim=[g.lire_flim(B, Y) for B, Y in H],
                          c_filtre=fm, filtre_bas=flo, filtre_haut=fhi, eta_filtre=ef, eta_filtre_sd=ef_sd))
    if lisser:
        pcl, pel, wl = moyenne_de_modeles(courses, causal=False)
        lm, llo, lhi = resume_c(pcl, gl.c); el, el_sd = resume_eta(pel, gl.eta)
        D = D.assign(c_lisse=lm, lisse_bas=llo, lisse_haut=lhi, eta_lisse=el, eta_lisse_sd=el_sd)
        D.attrs['poids_lisseur'] = dict(zip(f['s_eta_candidats'], wl[0]))
    duree = time.perf_counter() - t0
    for i, s in enumerate(f['s_eta_candidats']):                           # poids de chaque modèle de dynamique de η
        D[f'poids_filtre_s_eta_{s:g}'] = wf[:, i]
    D.attrs['poids_filtre_final'] = dict(zip(f['s_eta_candidats'], wf[-1]))
    D.attrs['ln_evidence'] = {s: float(r['ln_evidence'].sum()) for s, r in zip(f['s_eta_candidats'], courses)}
    return D, H, duree


def bilan_exemple(D):
    """Biais et erreur typique de chaque lecture présente, couverture des intervalles à 95 %, η, poids des modèles de η."""
    b = {}
    for nom, col in METHODES_EXEMPLE:
        if col in D:
            err = D[col] - D.c_vrai
            b[f'{col}_biais'] = err.mean(); b[f'{col}_erreur'] = np.sqrt((err ** 2).mean())
    b['couverture_filtre'] = np.mean((D.filtre_bas <= D.c_vrai) & (D.c_vrai <= D.filtre_haut))
    b['eta_filtre_final'] = D.eta_filtre.iloc[-1]
    b['erreur_eta_filtre_rms_pct'] = 100 * np.sqrt(np.mean(np.log(D.eta_filtre / D.eta_vrai) ** 2))
    if 'c_lisse' in D:
        b['couverture_lisse'] = np.mean((D.lisse_bas <= D.c_vrai) & (D.c_vrai <= D.lisse_haut))
        b['eta_lisse'] = D.eta_lisse.iloc[0]; b['eta_lisse_final'] = D.eta_lisse.iloc[-1]
        b['erreur_eta_lisse_rms_pct'] = 100 * np.sqrt(np.mean(np.log(D.eta_lisse / D.eta_vrai) ** 2))
    for nom in ('poids_filtre_final', 'poids_lisseur'):
        w = D.attrs.get(nom, {})
        for s, x in w.items():
            b[f'{nom}_s_eta_{s:g}'] = x
        if w:
            b[f's_eta_{nom}_moyen'] = sum(s * x for s, x in w.items())
    return b
