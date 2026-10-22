"""
exemple.py — une course de la section 15 : un soma, trames simulées le long de la trajectoire de l'exemple, lues par le
ratio (η = 1), la durée de vie seule, et le filtre hybride (temps réel et lissé) sur la grille (c, η) complète.
Sert à run.py (partie 3) et à exemples.py (la même course sous d'autres conditions).
"""
import time
import numpy as np
import pandas as pd

from .filtre import transition, filtre_lisseur, resume_c, resume_eta
from .simulation import trame_modele, trajectoire_exemple

METHODES_EXEMPLE = (('ratio avec la calibration (η = 1)', 'c_ratio_eta1'), ('durée de vie seule, trame par trame', 'c_flim'),
                    ('filtre hybride, temps réel', 'c_filtre'), ('filtre hybride, lissé après coup', 'c_lisse'))


def noyaux(g, f):
    return transition(np.log(g.c), f['s'], f['s_saut'], f['w']), transition(np.log(g.eta), f['s_eta_lisseur'])


def etirer(e, facteur):
    """La même trajectoire sur `facteur` fois plus de trames : montée, descente et constantes de temps à l'échelle."""
    e = dict(e)
    e['trames'] = facteur * e['trames']
    e['trame_montee'] = facteur * (e['trame_montee'] - 1) + 1
    e['trame_descente'] = facteur * (e['trame_descente'] - 1) + 1
    e['tau_montee_trames'] *= facteur; e['tau_descente_trames'] *= facteur
    return e


def course_exemple(S, g, f, e, rng, eta=None):
    """
    Une course : e est le bloc [exemple] (éventuellement étiré), f le bloc [filtre] ; eta, s'il est donné, est le η vrai
    de chaque trame (sinon e['eta'] partout). Renvoie (table par trame, histogrammes, durée du filtre et du lisseur).
    """
    T = e['trames']
    C = trajectoire_exemple(e)
    eta = np.full(T, e['eta']) if eta is None else np.asarray(eta, float)
    A = 1 + e['derive_A'] * np.arange(T) / (T - 1)            # photoblanchiment : le filtre ne le connaît pas
    H = [trame_modele(S, C[k], e['photons'] * A[k], rng, eta[k]) for k in range(T)]
    t0 = time.perf_counter()
    Ls = np.array([g.loglik(B, Y) for B, Y in H])
    Kc, Ke = noyaux(g, f)
    alpha, post = filtre_lisseur(Ls, Kc, Ke)
    duree = time.perf_counter() - t0
    out = []
    for k in range(T):
        B, Y = H[k]
        fm, flo, fhi = resume_c(alpha[k].sum(1), g.c); lm, llo, lhi = resume_c(post[k].sum(1), g.c)
        ef, ef_sd = resume_eta(alpha[k].sum(0), g.eta); el, el_sd = resume_eta(post[k].sum(0), g.eta)
        out.append(dict(trame=k + 1, c_vrai=C[k], eta_vrai=eta[k], N_B=B.sum(), N_Y=Y.sum(), R=Y.sum() / B.sum(),
                        c_ratio_eta1=g.lire_ratio(B, Y), c_flim=g.lire_flim(B, Y),
                        c_filtre=fm[0], filtre_bas=flo[0], filtre_haut=fhi[0], eta_filtre=ef, eta_filtre_sd=ef_sd,
                        c_lisse=lm[0], lisse_bas=llo[0], lisse_haut=lhi[0], eta_lisse=el, eta_lisse_sd=el_sd))
    return pd.DataFrame(out), H, duree


def bilan_exemple(D):
    b = {}
    for nom, col in METHODES_EXEMPLE:
        err = D[col] - D.c_vrai
        b[f'{col}_biais'] = err.mean(); b[f'{col}_erreur'] = np.sqrt((err ** 2).mean())
    b['couverture_filtre'] = np.mean((D.filtre_bas <= D.c_vrai) & (D.c_vrai <= D.filtre_haut))
    b['couverture_lisse'] = np.mean((D.lisse_bas <= D.c_vrai) & (D.c_vrai <= D.lisse_haut))
    b['eta_filtre_final'] = D.eta_filtre.iloc[-1]; b['eta_lisse'] = D.eta_lisse.iloc[0]
    return b
