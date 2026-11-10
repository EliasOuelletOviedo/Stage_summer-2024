"""
simulation.py — les données simulées, inspirées des mesures de titration du 2026-10-16 (SCLM_mCer3, SCLM).

Deux générateurs de trames (histogrammes bleu et jaune, cases de la fenêtre) :
  - « modèle »        : comptages de Poisson autour du modèle calibré (modele.Senseur), pour toute trajectoire
                        continue de [Cl⁻] ; η, l'amplitude A et ses dérives, le fond et les perturbations de la
                        section 16 se règlent ici, sans que le filtre les connaisse ;
  - « amincissement » : comptages de Poisson autour des histogrammes MESURÉS de la titration (mis à l'échelle de
                        la trame, jaune multiplié par η). Aucune hypothèse de modèle : fond, IRF réelle, décalages
                        d'horloge et défauts du modèle (déviance/ddl = 2,34) sont ceux des données. Seules les 8
                        concentrations de la titration sont disponibles.
"""
import os
import numpy as np
import pandas as pd

ETIQUETTES = ['0mM', '2p5mM', '5mM', '10mM', '20mM', '40mM', '80mM', '131mM']
CONC = np.array([0, 2.5, 5, 10, 20, 40, 80, 131], float)


def trame_modele(S, c, photons, rng, eta=1.0, fond=(0.0, 0.0), extra=None):
    """
    Une trame de Poisson autour du modèle. `photons` : photons du senseur attendus (bleu + jaune, fenêtre).
    fond : photons par case et par trame (bleu, jaune). extra : (μ_B, μ_Y) ajoutés (autofluorescence…).
    """
    b, y = S.mu(c, 1.0, eta)
    A = photons / (b.sum() + y.sum())
    muB, muY = A * b + fond[0], A * y + fond[1]
    if extra is not None:
        muB, muY = muB + extra[0], muY + extra[1]
    return rng.poisson(muB).astype(float), rng.poisson(muY).astype(float)


# ------------------------------------------------------------------------------------------- trajectoires
def trajectoire_exemple(e):
    """Section 15 : base, montée exponentielle vers c_haut à trame_montee, retour vers c_fin à trame_descente."""
    k = np.arange(e['trames'], dtype=float); c = np.full(e['trames'], e['c_base_mM'])
    on, off = e['trame_montee'] - 1, e['trame_descente'] - 1          # trames numérotées à partir de 1
    m = k >= on
    c[m] = e['c_haut_mM'] + (e['c_base_mM'] - e['c_haut_mM']) * np.exp(-(k[m] - on) / e['tau_montee_trames'])
    d = k >= off
    c[d] = e['c_fin_mM'] + (c[off - 1] - e['c_fin_mM']) * np.exp(-(k[d] - off) / e['tau_descente_trames'])
    return c


def trajectoires_roi(rng, n_roi, T):
    """Section 14 : [Cl⁻] propre à chaque ROI — base 6–12 mM, réponse de 8–16 mM, début et durée tirés au hasard."""
    k = np.arange(T, dtype=float); C = np.empty((n_roi, T))
    for m in range(n_roi):
        base, amp = rng.uniform(6, 12), rng.uniform(8, 16)
        on = int(rng.uniform(0.2, 0.4) * T); off = on + int(rng.uniform(0.2, 0.33) * T)
        c = np.full(T, base)
        c[on:] = base + amp * (1 - np.exp(-(k[on:] - on) / 4.0))
        c[off:] = base + (c[off - 1] - base) * np.exp(-(k[off:] - off) / 25.0)
        C[m] = c
    return C


# ------------------------------------------------------------------------------------------- amincissement
def _lire(dossier, fichier):
    """Un CSV d'histogrammes et ses deux préfixes de colonnes : voie1/voie2 (step17, mCer3) ou bleu/jaune (step16, SCLM)."""
    df = pd.read_csv(os.path.join(dossier, fichier)); df.columns = [c.strip() for c in df.columns]
    return df, (('voie1', 'voie2') if any(c.startswith('voie1_brut_') for c in df.columns) else ('bleu', 'jaune'))


class Mesures:
    """Les 16 histogrammes mesurés de la titration, restreints à la fenêtre d'analyse."""
    def __init__(self, dossier, fichier, cal):
        df, (vB, vY) = _lire(dossier, fichier)
        a, b = cal['fenetre']
        self.B = np.array([df[f'{vB}_brut_{e}'].values[a:b + 1] for e in ETIQUETTES], float)
        self.Y = np.array([df[f'{vY}_brut_{e}'].values[a:b + 1] for e in ETIQUETTES], float)
        self.c = CONC

    def indice(self, c):
        i = np.flatnonzero(np.isclose(self.c, c))
        if len(i) == 0:
            raise ValueError(f"{c} mM n'est pas une concentration de la titration ({self.c})")
        return int(i[0])

    def trame(self, c, photons, rng, eta=1.0):
        """Poisson autour de la forme mesurée : `photons` attendus au total (bleu + jaune) à η = 1, jaune × η."""
        i = self.indice(c); s = photons / (self.B[i].sum() + self.Y[i].sum())
        return rng.poisson(s * self.B[i]).astype(float), rng.poisson(eta * s * self.Y[i]).astype(float)


def inconnu(dossier, fichier, cal):
    """L'échantillon « inconnu » mesuré le même jour (histogrammes bleu et jaune, fenêtre)."""
    df, (vB, vY) = _lire(dossier, fichier)
    a, b = cal['fenetre']
    return df[f'{vB}_brut_inconnu'].values[a:b + 1].astype(float), df[f'{vY}_brut_inconnu'].values[a:b + 1].astype(float)


def fond_moyen(cal):
    """Part moyenne du fond dans la fenêtre (bleu, jaune), sur les 8 histogrammes de la titration."""
    t = cal['titration']
    return np.mean([x['part_fond_bleu'] for x in t]), np.mean([x['part_fond_jaune'] for x in t])
