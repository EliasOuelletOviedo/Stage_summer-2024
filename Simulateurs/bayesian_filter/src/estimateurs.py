"""
estimateurs.py — quatre lectures de [Cl⁻] sur un échantillon (histogrammes bleu et jaune de la fenêtre), calibrées une fois,
appliquées en bloc à des milliers d'échantillons (maximum de vraisemblance de Poisson sur une grille) :

  temps de vie unique : une exponentielle (⊗ IRF bleue, périodique) ajustée au canal bleu → τ, puis c par la courbe τ(c)
                        mesurée sur les histogrammes de calibration (interpolation monotone) ;
  FLIM du système     : le morceau « quand » de la section 11 — formes des deux canaux du modèle calibré, Σ n_i ln p_i(c) ;
  ratio               : le morceau « combien » à η = 1 — N_Y/N_B (fond soustrait) inversé sur R₁(c) ;
  hybride FRET-FLIM   : les deux morceaux ensemble, avec η fixé par le contexte de l'expérience : sur un échantillon de
                        contexte (beaucoup plus de photons, même état), le FLIM du système donne ĉ, et η̂ = R_ctx / R₁(ĉ) —
                        le maximum de la vraisemblance hybride à η libre ; l'échantillon est ensuite lu à η = η̂.

Le fond est une fraction connue f de chaque canal (calibration) : la forme attendue devient (1 − f) p_i + f/n, la même
pour tous les échantillons, et les totaux du ratio sont multipliés par (1 − f).
"""
import numpy as np
from scipy.interpolate import PchipInterpolator


class Estimateurs:
    def __init__(self, S, fond, c_grille, tau_grille, B_cal, Y_cal, c_cal, bloc=2000):
        """S : modele.Senseur calibré ; fond : (f_B, f_Y) ; B_cal, Y_cal : histogrammes de calibration (fenêtre), à c_cal."""
        self.c, self.tau, self.fB, self.fY, self.bloc = np.asarray(c_grille, float), np.asarray(tau_grille, float), *fond, bloc
        n = S.n_fen
        mB, mY = S.mu(self.c)
        self.MB, self.MY = mB.sum(1), mY.sum(1)
        melange = lambda m, f: np.log((1 - f) * m / m.sum(1, keepdims=True) + f / n)
        self.lB, self.lY = melange(mB, self.fB).T.copy(), melange(mY, self.fY).T.copy()          # (n, n_c)
        self.R1 = self.MY / self.MB
        mT = np.array([S.gB.declin(t)[S.fen] for t in self.tau])
        self.lT = melange(mT, self.fB).T.copy()                                                   # (n, n_τ)
        # courbe de calibration du temps de vie unique, sur les histogrammes complets de calibration
        self.c_cal = np.asarray(c_cal, float)
        self.tau_cal = self.temps_de_vie(np.asarray(B_cal, float))
        if np.any(np.diff(self.tau_cal) <= 0):
            raise ValueError(f"τ(c) non monotone sur la calibration : {self.tau_cal}")
        self._pchip = PchipInterpolator(self.tau_cal, self.c_cal, extrapolate=False)

    def _argmax(self, L, grille):
        """Maximum par ligne, affiné entre les points de la grille par une parabole sur les trois voisins."""
        i = np.argmax(L, axis=1); r = np.arange(len(L))
        j = np.clip(i, 1, L.shape[1] - 2)
        g, m, d = L[r, j - 1], L[r, j], L[r, j + 1]
        courbure = g - 2 * m + d
        with np.errstate(divide='ignore', invalid='ignore'):
            delta = np.where(courbure < 0, np.clip(0.5 * (g - d) / courbure, -0.5, 0.5), 0.0)
        x = np.where(i == j, i + delta, i)                             # pas d'affinage au bord de la grille
        return np.interp(x, np.arange(len(grille)), grille)

    def _par_blocs(self, f, *X):
        return np.concatenate([f(*(x[i:i + self.bloc] for x in X)) for i in range(0, len(X[0]), self.bloc)])

    # ------------------------------------------------------------------ temps de vie unique
    def temps_de_vie(self, B):
        return self._par_blocs(lambda b: self._argmax(b @ self.lT, self.tau), np.atleast_2d(B))

    def c_depuis_tau(self, t):
        """Inverse de τ(c) : PCHIP entre les points de calibration, droite des segments extrêmes au-delà, puis [0, c_max]."""
        t = np.asarray(t, float); c = self._pchip(t)
        tc, cc = self.tau_cal, self.c_cal
        bas, haut = t < tc[0], t > tc[-1]
        c[bas] = cc[0] + (t[bas] - tc[0]) * (cc[1] - cc[0]) / (tc[1] - tc[0])
        c[haut] = cc[-1] + (t[haut] - tc[-1]) * (cc[-1] - cc[-2]) / (tc[-1] - tc[-2])
        return np.clip(c, 0.0, self.c[-1])

    def c_temps_de_vie(self, B, Y=None):
        return self.c_depuis_tau(self.temps_de_vie(B))

    # ------------------------------------------------------------------ les trois lectures du système
    def _formes(self, B, Y):
        return B @ self.lB + Y @ self.lY

    def _ratio(self, B, Y, eta=1.0):
        NB = (1 - self.fB) * B.sum(1, keepdims=True); NY = (1 - self.fY) * Y.sum(1, keepdims=True)
        mB, mY = self.MB[None, :], self.MY[None, :] * np.broadcast_to(np.asarray(eta, float), (len(B),))[:, None]
        return NB * np.log(mB / (mB + mY)) + NY * np.log(mY / (mB + mY))

    def c_flim(self, B, Y):
        return self._par_blocs(lambda b, y: self._argmax(self._formes(b, y), self.c), B, Y)

    def _R(self, B, Y):
        return (1 - self.fY) * Y.sum(1) / ((1 - self.fB) * np.maximum(B.sum(1), 1))

    def c_ratio(self, B, Y):
        return np.interp(self._R(B, Y), self.R1[::-1], self.c[::-1])   # R₁ décroît : au-delà des bouts, 0 ou c_max

    def eta_contexte(self, B, Y):
        """η̂ = R / R₁(ĉ_FLIM) sur des échantillons de contexte : (η̂, ĉ)."""
        c = self.c_flim(B, Y)
        return self._R(B, Y) / np.interp(c, self.c, self.R1), c

    def c_hybride(self, B, Y, eta=1.0):
        eta = np.broadcast_to(np.asarray(eta, float), (len(B),))
        return self._par_blocs(lambda b, y, e: self._argmax(self._formes(b, y) + self._ratio(b, y, e), self.c), B, Y, eta)

    def toutes(self, B, Y, eta_hybride=1.0):
        return {'temps de vie unique': self.c_temps_de_vie(B), 'FLIM du système': self.c_flim(B, Y),
                'ratio': self.c_ratio(B, Y), 'hybride FRET-FLIM': self.c_hybride(B, Y, eta_hybride)}
