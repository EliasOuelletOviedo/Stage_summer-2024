"""Durées de vie : mono-exponentielle reconvoluée par l'IRF, maximum de vraisemblance.

Même modèle que FLIMApp (src/analysis/lifetime_analysis.jl, `mle_objective` à une
durée de vie) :

    modèle = N · [ (1 − b) · conv(IRF décalée, exp(−t/τ)) / Σ_fenêtre + b / n_fenêtre ]

convolution circulaire (FFT), normalisée sur la fenêtre utile, N le nombre de photons
de la fenêtre, b la part de fond uniforme. On minimise la déviance de Poisson
Σ y·ln(y/m) − y + m sur la fenêtre.

Deux différences assumées avec le direct : l'IRF est décalée par interpolation
linéaire de `s` canaux vers la droite (FLIMApp a son propre décalage d'origine,
absorbé par ce paramètre) ; et ce décalage, instrumental, est ajusté une fois par ROI
et par canal sur le déclin de toute l'acquisition (des millions de photons), puis fixé
pour les blocs de 3 s. Il ne reste alors que τ et b libres par bloc.

`temps_moyen` est l'estimateur sans ajustement : temps moyen d'arrivée moins celui de
l'IRF — le « Estimate » de FLIMApp.
"""
import numpy as np
from scipy.optimize import minimize

from config import COUPE_BAS, COUPE_HAUT, MIN_PHOTONS_FIT


class Modele:
    def __init__(self, irf, dt_ns, coupe_bas=COUPE_BAS, coupe_haut=COUPE_HAUT):
        self.n = len(irf)
        self.irf = np.asarray(irf, float) / np.sum(irf)
        self.dt = dt_ns
        self.x = np.arange(self.n) * dt_ns
        self.lo, self.hi = coupe_bas - 1, self.n - coupe_haut      # fenêtre [lo, hi)
        self.n_fen = self.hi - self.lo
        self.k = np.arange(self.n)
        self.t_irf = self._moyen(self.irf)

    # -------------------------------------------------------------- modèle
    def _irf_decalee(self, s):
        f = np.floor(s)
        w = s - f
        i = (self.k - int(f)) % self.n
        j = (self.k - int(f) - 1) % self.n
        return (1 - w) * self.irf[i] + w * self.irf[j]

    def forme(self, tau, s):
        """conv(IRF décalée de s, exp(−t/τ)), normalisée à 1 sur la fenêtre."""
        c = np.real(np.fft.ifft(np.fft.fft(self._irf_decalee(s)) * np.fft.fft(np.exp(-self.x / tau))))
        c = np.maximum(c, 0)
        return c / max(c[self.lo:self.hi].sum(), 1e-300)

    def attendu(self, y, tau, s, b):
        N = y[self.lo:self.hi].sum()
        return N * ((1 - b) * self.forme(tau, s) + b / self.n_fen)

    def _deviance(self, y, tau, s, b):
        m = np.maximum(self.attendu(y, tau, s, b)[self.lo:self.hi], 1e-300)
        yy = y[self.lo:self.hi]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(yy > 0, yy * np.log(yy / m), 0.0)
        return float(np.sum(t - yy + m))

    # -------------------------------------------------------------- estimateurs
    def _moyen(self, y):
        k = np.arange(self.lo, self.hi)
        v = y[self.lo:self.hi]
        return float((k * v).sum() / max(v.sum(), 1e-300)) * self.dt

    def temps_moyen(self, y):
        """Durée de vie « moyenne » sans ajustement : <t> − <t_IRF>, en ns."""
        y = np.asarray(y, float)
        if y[self.lo:self.hi].sum() < MIN_PHOTONS_FIT:
            return np.nan
        return self._moyen(y) - self.t_irf

    def ajuster(self, y, s=None, tau0=None):
        """(τ, s, b, déviance réduite). `s` donné : décalage fixé ; sinon ajusté."""
        y = np.asarray(y, float)
        if y[self.lo:self.hi].sum() < MIN_PHOTONS_FIT:
            return np.nan, np.nan, np.nan, np.nan
        tau0 = tau0 if tau0 is not None else float(np.clip(self.temps_moyen(y), 0.2, 8.0))
        libre = s is None
        # Déviance rapportée au nombre de photons : même minimum, mais une échelle qui ne
        # fait pas échouer la recherche linéaire de L-BFGS-B sur des millions de photons.
        N = max(y[self.lo:self.hi].sum(), 1.0)
        if libre:
            # Le décalage a des minima locaux : grille grossière, puis affinage libre.
            grille = [(self._deviance(y, tau0, si, 0.05), si) for si in np.arange(-20, 40.5, 1.0)]
            s0 = min(grille)[1]
            r = minimize(lambda p: self._deviance(y, p[0], p[1], p[2]) / N, [tau0, s0, 0.05],
                         method="L-BFGS-B",
                         bounds=[(2 * self.dt, 2 * 12.5), (-20, 40), (0, 1)])
            tau, s, b = r.x
        else:
            r = minimize(lambda p: self._deviance(y, p[0], s, p[1]) / N, [tau0, 0.05],
                         method="L-BFGS-B", bounds=[(2 * self.dt, 2 * 12.5), (0, 1)])
            tau, b = r.x
        # Un arrêt « ABNORMAL » de la recherche linéaire laisse un minimum valable : seul
        # un résultat non fini est rejeté.
        if not (np.isfinite(r.fun) and np.all(np.isfinite(r.x))):
            return np.nan, s, np.nan, np.nan
        dof = max(self.n_fen - (3 if libre else 2), 1)
        return float(tau), float(s), float(b), 2 * r.fun * N / dof
