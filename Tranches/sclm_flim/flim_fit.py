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

`ajuster2` : deux durées de vie, le modèle à 5 paramètres de FLIMApp (`conv_irf_data`)
qui a servi aux calibrations (tests/analyse_SCLM_3.jl) : a1·exp(−t/τ1) + (1 − a1)·exp(−t/τ2),
τ1 la plus longue et a1 son amplitude ; `tau_moyen` en donne la moyenne pondérée par
l'intensité. Trois paramètres libres demandent beaucoup de photons : ~6 000 (un soma,
un bloc de 3 s) ne suffisent pas, ~260 000 (34 somas réunis) oui.
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

    def _convoluer(self, declin, s):
        """conv(IRF décalée de s, declin), normalisée à 1 sur la fenêtre."""
        c = np.real(np.fft.ifft(np.fft.fft(self._irf_decalee(s)) * np.fft.fft(declin)))
        c = np.maximum(c, 0)
        return c / max(c[self.lo:self.hi].sum(), 1e-300)

    def forme(self, tau, s):
        """conv(IRF décalée de s, exp(−t/τ)), normalisée à 1 sur la fenêtre."""
        return self._convoluer(np.exp(-self.x / tau), s)

    def forme2(self, tau1, tau2, a1, s):
        """Comme `forme`, pour a1·exp(−t/τ1) + (1 − a1)·exp(−t/τ2)."""
        return self._convoluer(a1 * np.exp(-self.x / tau1) + (1 - a1) * np.exp(-self.x / tau2), s)

    def attendu(self, y, tau, s, b):
        N = y[self.lo:self.hi].sum()
        return N * ((1 - b) * self.forme(tau, s) + b / self.n_fen)

    def _deviance(self, y, tau, s, b):
        return self._deviance_modele(y, self.attendu(y, tau, s, b))

    def _deviance2(self, y, tau1, tau2, a1, s, b):
        N = y[self.lo:self.hi].sum()
        return self._deviance_modele(y, N * ((1 - b) * self.forme2(tau1, tau2, a1, s) + b / self.n_fen))

    def _deviance_modele(self, y, attendu):
        m = np.maximum(attendu[self.lo:self.hi], 1e-300)
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

    def ajuster2(self, y, s=None, depart=(2.5, 0.4, 0.4, 0.05)):
        """Deux durées de vie : (τ1, τ2, a1, s, b, déviance réduite), τ1 la plus longue.
        `s` donné : décalage fixé ; sinon ajusté, en partant de celui de la
        mono-exponentielle. `depart` : (τ1, τ2, a1, b) de départ."""
        y = np.asarray(y, float)
        if y[self.lo:self.hi].sum() < MIN_PHOTONS_FIT:
            return (np.nan,) * 6
        N = max(y[self.lo:self.hi].sum(), 1.0)
        bornes = [(2 * self.dt, 2 * 12.5)] * 2 + [(0, 1), (0, 1)]
        libre = s is None
        if libre:
            r = minimize(lambda p: self._deviance2(y, p[0], p[1], p[2], p[4], p[3]) / N,
                         list(depart) + [self.ajuster(y)[1]], method="L-BFGS-B",
                         bounds=bornes + [(-20, 40)])
            tau1, tau2, a1, b, s = r.x
        else:
            r = minimize(lambda p: self._deviance2(y, p[0], p[1], p[2], s, p[3]) / N,
                         list(depart), method="L-BFGS-B", bounds=bornes)
            tau1, tau2, a1, b = r.x
        if not (np.isfinite(r.fun) and np.all(np.isfinite(r.x))):
            return (np.nan,) * 6
        if tau2 > tau1:                                  # τ1 la plus longue, comme FLIMApp
            tau1, tau2, a1 = tau2, tau1, 1 - a1
        dof = max(self.n_fen - (5 if libre else 4), 1)
        return float(tau1), float(tau2), float(a1), float(s), float(b), 2 * r.fun * N / dof


def tau_moyen(tau1, tau2, a1):
    """Durée de vie moyenne pondérée par l'intensité, (a1·τ1² + a2·τ2²) / (a1·τ1 + a2·τ2) :
    le τ moyen de tests/analyse_SCLM_3.jl."""
    a2 = 1 - a1
    return (a1 * tau1 ** 2 + a2 * tau2 ** 2) / (a1 * tau1 + a2 * tau2)
