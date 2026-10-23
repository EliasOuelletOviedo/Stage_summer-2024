"""Durées de vie : exponentielle incomplète reconvoluée par l'IRF, maximum de vraisemblance.

Le laser tire toutes les P = 12,49 ns ; avec τ ~ 5 ns (faible Cl⁻), le déclin d'une
impulsion n'est pas fini quand arrive la suivante. Le modèle est donc celui du « déclin
incomplet » : la somme des déclins de toutes les impulsions précédentes,

    h(t) = Σ_n exp(−(t + nP)/τ) = exp(−t/τ) / (1 − exp(−P/τ)),   0 ≤ t < P,

reconvoluée par l'IRF mesurée de la carte et du canal (fichiers irf_16x_750nm_*.sdt,
sur la même grille de 256 canaux par période). h étant périodique, la convolution est
circulaire et exacte : l'IRF décalée de s canaux (interpolation linéaire) convolue h par
FFT. Avec un fond uniforme b :

    modèle = N · [ (1 − b) · (IRF_s ⊛ h_τ) / Σ_fenêtre + b / n_fenêtre ]

N = photons de la fenêtre de la carte (FENETRE) ; on minimise la déviance de Poisson
Σ y·ln(y/m) − y + m sur cette fenêtre. Paramètres libres : τ, b et, sur le déclin global
d'une acquisition, le décalage s ; pour ses tranches de taux, s est fixé à cette valeur
(l'empilement change la forme du déclin, pas l'instant de l'impulsion).
"""
import numpy as np
from scipy.optimize import minimize

from config import FENETRE, N_CANAUX, PERIODE_NS


class Modele:
    def __init__(self, irf, carte):
        self.n = N_CANAUX
        self.irf = np.asarray(irf, float) / np.sum(irf)
        self.dt = PERIODE_NS / N_CANAUX
        self.x = np.arange(self.n) * self.dt
        self.lo, self.hi = FENETRE[carte]
        self.n_fen = self.hi - self.lo
        self.k = np.arange(self.n)

    def _irf_decalee(self, s):
        f = np.floor(s)
        w = s - f
        i = (self.k - int(f)) % self.n
        j = (self.k - int(f) - 1) % self.n
        return (1 - w) * self.irf[i] + w * self.irf[j]

    def declin_incomplet(self, tau):
        """h(t) = exp(−t/τ) / (1 − exp(−P/τ)) : toutes les impulsions précédentes."""
        return np.exp(-self.x / tau) / -np.expm1(-PERIODE_NS / tau)

    def forme(self, tau, s):
        """IRF_s ⊛ h_τ (convolution circulaire), normalisée à 1 sur la fenêtre."""
        c = np.real(np.fft.ifft(np.fft.fft(self._irf_decalee(s)) * np.fft.fft(self.declin_incomplet(tau))))
        c = np.maximum(c, 0)
        return c / max(c[self.lo:self.hi].sum(), 1e-300)

    def attendu(self, y, tau, s, b):
        N = y[self.lo:self.hi].sum()
        return N * ((1 - b) * self.forme(tau, s) + b / self.n_fen)

    def deviance(self, y, m):
        m = np.maximum(m[self.lo:self.hi], 1e-300)
        yy = y[self.lo:self.hi]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(yy > 0, yy * np.log(yy / m), 0.0)
        return float(np.sum(t - yy + m))

    def ajuster(self, y, s=None, tau0=3.0):
        """(τ, s, b, déviance réduite, erreur-type de τ). `s` donné : décalage fixé."""
        y = np.asarray(y, float)
        N = max(y[self.lo:self.hi].sum(), 1.0)
        f = lambda tau, s_, b: self.deviance(y, self.attendu(y, tau, s_, b)) / N
        bornes_tau = (0.1, 60.0)
        if s is None:
            # Le décalage a des minima locaux : grille grossière, puis affinage libre.
            grille = [(f(tau0, si, 0.02), si) for si in np.arange(-20, 20.5, 1.0)]
            s0 = min(grille)[1]
            r = minimize(lambda p: f(p[0], p[1], p[2]), [tau0, s0, 0.02], method="L-BFGS-B",
                         bounds=[bornes_tau, (-25, 25), (0, 1)])
            tau, s, b = r.x
            libres = 3
        else:
            r = minimize(lambda p: f(p[0], s, p[1]), [tau0, 0.02], method="L-BFGS-B",
                         bounds=[bornes_tau, (0, 1)])
            tau, b = r.x
            libres = 2
        # Erreur-type de τ : courbure de −log-vraisemblance (= f · N) en τ, b fixé.
        h = 1e-3 * tau
        d2 = (f(tau + h, s, b) - 2 * f(tau, s, b) + f(tau - h, s, b)) * N / h ** 2
        err = float(np.sqrt(1 / d2)) if d2 > 0 else np.nan
        dof = max(self.n_fen - libres, 1)
        return float(tau), float(s), float(b), 2 * r.fun * N / dof, err
