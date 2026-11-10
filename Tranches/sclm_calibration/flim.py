"""Durées de vie : déclins reconvolués par l'IRF, maximum de vraisemblance de Poisson.

    modèle = N · [ (1 − b) · conv(IRF décalée de s, Σ aᵢ·exp(−t/τᵢ)) / Σ_fenêtre + b / n_fenêtre ]

Convolution circulaire sur une période (la somme des déclins de toutes les impulsions
précédentes a, sur une période, la forme de exp(−t/τ) : la convolution circulaire est donc
exacte), normalisée sur la fenêtre utile [17, 256) ; N = photons de la fenêtre, b = part de
fond uniforme. On minimise la déviance de Poisson Σ y·ln(y/m) − y + m. Mêmes choix que
sclm_flim/flim_fit.py (le modèle de FLIMApp), sauf le décalage de l'IRF : appliqué dans
l'espace de Fourier (sous-canal exact pour un signal périodique) au lieu d'une interpolation
linéaire, qui élargit l'IRF d'un demi-canal au plus — à 48,7 ps par canal, l'IRF (98 à
146 ps de largeur à mi-hauteur) n'occupe que 2 à 3 canaux.

- `ajuster_mono` : une durée de vie (le modèle de FLIMApp), décalage libre ou fixé ;
- `ajuster_multi` : n durées de vie (amplitudes ≥ 0), décalage libre ou fixé ;
- `amplitudes` : τ fixées, amplitudes et fond libres (populations à durées de vie connues) ;
- `temps_moyen` : temps moyen d'arrivée moins celui de l'IRF décalée, sans modèle.
"""
import numpy as np
from scipy.optimize import minimize

import config as C


class Modele:
    def __init__(self, irf, fenetre=C.FENETRE):
        self.n = C.N_CANAUX
        self.dt = C.DT_NS
        self.irf = np.asarray(irf, float) / np.sum(irf)
        self.F_irf = np.fft.fft(self.irf)
        self.freq = np.fft.fftfreq(self.n, d=1.0)            # cycles par canal
        self.x = np.arange(self.n) * self.dt
        self.lo, self.hi = fenetre
        self.n_fen = self.hi - self.lo
        self.k = np.arange(self.n)

    # -------------------------------------------------------------- modèle
    def _F_irf(self, s, sigma_b=0.0):
        F = self.F_irf * np.exp(-2j * np.pi * self.freq * s)
        if sigma_b:
            F = F * np.exp(-2 * (np.pi * self.freq * sigma_b / self.dt) ** 2)
        return F

    def irf_decalee(self, s):
        return np.maximum(np.real(np.fft.ifft(self._F_irf(s))), 0)

    def forme(self, taus, amps, s, sigma_b=0.0):
        """conv(IRF décalée de s, Σ aᵢ exp(−t/τᵢ)), normalisée à 1 sur la fenêtre."""
        dec = np.zeros(self.n)
        for a, t in zip(np.atleast_1d(amps), np.atleast_1d(taus)):
            dec += a * np.exp(-self.x / t)
        c = np.maximum(np.real(np.fft.ifft(self._F_irf(s, sigma_b) * np.fft.fft(dec))), 0)
        return c / max(c[self.lo:self.hi].sum(), 1e-300)

    def attendu(self, y, taus, amps, s, b, sigma_b=0.0):
        N = y[self.lo:self.hi].sum()
        return N * ((1 - b) * self.forme(taus, amps, s, sigma_b) + b / self.n_fen)

    def deviance(self, y, attendu):
        m = np.maximum(attendu[self.lo:self.hi], 1e-300)
        yy = y[self.lo:self.hi]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(yy > 0, yy * np.log(yy / m), 0.0)
        return 2.0 * float(np.sum(t - yy + m))

    # -------------------------------------------------------------- estimateurs
    def temps_moyen(self, y, s):
        """<t> − <t_IRF décalée>, dans la fenêtre, en ns (sans fond retiré)."""
        y = np.asarray(y, float)
        k = np.arange(self.lo, self.hi)
        v = y[self.lo:self.hi]
        if v.sum() <= 0:
            return np.nan
        g = self.irf_decalee(s)[self.lo:self.hi]
        return float((k * v).sum() / v.sum() - (k * g).sum() / g.sum()) * self.dt

    def ajuster_mono(self, y, s=None, tau0=2.0):
        """(τ, s, b, déviance/ddl). `s` donné : décalage fixé ; sinon ajusté (grille puis libre)."""
        y = np.asarray(y, float)
        N = y[self.lo:self.hi].sum()
        if N < C.MIN_PHOTONS_FIT:
            return np.nan, np.nan, np.nan, np.nan
        f = lambda p, s_: self.deviance(y, self.attendu(y, [p[0]], [1.0], s_, p[1])) / N
        if s is None:
            grille = min((f([tau0, 0.03], si), si) for si in np.arange(-15, 10.5, 0.5))
            r = minimize(lambda p: f(p[:2], p[2]), [tau0, 0.03, grille[1]], method="L-BFGS-B",
                         bounds=[(0.05, 20), (0, 0.9), (-20, 15)])
            tau, b, s = r.x
            ddl = self.n_fen - 3
        else:
            r = minimize(lambda p: f(p, s), [tau0, 0.03], method="L-BFGS-B",
                         bounds=[(0.05, 20), (0, 0.9)])
            tau, b = r.x
            ddl = self.n_fen - 2
        return float(tau), float(s), float(b), float(r.fun * N / ddl)

    def ajuster_multi(self, y, n_exp, s=None, depart=None, sigma_b_libre=False):
        """n exponentielles : dict(taus, amps (somme 1), s, b, sigma_b, dev_ddl).
        Amplitudes par softmax (≥ 0) ; plusieurs départs du décalage si `s` est libre."""
        y = np.asarray(y, float)
        N = y[self.lo:self.hi].sum()
        taus0 = depart or {1: [2.0], 2: [3.0, 0.4], 3: [3.0, 0.45, 0.06]}[n_exp]

        def unpack(p):
            k = 0
            s_ = p[k] if s is None else s
            k += s is None
            sb = p[k] if sigma_b_libre else 0.0
            k += sigma_b_libre
            b = p[k]
            taus = p[k + 1:k + 1 + n_exp]
            z = np.r_[p[k + 1 + n_exp:], 0.0]
            amps = np.exp(z - z.max())
            return s_, sb, b, taus, amps / amps.sum()

        def f(p):
            s_, sb, b, taus, amps = unpack(p)
            return self.deviance(y, self.attendu(y, taus, amps, s_, b, sb)) / N

        bnds = ([(-20, 15)] if s is None else []) + ([(0, 0.3)] if sigma_b_libre else []) \
            + [(0, 0.9)] + [(0.01, 20)] * n_exp + [(-10, 10)] * (n_exp - 1)
        departs = np.arange(-12, 4, 1.5) if s is None else [None]
        best = None
        for s0 in departs:
            p0 = ([s0] if s is None else []) + ([0.02] if sigma_b_libre else []) + [0.02] + list(taus0) + [0.0] * (n_exp - 1)
            r = minimize(f, p0, method="L-BFGS-B", bounds=bnds)
            if best is None or r.fun < best.fun:
                best = r
        s_, sb, b, taus, amps = unpack(best.x)
        o = np.argsort(taus)[::-1]                       # τ décroissantes
        return dict(taus=np.asarray(taus)[o], amps=amps[o], s=float(s_), b=float(b), sigma_b=float(sb),
                    dev_ddl=float(best.fun * N / (self.n_fen - len(best.x))), deviance=float(best.fun * N))

    def amplitudes(self, y, taus, s):
        """τ fixées : (amplitudes pré-exponentielles normalisées, fond b, déviance/ddl)."""
        y = np.asarray(y, float)
        N = y[self.lo:self.hi].sum()
        n_exp = len(taus)
        if N < C.MIN_PHOTONS_FIT:
            return np.full(n_exp, np.nan), np.nan, np.nan

        def f(p):
            z = np.r_[p[1:], 0.0]
            a = np.exp(z - z.max())
            return self.deviance(y, self.attendu(y, taus, a / a.sum(), s, p[0])) / N

        r = minimize(f, [0.02] + [0.0] * (n_exp - 1), method="L-BFGS-B",
                     bounds=[(0, 0.9)] + [(-12, 12)] * (n_exp - 1))
        z = np.r_[r.x[1:], 0.0]
        a = np.exp(z - z.max())
        return a / a.sum(), float(r.x[0]), float(r.fun * N / (self.n_fen - n_exp))


def tau_moyen_intensite(taus, amps):
    """Σ aᵢτᵢ² / Σ aᵢτᵢ : durée de vie moyenne pondérée par l'intensité."""
    taus, amps = np.asarray(taus, float), np.asarray(amps, float)
    return float((amps * taus ** 2).sum() / (amps * taus).sum())
