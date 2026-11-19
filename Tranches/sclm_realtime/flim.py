"""Durées de vie : déclins reconvolués par l'IRF, maximum de vraisemblance de Poisson.

    modèle = N · [ (1 − b) · bins(conv_P(IRF décalée de s, Σ aᵢ·exp(−t/τᵢ))) / Σ_fenêtre + b / n_fenêtre ]

Les canaux de FLIMApp font 12,5/256 ns, mais le laser se répète tous les P = 12,461 ns : la
convolution circulaire sur 256 canaux serait fausse de 0,8 canal. Le modèle est donc
calculé sur une grille fine d'une période exacte (2048 points), où la convolution
circulaire est exacte (la somme des déclins des impulsions précédentes a, sur une période,
la forme de exp(−t/τ)), puis intégré dans les canaux de FLIMApp par la fonction cumulée.
L'IRF est celle de FLIMApp (irf.csv, mêmes canaux), décalée de s (ns, sous-canal, par la
fonction cumulée) ; b = part de fond uniforme, N = photons de la fenêtre. On minimise la
déviance de Poisson.

- `ajuster_mono` : une durée de vie, décalage libre ou fixé, fond libre ;
- `ajuster_multi` : deux ou trois durées de vie (amplitudes ≥ 0), décalage libre ou fixé ;
- `temps_moyen` : temps moyen d'arrivée dans la fenêtre, sans modèle.
"""
import numpy as np
from scipy.optimize import minimize

import config as C


class Modele:
    def __init__(self, irf, fenetre=C.FENETRE, m=2048):
        self.P, self.m = C.PERIODE_NS, m
        self.h = self.P / m
        self.tf = np.arange(m) * self.h
        irf = np.asarray(irf, float)
        self.cum_irf = np.r_[0, np.cumsum(irf / irf.sum())]
        self.bords_irf = np.arange(len(irf) + 1) * C.DT_NS
        self.lo, self.hi = fenetre
        self.bords = np.arange(self.lo, self.hi + 1) * C.DT_NS          # canaux de la fenêtre (< P)
        self.n_fen = self.hi - self.lo
        self.tb = np.arange(m + 1) * self.h
        self.t_fen = (np.arange(self.lo, self.hi) + 0.5) * C.DT_NS
        self._cache = {}

    # -------------------------------------------------------------- modèle
    def F_irf(self, s):
        k = round(float(s), 6)
        F = self._cache.get(k)
        if F is None:
            g = self.tb - s
            c = np.interp(g, self.bords_irf, self.cum_irf, left=0.0, right=1.0)
            F = np.fft.rfft(np.diff(c))
            if len(self._cache) > 64:
                self._cache.clear()
            self._cache[k] = F
        return F

    def _binner(self, fine):
        cum = np.r_[0, np.cumsum(fine)]
        return np.diff(np.interp(self.bords, self.tb, cum))

    def forme(self, taus, amps, s):
        """Déclin modèle dans la fenêtre, somme 1."""
        dec = np.zeros(self.m)
        for a, t in zip(np.atleast_1d(amps), np.atleast_1d(taus)):
            dec += a * np.exp(-self.tf / t)
        c = np.maximum(np.fft.irfft(self.F_irf(s) * np.fft.rfft(dec), n=self.m), 0)
        y = self._binner(c)
        return y / max(y.sum(), 1e-300)

    def irf_fenetre(self, s):
        """IRF décalée, intégrée dans les canaux de la fenêtre (non normalisée)."""
        return self._binner(np.fft.irfft(self.F_irf(s), n=self.m))

    def attendu(self, yw, taus, amps, s, b):
        N = yw.sum()
        return N * ((1 - b) * self.forme(taus, amps, s) + b / self.n_fen)

    @staticmethod
    def deviance(yw, mu):
        mu = np.maximum(mu, 1e-300)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(yw > 0, yw * np.log(yw / mu), 0.0)
        return 2.0 * float(np.sum(t - yw + mu))

    def fen(self, y):
        return np.asarray(y, float)[self.lo:self.hi]

    # -------------------------------------------------------------- estimateurs
    def temps_moyen(self, y):
        yw = self.fen(y)
        return float((self.t_fen * yw).sum() / yw.sum()) if yw.sum() > 0 else np.nan

    def ajuster_mono(self, y, s=None, tau0=1.5, b0=0.01):
        """(τ, s, b, déviance/ddl). `s` donné : décalage fixé ; sinon ajusté."""
        yw = self.fen(y)
        if yw.sum() < 50:
            return np.nan, np.nan, np.nan, np.nan

        def dev(p):
            tau, ss, b = p[0], (p[1] if s is None else s), p[-1]
            if tau <= 0.02 or not 0 <= b < 0.9:
                return 1e30
            return self.deviance(yw, self.attendu(yw, [tau], [1.0], ss, b))

        x0 = [tau0, 0.0, b0] if s is None else [tau0, b0]
        if s is None:                       # départ du décalage : grille grossière
            best = min(((dev([tau0, ss, b0]), ss) for ss in np.arange(-0.5, 0.51, 0.05)))
            x0[1] = best[1]
        r = minimize(dev, x0, method="Nelder-Mead",
                     options=dict(xatol=1e-5, fatol=1e-4, maxiter=4000, maxfev=8000))
        p = r.x
        ddl = self.n_fen - len(p)
        return float(p[0]), float(p[1] if s is None else s), float(p[-1]), float(r.fun / ddl)

    def ajuster_multi(self, y, n_exp=2, s=None, depart=None):
        """Amplitudes (somme 1), τ (croissants), s, b, déviance/ddl. `s` donné : fixé."""
        yw = self.fen(y)
        if depart is None:
            depart = dict(taus=np.geomspace(0.3, 3.0, n_exp), amps=np.full(n_exp, 1 / n_exp), s=0.0, b=0.01)

        def unpack(p):
            taus = np.exp(p[:n_exp])
            amps = np.exp(p[n_exp:2 * n_exp])
            k = 2 * n_exp
            ss = p[k] if s is None else s
            b = 1 / (1 + np.exp(-p[-1]))
            return taus, amps / amps.sum(), ss, b

        def dev(p):
            taus, amps, ss, b = unpack(p)
            if np.any(taus > 30) or np.any(taus < 0.01):
                return 1e30
            return self.deviance(yw, self.attendu(yw, taus, amps, ss, b))

        x0 = list(np.log(depart["taus"])) + list(np.log(np.maximum(depart["amps"], 1e-3)))
        if s is None:
            x0.append(depart["s"])
        x0.append(np.log(max(depart["b"], 1e-4) / (1 - max(depart["b"], 1e-4))))
        r = minimize(dev, x0, method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-4, maxiter=20000, maxfev=40000, adaptive=True))
        r = minimize(dev, r.x, method="Nelder-Mead",
                     options=dict(xatol=1e-7, fatol=1e-5, maxiter=20000, maxfev=40000, adaptive=True))
        taus, amps, ss, b = unpack(r.x)
        o = np.argsort(taus)
        ddl = self.n_fen - len(r.x)
        return dict(taus=taus[o], amps=amps[o], s=float(ss), b=float(b), dev_ddl=float(r.fun / ddl))
