#!/usr/bin/env python3
"""
sclm_hybride.py — Calibration FLIM 2 canaux + ratio pour SuperClomeleon (mCerulean3 -> Topaz)
et filtre bayésien hybride pour les séries temporelles en tranche.

Usage :
    python sclm_hybride.py irf_bruts_mCer3.csv --out calibration_sclm.json

Étapes :
 1. IRF paramétrique (gaussienne ⊗ exponentielles) ajustée sur chaque IRF brute, sans écrêtage.
 2. Ajustement conjoint (vraisemblance de Poisson) des 16 histogrammes (8 [Cl-] × 2 canaux) :
      bleu  : A [ f1 e^{-t/τ1} + f2 e^{-t/τ2} + (1-f1-f2) e^{-t/τD} ]
      jaune : η { λ · donneur(t)                                     (fuite spectrale du donneur)
                + G Σ_i g_i f_i k_i (e^{-t/τi} ⊗ e^{-t/τA})          (émission sensibilisée, k_i = 1/τi - 1/τD)
                + κ (f1+f2) e^{-t/τA} }                              (excitation directe des accepteurs anioniques)
      g_1 = 1, g_2 = G2r ; η = 1 par définition dans la calibration ; décalage temporel propre à chaque acquisition.
      Modèle périodique (queue des impulsions précédentes incluse), convolution sur grille fine.
 3. Isothermes f_i([Cl]) = f_min + (f_max - f_min) / (1 + [Cl]/K).
 4. Précision attendue (Cramér-Rao) : FLIM bleu, FLIM jaune, FLIM 2 canaux (η inconnu), ratio/hybride (η connu).
 5. HybridGrid + forward_backward : filtre/lisseur bayésien hors ligne sur la grille ([Cl], η).
 6. RealTimeHybrid : filtre causal temps réel, exact pour un η statique partagé par un groupe de ROI,
    avec portes temporelles (log_gates) et bascule vers un filtre 1D sur le ratio quand η est connu.
    Démo : python sclm_hybride.py irf_bruts_mCer3.csv --demo-rt --photons 3000

Dépendances : numpy, scipy, pandas.

Exemple d'utilisation sur une série temporelle en tranche (B[k], Y[k] : histogrammes d'une ROI pour la trame k,
restreints aux canaux de la fenêtre, acquis avec les mêmes réglages que la calibration) :
    d = load('irf_bruts_mCer3.csv'); irfp = {1: fit_irf(d, 1, 2), 2: fit_irf(d, 2, 3)}
    R = fit_joint(d, irfp); iso = fit_isotherms(R)
    consts = {k: R['p'][k] for k in ['t1', 't2', 'tD', 'tA', 'lam', 'G', 'G2r', 'k1', 'k0']}
    S = Sensor(d, consts, iso, irfp, shifts=(s_bleu, s_jaune))      # décalages mesurés pour la session de tranche
    H = HybridGrid(S, np.geomspace(0.5, 300, 260), np.geomspace(0.6, 1.8, 61))
    Ls = [H.frame_loglik(B[k], Y[k]) for k in range(len(B))]
    Kc = transition(np.log(H.c), 0.04, jump_sd=0.5, jump_w=0.03)    # [Cl] : marche aléatoire log + sauts rares
    Ke = transition(np.log(H.eta), 0.002)                            # eta : quasi constant pour une ROI
    alpha, post = forward_backward(Ls, Kc, Ke)                       # alpha : filtre causal ; post : lissage
    c_moyen = (post.sum(2) * H.c).sum(1)
À ajouter pour des données réelles : un terme d'autofluorescence par canal (gabarit de déclin mesuré dans du tissu
non exprimant, amplitude traitée comme une nuisance lente), le fond, et la ré-estimation des tau sur l'ensemble
des pixels de la tranche (a priori centré sur la calibration).
"""
import argparse, json, time
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, curve_fit
from scipy.optimize._numdiff import approx_derivative
from scipy.special import erfcx, erfc

CONC_LABELS = ['0mM', '2p5mM', '5mM', '10mM', '20mM', '40mM', '80mM', '131mM']
CONC = np.array([0, 2.5, 5, 10, 20, 40, 80, 131], dtype=float)
M_FINE = 16


# ============================================================================ données
def load(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    for c in ['dans_fenetre_voie1', 'dans_fenetre_voie2']:
        df[c] = df[c].astype(str).str.strip().str.lower().map({'true': True, 'false': False})
    T = float(df.t_fin_ns.iloc[-1]); NB = len(df)
    return dict(T=T, NB=NB, DT=T / NB, t=df.t_debut_ns.values + 0.5 * T / NB,
                win={1: df.dans_fenetre_voie1.values.astype(bool), 2: df.dans_fenetre_voie2.values.astype(bool)},
                irf_raw={1: df.irf_voie1_brute.values.astype(float), 2: df.irf_voie2_brute.values.astype(float)},
                B=np.array([df[f'voie1_brut_{c}'].values for c in CONC_LABELS], float),
                Y=np.array([df[f'voie2_brut_{c}'].values for c in CONC_LABELS], float))


def dev_resid(y, mu):
    """Résidus de déviance de Poisson (somme des carrés = déviance)."""
    mu = np.maximum(mu, 1e-12)
    with np.errstate(divide='ignore', invalid='ignore'):
        term = np.where(y > 0, y * np.log(y / mu), 0.0)
    return np.sign(y - mu) * np.sqrt(np.maximum(2.0 * (mu - y + term), 0.0))


class Params:
    def __init__(self): self.names, self.x0, self.lo, self.hi = [], [], [], []
    def add(self, n, x0, lo=-np.inf, hi=np.inf): self.names.append(n); self.x0.append(x0); self.lo.append(lo); self.hi.append(hi)
    def unpack(self, v): return dict(zip(self.names, v))
    @property
    def bounds(self): return np.array(self.lo, float), np.array(self.hi, float)
    @property
    def start(self):
        x = np.array(self.x0, float); lo, hi = self.bounds
        with np.errstate(invalid='ignore'):
            lo_m = np.where(np.isfinite(lo), lo + 1e-6 * (np.abs(lo) + 1), -np.inf)
            hi_m = np.where(np.isfinite(hi), hi - 1e-6 * (np.abs(hi) + 1), np.inf)
        return np.clip(x, lo_m, hi_m)


# ============================================================================ IRF paramétrique et moteur de convolution
def emg_density(t, mu, sig, a):
    if a <= 0:
        return np.exp(-0.5 * ((t - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    z = (sig / a - (t - mu) / sig) / np.sqrt(2)
    out = np.empty_like(t, dtype=float); pos = z >= 0; neg = ~pos
    out[pos] = 0.5 / a * np.exp(-0.5 * ((t[pos] - mu) / sig) ** 2) * erfcx(z[pos])
    out[neg] = 0.5 / a * np.exp(0.5 * (sig / a) ** 2 - (t[neg] - mu) / a) * erfc(z[neg])
    return out


class Fine:
    """Grille fine périodique : IRF paramétrique, exponentielles périodiques intégrées par canal, convolution FFT."""
    def __init__(self, T, NB, M=M_FINE):
        self.T, self.NB, self.M = T, NB, M
        self.NF = NB * M; self.dt = T / self.NF
        self.tf = (np.arange(self.NF) + 0.5) * self.dt; self.i = np.arange(self.NF)
        self.F = None

    def irf(self, mu, sig, weights, taus):
        out = np.zeros(self.NF)
        for wk, ak in zip(weights, [0.0] + list(taus)):
            for sh in (-self.T, 0.0, self.T):
                out += wk * emg_density(self.tf + sh, mu, sig, ak)
        return out / out.sum()

    def set_irf(self, irf_fine): self.F = np.fft.rfft(irf_fine)

    def conv_exp(self, tau):
        e = tau * (-np.expm1(-self.dt / tau)) * np.exp(-self.i * self.dt / tau) / (-np.expm1(-self.T / tau))
        return np.fft.irfft(np.fft.rfft(e) * self.F, n=self.NF).reshape(self.NB, self.M).sum(1)

    def coarse(self, v): return v.reshape(self.NB, self.M).sum(1)


def fit_irf(d, ch, ntail=2):
    """IRF brute = S * [gauss ⊗ (w0 δ + Σ w_k exp(a_k))] + fond. Ajustement Poisson, sans écrêtage."""
    y = d['irf_raw'][ch]; w = d['win'][ch]; pk = int(np.argmax(y)); E = Fine(d['T'], d['NB'])
    P = Params(); P.add('mu', d['t'][pk] - 0.03, d['t'][pk] - 0.3, d['t'][pk] + 0.3); P.add('sig', 0.02, 0.003, 0.2)
    for k, a0 in enumerate([0.08, 2.0, 0.5][:ntail]): P.add(f'a{k}', a0, 0.005, 8.0)
    for k in range(ntail + 1): P.add(f'w{k}', 1.0 if k == 0 else 0.3, 0, np.inf)
    P.add('S', (y[w] - np.median(y[w])).sum(), 0, np.inf); P.add('bg', np.median(y[w][:12]), 0, np.inf)
    def model(v):
        p = P.unpack(v); ws = np.array([p[f'w{k}'] for k in range(ntail + 1)]); ws /= ws.sum()
        return p['S'] * E.coarse(E.irf(p['mu'], p['sig'], ws, [p[f'a{k}'] for k in range(ntail)])) + p['bg']
    r = least_squares(lambda v: dev_resid(y[w], model(v)[w]), P.start, bounds=P.bounds, x_scale='jac', max_nfev=4000, xtol=1e-13, ftol=1e-13)
    p = P.unpack(r.x); ws = np.array([p[f'w{k}'] for k in range(ntail + 1)]); ws /= ws.sum()
    return dict(mu=p['mu'], sig=p['sig'], weights=ws.tolist(), taus=[p[f'a{k}'] for k in range(ntail)],
                dev_per_dof=float(np.sum(r.fun ** 2) / (w.sum() - len(r.x))))


# ============================================================================ modèle conjoint de calibration
def build_joint(d, irfp, G2=True):
    E = Fine(d['T'], d['NB']); NC = len(CONC)
    P = Params()
    P.add('t1', 0.09, 0.02, 0.5); P.add('t2', 0.9, 0.2, 3.0); P.add('tD', 3.75, 2.5, 6.0); P.add('tA', 3.45, 1.5, 6.0)
    P.add('lam', 0.75, 0, 5); P.add('G', 2.0, 0, 50)
    if G2: P.add('G2r', 1.0, 0, 5)
    P.add('k1', 0.8, 0, 50); P.add('k0', 0.01, 0, 50)
    for c in range(NC):
        P.add(f's1_{c}', 0.0, -0.3, 0.3); P.add(f's2_{c}', 0.0, -0.3, 0.3)
        P.add(f'A{c}', d['B'][c].max() / d['DT'], 0, np.inf)
        P.add(f'f1_{c}', 0.3, 0, 1); P.add(f'f2_{c}', 0.1, 0, 1)
        P.add(f'bgB{c}', 300.0, 0, np.inf); P.add(f'bgY{c}', 300.0, 0, np.inf)

    def model_fn(v, parts=False):
        p = P.unpack(v); taus = [p['t1'], p['t2']]
        Gs = [p['G'], p['G'] * p.get('G2r', 1.0)]
        kT = [1 / t - 1 / p['tD'] for t in taus]
        outB, outY, comp = [], [], []
        for c in range(NC):
            fr = [p[f'f1_{c}'], p[f'f2_{c}']]; fD = max(1 - sum(fr), 0); A = p[f'A{c}']
            E.set_irf(E.irf(irfp[1]['mu'] + p[f's1_{c}'], irfp[1]['sig'], irfp[1]['weights'], irfp[1]['taus']))
            bB = [E.conv_exp(t) for t in taus]; bD = E.conv_exp(p['tD'])
            muB = A * (fr[0] * bB[0] + fr[1] * bB[1] + fD * bD) + p[f'bgB{c}']
            E.set_irf(E.irf(irfp[2]['mu'] + p[f's2_{c}'], irfp[2]['sig'], irfp[2]['weights'], irfp[2]['taus']))
            yB = [E.conv_exp(t) for t in taus]; yD = E.conv_exp(p['tD']); yA = E.conv_exp(p['tA'])
            bleed = p['lam'] * A * (fr[0] * yB[0] + fr[1] * yB[1] + fD * yD)
            rises = [t * p['tA'] / (p['tA'] - t) * (yA - b) for t, b in zip(taus, yB)]
            sens = A * sum(G_ * f * k * r for G_, f, k, r in zip(Gs, fr, kT, rises))
            dirc = A * (p['k1'] * sum(fr) + p['k0']) * yA
            outB.append(muB); outY.append(bleed + sens + dirc + p[f'bgY{c}'])
            if parts: comp.append(dict(bleed=bleed, sens=sens, direct=dirc))
        return (np.array(outB), np.array(outY), comp) if parts else (np.array(outB), np.array(outY))
    return P, model_fn


def fit_joint(d, irfp, G2=True):
    P, model_fn = build_joint(d, irfp, G2)
    wB, wY = d['win'][1], d['win'][2]; NC = len(CONC)
    def resid(v):
        mB, mY = model_fn(v)
        return np.concatenate([dev_resid(d['B'][c][wB], mB[c][wB]) for c in range(NC)] +
                              [dev_resid(d['Y'][c][wY], mY[c][wY]) for c in range(NC)])
    r = least_squares(resid, P.start, bounds=P.bounds, x_scale='jac', max_nfev=20000, xtol=1e-14, ftol=1e-14, gtol=1e-12)
    dev = float(np.sum(r.fun ** 2)); dof = (wB.sum() + wY.sum()) * NC - len(r.x)
    J = approx_derivative(resid, r.x, method='3-point', bounds=P.bounds)
    cov = np.linalg.pinv(J.T @ J) * dev / dof          # erreurs-types gonflées par la surdispersion
    return dict(P=P, model_fn=model_fn, x=r.x, p=P.unpack(r.x), cov=cov, dev=dev, dof=int(dof))


def isotherm(c, fmax, fmin, K):
    return fmin + (fmax - fmin) / (1.0 + np.maximum(np.asarray(c, float), 0) / K)


def fit_isotherms(R):
    p, cov, idx = R['p'], R['cov'], {n: i for i, n in enumerate(R['P'].names)}
    NC = len(CONC); out = {}
    f1 = np.array([p[f'f1_{c}'] for c in range(NC)]); f2 = np.array([p[f'f2_{c}'] for c in range(NC)])
    s1 = np.array([np.sqrt(cov[idx[f'f1_{c}'], idx[f'f1_{c}']]) for c in range(NC)])
    s2 = np.array([np.sqrt(cov[idx[f'f2_{c}'], idx[f'f2_{c}']]) for c in range(NC)])
    sx = np.array([np.sqrt(cov[idx[f'f1_{c}'], idx[f'f1_{c}']] + cov[idx[f'f2_{c}'], idx[f'f2_{c}']] + 2 * cov[idx[f'f1_{c}'], idx[f'f2_{c}']]) for c in range(NC)])
    for name, y, sy in [('f1', f1, s1), ('f2', f2, s2), ('x', f1 + f2, sx)]:
        pp, cc = curve_fit(isotherm, CONC, y, p0=[y[0], 0.0, 60], sigma=sy, bounds=([0, -0.5, 1], [1, 1, 1e4]))
        out[name] = dict(values=y.tolist(), se=sy.tolist(), fmax=pp[0], fmin=pp[1], K=pp[2], se_params=np.sqrt(np.diag(cc)).tolist())
    return out


# ============================================================================ modèle direct du senseur (après calibration)
class Sensor:
    """Histogrammes attendus (bleu, jaune) en fonction de [Cl], de A (expression×excitation) et de η."""
    def __init__(self, d, consts, iso, irfp, shifts=(0.0, 0.0)):
        self.k, self.iso, self.win = consts, iso, d['win'][1]
        self.EB = Fine(d['T'], d['NB']); self.EY = Fine(d['T'], d['NB'])
        self.EB.set_irf(self.EB.irf(irfp[1]['mu'] + shifts[0], irfp[1]['sig'], irfp[1]['weights'], irfp[1]['taus']))
        self.EY.set_irf(self.EY.irf(irfp[2]['mu'] + shifts[1], irfp[2]['sig'], irfp[2]['weights'], irfp[2]['taus']))
        k = consts; self.taus = [k['t1'], k['t2']]
        self.bB = [self.EB.conv_exp(t) for t in self.taus]; self.bBD = self.EB.conv_exp(k['tD'])
        self.bY = [self.EY.conv_exp(t) for t in self.taus]; self.bYD = self.EY.conv_exp(k['tD']); self.bYA = self.EY.conv_exp(k['tA'])
        self.rise = [t * k['tA'] / (k['tA'] - t) * (self.bYA - b) for t, b in zip(self.taus, self.bY)]
        self.kT = [1 / t - 1 / k['tD'] for t in self.taus]
        self.Gs = [k['G'], k['G'] * k.get('G2r', 1.0)]

    def mu(self, c, A=1.0, eta=1.0):
        fr = [isotherm(c, self.iso['f1']['fmax'], self.iso['f1']['fmin'], self.iso['f1']['K']),
              isotherm(c, self.iso['f2']['fmax'], self.iso['f2']['fmin'], self.iso['f2']['K'])]
        fD = 1 - sum(fr); k = self.k
        muB = A * (fr[0] * self.bB[0] + fr[1] * self.bB[1] + fD * self.bBD)
        bleed = k['lam'] * A * (fr[0] * self.bY[0] + fr[1] * self.bY[1] + fD * self.bYD)
        sens = A * sum(G_ * f * kt * r for G_, f, kt, r in zip(self.Gs, fr, self.kT, self.rise))
        dirc = A * (k['k1'] * sum(fr) + k['k0']) * self.bYA
        return muB, eta * (bleed + sens + dirc)


def cramer_rao(S, c, N=1e4):
    """Écart-type minimal de [Cl] (mM) pour N photons (bleu+jaune), selon l'information utilisée."""
    w = S.win
    mB0, mY0 = S.mu(c); A = N / (mB0[w].sum() + mY0[w].sum())
    def f(th):
        b, y = S.mu(th[0], np.exp(th[1]), np.exp(th[2])); return b[w], y[w]
    th = np.array([c, np.log(A), 0.0]); h = [1e-3 * max(c, 1), 1e-4, 1e-4]
    mB, mY = f(th); dB = []; dY = []
    for i in range(3):
        e = np.zeros(3); e[i] = h[i]; bp, yp = f(th + e); bm, ym = f(th - e)
        dB.append((bp - bm) / (2 * h[i])); dY.append((yp - ym) / (2 * h[i]))
    dB, dY = np.array(dB), np.array(dY)
    FB = (dB / mB) @ dB.T; FY = (dY / mY) @ dY.T
    def schur(F, nuis):
        idx = [0] + nuis; Fs = F[np.ix_(idx, idx)]
        return Fs[0, 0] - Fs[0, 1:] @ np.linalg.pinv(Fs[1:, 1:]) @ Fs[1:, 0]
    return {'FLIM bleu seul': 1 / np.sqrt(schur(FB, [1])), 'FLIM jaune seul': 1 / np.sqrt(schur(FY, [1])),
            'FLIM 2 canaux (eta inconnu)': 1 / np.sqrt(schur(FB + FY, [1, 2])),
            'ratio + FLIM (eta connu)': 1 / np.sqrt(schur(FB + FY, [1]))}


# ============================================================================ filtre bayésien hybride (grille [Cl] × η)
class HybridGrid:
    """log L(c, η) par trame = formes (FLIM, indépendantes de A et η) + totaux (ratio, dépendent de η) ; A profilé."""
    def __init__(self, S, c_grid, eta_grid):
        self.c, self.eta, w = np.asarray(c_grid, float), np.asarray(eta_grid, float), S.win
        mB, mY = map(np.array, zip(*[(b[w], y[w]) for b, y in (S.mu(c) for c in self.c)]))
        self.mBtot, self.mYtot = mB.sum(1), mY.sum(1)
        self.logpB = np.log(mB / self.mBtot[:, None]); self.logpY = np.log(mY / self.mYtot[:, None])

    def frame_loglik(self, B, Y):
        """B, Y : histogrammes d'une trame (canaux de la fenêtre seulement)."""
        L = (self.logpB @ B + self.logpY @ Y)[:, None]
        mB = self.mBtot[:, None]; mY = self.mYtot[:, None] * self.eta[None, :]
        return L + B.sum() * np.log(mB / (mB + mY)) + Y.sum() * np.log(mY / (mB + mY))


def transition(logx, sd, jump_sd=None, jump_w=0.0):
    dd = logx[:, None] - logx[None, :]
    K = np.exp(-0.5 * (dd / sd) ** 2); K /= K.sum(0, keepdims=True)
    if jump_sd:
        K2 = np.exp(-0.5 * (dd / jump_sd) ** 2); K2 /= K2.sum(0, keepdims=True); K = (1 - jump_w) * K + jump_w * K2
    return K


def forward_backward(Ls, Kc, Ke):
    """Filtre (causal) puis lissage avant-arrière. Renvoie (alpha, post) de forme (T, Nc, Neta)."""
    T = len(Ls); Nc, Ne = Ls[0].shape; p = np.ones((Nc, Ne)) / (Nc * Ne)
    alpha = np.zeros((T, Nc, Ne)); preds = []
    for t in range(T):
        pred = Kc @ p @ Ke.T if t > 0 else p; preds.append(pred)
        p = pred * np.exp(Ls[t] - Ls[t].max()); p /= p.sum(); alpha[t] = p
    post = np.zeros_like(alpha); post[-1] = alpha[-1]
    for t in range(T - 2, -1, -1):
        post[t] = alpha[t] * (Kc.T @ (post[t + 1] / np.maximum(preds[t + 1], 1e-300)) @ Ke); post[t] /= post[t].sum()
    return alpha, post


# ============================================================================ version temps réel
def log_gates(n_bins, n_gates=16):
    """Matrice de portes contiguës, plus fines au début de la fenêtre (16 portes ≈ 93 % de l'info FLIM)."""
    edges = np.unique(np.round(np.geomspace(1, n_bins + 1, n_gates + 1) - 1).astype(int)); edges[0], edges[-1] = 0, n_bins
    G = np.zeros((len(edges) - 1, n_bins))
    for g in range(len(edges) - 1): G[g, edges[g]:edges[g + 1]] = 1
    return G


class RealTimeHybrid:
    """Filtre bayésien causal, exact quand eta est statique et partagé par un groupe de ROI (même profondeur) :
        p(c_1..c_M, eta | données) = p(eta | données) · Π_m p(c_m | eta, données_m)
    Chaque ROI garde une grille p(c | eta) (colonnes normalisées) ; p(eta) est commun au groupe.
    Quand sd(log eta) < collapse_sd, eta est figé et chaque ROI passe à un filtre 1D sur le ratio (≈ 15× plus rapide).
    """
    def __init__(self, S, n_roi, gates=None, c_grid=None, eta_grid=None, step_sd=0.04, jump_sd=0.5, jump_w=0.03,
                 eta_prior=None, collapse_sd=None):
        w = S.win
        self.c = np.geomspace(0.5, 300, 260) if c_grid is None else np.asarray(c_grid, float)
        self.eta = np.geomspace(0.6, 1.8, 61) if eta_grid is None else np.asarray(eta_grid, float)
        self.le = np.log(self.eta); self.G = gates
        mB, mY = map(np.array, zip(*[(b[w], y[w]) for b, y in (S.mu(c) for c in self.c)]))
        if gates is not None: mB, mY = mB @ gates.T, mY @ gates.T
        self.mBt, self.mYt = mB.sum(1), mY.sum(1)
        self.lpB, self.lpY = np.log(mB / self.mBt[:, None]), np.log(mY / self.mYt[:, None])
        self.Kc = transition(np.log(self.c), step_sd, jump_sd, jump_w)
        self.P = np.ones((n_roi, len(self.c), len(self.eta))) / len(self.c)
        self.log_pe = np.zeros(len(self.eta)) if eta_prior is None else \
            -0.5 * ((self.le - np.log(eta_prior[0])) / eta_prior[1]) ** 2       # eta_prior = (centre, sd en log)
        self.collapse_sd, self.eta_fix, self.pc, self.k = collapse_sd, None, None, 0

    def _totals(self, NB, NY, eta):
        mB = self.mBt[:, None]; mY = self.mYt[:, None] * np.atleast_1d(eta)[None, :]
        return NB * np.log(mB / (mB + mY)) + NY * np.log(mY / (mB + mY))

    def update(self, B, Y):
        """B, Y : (n_roi, n_canaux) — histogrammes de la trame (canaux de la fenêtre ; ou portes si gates est donné).
        Renvoie dict(c=moyenne, lo=2,5 %, hi=97,5 % par ROI, eta=moyenne, eta_sd=sd(log eta))."""
        B = np.atleast_2d(np.asarray(B, float)); Y = np.atleast_2d(np.asarray(Y, float))
        if self.G is not None and B.shape[1] != self.G.shape[0]: B, Y = B @ self.G.T, Y @ self.G.T
        M = B.shape[0]; out_pc = []
        for m in range(M):
            NB, NY = B[m].sum(), Y[m].sum()
            if self.eta_fix is None:
                if self.k > 0: self.P[m] = self.Kc @ self.P[m]
                L = (self.lpB @ B[m] + self.lpY @ Y[m])[:, None] + self._totals(NB, NY, self.eta)
                Q = self.P[m] * np.exp(L - L.max()); lik = Q.sum(0)
                self.P[m] = Q / np.maximum(lik, 1e-300); self.log_pe += np.log(np.maximum(lik, 1e-300))
            else:
                if self.k > 0: self.pc[m] = self.Kc @ self.pc[m]
                L = self._totals(NB, NY, self.eta_fix)[:, 0]
                self.pc[m] = self.pc[m] * np.exp(L - L.max()); self.pc[m] /= self.pc[m].sum()
        self.log_pe -= self.log_pe.max(); pe = np.exp(self.log_pe); pe /= pe.sum()
        mle = (pe * self.le).sum(); sle = np.sqrt((pe * (self.le - mle) ** 2).sum())
        if self.eta_fix is None:
            pcs = np.array([(self.P[m] * pe[None, :]).sum(1) for m in range(M)])
            if self.collapse_sd is not None and sle < self.collapse_sd:
                self.eta_fix, self.pc = np.exp(mle), pcs.copy()
        else:
            pcs = self.pc
        self.k += 1
        cdf = np.cumsum(pcs, 1)
        return dict(c=(pcs * self.c).sum(1), lo=self.c[np.argmax(cdf >= 0.025, 1)], hi=self.c[np.argmax(cdf >= 0.975, 1)],
                    eta=float(np.exp(mle)), eta_sd=float(sle))


def demo_temps_reel(d, S, n_roi=10, photons=3000, T=300, eta_true=1.25, seed=0):
    """Simulation : n_roi ROI, eta inconnu partagé, [Cl] propre à chaque ROI ; compare à l'oracle (eta connu)."""
    rng = np.random.default_rng(seed); w = S.win; t = np.arange(T); C = []
    for m in range(n_roi):
        base, amp, on = rng.uniform(6, 12), rng.uniform(8, 16), int(rng.uniform(60, 120)); off = on + int(rng.uniform(60, 100))
        c = np.full(T, base); c[on:] = base + amp * (1 - np.exp(-(t[on:] - on) / 4.0))
        c[off:] = base + (c[off - 1] - base) * np.exp(-(t[off:] - off) / 25.0); C.append(c)
    C = np.array(C); G = log_gates(int(w.sum()), 16)
    rt = RealTimeHybrid(S, n_roi, gates=G, collapse_sd=0.01)
    orc = RealTimeHybrid(S, n_roi, eta_grid=[eta_true])
    est = np.zeros((n_roi, T)); eo = np.zeros((n_roi, T)); tt = 0.0
    for k in range(T):
        B, Y = [], []
        for m in range(n_roi):
            b, y = S.mu(C[m, k], 1.0, eta_true); A = photons / (b[w].sum() + y[w].sum())
            B.append(rng.poisson(A * b[w])); Y.append(rng.poisson(A * y[w]))
        B, Y = np.array(B, float), np.array(Y, float)
        t0 = time.perf_counter(); r = rt.update(B @ G.T, Y @ G.T); tt += time.perf_counter() - t0
        est[:, k] = r['c']; eo[:, k] = orc.update(B, Y)['c']
    sl = slice(100, None)
    print(f"\nDémo temps réel : {n_roi} ROI, {photons} photons/trame, eta vrai {eta_true} (inconnu du filtre)")
    print(f"   eta estimé {r['eta']:.3f} ± {100*r['eta_sd']:.1f} % | RMSE trames 100+ : {np.sqrt(np.mean((est[:, sl]-C[:, sl])**2)):.2f} mM "
          f"(oracle eta connu {np.sqrt(np.mean((eo[:, sl]-C[:, sl])**2)):.2f} mM) | {1e3*tt/T/n_roi:.2f} ms par ROI et par trame")


# ============================================================================ programme principal
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('csv'); ap.add_argument('--out', default='calibration_sclm.json')
    ap.add_argument('--no-G2', action='store_true', help='un seul G pour les deux populations FRET')
    ap.add_argument('--demo-rt', action='store_true', help='simulation du filtre temps réel après la calibration')
    ap.add_argument('--photons', type=int, default=3000, help='photons par ROI et par trame pour --demo-rt')
    a = ap.parse_args()
    t0 = time.time()
    d = load(a.csv)
    irfp = {ch: fit_irf(d, ch, ntail=2 if ch == 1 else 3) for ch in (1, 2)}
    for ch in (1, 2):
        print(f"IRF voie{ch}: sigma={irfp[ch]['sig']*1e3:.1f} ps, queues={np.round(irfp[ch]['taus'], 3)} ns, poids={np.round(irfp[ch]['weights'], 3)}, dev/dof={irfp[ch]['dev_per_dof']:.2f}")
    R = fit_joint(d, irfp, G2=not a.no_G2)
    p = R['p']; se = dict(zip(R['P'].names, np.sqrt(np.maximum(np.diag(R['cov']), 0))))
    keys = ['t1', 't2', 'tD', 'tA', 'lam', 'G', 'G2r', 'k1', 'k0']
    consts = {k: float(p[k]) for k in keys if k in p}
    print(f"\nAjustement conjoint : déviance/ddl = {R['dev']:.0f}/{R['dof']} = {R['dev']/R['dof']:.3f}")
    for k in consts: print(f"   {k:4s} = {consts[k]:.4f} ± {se[k]:.4f}")
    NC = len(CONC)
    print("   décalages bleu (ps) :", np.round([p[f's1_{c}'] * 1e3 for c in range(NC)], 1))
    print("   décalages jaune (ps):", np.round([p[f's2_{c}'] * 1e3 for c in range(NC)], 1))
    iso = fit_isotherms(R)
    for k, v in iso.items():
        print(f"   {k:2s} = {np.round(v['values'], 4)}  -> fmax={v['fmax']:.4f} fmin={v['fmin']:.4f} K={v['K']:.1f} mM")
    S = Sensor(d, consts, iso, irfp, shifts=(np.mean([p[f's1_{c}'] for c in range(NC)]), np.mean([p[f's2_{c}'] for c in range(NC)])))
    print("\nÉcart-type minimal de [Cl] (mM) pour 10 000 photons :")
    for c in (5.0, 10.0, 20.0, 40.0):
        cr = cramer_rao(S, c)
        print(f"   [Cl]={c:5.1f} mM : " + " | ".join(f"{k} {v:6.2f}" for k, v in cr.items()))
    out = dict(constantes=consts, erreurs_types={k: float(se[k]) for k in consts}, isothermes=iso, irf=irfp,
               deviance=R['dev'], ddl=R['dof'],
               decalages_ps={'bleu': [float(p[f's1_{c}'] * 1e3) for c in range(NC)], 'jaune': [float(p[f's2_{c}'] * 1e3) for c in range(NC)]},
               note="eta = 1 par définition dans les conditions de la calibration ; lambda, G et kappa sont relatifs à ces conditions.")
    json.dump(out, open(a.out, 'w'), indent=2, default=float)
    print(f"\nCalibration écrite dans {a.out}  ({time.time()-t0:.0f} s)")
    if a.demo_rt:
        demo_temps_reel(d, S, photons=a.photons)


if __name__ == '__main__':
    main()