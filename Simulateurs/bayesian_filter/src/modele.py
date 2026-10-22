"""
modele.py — le modèle direct du senseur (sections 3 à 7 du document), construit à partir de calibration.json.

    bleu  : μ_B,i = A [IRF_B ⊗ d_B]_i,     d_B = f1 e^{-t/τ1} + f2 e^{-t/τ2} + (1 − f1 − f2) e^{-t/τD}
    jaune : μ_Y,i = η A [IRF_Y ⊗ (λ d_B + G (f1 k1 s1 + g2 f2 k2 s2) + κ x e^{-t/τA})]_i
            s_j(t) = τj τA/(τA − τj) (e^{-t/τA} − e^{-t/τj}),  k_j = 1/τj − 1/τD,  x = f1 + f2
    f_j(c) = f_j,min + (f_j,max − f_j,min) / (1 + c/K_j)                                (isothermes, section 6)

Modèle périodique (traîne des impulsions précédentes, section 3), IRF paramétrique (gaussienne ⊗ traînes
exponentielles) décalée de s, convolution sur une grille fine de 16 sous-cases par case. Les μ rendus sont
restreints à la fenêtre d'analyse (cases 13 à 242). Le fond n'est pas dans le modèle : le simulateur l'ajoute.
"""
import json, os, tomllib
import numpy as np
from scipy.special import erfc, erfcx

ICI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOUS_CASES = 16


def charger_params(chemin=None):
    with open(chemin or os.path.join(ICI, 'params.toml'), 'rb') as f:
        return tomllib.load(f)


def charger_calibration(chemin=None):
    chemin = chemin or os.path.join(ICI, 'calibration.json')
    if not os.path.exists(chemin):
        raise FileNotFoundError(f"{chemin} absent : lancer d'abord python3 step0_calibration.py")
    with open(chemin) as f:
        return json.load(f)


def isotherme(c, fmax, fmin, K):
    return fmin + (fmax - fmin) / (1.0 + np.maximum(np.asarray(c, float), 0.0) / K)


def x_vers_c(x, iso):
    """Section 6 : c = K (x_max − x) / (x − x_min), sur l'isotherme de x = f1 + f2."""
    return iso['K'] * (iso['fmax'] - x) / (x - iso['fmin'])


def _emg(t, mu, sig, a):
    """Densité gaussienne (a = 0) ou gaussienne ⊗ exponentielle de constante a (ns)."""
    if a <= 0:
        return np.exp(-0.5 * ((t - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    z = (sig / a - (t - mu) / sig) / np.sqrt(2)
    out = np.empty_like(t); pos = z >= 0; neg = ~pos
    out[pos] = 0.5 / a * np.exp(-0.5 * ((t[pos] - mu) / sig) ** 2) * erfcx(z[pos])
    out[neg] = 0.5 / a * np.exp(0.5 * (sig / a) ** 2 - (t[neg] - mu) / a) * erfc(z[neg])
    return out


class GrilleFine:
    """Grille fine sur une période ; IRF paramétrique et déclins périodiques convolués, intégrés par case."""
    def __init__(self, T, n_cases, irf, decalage_ns=0.0):
        self.T, self.n, self.NF = T, n_cases, n_cases * SOUS_CASES
        self.dt = T / self.NF
        tf = (np.arange(self.NF) + 0.5) * self.dt
        h = np.zeros(self.NF)
        for wk, ak in zip(irf['poids'], [0.0] + list(irf['taus'])):
            for per in (-T, 0.0, T):
                h += wk * _emg(tf + per, irf['mu'] + decalage_ns, irf['sig'], ak)
        self.F = np.fft.rfft(h / h.sum())

    def declin(self, tau):
        """[IRF ⊗ d_pér]_i pour d(t) = e^{-t/τ} répété toutes les T ns (section 3), intégré sur chaque case."""
        i = np.arange(self.NF)
        e = tau * (-np.expm1(-self.dt / tau)) * np.exp(-i * self.dt / tau) / (-np.expm1(-self.T / tau))
        return np.fft.irfft(np.fft.rfft(e) * self.F, n=self.NF).reshape(self.n, SOUS_CASES).sum(1)


class Senseur:
    """
    Histogrammes attendus (bleu, jaune) dans la fenêtre, en fonction de c, A et η.

    decalage_ps         : décalage d'horloge ajouté aux décalages moyens de la calibration (section 3) ;
                          un nombre (les deux canaux) ou (bleu, jaune)
    echelle_tau         : durées de vie multipliées par ce facteur (perturbation de la section 16)
    part_active         : fraction des accepteurs fonctionnels relativement à la calibration (section 16)
    """
    def __init__(self, cal, decalage_ps=0.0, echelle_tau=1.0, part_active=1.0):
        self.cal, self.k, self.iso = cal, dict(cal['constantes']), cal['isothermes']
        self.part_active = part_active
        T, n = cal['periode_ns'], cal['n_cases']
        a, b = cal['fenetre']; self.fen = slice(a, b + 1); self.n_fen = b + 1 - a
        dm = cal['decalage_moyen_ps']; sB, sY = np.broadcast_to(np.asarray(decalage_ps, float), (2,))
        self.gB = GrilleFine(T, n, cal['irf']['1'], 1e-3 * (dm['bleu'] + sB))
        self.gY = GrilleFine(T, n, cal['irf']['2'], 1e-3 * (dm['jaune'] + sY))
        k = self.k
        self.taus = [echelle_tau * k['t1'], echelle_tau * k['t2']]
        tD, tA = echelle_tau * k['tD'], echelle_tau * k['tA']
        f = self.fen
        self.bB = np.array([self.gB.declin(t)[f] for t in self.taus]); self.bBD = self.gB.declin(tD)[f]
        bY = [self.gY.declin(t)[f] for t in self.taus]; self.bYD = self.gY.declin(tD)[f]
        self.bYA = self.gY.declin(tA)[f]; self.bY = np.array(bY)
        # émission sensibilisée : G g_j k_j s_j, avec k_j = 1/τj − 1/τD (section 5)
        Gs = [k['G'], k['G'] * k['G2r']]
        self.sens = np.array([G_ * (1 / t - 1 / tD) * t * tA / (tA - t) * (self.bYA - b)
                              for G_, t, b in zip(Gs, self.taus, bY)])

    def fractions(self, c):
        a = self.part_active
        f1 = a * isotherme(c, self.iso['f1']['fmax'], self.iso['f1']['fmin'], self.iso['f1']['K'])
        f2 = a * isotherme(c, self.iso['f2']['fmax'], self.iso['f2']['fmin'], self.iso['f2']['K'])
        return f1, f2

    def mu(self, c, A=1.0, eta=1.0, composantes=False):
        """c scalaire ou vecteur (n,) → μ_B, μ_Y de forme (n_fen,) ou (n, n_fen)."""
        c = np.asarray(c, float); f1, f2 = self.fractions(c)
        f1, f2 = f1[..., None], f2[..., None]; fD = 1 - f1 - f2; k = self.k
        dB = f1 * self.bB[0] + f2 * self.bB[1] + fD * self.bBD
        fuite = k['lam'] * (f1 * self.bY[0] + f2 * self.bY[1] + fD * self.bYD)
        sens = f1 * self.sens[0] + f2 * self.sens[1]
        directe = (k['k1'] * (f1 + f2) + k['k0']) * self.bYA
        muB, muY = A * dB, eta * A * (fuite + sens + directe)
        if composantes:
            return muB, muY, dict(fuite=eta * A * fuite, sensibilisee=eta * A * sens, directe=eta * A * directe)
        return muB, muY

    def ratio_calibration(self, c):
        """R₁(c) = M_Y(c) / M_B(c) (section 7)."""
        muB, muY = self.mu(c)
        return muY.sum(-1) / muB.sum(-1)

    def decroissance(self, poids_taus, canal):
        """Forme normalisée d'une somme d'exponentielles [[poids, τ], …] vue par l'IRF d'un canal (autofluorescence)."""
        g = self.gB if canal == 'bleu' else self.gY
        h = sum(wt * g.declin(t)[self.fen] for wt, t in poids_taus)
        return h / h.sum()
