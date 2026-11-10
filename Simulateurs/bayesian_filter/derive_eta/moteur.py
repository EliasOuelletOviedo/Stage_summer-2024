"""
Filtre temps réel à M ROI sur une grille (c, θ), θ = état de η partagé par les ROI.
- η « marche aléatoire » : θ = η sur une grille, transition Ke à chaque trame ; pour M > 1, après la transition de η, la
  loi est reprojetée sur p(η) Π_m p(c_m | η) (marginales exactes par ROI ; exact pour M = 1).
- η « dérive paramétrique » : θ = (η0, D) statiques, η_k(θ) = η0 · exp(D · k/(T−1)) ; pas de transition (exact pour tout M).
Le ratio : NB·lB(c, η) + NY·lY(c, η), tables précalculées sur une grille fine de η (plus proche voisin).
Renvoie, par ROI et par trame, la marginale de c, et les moments de ln η, et ln p(trame k | 1..k−1) (évidence).
"""
import os, sys, numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # bayesian_filter
from src.filtre import transition

class Ratio:
    def __init__(self, g, n_fin=4001, eta_min=0.5, eta_max=2.0):
        self.le = np.linspace(np.log(eta_min), np.log(eta_max), n_fin)
        mB = g.MB[:, None]; mY = g.MY[:, None] * np.exp(self.le)[None, :]
        self.lB = np.log(mB / (mB + mY)); self.lY = np.log(mY / (mB + mY))
    def idx(self, ln_eta):
        return np.clip(np.rint((ln_eta - self.le[0]) / (self.le[1] - self.le[0])).astype(int), 0, len(self.le) - 1)

def filtrer(formes, NB, NY, ratio, Kc, ln_eta_k, Ke=None):
    """
    formes : (M, T, n_c) ; NB, NY : (M, T) ; ln_eta_k(k) -> (n_θ,) valeurs de ln η à la trame k ; Ke : (n_θ, n_θ) ou None.
    Renvoie pc (M, T, n_c) float32, m_eta (T,), v_eta (T,), ln_ev (T,).
    """
    M, T, nc = formes.shape; le0 = ln_eta_k(0); nt = len(le0)
    Kc = Kc.astype(np.float32); Ke = None if Ke is None else Ke.astype(np.float32)
    lBf, lYf = ratio.lB.astype(np.float32), ratio.lY.astype(np.float32)
    P = np.full((M, nc, nt), 1.0 / nc, np.float32); pe = np.full(nt, 1.0 / nt)
    pc = np.empty((M, T, nc), np.float32); m_eta = np.empty(T); v_eta = np.empty(T); ln_ev = np.empty(T)
    for k in range(T):
        lek = ln_eta_k(k); ix = ratio.idx(lek)
        lB, lY = lBf[:, ix], lYf[:, ix]
        if k > 0:
            for m in range(M):
                P[m] = Kc @ P[m]
            if Ke is not None:
                pe_n = Ke.astype(np.float64) @ pe
                for m in range(M):
                    P[m] = ((P[m] * pe[None, :].astype(np.float32)) @ Ke.T) / np.maximum(pe_n[None, :], 1e-30).astype(np.float32)
                pe = pe_n / pe_n.sum()
        lse = np.log(np.maximum(pe, 1e-300)); tot = 0.0
        for m in range(M):
            L = formes[m, k].astype(np.float32)[:, None] + np.float32(NB[m, k]) * lB + np.float32(NY[m, k]) * lY
            mx = L.max(); Q = P[m] * np.exp(L - mx); ell = np.maximum(Q.sum(0, dtype=np.float64), 1e-300)
            P[m] = Q / np.maximum(ell, 1e-30).astype(np.float32); lse = lse + np.log(ell); tot += float(mx)   # colonne vide : reste à 0
        mx2 = lse.max(); w = np.exp(lse - mx2); z = w.sum(); pe = w / z
        ln_ev[k] = np.log(z) + mx2 + tot
        for m in range(M):
            pc[m, k] = P[m] @ pe.astype(np.float32)
        m_eta[k] = (pe * lek).sum(); v_eta[k] = (pe * (lek - m_eta[k]) ** 2).sum()
    return pc, m_eta, v_eta, ln_ev

def moyenne_modeles(sorties):
    """Moyenne bayésienne causale sur des modèles globaux : poids de la trame k ∝ p(trames 1..k | modèle)."""
    cum = np.array([np.cumsum(s[3]) for s in sorties]).T
    w = np.exp(cum - cum.max(1, keepdims=True)); w /= w.sum(1, keepdims=True)
    pc = sum(w[None, :, i, None] * s[0] for i, s in enumerate(sorties))
    m = sum(w[:, i] * s[1] for i, s in enumerate(sorties))
    v = sum(w[:, i] * (s[2] + s[1] ** 2) for i, s in enumerate(sorties)) - m ** 2
    return pc, m, np.maximum(v, 0), w
