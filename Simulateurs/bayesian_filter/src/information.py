"""
information.py — combien d'information contient une trame (section 10 du document).

Information de Fisher de Poisson sur θ = (c, ln A, ln η) :  I = Σ_i (1/μ_i) ∂μ_i/∂θ ∂μ_i/∂θᵀ, par canal.
Les nuisances se retirent par complément de Schur : I_eff = I_cc − I_cν I_νν⁻¹ I_νc. Borne : σ_c ≥ 1/√I_eff.
"""
import numpy as np


def fisher(S, c, photons):
    """Matrices de Fisher (bleu, jaune) 3×3 sur (c, ln A, ln η), pour `photons` attendus (bleu + jaune) à η = 1."""
    b0, y0 = S.mu(c); A = photons / (b0.sum() + y0.sum())
    th = np.array([c, np.log(A), 0.0]); h = np.array([1e-3 * max(c, 1.0), 1e-4, 1e-4])
    f = lambda t: S.mu(t[0], np.exp(t[1]), np.exp(t[2]))
    mB, mY = f(th); dB, dY = [], []
    for i in range(3):
        e = np.zeros(3); e[i] = h[i]
        (bp, yp), (bm, ym) = f(th + e), f(th - e)
        dB.append((bp - bm) / (2 * h[i])); dY.append((yp - ym) / (2 * h[i]))
    dB, dY = np.array(dB), np.array(dY)
    return (dB / mB) @ dB.T, (dY / mY) @ dY.T


def efficace(F, garder, nuisances):
    """Information effective sur les paramètres `garder`, les `nuisances` retirées (complément de Schur)."""
    g, n = list(garder), list(nuisances)
    Fgg = F[np.ix_(g, g)]
    if not n:
        return Fgg
    return Fgg - F[np.ix_(g, n)] @ np.linalg.pinv(F[np.ix_(n, n)]) @ F[np.ix_(n, g)]


def bornes(S, c, photons):
    """σ_c minimal (mM) selon l'information utilisée (tableau « bilan pour 10 000 photons » de la section 10)."""
    FB, FY = fisher(S, c, photons)
    sd = lambda F, nu: float(1 / np.sqrt(efficace(F, [0], nu)[0, 0]))
    # ratio seul : information d'un partage entre deux couleurs, N (∂p_Y/∂c)² / (p_Y (1 − p_Y))
    R = lambda x: float(S.ratio_calibration(x)); hc = 1e-3 * max(c, 1.0)
    pY = R(c) / (1 + R(c)); dp = (R(c + hc) / (1 + R(c + hc)) - R(c - hc) / (1 + R(c - hc))) / (2 * hc)
    return {'forme bleu seule': sd(FB, [1]),
            'forme jaune seule': sd(FY, [1]),
            'formes 2 canaux (eta inconnu)': sd(FB + FY, [1, 2]),
            'ratio (eta connu)': float(np.sqrt(pY * (1 - pY) / photons) / abs(dp)),
            'ratio + formes (eta connu)': sd(FB + FY, [1])}


def sigma_eta(S, c, photons):
    """sd(ln η) minimal quand c et A sont inconnus : ce que les photons mis en commun disent de η (section 14)."""
    FB, FY = fisher(S, c, photons)
    return float(1 / np.sqrt(efficace(FB + FY, [2], [0, 1])[0, 0]))
