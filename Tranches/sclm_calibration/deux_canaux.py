"""Modèle FLIM à deux canaux de Simulateurs/bayesian_filter (sclm_hybride.py), appliqué à la tranche.

Le modèle mathématique est celui du simulateur, sans changement (« Mesurer le chlore avec SuperClomeleon »,
sections 3 à 7) :

    bleu  : μ_B = A [IRF_B ⊗ (f1 e^{-t/τ1} + f2 e^{-t/τ2} + (1 − f1 − f2) e^{-t/τD})] + fond_B
    jaune : μ_Y = η A [IRF_Y ⊗ (λ·donneur + G Σ g_j f_j k_j (e^{-t/τj} ⊗ e^{-t/τA}) + (κ (f1 + f2) + k0) e^{-t/τA})] + fond_Y
            k_j = 1/τj − 1/τD, g_1 = 1, g_2 = G2r ;  f_j(c) = f_j,min + (f_j,max − f_j,min) / (1 + c/K_j)

modèle périodique, IRF paramétrique (gaussienne ⊗ traînes exponentielles) ajustée sur l'IRF mesurée, convolution
sur une grille fine de 16 sous-cases. On réutilise les fonctions de sclm_hybride.py (IRF, grille fine, `Sensor`,
isothermes) ; seul l'ajustement conjoint est recopié ici, avec des bornes de durées de vie élargies pour le tissu
(τD 1–6 ns au lieu de 2,5–6 : le donneur en tranche décline plus vite qu'en solution, sclm_flim/step8).

Lecture d'un soma (`LecteurSoma`), par maximum de vraisemblance sur une grille de c, sans filtre temporel :
pour chaque c, chaque canal reçoit son amplitude a et son fond β libres (Poisson, Newton sur les deux paramètres
linéaires), d'où la vraisemblance des **formes** ℓ_B(c) + ℓ_Y(c). La **quantité** (le ratio) est le rapport des
amplitudes : ln(â_Y/â_B) ≈ ln η + ln R₁(c), d'écart-type tiré de l'information de Fisher des deux ajustements —
le découpage « formes + totaux » de la section 11, avec un fond libre.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.optimize._numdiff import approx_derivative

import config as C

SCLM_HYBRIDE = (Path(__file__).resolve().parents[2] / "Simulateurs" / "FLIM-FRET_2-channels" / "sorties"
                / "declins_purs_2026-10-16_SCLM" / "sclm_hybride.py")


def module_hybride():
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("sclm_hybride", SCLM_HYBRIDE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


SH = module_hybride()


def _isotherme(c, fmax, fmin, K):
    """Isotherme de sclm_hybride **sans** le plancher c ≥ 0 : une [Cl⁻] lue négative est gardée (la tronquer
    biaiserait les moyennes vers le haut). Définie pour c > −K."""
    return fmin + (fmax - fmin) / (1.0 + np.asarray(c, float) / K)


SH.isotherm = _isotherme          # utilisée par Sensor.mu au moment de l'appel


def c_min_modele(iso, marge=1e-3):
    """Plus petite [Cl⁻] (négative) où le modèle garde un sens : f1, f2 ≥ 0 et f1 + f2 ≤ 1 (fraction du donneur
    sans FRET ≥ 0), et c > −K de chaque isotherme. Au-delà, le prolongement de l'isotherme n'est plus un déclin."""
    K = [v["K"] for v in iso.values() if v["fmax"] != v["fmin"]]      # une isotherme constante (f2 = 0) ne borne rien
    c = -np.linspace(0, 0.999 * min(K), 20001)
    f1 = _isotherme(c, iso["f1"]["fmax"], iso["f1"]["fmin"], iso["f1"]["K"])
    f2 = _isotherme(c, iso["f2"]["fmax"], iso["f2"]["fmin"], iso["f2"]["K"])
    ok = (f1 >= 0) & (f2 >= 0) & (f1 + f2 <= 1) & (f1 <= 1) & (f2 <= 1)
    k = np.argmin(ok) if not ok.all() else len(c)          # premier c qui sort du domaine
    return float(c[max(k - 1, 0)] * (1 - marge))


# ------------------------------------------------------------------ instrument
def irf_brute(root, canal):
    """IRF mesurée (comptes, fond compris), ramenée sur les 256 canaux de la période par sa fonction cumulée."""
    import sdtfile
    s = sdtfile.SdtFile(C.irf_sdt(root))
    y = s.data[0][C.IRF_LIGNE[canal]].astype(float)
    t = s.times[0] * 1e9
    dt = t[1] - t[0]
    bords = np.r_[t, t[-1] + dt]
    cum = np.r_[0, np.cumsum(y)]
    g = np.arange(C.N_CANAUX + 1) * C.DT_NS
    return np.diff(np.interp(g, bords, cum)) + np.diff(np.interp(g + C.PERIODE_NS, bords, cum, right=cum[-1]))


def donnees(root, B, Y, conc, fenetre):
    """Le dictionnaire `d` de sclm_hybride (voie 1 = bleu = CFP, voie 2 = jaune = YFP)."""
    win = np.zeros(C.N_CANAUX, bool)
    win[fenetre[0]:fenetre[1]] = True
    return dict(T=C.PERIODE_NS, NB=C.N_CANAUX, DT=C.DT_NS, t=C.temps_canaux(), win={1: win, 2: win.copy()},
                irf_raw={1: irf_brute(root, "cfp"), 2: irf_brute(root, "yfp")},
                B=np.asarray(B, float), Y=np.asarray(Y, float), conc=np.asarray(conc, float))


def ajuster_irf(d):
    return {1: SH.fit_irf(d, 1, ntail=2), 2: SH.fit_irf(d, 2, ntail=3)}


# ------------------------------------------------------------------ calibration (ajustement conjoint)
def construire(d, irfp, n_transfert=2):
    """build_joint de sclm_hybride, bornes des τ élargies ; une « acquisition » par [Cl⁻] imposée.
    n_transfert = 1 : une seule population en FRET (τ1) — f2, τ2 et g2 retirés du modèle."""
    E = SH.Fine(d["T"], d["NB"])
    NC = len(d["conc"])
    P = SH.Params()
    deux = n_transfert == 2
    P.add("t1", 0.06 if deux else 0.4, 0.005, 3.0 if not deux else 0.5)
    if deux:
        P.add("t2", 0.5, 0.15, 3.0)
    P.add("tD", 2.5, 1.0, 6.0); P.add("tA", 3.2, 1.0, 6.0)
    P.add("lam", 0.75, 0, 5); P.add("G", 2.0, 0, 50)
    if deux:
        P.add("G2r", 1.0, 0, 5)
    P.add("k1", 0.8, 0, 50); P.add("k0", 0.01, 0, 50)
    for c in range(NC):
        P.add(f"s1_{c}", 0.0, -0.3, 0.3); P.add(f"s2_{c}", 0.0, -0.3, 0.3)
        P.add(f"A{c}", d["B"][c].max() / d["DT"], 0, np.inf)
        P.add(f"f1_{c}", 0.3 if deux else 0.6, 0, 1)
        if deux:
            P.add(f"f2_{c}", 0.2, 0, 1)
        P.add(f"bgB{c}", 0.02 * d["B"][c].mean(), 0, np.inf); P.add(f"bgY{c}", 0.02 * d["Y"][c].mean(), 0, np.inf)

    def model_fn(v, parts=False):
        p = P.unpack(v)
        taus = [p["t1"], p.get("t2", 1.0)]
        Gs = [p["G"], p["G"] * p.get("G2r", 0.0)]
        kT = [1 / t - 1 / p["tD"] for t in taus]
        outB, outY, comp = [], [], []
        for c in range(NC):
            fr = [p[f"f1_{c}"], p.get(f"f2_{c}", 0.0)]
            fD = max(1 - sum(fr), 0)
            A = p[f"A{c}"]
            E.set_irf(E.irf(irfp[1]["mu"] + p[f"s1_{c}"], irfp[1]["sig"], irfp[1]["weights"], irfp[1]["taus"]))
            bB = [E.conv_exp(t) for t in taus]
            bD = E.conv_exp(p["tD"])
            muB = A * (fr[0] * bB[0] + fr[1] * bB[1] + fD * bD) + p[f"bgB{c}"]
            E.set_irf(E.irf(irfp[2]["mu"] + p[f"s2_{c}"], irfp[2]["sig"], irfp[2]["weights"], irfp[2]["taus"]))
            yB = [E.conv_exp(t) for t in taus]
            yD = E.conv_exp(p["tD"])
            yA = E.conv_exp(p["tA"])
            fuite = p["lam"] * A * (fr[0] * yB[0] + fr[1] * yB[1] + fD * yD)
            montees = [t * p["tA"] / (p["tA"] - t) * (yA - b) for t, b in zip(taus, yB)]
            sens = A * sum(G_ * f * k * r for G_, f, k, r in zip(Gs, fr, kT, montees))
            dirc = A * (p["k1"] * sum(fr) + p["k0"]) * yA
            outB.append(muB)
            outY.append(fuite + sens + dirc + p[f"bgY{c}"])
            if parts:
                comp.append(dict(donneur=muB - p[f"bgB{c}"], fuite=fuite, sensibilisee=sens, directe=dirc))
        return (np.array(outB), np.array(outY), comp) if parts else (np.array(outB), np.array(outY))
    return P, model_fn


def calibrer(d, irfp, n_transfert=2):
    """fit_joint de sclm_hybride : ajustement de Poisson conjoint des histogrammes bleus et jaunes."""
    P, model_fn = construire(d, irfp, n_transfert)
    w = d["win"][1]
    NC = len(d["conc"])

    def resid(v):
        mB, mY = model_fn(v)
        return np.concatenate([SH.dev_resid(d["B"][c][w], mB[c][w]) for c in range(NC)]
                              + [SH.dev_resid(d["Y"][c][w], mY[c][w]) for c in range(NC)])
    r = least_squares(resid, P.start, bounds=P.bounds, x_scale="jac", max_nfev=20000, xtol=1e-14, ftol=1e-14, gtol=1e-12)
    dev = float(np.sum(r.fun ** 2))
    dof = 2 * w.sum() * NC - len(r.x)
    J = approx_derivative(resid, r.x, method="3-point", bounds=P.bounds)
    cov = np.linalg.pinv(J.T @ J) * dev / dof
    p = P.unpack(r.x)
    if n_transfert == 1:                       # complète p pour Sensor : population 2 absente
        p.update(t2=1.0, G2r=0.0, **{f"f2_{c}": 0.0 for c in range(NC)})
    return dict(P=P, model_fn=model_fn, x=r.x, p=p, cov=cov, dev=dev, dof=int(dof), resid=r.fun,
                n_transfert=n_transfert)


def isothermes(R, conc):
    """fit_isotherms de sclm_hybride, sur les concentrations du jeu ; avec une seule population en FRET,
    l'isotherme de f1 seule (f2 = 0)."""
    SH.CONC = np.asarray(conc, float)
    if R.get("n_transfert", 2) == 2:
        return SH.fit_isotherms(R)
    from scipy.optimize import curve_fit
    y = np.array([R["p"][f"f1_{c}"] for c in range(len(conc))])
    pp, _ = curve_fit(SH.isotherm, np.asarray(conc, float), y, p0=[y[0], y[-1] * 0.8, 60],
                      bounds=([0, -0.5, 0.1], [1, 1, 1e4]), maxfev=20000)
    f1 = dict(values=y.tolist(), fmax=pp[0], fmin=pp[1], K=pp[2])
    return {"f1": f1, "f2": dict(values=[0.0] * len(conc), fmax=0.0, fmin=0.0, K=1.0), "x": dict(f1)}


CLES = ("t1", "t2", "tD", "tA", "lam", "G", "G2r", "k1", "k0")


def senseur(d, R, iso):
    """Modèle direct calibré (`Sensor` de sclm_hybride), décalages d'horloge moyens de la calibration."""
    NC = len(d["conc"])
    s1 = np.mean([R["p"][f"s1_{c}"] for c in range(NC)])
    s2 = np.mean([R["p"][f"s2_{c}"] for c in range(NC)])
    return SH.Sensor(d, {k: R["p"][k] for k in CLES}, iso, d["irfp"], shifts=(s1, s2))


# ------------------------------------------------------------------ lecture d'un soma
def ajuster_fond(y, P, iters=25):
    """Pour chaque forme P[k] (normalisée, n canaux) : μ = a·P[k] + β/n, a et β ≥ 0 libres (Poisson, Newton).
    Rend (log-vraisemblance (K,), a (K,), var(ln a) (K,), β (K,)) ; constantes de y omises."""
    y = np.asarray(y, float)
    K, n = P.shape
    u = 1.0 / n
    N = y.sum()
    a = np.full(K, 0.97 * N)
    b = np.full(K, 0.03 * N)
    for _ in range(iters):
        mu = np.maximum(a[:, None] * P + b[:, None] * u, 1e-12)
        r = y[None, :] / mu
        ga = (r * P).sum(1) - 1.0                     # Σ p = 1
        gb = (r * u).sum(1) - 1.0                     # Σ u = 1
        w = y[None, :] / mu ** 2
        haa = (w * P * P).sum(1)
        hab = (w * P * u).sum(1)
        hbb = (w * u * u).sum(1)
        det = haa * hbb - hab ** 2
        da = (hbb * ga - hab * gb) / det
        db = (haa * gb - hab * ga) / det
        # pas amorti : garder a, b > 0
        t = np.ones(K)
        for x, dx in ((a, da), (b, db)):
            neg = dx < 0
            t[neg] = np.minimum(t[neg], 0.9 * x[neg] / -dx[neg])
        a = np.maximum(a + t * da, 1e-9)
        b = np.maximum(b + t * db, 1e-9)
    mu = np.maximum(a[:, None] * P + b[:, None] * u, 1e-12)
    ll = (y[None, :] * np.log(mu) - mu).sum(1)
    # information de Fisher attendue → var(a)
    f = 1.0 / mu
    iaa, iab, ibb = (f * P * P).sum(1), (f * P * u).sum(1), (f * u * u).sum(1)
    var_a = ibb / (iaa * ibb - iab ** 2)
    return ll, a, var_a / a ** 2, b


class LecteurSoma:
    """Formes calibrées sur la grille de c (et de τ pour le temps de vie unique), dans la fenêtre."""

    def __init__(self, S, c_grille, tau_grille):
        w = S.win
        self.c, self.tau = np.asarray(c_grille, float), np.asarray(tau_grille, float)
        mB, mY = map(np.array, zip(*[(b[w], y[w]) for b, y in (S.mu(c) for c in self.c)]))
        self.MB, self.MY = mB.sum(1), mY.sum(1)
        self.PB, self.PY = mB / self.MB[:, None], mY / self.MY[:, None]
        self.lnR1 = np.log(self.MY / self.MB)                   # R₁(c) à η = 1
        T = np.array([S.EB.conv_exp(t)[w] for t in self.tau])
        self.PT = T / T.sum(1, keepdims=True)

    def lire(self, B, Y):
        """Profils sur la grille de c pour un soma : formes ℓ_B + ℓ_Y, ln(â_Y/â_B), sa variance, fonds."""
        lB, aB, vB, bB = ajuster_fond(B, self.PB)
        lY, aY, vY, bY = ajuster_fond(Y, self.PY)
        return dict(formes=lB + lY, ln_q=np.log(aY / aB), v=vB + vY, aB=aB, aY=aY, bB=bB, bY=bY)

    def hybride(self, prof, ln_eta):
        """log L(c) à η donné : formes + quantité (ln â_Y/â_B autour de ln η + ln R₁(c))."""
        return prof["formes"] - 0.5 * (prof["ln_q"] - ln_eta - self.lnR1) ** 2 / prof["v"]

    def tau_unique(self, B):
        """Une exponentielle (⊗ IRF bleue paramétrique, périodique) et un fond libre : τ̂ (ns)."""
        ll, _, _, _ = ajuster_fond(B, self.PT)
        return argmax_affine(ll, self.tau)


def argmax_affine(L, grille):
    """Maximum d'un profil, affiné par une parabole sur les trois voisins (comme estimateurs.py)."""
    i = int(np.argmax(L))
    if 0 < i < len(L) - 1:
        g, m, d = L[i - 1], L[i], L[i + 1]
        k = g - 2 * m + d
        delta = float(np.clip(0.5 * (g - d) / k, -0.5, 0.5)) if k < 0 else 0.0
        return float(np.interp(i + delta, np.arange(len(grille)), grille))
    return float(grille[i])


def intervalle(L, grille, seuil=1.92):
    """Intervalle de vraisemblance profilée (ΔlnL = 1,92 : ~95 %) sur la grille."""
    ok = L >= L.max() - seuil
    idx = np.nonzero(ok)[0]
    return float(grille[idx[0]]), float(grille[idx[-1]])
