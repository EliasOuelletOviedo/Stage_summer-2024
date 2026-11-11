"""Lectures candidates de la [Cl⁻] à partir d'une paire de déclins (CFP, YFP) — banc d'essai de l'étape 12.

Chaque lecture prend des déclins complets (M, 256) et rend un nombre par ligne, monotone avec la [Cl⁻] si elle
marche. Familles (les noms sont ceux des colonnes de results/banc_lectures.csv) :

« combien » — rapports de photons
    R            YFP / CFP, tous les photons de la fenêtre
    R_tard       le même, photons arrivés plus de 1 ns après le pic (écarte ce qui décline très vite)
« quand », canal CFP — forme du déclin du donneur
    tm_B         temps moyen d'arrivée (aucun ajustement, aucune IRF)
    tm_queue_B   temps moyen des photons arrivés plus de 0,6 ns après le pic
    tm_B_sf      temps moyen, fond plat retiré (le fond est la composante uniforme du mélange à trois populations)
    porte_B      part des photons arrivés plus de 2 ns après le pic
    tau_rld_B    deux portes de 3 ns : Δt / ln(N1 / N2)
    tphi_B, tmod_B, tphi2_B   phaseur rapporté à celui de l'IRF : durée de vie de phase (harmoniques 1 et 2), de modulation
    tau_mono     une exponentielle ⊗ IRF mesurée, décalage commun fixé, fond libre (Poisson)
    tau_mono_s   la même, décalage libre pour chaque déclin (le modèle de FLIMApp)
    tau_bi, tau_bi_amp   deux exponentielles libres : moyenne pondérée par l'intensité, par l'amplitude
    tau_3pop, f_long     trois durées de vie communes (étape 5) : τ moyen, part des photons du donneur sans FRET
« quand », canal YFP
    tm_Y, tphi_Y, tau_mono_Y
gabarits — mélange des motifs mesurés à 0 et à 131 mM (ni IRF ni modèle de durée de vie) ; f = part des photons
CFP qui viennent du motif « 131 mM » ; pour un capteur à deux états, R = R0·(1 − f) + R1·f
    gab_B, gab_Y   forme d'un seul canal
    gab_formes     formes des deux canaux, amplitudes libres (sans le rapport)
    gab_joint      formes et rapport ensemble (η = 1)
    gab_formes_fond  comme gab_formes, avec un fond plat libre dans chaque canal
    phas_axe_B     projection du phaseur CFP sur l'axe 0 → 131 mM
    R_eta          R divisé par le η du champ, η = R du champ / R prédit par la forme des photons réunis du champ
apprentissage — 20 portes temporelles (10 par canal), parts de photons, racine carrée
    lda_portes     discriminant linéaire 0 contre 20 mM ;  lda_3 : le même sur (ln R, tm_B, tm_Y)
    ridge_portes   régression ridge sur c / (c + 15) ;  gbr_portes : arbres boostés ;  pc1 : 1re composante principale
modèle — celui de l'étape 9, formes seules (FLIM du système), fond libre
    mod2_flim, mod1_flim   deux populations en FRET, une seule
"""
import numpy as np

import config as C
import deux_canaux as DC
import flim

C.MIN_PHOTONS_FIT = 1          # le banc ajuste aussi des déclins de quelques dizaines de photons

FEN = (17, 254)                # canaux pleins (comme l'étape 9)
W = slice(*FEN)
NW = FEN[1] - FEN[0]
T = C.temps_canaux()[W]
OMEGA = 2 * np.pi / C.PERIODE_NS
CONDITIONS = ("0 mM", "20 mM", "131 mM")
K_CIBLE = 15.0                 # cible des régressions : c / (c + K_CIBLE)
APPRISES = ("lda_portes", "lda_3", "ridge_portes", "gbr_portes")      # lues hors de leurs champs d'entraînement


class Instrument:
    """IRF mesurée de chaque canal ; décalages (canaux) communs à toutes les acquisitions ; pic du déclin CFP."""

    def __init__(self, root, sB, sY, pic, taus3):
        self.MB = flim.Modele(C.charger_irf("cfp", root), fenetre=FEN)
        self.MY = flim.Modele(C.charger_irf("yfp", root), fenetre=FEN)
        self.sB, self.sY, self.pic = float(sB), float(sY), int(pic)
        self.taus3 = np.asarray(taus3, float)
        P3 = np.array([self.MB.forme([t], [1.0], self.sB)[W] for t in self.taus3])
        self.P3 = np.vstack([P3 / P3.sum(1, keepdims=True), np.full((1, NW), 1.0 / NW)])
        tt = C.temps_canaux()
        self.z_irf = {}
        for c, M, s in (("B", self.MB, self.sB), ("Y", self.MY, self.sY)):
            g = M.irf_decalee(s)
            self.z_irf[c] = {h: complex((g @ np.exp(1j * h * OMEGA * tt)) / g.sum()) for h in (1, 2)}


# ------------------------------------------------------------------ lectures fermées
def simples(B, Y, ins):
    Bw, Yw = B[:, W].astype(float), Y[:, W].astype(float)
    NB, NY = Bw.sum(1), Yw.sum(1)
    pk = ins.pic
    o = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        o["R"] = NY / NB
        o["R_tard"] = Yw[:, pk + 21:].sum(1) / Bw[:, pk + 21:].sum(1)
        o["tm_B"] = Bw @ T / NB - T[pk]
        o["tm_Y"] = Yw @ T / NY - T[pk]
        q = pk + 12
        o["tm_queue_B"] = Bw[:, q:] @ T[q:] / Bw[:, q:].sum(1) - T[q]
        o["porte_B"] = Bw[:, pk + 41:].sum(1) / NB
        g1, g2 = Bw[:, pk + 6:pk + 68].sum(1), Bw[:, pk + 68:pk + 130].sum(1)
        o["tau_rld_B"] = 62 * C.DT_NS / np.log(g1 / g2)
        e1, e2 = np.exp(1j * OMEGA * T), np.exp(2j * OMEGA * T)
        zB1 = (Bw @ e1) / NB / ins.z_irf["B"][1]
        zB2 = (Bw @ e2) / NB / ins.z_irf["B"][2]
        zY1 = (Yw @ e1) / NY / ins.z_irf["Y"][1]
        o["tphi_B"] = zB1.imag / zB1.real / OMEGA
        o["tmod_B"] = np.sqrt(np.maximum(1 / np.abs(zB1) ** 2 - 1, 0)) / OMEGA
        o["tphi2_B"] = zB2.imag / zB2.real / (2 * OMEGA)
        o["tphi_Y"] = zY1.imag / zY1.real / OMEGA
        o["_gB"], o["_sB"], o["_NB"], o["_NY"] = zB1.real, zB1.imag, NB, NY
    return o


# ------------------------------------------------------------------ ajustements (un déclin à la fois)
def fits(B, Y, ins):
    M = len(B)
    o = {k: np.full(M, np.nan) for k in ("tau_mono", "tau_mono_s", "tau_bi", "tau_bi_amp", "tau_mono_Y")}
    for i in range(M):
        b, y = B[i].astype(float), Y[i].astype(float)
        try:
            if b[W].sum() >= 20:
                o["tau_mono"][i] = ins.MB.ajuster_mono(b, s=ins.sB)[0]
                o["tau_mono_s"][i] = ins.MB.ajuster_mono(b)[0]
                r = ins.MB.ajuster_multi(b, 2, s=ins.sB)
                ta, am = r["taus"], r["amps"]
                o["tau_bi"][i] = (am * ta ** 2).sum() / (am * ta).sum()
                o["tau_bi_amp"][i] = (am * ta).sum()
            if y[W].sum() >= 20:
                o["tau_mono_Y"][i] = ins.MY.ajuster_mono(y, s=ins.sY)[0]
        except Exception:
            pass
    return o


def em_melange(Yw, P, iters=400):
    """Mélange de formes fixes P (K, n), normalisées : photons de chaque forme (M, K), maximum de vraisemblance de
    Poisson par EM (poids ≥ 0)."""
    a = np.repeat(Yw.sum(1, keepdims=True) / len(P), len(P), 1)
    for _ in range(iters):
        mu = np.maximum(a @ P, 1e-12)
        a *= (Yw / mu) @ P.T
    return a


def trois_pop(B, ins):
    Bw = B[:, W].astype(float)
    a = em_melange(Bw, ins.P3)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = a[:, :3].sum(1)
        return {"tau_3pop": a[:, :3] @ ins.taus3 / s, "f_long": a[:, 0] / s,
                "tm_B_sf": (Bw @ T - a[:, 3] * T.mean()) / (Bw.sum(1) - a[:, 3]) - T[ins.pic]}


# ------------------------------------------------------------------ gabarits mesurés
def argmax_lignes(L, grille):
    """Maximum de chaque ligne, affiné par une parabole ; rend aussi les lignes dont le maximum est au bord."""
    i = np.argmax(L, 1)
    r = np.arange(len(L))
    j = np.clip(i, 1, L.shape[1] - 2)
    g, m, d = L[r, j - 1], L[r, j], L[r, j + 1]
    k = g - 2 * m + d
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = np.where(k < 0, np.clip(0.5 * (g - d) / k, -0.5, 0.5), 0.0)
    x = np.where(i == j, i + delta, i).astype(float)
    return np.interp(x, np.arange(len(grille)), grille), (i == 0) | (i == L.shape[1] - 1)


def _bornes(p0, d):
    """Intervalle de f où p0 + f·d reste positif partout."""
    return np.max(-p0[d > 0] / d[d > 0]), np.min(p0[d < 0] / -d[d < 0])


class Gabarits:
    """Motifs à 0 mM (indice 0) et à 131 mM (indice 1) : moyenne des champs fournis, chaque champ normalisé sur la
    somme de ses deux canaux. Un déclin est lu comme (1 − f)·motif₀ + f·motif₁ dans le canal CFP ; dans le YFP la part
    correspondante est f·ρ / (1 − f·(1 − ρ)), ρ = R1 / R0 (capteur à deux états)."""

    def __init__(self, B0, Y0, B1, Y1, n_f=401):
        def motif(Bs, Ys):
            j = [b[W].sum() + y[W].sum() for b, y in zip(Bs, Ys)]
            return (np.mean([b[W] / t for b, t in zip(Bs, j)], 0), np.mean([y[W] / t for y, t in zip(Ys, j)], 0))
        P0B, P0Y = motif(B0, Y0)
        P1B, P1Y = motif(B1, Y1)
        self.R0, self.R1 = float(P0Y.sum() / P0B.sum()), float(P1Y.sum() / P1B.sum())
        self.rho = self.R1 / self.R0
        self.p0B, self.p1B = P0B / P0B.sum(), P1B / P1B.sum()
        self.p0Y, self.p1Y = P0Y / P0Y.sum(), P1Y / P1Y.sum()
        dB, dY = self.p1B - self.p0B, self.p1Y - self.p0Y
        loB, hiB = _bornes(self.p0B, dB)
        loY, hiY = _bornes(self.p0Y, dY)
        vers_fB = lambda fy: fy / (self.rho + fy * (1 - self.rho))
        lo = max(loB, vers_fB(loY) if loY > -self.rho / (1 - self.rho) else -np.inf, -1.0)
        hi = min(hiB, vers_fB(hiY), 1 / (1 - self.rho), 3.0)
        m = 0.02 * (hi - lo)
        self.f = np.linspace(lo + m, hi - m, n_f)
        fY = self.f * self.rho / (1 - self.f * (1 - self.rho))
        self.lB = np.log(self.p0B[None, :] + self.f[:, None] * dB[None, :])
        self.lY = np.log(self.p0Y[None, :] + fY[:, None] * dY[None, :])
        R = self.R0 * (1 - self.f) + self.R1 * self.f
        self.lrB, self.lrY = -np.log1p(R), np.log(R) - np.log1p(R)
        pas = (n_f - 1) // 80
        self.f_g = self.f[::pas]
        self.PB_g = self.p0B[None, :] + self.f_g[:, None] * dB[None, :]
        self.PY_g = self.p0Y[None, :] + fY[::pas][:, None] * dY[None, :]

    def lire(self, B, Y, lot=4000):
        noms = ("gab_B", "gab_Y", "gab_formes", "gab_joint")
        o = {k: np.empty(len(B)) for k in noms}
        o.update({"_bord_" + k: np.zeros(len(B), bool) for k in noms})
        for i in range(0, len(B), lot):
            Bw, Yw = B[i:i + lot, W].astype(float), Y[i:i + lot, W].astype(float)
            LB, LY = Bw @ self.lB.T, Yw @ self.lY.T
            J = LB + LY + Bw.sum(1)[:, None] * self.lrB[None, :] + Yw.sum(1)[:, None] * self.lrY[None, :]
            for k, L in zip(noms, (LB, LY, LB + LY, J)):
                o[k][i:i + lot], o["_bord_" + k][i:i + lot] = argmax_lignes(L, self.f)
        return o

    def lire_fond(self, B, Y):
        """Formes des deux canaux avec un fond plat libre dans chacun (Newton sur amplitude et fond, pour chaque f)."""
        out = np.full(len(B), np.nan)
        for i in range(len(B)):
            b, y = B[i, W].astype(float), Y[i, W].astype(float)
            if b.sum() < 20 or y.sum() < 20:
                continue
            ll = DC.ajuster_fond(b, self.PB_g)[0] + DC.ajuster_fond(y, self.PY_g)[0]
            out[i] = DC.argmax_affine(ll, self.f_g)
        return out


# ------------------------------------------------------------------ apprentissage
def traits(B, Y, pic):
    """20 portes temporelles (10 par canal, de plus en plus larges après le pic) : racine de la part des photons."""
    e = [0, pic - 2, pic + 2, pic + 6, pic + 12, pic + 21, pic + 35, pic + 55, pic + 85, pic + 130]
    g = np.concatenate([np.add.reduceat(Z[:, W].astype(float), e, axis=1) for Z in (B, Y)], 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(g / g.sum(1, keepdims=True))


def traits3(o):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.column_stack([np.log(o["R"]), o["tm_B"], o["tm_Y"]])


def _lda(F0, F1):
    F0, F1 = F0[np.isfinite(F0).all(1)], F1[np.isfinite(F1).all(1)]
    mu0, mu1 = F0.mean(0), F1.mean(0)
    Sw = (np.cov(F0.T) * (len(F0) - 1) + np.cov(F1.T) * (len(F1) - 1)) / (len(F0) + len(F1) - 2)
    Sw += 1e-6 * np.trace(Sw) / len(Sw) * np.eye(len(Sw))
    return np.linalg.solve(Sw, mu1 - mu0), (mu0 + mu1) / 2


class Appris:
    """Lectures apprises sur des somas d'entraînement (déclins réunis de chaque soma, [Cl⁻] imposée)."""

    def __init__(self, B, Y, cond, acq, ins):
        from sklearn.ensemble import HistGradientBoostingRegressor
        self.pic = ins.pic
        o = simples(B, Y, ins)
        F, F3 = traits(B, Y, self.pic), traits3(o)
        m = {c: (cond == c) for c in CONDITIONS}
        self.lda = _lda(F[m["0 mM"]], F[m["20 mM"]])
        self.lda3 = _lda(F3[m["0 mM"]], F3[m["20 mM"]])
        cl = np.select([m[c] for c in CONDITIONS], [0.0, 20.0, 131.0], np.nan)
        ok = np.isfinite(cl) & np.isfinite(F).all(1)
        t = cl[ok] / (cl[ok] + K_CIBLE)
        self.mu, self.sd, self.t0 = F[ok].mean(0), F[ok].std(0), t.mean()
        Z = (F[ok] - self.mu) / self.sd
        self.w = np.linalg.solve(Z.T @ Z + 1.0 * np.eye(Z.shape[1]), Z.T @ (t - self.t0))
        self.gbr = HistGradientBoostingRegressor(max_iter=200, max_depth=3, learning_rate=0.05,
                                                 random_state=0).fit(F[ok], t)
        self.pc_mu = F[ok].mean(0)
        self.pc = np.linalg.svd(F[ok] - self.pc_mu, full_matrices=False)[2][0]
        if ((F[ok] - self.pc_mu) @ self.pc)[cl[ok] == 131].mean() < ((F[ok] - self.pc_mu) @ self.pc)[cl[ok] == 0].mean():
            self.pc = -self.pc
        # axe du phaseur CFP : moyenne des champs (photons réunis de chaque champ) à 0 et à 131 mM
        z = o["_gB"] + 1j * o["_sB"]
        ref = {}
        for c in ("0 mM", "131 mM"):
            zs = [np.sum(o["_NB"][m[c] & (acq == a)] * z[m[c] & (acq == a)]) / np.sum(o["_NB"][m[c] & (acq == a)])
                  for a in np.unique(acq[m[c]])]
            ref[c] = np.mean(zs)
        self.z0, self.z1 = ref["0 mM"], ref["131 mM"]

    def lire(self, B, Y, o):
        F = traits(B, Y, self.pic)
        d = self.z1 - self.z0
        z = o["_gB"] + 1j * o["_sB"]
        return {"lda_portes": (F - self.lda[1]) @ self.lda[0],
                "lda_3": (traits3(o) - self.lda3[1]) @ self.lda3[0],
                "ridge_portes": self.t0 + ((F - self.mu) / self.sd) @ self.w,
                "gbr_portes": self.gbr.predict(F),
                "pc1": (F - self.pc_mu) @ self.pc,
                "phas_axe_B": ((z - self.z0) * np.conj(d)).real / abs(d) ** 2}


# ------------------------------------------------------------------ modèle de l'étape 9
def lecteur_modele(root, Bc, Yc, conc, irfp, nt):
    """Calibration conjointe (étape 9) sur les déclins réunis de chaque [Cl⁻], lecteur sur une grille de c réduite."""
    d = DC.donnees(root, Bc, Yc, conc, FEN)
    d["irfp"] = irfp
    R = DC.calibrer(d, irfp, nt)
    iso = DC.isothermes(R, conc)
    S = DC.senseur(d, R, iso)
    cg = np.r_[np.linspace(DC.c_min_modele(iso), 0, 40, endpoint=False), 0.0, np.geomspace(0.05, 1e5, 160)]
    return DC.LecteurSoma(S, cg, np.array([1.0, 2.0])), d["win"][1]


def lire_modele(L, w, B, Y):
    out = np.full(len(B), np.nan)
    haut = L.c[-1] * (1 - 1e-9)
    for i in range(len(B)):
        b, y = B[i].astype(float)[w], Y[i].astype(float)[w]
        if b.sum() < 20 or y.sum() < 20:
            continue
        c = DC.argmax_affine(L.lire(b, y)["formes"], L.c)
        out[i] = np.inf if c >= haut else c
    return out


# ------------------------------------------------------------------ lectures lentes, en parallèle
_CTX = {}


def init_lent(ins, gab, lecteurs):
    _CTX.update(ins=ins, gab=gab, lecteurs=lecteurs)


def lent(bloc):
    B, Y = bloc
    o = fits(B, Y, _CTX["ins"])
    o["gab_formes_fond"] = _CTX["gab"].lire_fond(B, Y)
    for nom, (L, w) in _CTX["lecteurs"].items():
        o[nom] = lire_modele(L, w, B, Y)
    return o
