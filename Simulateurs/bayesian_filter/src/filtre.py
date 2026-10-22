"""
filtre.py — le filtre bayésien hybride (sections 9 à 14 du document) et les lectures de référence.

Vraisemblance d'une trame (section 11), A profilé, découpage exact « combien » / « quand » :
    ln L(c, η) = N_B ln[M_B/(M_B + η M_Y)] + N_Y ln[η M_Y/(M_B + η M_Y)]      ratio : dépend de c et de η
               + Σ_i n_B,i ln p_B,i(c) + Σ_i n_Y,i ln p_Y,i(c)                   FLIM  : dépend de c seulement
Filtre dans le temps (section 13) : prédiction par le noyau K (marche aléatoire sur ln c avec sauts rares), puis
mise à jour bayésienne ; lissage avant-arrière hors ligne. Version temps réel (section 14) : η statique partagé par
M ROI, p(c_1..c_M, η) = p(η) Π_m p(c_m | η), portes temporelles, bascule vers un filtre sur le ratio seul.

Fond (absent du document, présent dans les données : ~2,3 % du bleu) : facultatif, en photons par case et par trame,
supposé connu. Il est soustrait des totaux (« fond soustrait », section 7) et ajouté aux formes :
p_i = (Ŝ p_i(c) + b) / N, avec Ŝ = N − n b le signal estimé de la trame.
"""
import numpy as np


def grille_c(f):
    return np.geomspace(f['c_min_mM'], f['c_max_mM'], f['n_c'])


def grille_eta(f):
    return np.geomspace(f['eta_min'], f['eta_max'], f['n_eta'])


def portes_log(n_cases, n_portes=16):
    """Matrice (n_portes, n_cases) de portes contiguës, plus fines au début du déclin (section 14)."""
    bords = np.round(np.geomspace(1, n_cases + 1, n_portes + 1) - 1).astype(int)
    bords[0], bords[-1] = 0, n_cases
    for i in range(1, n_portes):                      # au moins une case par porte
        bords[i] = max(bords[i], bords[i - 1] + 1)
    G = np.zeros((n_portes, n_cases))
    for g in range(len(bords) - 1):
        G[g, bords[g]:bords[g + 1]] = 1
    return G


class Grille:
    """Formes p_i(c) et totaux M(c) du modèle calibré sur la grille de c (A = 1, η = 1), en cases ou en portes."""
    def __init__(self, S, c, eta, portes=None):
        self.c, self.eta, self.G = np.asarray(c, float), np.asarray(eta, float), portes
        mB, mY = S.mu(self.c)
        self.largeur = np.ones(mB.shape[1])
        if portes is not None:
            mB, mY, self.largeur = mB @ portes.T, mY @ portes.T, portes.sum(1)
        self.MB, self.MY = mB.sum(1), mY.sum(1)
        self.pB, self.pY = mB / self.MB[:, None], mY / self.MY[:, None]
        self.lpB, self.lpY = np.log(self.pB), np.log(self.pY)
        self.R1 = self.MY / self.MB
        self.R1_zero = float(S.ratio_calibration(0.0))

    def regrouper(self, h):
        """Histogramme en cases de la fenêtre → portes (sans effet si la grille est en cases)."""
        h = np.asarray(h, float)
        return h @ self.G.T if self.G is not None and h.shape[-1] == self.G.shape[1] else h

    def _sous_fond(self, B, Y, fond):
        NB, NY = B.sum(), Y.sum()
        if fond is None:
            return NB, NY
        return max(NB - fond[0] * self.largeur.sum(), 1.0), max(NY - fond[1] * self.largeur.sum(), 1.0)

    def formes(self, B, Y, fond=None, canaux=('B', 'Y')):
        """Morceau « quand » : Σ n_i ln p_i(c), vecteur sur la grille de c."""
        out = np.zeros(len(self.c))
        for nom, h, lp, p, k in (('B', B, self.lpB, self.pB, 0), ('Y', Y, self.lpY, self.pY, 1)):
            if nom not in canaux:
                continue
            if fond is None:
                out += lp @ h
            else:
                b = fond[k] * self.largeur; N = h.sum(); S_ = max(N - b.sum(), 1.0)
                out += np.log((S_ * p + b) / (S_ + b.sum())) @ h
        return out

    def ratio(self, NB, NY, eta=None):
        """Morceau « combien » : partage binomial des photons entre les couleurs, matrice (n_c, n_η)."""
        eta = self.eta if eta is None else np.atleast_1d(eta)
        mB = self.MB[:, None]; mY = self.MY[:, None] * eta[None, :]
        return NB * np.log(mB / (mB + mY)) + NY * np.log(mY / (mB + mY))

    def loglik(self, B, Y, fond=None, eta=None):
        B, Y = self.regrouper(B), self.regrouper(Y)
        NB, NY = self._sous_fond(B, Y, fond)
        return self.formes(B, Y, fond)[:, None] + self.ratio(NB, NY, eta)

    # ---------------------------------------------------------------- lectures de référence
    def lire_ratio(self, B, Y, eta=1.0, fond=None):
        """Lecture classique : R observé / η inversé sur R₁(c) ; R > R₁(0) donne 0 mM (section 15)."""
        B, Y = self.regrouper(B), self.regrouper(Y)
        NB, NY = self._sous_fond(B, Y, fond)
        R = NY / NB / eta
        cs = np.concatenate([[0.0], self.c]); Rs = np.concatenate([[self.R1_zero], self.R1])
        return float(np.interp(R, Rs[::-1], cs[::-1]))

    def lire_flim(self, B, Y, fond=None, canaux=('B', 'Y')):
        """Durée de vie seule : maximum de vraisemblance des formes, trame par trame."""
        B, Y = self.regrouper(B), self.regrouper(Y)
        return float(self.c[np.argmax(self.formes(B, Y, fond, canaux))])


def transition(logx, s, s_saut=None, w=0.0):
    """K[i, j] = P(x_i | x_j) : marche aléatoire gaussienne sur ln x, plus des sauts rares (section 13)."""
    dd = logx[:, None] - logx[None, :]
    K = np.exp(-0.5 * (dd / s) ** 2); K /= K.sum(0, keepdims=True)
    if s_saut:
        K2 = np.exp(-0.5 * (dd / s_saut) ** 2); K2 /= K2.sum(0, keepdims=True)
        K = (1 - w) * K + w * K2
    return K


def resume_c(pc, c):
    """Moyenne et intervalle crédible à 95 % de distributions sur la grille de c (dernière dimension)."""
    pc = np.atleast_2d(pc); cdf = np.cumsum(pc, -1)
    return (pc * c).sum(-1), c[np.argmax(cdf >= 0.025, -1)], c[np.argmax(cdf >= 0.975, -1)]


def resume_eta(pe, eta):
    """η moyen (géométrique) et sd(ln η), qui vaut l'incertitude relative sur η."""
    le = np.log(eta); m = (pe * le).sum()
    return float(np.exp(m)), float(np.sqrt((pe * (le - m) ** 2).sum()))


def filtre_lisseur(Ls, Kc, Ke):
    """
    Filtre causal puis lissage avant-arrière sur la grille (c, η) (sections 12 et 13).
    Ls : (T, n_c, n_η) log-vraisemblances par trame. Renvoie (filtre, lisse), chacun (T, n_c, n_η), normalisés.
    """
    T, nc, ne = Ls.shape
    p = np.full((nc, ne), 1.0 / (nc * ne))
    alpha = np.empty_like(Ls); pred = np.empty_like(Ls)
    for k in range(T):
        pred[k] = Kc @ p @ Ke.T if k > 0 else p
        p = pred[k] * np.exp(Ls[k] - Ls[k].max()); p /= p.sum(); alpha[k] = p
    post = np.empty_like(alpha); post[-1] = alpha[-1]
    for k in range(T - 2, -1, -1):
        post[k] = alpha[k] * (Kc.T @ (post[k + 1] / np.maximum(pred[k + 1], 1e-300)) @ Ke)
        post[k] /= post[k].sum()
    return alpha, post


class FiltreTempsReel:
    """
    Filtre causal exact pour un η statique partagé par M ROI (même profondeur) (section 14).
    Chaque ROI garde p(c | η, trames 1..k) (colonnes normalisées) ; le groupe garde ln p(η).
    ℓ_m,k(η) = Σ_c L_m,k(c, η) p(c | η, 1..k−1) ;   p(η | 1..k) ∝ p(η | 1..k−1) Π_m ℓ_m,k(η).
    Quand sd(ln η) < bascule, η est figé et chaque ROI passe à un filtre sur c seul, nourri par N_B et N_Y.
    eta_fixe : η connu d'avance (référence « η connu ») ; la grille de η est alors réduite à cette valeur.
    """
    def __init__(self, grille, n_roi, Kc, bascule=None, eta_fixe=None):
        self.g, self.Kc, self.bascule = grille, Kc, bascule
        self.eta = grille.eta if eta_fixe is None else np.array([float(eta_fixe)])
        self.P = np.full((n_roi, len(grille.c), len(self.eta)), 1.0 / len(grille.c))
        self.log_pe = np.zeros(len(self.eta))
        self.eta_fige = float(eta_fixe) if eta_fixe is not None else None
        self.pc = self.P[:, :, 0].copy() if eta_fixe is not None else None
        self.k = 0; self.trame_bascule = 0 if eta_fixe is not None else None

    def pe(self):
        p = np.exp(self.log_pe - self.log_pe.max()); return p / p.sum()

    def mettre_a_jour(self, B, Y, fond=None):
        """B, Y : (M, n_cases) ou (M, n_portes) pour la trame. Renvoie un dict de résumés par ROI."""
        g = self.g
        B, Y = g.regrouper(np.atleast_2d(B)), g.regrouper(np.atleast_2d(Y))
        for m in range(B.shape[0]):
            NB, NY = g._sous_fond(B[m], Y[m], fond)
            if self.eta_fige is None:
                if self.k > 0:
                    self.P[m] = self.Kc @ self.P[m]
                L = g.formes(B[m], Y[m], fond)[:, None] + g.ratio(NB, NY, self.eta)
                Q = self.P[m] * np.exp(L - L.max())
                ell = np.maximum(Q.sum(0), 1e-300)
                self.P[m] = Q / ell; self.log_pe += np.log(ell)
            else:
                if self.k > 0:
                    self.pc[m] = self.Kc @ self.pc[m]
                L = g.ratio(NB, NY, self.eta_fige)[:, 0]
                q = self.pc[m] * np.exp(L - L.max()); self.pc[m] = q / q.sum()
        self.log_pe -= self.log_pe.max()
        pe = self.pe(); eta_m, eta_sd = resume_eta(pe, self.eta)
        if self.eta_fige is None:
            pc = np.einsum('mce,e->mc', self.P, pe)
            if self.bascule is not None and eta_sd < self.bascule:
                self.eta_fige, self.pc, self.trame_bascule = float(eta_m), pc.copy(), self.k
        else:
            pc = self.pc
        self.k += 1
        m_, lo, hi = resume_c(pc, g.c)
        return dict(c=m_, bas=lo, haut=hi, eta=float(eta_m), eta_sd=float(eta_sd))
