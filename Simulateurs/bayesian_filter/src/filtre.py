"""
filtre.py — le filtre bayésien hybride (sections 9 à 14 du document) et les lectures de référence.

Vraisemblance d'une trame (section 11), A profilé, découpage exact « combien » / « quand » :
    ln L(c, η) = N_B ln[M_B/(M_B + η M_Y)] + N_Y ln[η M_Y/(M_B + η M_Y)]      ratio : dépend de c et de η
               + Σ_i n_B,i ln p_B,i(c) + Σ_i n_Y,i ln p_Y,i(c)                   FLIM  : dépend de c seulement
Filtre dans le temps (section 13) : prédiction par le noyau K (marche aléatoire sur ln c avec sauts rares), puis
mise à jour bayésienne ; lissage avant-arrière hors ligne. Version temps réel (section 14) : η statique partagé par
M ROI, p(c_1..c_M, η) = p(η) Π_m p(c_m | η), portes temporelles. η n'est jamais figé : la bascule de la section 14
(filtre léger sur le ratio seul une fois sd(ln η) < 1 %) a été retirée, elle figeait η sur un biais.

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

    def logs_ratio(self):
        """(ln M_B/(M_B + η M_Y), ln η M_Y/(M_B + η M_Y)) sur la grille (n_c, n_η) : ratio(NB, NY) = NB·[0] + NY·[1]."""
        mB = self.MB[:, None]; mY = self.MY[:, None] * self.eta[None, :]
        return np.log(mB / (mB + mY)), np.log(mY / (mB + mY))

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
    """
    K[i, j] = P(x_i | x_j) : marche aléatoire de pas s (sd par trame) sur ln x, plus des sauts rares (section 13).

    Sur une grille de pas Δ (uniforme en ln x), une gaussienne échantillonnée ne représente la marche que si s ≳ Δ : pour
    s ≪ Δ, le poids vers le voisin, exp(−Δ²/2s²), s'effondre (~ 1e-19 pour s = 0,2 %, Δ = 1,83 %) et x ne bouge plus.
    Pour s < Δ, on prend donc la discrétisation standard de la diffusion : vers chaque voisin avec la probabilité
    s²/2Δ², sur place sinon ; la variance ajoutée par trame vaut alors exactement s² (bords réfléchissants).
    s = 0 : x constant (identité).
    """
    n = len(logx)
    if s == 0:
        K = np.eye(n)
    elif n > 1 and s < (pas := float(np.mean(np.diff(logx)))):
        q = s ** 2 / (2 * pas ** 2)
        K = np.diag(np.full(n, 1 - 2 * q)) + np.diag(np.full(n - 1, q), 1) + np.diag(np.full(n - 1, q), -1)
        K[0, 0] += q; K[-1, -1] += q                                  # bords : la part qui sortirait reste sur place
    else:
        dd = logx[:, None] - logx[None, :]
        K = np.exp(-0.5 * (dd / s) ** 2); K /= K.sum(0, keepdims=True)
    if s_saut:
        dd = logx[:, None] - logx[None, :]
        K2 = np.exp(-0.5 * (dd / s_saut) ** 2); K2 /= K2.sum(0, keepdims=True)
        K = (1 - w) * K + w * K2
    return K


def resume_c(pc, c):
    """Moyenne et intervalle crédible à 95 % de distributions sur la grille de c (dernière dimension)."""
    pc = np.atleast_2d(pc); cdf = np.cumsum(pc, -1)
    return (pc * c).sum(-1), c[np.argmax(cdf >= 0.025, -1)], c[np.argmax(cdf >= 0.975, -1)]


def resume_eta(pe, eta):
    """η moyen (géométrique) et sd(ln η), qui vaut l'incertitude relative sur η ; pe : (n_η,) ou (T, n_η)."""
    le = np.log(eta); m = (pe * le).sum(-1)
    sd = np.sqrt((pe * (le - np.expand_dims(m, -1)) ** 2).sum(-1))
    return (float(np.exp(m)), float(sd)) if np.ndim(pe) == 1 else (np.exp(m), sd)


def filtre_lisseur(loglik, T, Kc, Ke, lisser=True):
    """
    Filtre causal puis lissage avant-arrière sur la grille (c, η) (sections 12 et 13), η pouvant évoluer selon Ke.
    loglik(k) : ln L de la trame k, (n_c, n_η). Renvoie un dict de marginales par trame — filtre_c (T, n_c),
    filtre_eta (T, n_η), lisse_c, lisse_eta — et ln_evidence (T,) : ln p(trame k | trames 1..k−1), dont la somme
    cumulée est la vraisemblance marginale des k premières trames sous ce modèle de dynamique.
    Mémoire : seul le filtre est gardé (float32) ; les prédictions sont recalculées pendant la passe arrière.
    lisser=False : filtre causal seul (pas de passe arrière, rien d'autre que les marginales gardé).
    """
    p = None; alpha = None; ln_ev = np.empty(T); fc = fe = None
    for k in range(T):
        L = loglik(k)
        if fc is None:
            fc, fe = np.empty((T, L.shape[0])), np.empty((T, L.shape[1]))
            alpha = np.empty((T,) + L.shape, np.float32) if lisser else None
            pred = np.full(L.shape, 1.0 / L.size)
        else:
            pred = Kc @ p @ Ke.T
        m = L.max(); q = pred * np.exp(L - m); z = q.sum()
        p = q / z; ln_ev[k] = np.log(z) + m
        fc[k], fe[k] = p.sum(1), p.sum(0)
        if lisser:
            alpha[k] = p
    out = dict(filtre_c=fc, filtre_eta=fe, ln_evidence=ln_ev)
    if not lisser:
        return out
    lc = np.empty((T, alpha.shape[1])); le = np.empty((T, alpha.shape[2]))
    post = alpha[-1].astype(np.float64); lc[-1], le[-1] = post.sum(1), post.sum(0)
    for k in range(T - 2, -1, -1):
        a = alpha[k].astype(np.float64)
        pred = Kc @ a @ Ke.T
        post = a * (Kc.T @ (post / np.maximum(pred, 1e-300)) @ Ke); post /= post.sum()
        lc[k], le[k] = post.sum(1), post.sum(0)
    out.update(lisse_c=lc, lisse_eta=le)
    return out


def moyenne_de_modeles(courses, causal):
    """
    Moyenne bayésienne sur les modèles de dynamique de η (a priori uniforme), à partir des sorties de filtre_lisseur.
    causal=True  : poids de la trame k ∝ p(trames 1..k | modèle) — pour le filtre temps réel ;
    causal=False : poids ∝ p(toutes les trames | modèle), les mêmes à toutes les trames — pour le lisseur.
    Renvoie (marginale de c, marginale de η, poids (T, n_modèles)).
    """
    cum = np.array([np.cumsum(r['ln_evidence']) for r in courses]).T          # (T, n_modèles)
    if not causal:
        cum = np.repeat(cum[-1:], len(cum), 0)
    w = np.exp(cum - cum.max(1, keepdims=True)); w /= w.sum(1, keepdims=True)
    cle = 'filtre' if causal else 'lisse'
    pc = sum(w[:, i, None] * r[f'{cle}_c'] for i, r in enumerate(courses))
    pe = sum(w[:, i, None] * r[f'{cle}_eta'] for i, r in enumerate(courses))
    return pc, pe, w


class FiltreTempsReel:
    """
    Filtre causal exact pour un η statique partagé par M ROI (même profondeur) (section 14).
    Chaque ROI garde p(c | η, trames 1..k) (colonnes normalisées) ; le groupe garde ln p(η).
    ℓ_m,k(η) = Σ_c L_m,k(c, η) p(c | η, 1..k−1) ;   p(η | 1..k) ∝ p(η | 1..k−1) Π_m ℓ_m,k(η).
    η est mis à jour à chaque trame, jusqu'à la fin : il n'est jamais figé.
    eta_fixe : η connu d'avance (référence « η connu ») ; la grille de η est réduite à cette valeur, la mise à jour
    (formes et ratio) est la même.
    """
    def __init__(self, grille, n_roi, Kc, eta_fixe=None):
        self.g, self.Kc = grille, Kc
        self.eta = grille.eta if eta_fixe is None else np.array([float(eta_fixe)])
        self.P = np.full((n_roi, len(grille.c), len(self.eta)), 1.0 / len(grille.c))
        self.log_pe = np.zeros(len(self.eta))
        self.k = 0

    def pe(self):
        p = np.exp(self.log_pe - self.log_pe.max()); return p / p.sum()

    def mettre_a_jour(self, B, Y, fond=None):
        """B, Y : (M, n_cases) ou (M, n_portes) pour la trame. Renvoie un dict de résumés par ROI."""
        g = self.g
        B, Y = g.regrouper(np.atleast_2d(B)), g.regrouper(np.atleast_2d(Y))
        for m in range(B.shape[0]):
            NB, NY = g._sous_fond(B[m], Y[m], fond)
            if self.k > 0:
                self.P[m] = self.Kc @ self.P[m]
            L = g.formes(B[m], Y[m], fond)[:, None] + g.ratio(NB, NY, self.eta)
            Q = self.P[m] * np.exp(L - L.max())
            ell = np.maximum(Q.sum(0), 1e-300)
            self.P[m] = Q / ell; self.log_pe += np.log(ell)
        self.log_pe -= self.log_pe.max()
        pe = self.pe(); eta_m, eta_sd = resume_eta(pe, self.eta)
        pc = np.einsum('mce,e->mc', self.P, pe)
        self.k += 1
        m_, lo, hi = resume_c(pc, g.c)
        return dict(c=m_, bas=lo, haut=hi, eta=float(eta_m), eta_sd=float(eta_sd))
