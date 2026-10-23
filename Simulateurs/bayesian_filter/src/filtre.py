"""
filtre.py — le filtre bayésien hybride (sections 9 à 14 du document) et les lectures de référence.

Vraisemblance d'une trame (section 11), A profilé, découpage exact « combien » / « quand » :
    ln L(c, η) = N_B ln[M_B/(M_B + η M_Y)] + N_Y ln[η M_Y/(M_B + η M_Y)]      ratio : dépend de c et de η
               + Σ_i n_B,i ln p_B,i(c) + Σ_i n_Y,i ln p_Y,i(c)                   FLIM  : dépend de c seulement
Filtre dans le temps (section 13) : prédiction par le noyau K (marche aléatoire sur ln c avec sauts rares), puis
mise à jour bayésienne ; lissage avant-arrière hors ligne. η n'est jamais figé : la bascule de la section 14
(filtre léger sur le ratio seul une fois sd(ln η) < 1 %) a été retirée, elle figeait η sur un biais.

Modèle par défaut, « 1-2-3 » (FiltreTempsReel ; essais et figures : derive_eta/README.md) :
  1. η dérive : ln η_k = ln η0 + β·k/1000, θ = (η0, β) statique et inconnu, sur une grille ;
  2. θ partagé par les M ROI d'un même champ : p(θ, c_1..c_M) = p(θ) Π_m p(c_m | θ), exact (section 14) ;
  3. vitesse de c inconnue : un filtre par vitesse s_c candidate, moyennés selon leur vraisemblance marginale.
L'ancien modèle d'une ROI (η en marche aléatoire de pas s_η inconnu, s fixe : filtre_lisseur, moyenne_de_modeles)
reste disponible ([filtre] dynamique_eta = "marche").

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


def grille_theta(f, derive=True):
    """
    (1) θ = (ln η0, β) de la dérive ln η_k = ln η0 + β·k/1000 : η0 de eta0_min à eta0_max (n_eta0 valeurs, pas uniforme
    en ln η), β de −derive_max à +derive_max (n_derive valeurs), en ln η par 1 000 trames. derive=False : β = 0.
    Renvoie (ln η0, β), à plat ; a priori uniforme sur la grille.
    """
    le0 = np.log(np.geomspace(f['eta0_min'], f['eta0_max'], f['n_eta0']))
    beta = np.linspace(-f['derive_max'], f['derive_max'], f['n_derive']) if derive else np.zeros(1)
    L0, Bt = np.meshgrid(le0, beta, indexing='ij')
    return L0.ravel(), Bt.ravel()


def noyaux_c(c, f):
    """(3) Un noyau de transition de c par vitesse candidate ([filtre] s_c_candidats), sauts du document compris."""
    return [transition(np.log(c), s, f['s_saut'], f['w']) for s in f['s_c_candidats']]


def logs_ratio_eta(g, ln_eta):
    """Le partage pour des valeurs quelconques de ln η : (ln M_B/(M_B + η M_Y), ln η M_Y/(M_B + η M_Y)), (n_c, n)."""
    r = (g.MY / g.MB)[:, None] * np.exp(np.asarray(ln_eta, float))[None, :]
    lB = -np.log1p(r)
    return lB, np.log(r) + lB


# En float32, les valeurs sous ~1e-38 sont « dénormalisées » et les calculer est très lent (produits matriciels ×4) :
# les probabilités, les noyaux et les rapports de vraisemblance sous ce plancher sont mis à 0 (ils pèsent moins de 1e-18
# relativement au plus grand de leur colonne ; un produit de deux valeurs gardées reste au-dessus de 1e-36).
PLANCHER = 1e-18


def _lse(x, axis=-1):
    m = x.max(axis, keepdims=True)
    return np.squeeze(m, axis) + np.log(np.exp(x - m).sum(axis))


class NoyauStructure:
    """
    (6) K @ P pour un noyau de c « bande + rang faible » : la bande |i − j| ≤ largeur (la marche) en matrice creuse, le
    reste (les sauts, lisses) par une SVD tronquée à `tol` de sa plus grande valeur singulière. Approché à ~tol près.
    """
    def __init__(self, K, largeur=15, tol=1e-7):
        import scipy.sparse as sp
        n = K.shape[0]; i, j = np.indices(K.shape); bande = np.where(np.abs(i - j) <= largeur, K, 0.0)
        U, s, Vt = np.linalg.svd(K - bande); r = max(1, int((s > tol * s[0]).sum()))
        self.B = sp.csr_matrix(bande.astype(np.float32)); self.U = (U[:, :r] * s[:r]).astype(np.float32)
        self.V = Vt[:r].astype(np.float32); self.rang = r

    def __matmul__(self, P):
        out = np.matmul(self.U, np.matmul(self.V, P))
        for m in range(P.shape[0]):
            out[m] += self.B @ P[m]
        return out


try:                                                       # (5) mise à jour fusionnée, si numba est installé
    from numba import njit, prange

    @njit(parallel=True, cache=False)
    def _maj_numba(P, E, plancher):
        """P ← P·E normalisé par colonne (sur c), valeurs sous `plancher` mises à 0 ; renvoie les sommes ℓ (M, n_θ)."""
        M, nc, nt = P.shape; ell = np.empty((M, nt))
        for m in prange(M):
            acc = np.zeros(nt)
            for c in range(nc):
                for t in range(nt):
                    q = P[m, c, t] * E[m, c, t]; P[m, c, t] = q; acc[t] += q
            for t in range(nt):
                ell[m, t] = max(acc[t], 1e-300)
            for c in range(nc):
                for t in range(nt):
                    v = P[m, c, t] / max(ell[m, t], 1e-30)
                    P[m, c, t] = v if v >= plancher else 0.0
        return ell
    NUMBA = True
except ImportError:                                        # pragma: no cover
    NUMBA = False


def _maj_numpy(P, E, plancher):
    Q = P * E; ell = np.maximum(Q.sum(1, dtype=np.float64), 1e-300)
    Q /= np.maximum(ell, 1e-30)[:, None, :].astype(np.float32); Q[Q < plancher] = 0; P[...] = Q
    return ell


class FiltreTempsReel:
    """
    Filtre causal exact pour M ROI d'un même champ ; par défaut le modèle « 1-2-3 » :
      1. η dérive : ln η_k = ln η0 + β·k/1000 (k : trames depuis le début), θ = (η0, β) statique, grille de grille_theta ;
      2. θ partagé par les M ROI (section 14) : chaque ROI garde p(c | θ, trames 1..k) (colonnes normalisées) et
         ℓ_m,k(θ) = Σ_c L_m,k(c, θ) p(c | θ, 1..k−1) ;  p(θ | 1..k) ∝ p(θ | 1..k−1) Π_m ℓ_m,k(θ) — exact ;
      3. un filtre par noyau de c (vitesses [filtre] s_c_candidats), moyennés selon leur vraisemblance marginale causale
         p(trames 1..k | noyau) (a priori uniforme sur les noyaux).
    θ est mis à jour à chaque trame, jusqu'à la fin : il n'est jamais figé.
    Options : derive=False → β = 0 (η constant) ; eta_fixe → η connu d'avance (θ réduit à ce point, β = 0) ;
    noyaux → liste de noyaux de c imposée ; partage=False → chaque ROI a son propre θ et sa propre moyenne de modèles
    (M filtres indépendants, calculés ensemble).
    Accélérations ([filtre] ou arguments ; README, « Accélérations ») :
      elaguer_theta    θ dont la probabilité relative tombe sous seuil_theta (pour tous les noyaux et groupes) retirés ;
      elaguer_noyaux   noyaux de c dont le poids tombe sous seuil_noyaux retirés ;
      numba            mise à jour fusionnée (une passe au lieu de cinq) ;
      noyau_structure  K @ P par bande creuse + rang faible ;
      adaptatif        grille de θ recentrée à ±6 sd (θ partagé seulement) ;
      trames           durée prévue : la grille de β est réduite à n_derive·trames/3000 valeurs (au moins 5).
    """
    def __init__(self, grille, n_roi, f, eta_fixe=None, derive=True, noyaux=None, partage=True, trames=None, **acc):
        self.g, self.M, self.partage, self.f = grille, n_roi, partage, f
        o = lambda k, d: acc[k] if k in acc else f.get(k, d)
        self.elaguer_theta, self.seuil_theta = o('elaguer_theta', False), o('seuil_theta', PLANCHER)
        self.elaguer_noyaux, self.seuil_noyaux = o('elaguer_noyaux', False), o('seuil_noyaux', 1e-12)
        self.numba = o('numba', False) and NUMBA
        self.adaptatif = o('adaptatif', False) and partage and eta_fixe is None and derive
        Ks = noyaux if noyaux is not None else noyaux_c(grille.c, f)
        self.K = [np.where(np.asarray(K) < PLANCHER, 0, K).astype(np.float32) for K in Ks]
        self.op = [NoyauStructure(K) for K in self.K] if o('noyau_structure', False) else list(self.K)
        if eta_fixe is not None:
            self.le0, self.beta = np.array([np.log(float(eta_fixe))]), np.zeros(1)
        else:
            fd = dict(f)
            if trames and derive:                           # (8) moins de dérives pour une course courte
                fd['n_derive'] = int(max(5, min(f['n_derive'], 2 * np.ceil(f['n_derive'] * trames / 6000) + 1)))
            self.le0, self.beta = grille_theta(fd, derive)
            self.a_vals = np.unique(self.le0); self.b_vals = np.unique(self.beta); self.tau = 0.0
        nk, nc, nt, ng = len(self.K), len(grille.c), len(self.le0), (1 if partage else n_roi); self.n_noyaux = nk
        self.P = np.full((nk, n_roi, nc, nt), 1.0 / nc, np.float32)
        self.lpe = np.full((nk, ng, nt), -np.log(nt))           # ln p(θ | trames) par noyau et par groupe
        self.cum = np.zeros((nk, ng))                           # ln p(trames 1..k | noyau)
        self.k = 0; self.n_regrilles = 0; self.ids = np.arange(nk)        # ids : noyaux d'origine encore gardés

    def ln_eta(self, k):
        """ln η à la trame k (0 = la première) en chaque point de θ."""
        return self.le0 + self.beta * k / 1000.0

    def poids(self):
        """Poids des noyaux de c gardés, (n_noyaux, groupes)."""
        w = np.exp(self.cum - self.cum.max(0)); return w / w.sum(0)

    def poids_complets(self):
        """Poids de tous les noyaux d'origine (0 pour ceux qu'on a retirés)."""
        w = np.zeros((self.n_noyaux, self.cum.shape[1])); w[self.ids] = self.poids()
        return w

    def pe(self):
        """Loi de θ, moyennée sur les noyaux : (n_θ,) si θ est partagé, (M, n_θ) sinon."""
        p = np.einsum('kg,kgt->gt', self.poids(), np.exp(self.lpe))
        return p[0] if self.partage else p

    def mettre_a_jour(self, B, Y, fond=None):
        """B, Y : (M, n_cases) ou (M, n_portes) pour la trame. Renvoie un dict de résumés par ROI."""
        g = self.g
        if self.k >= 10 and self.k % 10 == 0:
            self._entretien()
        B, Y = g.regrouper(np.atleast_2d(B)), g.regrouper(np.atleast_2d(Y))
        le = self.ln_eta(self.k); lB, lY = logs_ratio_eta(g, le)
        E = np.empty(self.P.shape[1:], np.float32); mx = np.empty(self.M)
        for m in range(self.M):
            NB, NY = g._sous_fond(B[m], Y[m], fond)
            L = g.formes(B[m], Y[m], fond)[:, None] + NB * lB + NY * lY          # float64 : |ln L| ~ 1e5 à 1e5 photons
            mx[m] = L.max(); L -= mx[m]; L[L < np.log(PLANCHER)] = -np.inf; E[m] = np.exp(L)
        for i, K in enumerate(self.op):
            P = (K @ self.P[i]) if self.k > 0 else self.P[i].copy()
            ell = _maj_numba(P, E, np.float32(PLANCHER)) if self.numba else _maj_numpy(P, E, PLANCHER)
            self.P[i] = P                                                         # colonne vide : reste à 0
            if self.partage:
                x = self.lpe[i, 0] + np.log(ell).sum(0); z = _lse(x)
                self.lpe[i, 0] = x - z; self.cum[i, 0] += z + mx.sum()
            else:
                x = self.lpe[i] + np.log(ell); z = _lse(x)
                self.lpe[i] = x - z[:, None]; self.cum[i] += z + mx
        self.k += 1
        return self.resume(le)

    # ------------------------------------------------------------------------------------------------- accélérations
    def _entretien(self):
        if self.elaguer_noyaux and len(self.K) > 1:
            garde = np.flatnonzero((self.poids() >= self.seuil_noyaux).any(1))
            if len(garde) < len(self.K):
                self.K = [self.K[i] for i in garde]; self.op = [self.op[i] for i in garde]; self.ids = self.ids[garde]
                self.P, self.lpe, self.cum = self.P[garde], self.lpe[garde], self.cum[garde]
        if self.adaptatif:
            self._regriller_si_besoin()
        elif self.elaguer_theta and len(self.le0) > 1:
            vivant = (self.lpe - self.lpe.max(-1, keepdims=True) >= np.log(self.seuil_theta)).any((0, 1))
            if vivant.sum() < 0.8 * len(vivant):
                self.le0, self.beta = self.le0[vivant], self.beta[vivant]
                self.P, self.lpe = np.ascontiguousarray(self.P[..., vivant]), self.lpe[..., vivant]

    def _regriller_si_besoin(self, NA=41, NB=31):
        """(2) θ = (a, b), ln η_k = a + b·(k − τ)/1000 ; nouvelle grille à ±6 sd quand le pas dépasse 0,6 sd, que (a, b)
        sont corrélés à plus de 0,5, ou que la masse touche un bord ; τ choisi pour décorréler ; interpolation bilinéaire
        (log-densité de θ, et p(c | θ) de chaque ROI)."""
        import scipy.sparse as sp
        a, b, tau, k = self.a_vals, self.b_vals, self.tau, self.k
        A = self.le0 + self.beta * tau / 1000.0; Bv = self.beta
        pe = self.pe(); mA, mB = pe @ A, pe @ Bv; dA, dB = A - mA, Bv - mB
        vA, vB, cAB = pe @ dA ** 2, pe @ dB ** 2, pe @ (dA * dB)
        spA, spB = a[1] - a[0], (b[1] - b[0]) if len(b) > 1 else np.inf
        corr = cAB / np.sqrt(vA * vB + 1e-300)
        Pm = pe.reshape(len(a), len(b)); bmax = self.f['derive_max']
        bord = Pm[:2].sum() + Pm[-2:].sum() + (Pm[:, :2].sum() if b[0] > -bmax + 1e-12 else 0) \
            + (Pm[:, -2:].sum() if b[-1] < bmax - 1e-12 else 0)
        if not (spA > 0.6 * np.sqrt(vA) or spB > 0.6 * np.sqrt(vB) or abs(corr) > 0.5 or bord > 1e-4):
            return
        tau2 = float(np.clip(tau - 1000.0 * cAB / max(vB, 1e-300), 0.0, k)); dt = (tau2 - tau) / 1000.0
        mA2 = mA + dt * mB; sdA = np.sqrt(max(vA + 2 * dt * cAB + dt * dt * vB, 1e-14)); sdB = np.sqrt(max(vB, 1e-14))
        a2 = np.linspace(mA2 - 6 * sdA, mA2 + 6 * sdA, NA)
        b2 = np.linspace(max(mB - 6 * sdB, -bmax), min(mB + 6 * sdB, bmax), NB)
        A2, B2 = (x.ravel() for x in np.meshgrid(a2, b2, indexing='ij')); n2 = A2.size
        fa = (A2 - B2 * dt - a[0]) / spA; fb = (B2 - b[0]) / spB
        na, nb = len(a), len(b)
        dedans = (fa >= 0) & (fa <= na - 1) & (fb >= 0) & (fb <= nb - 1)
        ia = np.clip(np.floor(fa).astype(int), 0, na - 2); ib = np.clip(np.floor(fb).astype(int), 0, nb - 2)
        ua, ub = np.clip(fa - ia, 0, 1), np.clip(fb - ib, 0, 1)
        r, cc, v = [], [], []
        for da, db, wt in ((0, 0, (1 - ua) * (1 - ub)), (1, 0, ua * (1 - ub)), (0, 1, (1 - ua) * ub), (1, 1, ua * ub)):
            r.append((ia + da) * nb + ib + db); cc.append(np.arange(n2)); v.append(wt * dedans)
        W = sp.csr_matrix((np.concatenate(v), (np.concatenate(cc), np.concatenate(r))), shape=(n2, na * nb))
        lpe2 = np.empty((len(self.K), 1, n2)); P2 = np.empty(self.P.shape[:3] + (n2,), np.float32)
        for i in range(len(self.K)):
            ld = np.where(dedans, W @ np.maximum(self.lpe[i, 0], -1e4), -1e4); lpe2[i, 0] = ld - _lse(ld)
            for m in range(self.M):
                Pn = np.asarray(W @ self.P[i, m].T.astype(np.float64)).T
                s = Pn.sum(0); Pn[:, s <= 0] = 1.0 / Pn.shape[0]; s[s <= 0] = 1.0
                P2[i, m] = Pn / s
        self.a_vals, self.b_vals, self.tau = a2, b2, tau2
        self.le0, self.beta = A2 - B2 * tau2 / 1000.0, B2
        self.P, self.lpe = P2, lpe2; self.n_regrilles += 1

    def resume(self, le):
        """c (moyenne, intervalle à 95 %) par ROI ; η de la trame (moyenne géométrique, sd(ln η)) ; poids des noyaux."""
        w = self.poids(); pes = np.exp(self.lpe)                                   # (n_noyaux, groupes), (…, n_θ)
        wc = np.zeros((self.n_noyaux, w.shape[1])); wc[self.ids] = w
        grp = np.zeros(self.M, int) if self.partage else np.arange(self.M)
        pc = sum(w[i, grp, None] * np.einsum('mct,mt->mc', self.P[i], pes[i, grp].astype(np.float32))
                 for i in range(len(self.K)))
        m_i = pes @ le; v_i = np.einsum('kgt,kgt->kg', pes, (le[None, None, :] - m_i[..., None]) ** 2)
        m = (w * m_i).sum(0); sd = np.sqrt(np.maximum((w * (v_i + m_i ** 2)).sum(0) - m ** 2, 0))
        c_, lo, hi = resume_c(pc, self.g.c)
        if self.partage:
            return dict(c=c_, bas=lo, haut=hi, eta=float(np.exp(m[0])), eta_sd=float(sd[0]), poids=wc[:, 0])
        return dict(c=c_, bas=lo, haut=hi, eta=np.exp(m), eta_sd=sd, poids=wc)


def lisseur_derive(fl, formes, NB, NY, masse=1 - 1e-4, n_max=None, memoire=4e8):
    """
    Lissage hors ligne du modèle 1-2-3 pour une ROI (fl : FiltreTempsReel à 1 ROI, passé sur toutes les trames).
    θ étant statique, p(c_k | tout) = Σ_noyau w Σ_θ p(θ | tout) p(c_k | θ, tout) : pour chaque noyau, on garde les θ
    les plus probables jusqu'à la masse `masse` (au plus n_max, et ce que permet `memoire` octets), puis lissage
    avant-arrière exact de c sachant θ. Les poids des noyaux sont ceux de toutes les trames. η : ln η_k(θ) est connu
    pour chaque θ, sa loi lissée vient de p(θ | tout), sans troncature.
    formes : (T, n_c) morceau « quand » ; NB, NY : (T,) totaux (fond déjà soustrait). Renvoie (p(c_k), moyenne de
    ln η_k, sd, masse de θ gardée par noyau).
    """
    T, nc = formes.shape; w = fl.poids()[:, 0]; pes = np.exp(fl.lpe[:, 0])
    n_max = n_max or max(1, int(memoire / (4 * T * nc)))
    pc = np.zeros((T, nc)); garde = []
    for i, K in enumerate(fl.K):
        o = np.argsort(pes[i])[::-1]; cum = np.cumsum(pes[i][o])
        sel = o[:min(n_max, int(np.searchsorted(cum, masse)) + 1)]; p_sel = pes[i][sel] / pes[i][sel].sum()
        garde.append(float(pes[i][sel].sum()))
        K64 = K.astype(np.float64); alpha = np.empty((T, nc, len(sel)), np.float32); a = None
        for k in range(T):
            lB, lY = logs_ratio_eta(fl.g, fl.le0[sel] + fl.beta[sel] * k / 1000.0)
            L = formes[k][:, None] + NB[k] * lB + NY[k] * lY
            pred = np.full((nc, len(sel)), 1.0 / nc) if a is None else K64 @ a
            q = pred * np.exp(L - L.max()); a = q / np.maximum(q.sum(0), 1e-300); alpha[k] = a
        post = a; pc[-1] += w[i] * (post @ p_sel)
        for k in range(T - 2, -1, -1):
            a = alpha[k].astype(np.float64); pred = K64 @ a
            post = a * (K64.T @ (post / np.maximum(pred, 1e-300))); post /= np.maximum(post.sum(0), 1e-300)
            pc[k] += w[i] * (post @ p_sel)
    pe = np.einsum('k,kt->t', w, pes); le = fl.le0[None, :] + fl.beta[None, :] * np.arange(T)[:, None] / 1000.0
    m = le @ pe; sd = np.sqrt(np.maximum(((le - m[:, None]) ** 2) @ pe, 0))
    return pc, m, sd, garde
