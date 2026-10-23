"""
Variante « tout » du test de dérive de η (essai.py) : 1-2-3 plus toutes les modifications applicables à la simulation.
  1-2-3  dérive paramétrique de η partagée par les 10 ROI ; dynamique de c apprise (moyenne de modèles sur s_c)
  [8]    composante de saut uniforme (p = 0,005) dans la transition de c
  [7]    grille de θ adaptative : θ = (a, D), ln η_k = a + D·(k/(T−1) − τ) ; recentrée à ±6 sd quand le pas dépasse
         0,6 sd ou que la corrélation (a, D) dépasse 0,5 (τ choisi pour décorréler) ; interpolation bilinéaire
  [2]    décalage d'horloge δ partagé, appris
  [3]    fond additif par ROI et par canal, appris
  [1]    ancre sur la forme bleue seule : la forme jaune (sensible à un changement de l'accepteur seul, photoblanchiment
         du YFP) est retirée de la vraisemblance ; le jaune n'entre que par le partage (η)
         (essai précédent : ρ_A appris en ligne n'est pas identifiable avec c fixé à sa moyenne filtrée — il partait à 1,2-1,26
         pour 0,98 vrai : dans la forme jaune, ρ_A et c agissent presque tous deux par la part de fuite)
  [2][3] : paramètres lents estimés en ligne par EM (vraisemblance des formes, amplitude libre par trame ; l'étape E prend
  c aux quantiles 10/30/50/70/90 % de sa loi filtrée de chaque trame ; δ et fond bleu sur la forme bleue, fond jaune sur
  la forme jaune à ρ_A = 1). À la trame 100 : 4 itérations EM, chacune retraitant les trames 0-99 avec les paramètres du moment (les
  sorties temps réel de ces trames restent celles du premier passage, sans correction) ; ensuite, toutes les 100 trames,
  une estimation sur les 600 dernières trames, sans retraitement.
Non inclus : [4] incertitude de calibration (la simulation utilise la calibration exacte), [5] point de référence in
situ, [6] diagnostics, [9] (remplacé par la moyenne de modèles sur s_c), [10] empilement (absent de la simulation).

Scénarios : A = données du test précédent (même graine, même code) ; B = mêmes trajectoires de c et de A, η de
détection 1 → 1,1, plus : fond additif constant (4 % du bleu et 2 % du jaune d'une trame typique, × U(0,5 ; 1,5) par ROI),
accepteur blanchi linéairement 1 → 0,95. (Un décalage d'horloge de 25 ps a été essayé : sans correction il multiplie c
par 3 à 4,6 ; appris en ligne, δ reste biaisé de +2 à +4 ps, ce qui suffit à biaiser c de 12 à 27 % — retiré de B.)
Usage : python3 tout.py T A|B variante... (variantes : tout, tout_sans_delta, 123, 123moteur, ref)
"""
import sys, os, time, json, numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize_scalar
sys.dont_write_bytecode = True
ICI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ICI)
from moteur import Ratio, filtrer, moyenne_modeles
from src.modele import charger_params, charger_calibration, Senseur
from src.filtre import Grille, grille_c, grille_eta, resume_c, transition
from src.simulation import trame_modele

P_ = charger_params(); f = P_['filtre']; CAL = charger_calibration(); S = Senseur(CAL)
g = Grille(S, grille_c(f), grille_eta(dict(f, n_eta=f['n_eta_lisseur'])))
T = int(sys.argv[1]); SCEN = sys.argv[2]; VARIANTES = sys.argv[3:] or ['tout']
M = 10; N = 1000; N_FEN = S.n_fen; NC = len(g.c); LC = np.log(g.c)
P_U = 0.005; S_C = (0.004, 0.013, 0.04); CANAUX = ('B',); QUANTILES_E = (0.1, 0.3, 0.5, 0.7, 0.9)
DELTA_V, RHO_FIN_V, FOND_V = 0.0, 0.95, (0.04, 0.02)

# ------------------------------------------------------------------ données (A : identique à essai.py)
eta_vrai = np.linspace(1.0, 1.1, T)
rng = np.random.default_rng(20261023); kk = np.arange(T, dtype=float); C = np.empty((M, T))
for m in range(M):
    base, amp = rng.uniform(6, 12), rng.uniform(8, 16); on = int(rng.uniform(0.2, 0.4) * T); off = on + int(rng.uniform(0.2, 0.33) * T)
    c = np.full(T, base); c[on:] = base + amp * (1 - np.exp(-(kk[on:] - on) / (T / 75)))
    c[off:] = base + (c[off - 1] - base) * np.exp(-(kk[off:] - off) / (T / 12)); C[m] = c
A = 1 - 0.3 * kk / (T - 1)
H = np.empty((M, T, 2, N_FEN), np.float32)
rho_vrai = np.linspace(1.0, RHO_FIN_V, T) if SCEN == 'B' else np.ones(T)
if SCEN == 'A':
    for m in range(M):
        for j in range(T):
            H[m, j, 0], H[m, j, 1] = trame_modele(S, C[m, j], N * A[j], rng, eta_vrai[j])
    fond_vrai = (np.zeros(M), np.zeros(M)); delta_vrai = 0.0
else:
    Sv = Senseur(CAL, DELTA_V); b0, y0 = S.mu(10.0); fb = b0.sum() / (b0.sum() + y0.sum())
    echelle = rng.uniform(0.5, 1.5, M)
    fond_vrai = (FOND_V[0] * N * fb / N_FEN * echelle, FOND_V[1] * N * (1 - fb) / N_FEN * echelle); delta_vrai = DELTA_V
    for m in range(M):
        for j in range(T):
            mB, _, p = Sv.mu(C[m, j], 1.0, 1.0, composantes=True)
            y = eta_vrai[j] * (p['fuite'] + rho_vrai[j] * (p['sensibilisee'] + p['directe']))
            a = N * A[j] / (mB.sum() + y.sum())
            H[m, j, 0] = rng.poisson(a * mB + fond_vrai[0][m]); H[m, j, 1] = rng.poisson(a * y + fond_vrai[1][m])

# ------------------------------------------------------------------ tables du senseur
LE_FIN = np.linspace(np.log(0.5), np.log(2.0), 4001); DLE = LE_FIN[1] - LE_FIN[0]
idx_fin = lambda le: np.clip(np.rint((le - LE_FIN[0]) / DLE).astype(int), 0, len(LE_FIN) - 1)


class Banque:
    """μ(c) sur la grille de c pour un décalage δ (banque au pas de 5 ps, interpolée) et un facteur d'accepteur ρ."""
    def __init__(self, deltas=np.arange(-100.0, 100.1, 5.0)):
        self.d = deltas; B, F, Ac = [], [], []
        for dl in deltas:
            mB, _, p = Senseur(CAL, dl).mu(g.c, 1.0, 1.0, composantes=True)
            B.append(mB); F.append(p['fuite']); Ac.append(p['sensibilisee'] + p['directe'])
        self.B, self.F, self.Ac = (np.array(x) for x in (B, F, Ac))

    def mu(self, delta, rho):
        x = float(np.clip((delta - self.d[0]) / (self.d[1] - self.d[0]), 0, len(self.d) - 1.000001)); i = int(x); w = x - i
        mB = (1 - w) * self.B[i] + w * self.B[i + 1]
        mY = (1 - w) * (self.F[i] + rho * self.Ac[i]) + w * (self.F[i + 1] + rho * self.Ac[i + 1])
        return mB, mY


class TablesFiltre:
    """Formes et ratio du filtre pour les paramètres lents courants (δ, ρ_A, fonds)."""
    def __init__(self, banque, nuis, canaux=('B', 'Y')):
        self.canaux = canaux
        mB, mY = banque.mu(nuis['delta'], nuis['rho'])
        self.MB, self.MY = mB.sum(1), mY.sum(1)
        self.pB, self.pY = mB / self.MB[:, None], mY / self.MY[:, None]
        self.lpB, self.lpY = np.log(self.pB), np.log(self.pY)
        MB_, MY_ = self.MB[:, None], self.MY[:, None] * np.exp(LE_FIN)[None, :]
        self.lB = np.log(MB_ / (MB_ + MY_)).astype(np.float32); self.lY = np.log(MY_ / (MB_ + MY_)).astype(np.float32)
        self.bB, self.bY = np.asarray(nuis['bB'], float), np.asarray(nuis['bY'], float)

    def formes(self, m, hB, hY):
        out = np.zeros(NC)
        for nom, h, p, lp, b in (('B', hB, self.pB, self.lpB, self.bB[m]), ('Y', hY, self.pY, self.lpY, self.bY[m])):
            if nom not in self.canaux:
                continue
            if b <= 0:
                out += lp @ h
            else:
                S_ = max(h.sum() - N_FEN * b, 1.0); out += np.log((S_ * p + b) / (S_ + N_FEN * b)) @ h
        return out

    def totaux(self, m, hB, hY):
        return max(hB.sum() - N_FEN * self.bB[m], 1.0), max(hY.sum() - N_FEN * self.bY[m], 1.0)


# ------------------------------------------------------------------ grille de θ
LN0_MIN, LN0_MAX, D_MIN, D_MAX = np.log(0.7), np.log(1.5), -0.15, 0.15
NA, ND = 41, 31


class GrilleTheta:
    def __init__(self, a, D, tau=0.0):
        self.a, self.D, self.tau = np.asarray(a, float), np.asarray(D, float), float(tau)
        Ag, Dg = np.meshgrid(self.a, self.D, indexing='ij'); self.A, self.Dv = Ag.ravel(), Dg.ravel(); self.n = self.A.size

    def ln_eta(self, k):
        return self.A + self.Dv * (k / (T - 1) - self.tau)


def grille_initiale():
    return GrilleTheta(np.log(np.geomspace(0.7, 1.5, 81)), np.linspace(D_MIN, D_MAX, 31))


def moments(gr, pe):
    mA, mD = pe @ gr.A, pe @ gr.Dv; dA, dD = gr.A - mA, gr.Dv - mD
    return mA, mD, pe @ dA ** 2, pe @ dD ** 2, pe @ (dA * dD)


class Modele:
    def __init__(self, Kc, gr):
        self.Kc = Kc.astype(np.float32); self.P = np.full((M, NC, gr.n), 1.0 / NC, np.float32)
        self.pe = np.full(gr.n, 1.0 / gr.n); self.cum = 0.0; self.vu = False


class Etat:
    def __init__(self, Kcs, tables):
        self.gr = grille_initiale(); self.modeles = [Modele(K, self.gr) for K in Kcs]; self.tab = tables; self.n_regrilles = 0

    def poids(self):
        cum = np.array([mo.cum for mo in self.modeles]); w = np.exp(cum - cum.max()); return w / w.sum()

    def regriller_si_besoin(self, k):
        gr = self.gr; w = self.poids(); pe = sum(wi * mo.pe for wi, mo in zip(w, self.modeles))
        mA, mD, vA, vD, cAD = moments(gr, pe)
        spA, spD = gr.a[1] - gr.a[0], gr.D[1] - gr.D[0]
        corr = cAD / np.sqrt(vA * vD + 1e-300)
        Pm = pe.reshape(len(gr.a), len(gr.D)); bord = Pm[:2].sum() + Pm[-2:].sum()
        if gr.D[0] > D_MIN + 1e-9: bord += Pm[:, :2].sum()
        if gr.D[-1] < D_MAX - 1e-9: bord += Pm[:, -2:].sum()
        if not (spA > 0.6 * np.sqrt(vA) or spD > 0.6 * np.sqrt(vD) or abs(corr) > 0.5 or bord > 1e-4):
            return
        tau2 = float(np.clip(gr.tau - cAD / max(vD, 1e-300), 0.0, k / (T - 1))); dt = tau2 - gr.tau
        mA2 = mA + dt * mD; sdA = np.sqrt(max(vA + 2 * dt * cAD + dt * dt * vD, 1e-14)); sdD = np.sqrt(max(vD, 1e-14))
        lo, hi = max(mD - 6 * sdD, D_MIN), min(mD + 6 * sdD, D_MAX)
        ng = GrilleTheta(np.linspace(mA2 - 6 * sdA, mA2 + 6 * sdA, NA), np.linspace(lo, hi, ND), tau2)
        fa = (ng.A - ng.Dv * dt - gr.a[0]) / spA; fD = (ng.Dv - gr.D[0]) / spD
        na, nd = len(gr.a), len(gr.D)
        dedans = (fa >= 0) & (fa <= na - 1) & (fD >= 0) & (fD <= nd - 1)
        ia = np.clip(np.floor(fa).astype(int), 0, na - 2); iD = np.clip(np.floor(fD).astype(int), 0, nd - 2)
        ua, uD = np.clip(fa - ia, 0, 1), np.clip(fD - iD, 0, 1)
        r, cidx, v = [], [], []
        for da, dd, wt in ((0, 0, (1 - ua) * (1 - uD)), (1, 0, ua * (1 - uD)), (0, 1, (1 - ua) * uD), (1, 1, ua * uD)):
            r.append((ia + da) * nd + iD + dd); cidx.append(np.arange(ng.n)); v.append(wt * dedans)
        Wt = sp.csr_matrix((np.concatenate(v), (np.concatenate(cidx), np.concatenate(r))), shape=(ng.n, gr.n))
        for mo in self.modeles:
            ld = Wt @ np.log(np.maximum(mo.pe, 1e-300))          # log-densité interpolée (même volume de case partout)
            ld = np.where(dedans, ld, -np.inf); pe = np.exp(ld - ld[dedans].max()); mo.pe = pe / pe.sum()
            P2 = np.empty((M, NC, ng.n), np.float32)
            for m in range(M):
                Pn = np.asarray(Wt @ mo.P[m].T.astype(np.float64)).T
                s = Pn.sum(0); Pn[:, s <= 0] = 1.0 / NC; s[s <= 0] = 1.0
                P2[m] = Pn / s
            mo.P = P2
        self.gr = ng; self.n_regrilles += 1

    def trame(self, k, chat, sortie=None, adaptatif=True):
        if adaptatif and k >= 20 and k % 10 == 0:
            self.regriller_si_besoin(k)
        lek = self.gr.ln_eta(k); ix = idx_fin(lek)
        lB, lY = self.tab.lB[:, ix], self.tab.lY[:, ix]
        E, tot = [], 0.0
        for m in range(M):
            hB, hY = H[m, k, 0].astype(float), H[m, k, 1].astype(float)
            fo = self.tab.formes(m, hB, hY).astype(np.float32); NBp, NYp = self.tab.totaux(m, hB, hY)
            L = fo[:, None] + np.float32(NBp) * lB + np.float32(NYp) * lY
            mx = L.max(); E.append(np.exp(L - mx)); tot += float(mx)
        for mo in self.modeles:
            if mo.vu:
                for m in range(M):
                    mo.P[m] = mo.Kc @ mo.P[m]
            mo.vu = True
            lse = np.log(np.maximum(mo.pe, 1e-300))
            for m in range(M):
                Q = mo.P[m] * E[m]; ell = np.maximum(Q.sum(0, dtype=np.float64), 1e-300)
                mo.P[m] = Q / np.maximum(ell, 1e-30).astype(np.float32); lse = lse + np.log(ell)
            mx2 = lse.max(); w = np.exp(lse - mx2); z = w.sum(); mo.pe = w / z; mo.cum += np.log(z) + mx2 + tot
        w = self.poids()
        pc = sum(wi * np.stack([mo.P[m] @ mo.pe.astype(np.float32) for m in range(M)]) for wi, mo in zip(w, self.modeles))
        cdf = np.cumsum(pc, 1)
        chat[:, k] = np.stack([np.minimum(np.argmax(cdf >= q, 1), NC - 1) for q in QUANTILES_E], 1)
        if sortie is not None:
            me = np.array([mo.pe @ lek for mo in self.modeles]); ve = np.array([mo.pe @ (lek - x) ** 2 for mo, x in zip(self.modeles, me)])
            sortie['pc'][:, k] = pc; sortie['m_eta'][k] = w @ me; sortie['v_eta'][k] = max(w @ (ve + me ** 2) - (w @ me) ** 2, 0)
            sortie['w'][k] = w


def mstep(banque, nuis, chat, k_hi, fen, apprendre_delta=True):
    """Paramètres lents : maximum de vraisemblance des formes (amplitude libre par trame et par canal), c fixé."""
    lo = max(0, k_hi - fen); n = k_hi - lo
    hB = H[:, lo:k_hi, 0].reshape(-1, N_FEN).astype(float); hY = H[:, lo:k_hi, 1].reshape(-1, N_FEN).astype(float)
    nq = len(QUANTILES_E)                          # chaque trame : nq pseudo-trames de poids 1/nq (quadrature de E_c)
    ci = chat[:, lo:k_hi].ravel(); roi = np.repeat(np.arange(M), n * nq)
    hB, hY = np.repeat(hB, nq, 0), np.repeat(hY, nq, 0); NBr, NYr = hB.sum(1), hY.sum(1)

    def ll(h, Nr, prof, b):
        S_ = np.maximum(Nr - N_FEN * b, 1.0)
        return (h * np.log((S_[:, None] * prof + b[:, None]) / (S_ + N_FEN * b)[:, None])).sum() / nq

    def profils(delta, rho):
        mB, mY = banque.mu(delta, rho)
        return (mB / mB.sum(1, keepdims=True))[ci], (mY / mY.sum(1, keepdims=True))[ci]

    d, r = nuis['delta'], 1.0; bB, bY = np.array(nuis['bB'], float), np.array(nuis['bY'], float)
    for it in range(2):
        pB, pY = profils(d, r)
        for m in range(M):
            s = roi == m
            for h, Nr, p, b in ((hB, NBr, pB, bB), (hY, NYr, pY, bY)) if it == 1 else ((hB, NBr, pB, bB),):
                bmax = 0.25 * Nr[s].mean() / N_FEN
                b[m] = minimize_scalar(lambda x: -ll(h[s], Nr[s], p[s], np.full(s.sum(), x)), bounds=(0.0, bmax),
                                       method='bounded', options=dict(xatol=1e-4 * bmax)).x
        bBr = bB[roi]
        if apprendre_delta:
            d = minimize_scalar(lambda x: -ll(hB, NBr, profils(x, r)[0], bBr), bounds=(-100.0, 100.0), method='bounded',
                                options=dict(xatol=0.05)).x
    return dict(delta=float(d), rho=1.0, bB=bB, bY=bY)


def moteur_tout(Kcs, adaptatif=True, nuisances=True, journal=None, apprendre_delta=True, canaux=CANAUX):
    banque = Banque() if nuisances else None
    nuis = dict(delta=0.0, rho=1.0, bB=np.zeros(M), bY=np.zeros(M))
    tab0 = TablesFiltre(banque, nuis, canaux) if nuisances else None
    if not nuisances:                                         # tables par défaut, identiques à g et Ratio(g)
        tab0 = TablesFiltre.__new__(TablesFiltre); tab0.canaux = canaux
        tab0.MB, tab0.MY, tab0.pB, tab0.pY = g.MB, g.MY, g.pB, g.pY; tab0.lpB, tab0.lpY = g.lpB, g.lpY
        r_ = Ratio(g); tab0.lB, tab0.lY = r_.lB.astype(np.float32), r_.lY.astype(np.float32); tab0.bB = tab0.bY = np.zeros(M)
    sortie = dict(pc=np.empty((M, T, NC), np.float32), m_eta=np.empty(T), v_eta=np.empty(T), w=np.empty((T, len(Kcs))), hist=[])
    chat = np.zeros((M, T, len(QUANTILES_E)), int); et = Etat(Kcs, tab0); t0 = time.time()
    for k in range(T):
        if nuisances and k >= 100 and k % 100 == 0:
            if k == 100:
                for _ in range(4):
                    nuis = mstep(banque, nuis, chat, k, k, apprendre_delta)
                    et = Etat(Kcs, TablesFiltre(banque, nuis, canaux))
                    for k2 in range(k):
                        et.trame(k2, chat, None, adaptatif)
            else:
                nuis = mstep(banque, nuis, chat, k, 600, apprendre_delta); et.tab = TablesFiltre(banque, nuis, canaux)
            sortie['hist'].append(dict(trame=k, delta=nuis['delta'], rho=nuis['rho'], bB=nuis['bB'].tolist(), bY=nuis['bY'].tolist()))
            if journal:
                journal(f"  trame {k}: δ = {nuis['delta']:+.2f} ps (vrai {delta_vrai:+.0f}), "
                        f"fond bleu {np.mean(nuis['bB']):.4f} (vrai {np.mean(fond_vrai[0]):.4f}), jaune {np.mean(nuis['bY']):.4f} "
                        f"(vrai {np.mean(fond_vrai[1]):.4f}) photons/case/trame | {et.n_regrilles} regrilles | {time.time() - t0:.0f} s")
        et.trame(k, chat, sortie, adaptatif)
    sortie['n_regrilles'] = et.n_regrilles; sortie['nuis'] = nuis
    return sortie


def resultat(nom, pc, m_eta, v_eta, ROIs):
    moy, lo, hi = resume_c(pc.reshape(-1, NC), g.c)
    moy, lo, hi = (x.reshape(len(ROIs), T) for x in (moy, lo, hi)); Cv = C[ROIs]
    err = moy - Cv; e_eta = m_eta - np.log(eta_vrai); sl = slice(T // 3, None)
    return dict(variante=nom, rmse_c=float(np.sqrt(np.mean(err ** 2))), rmse_c_apres_1tiers=float(np.sqrt(np.mean(err[:, sl] ** 2))),
                biais_c=float(err.mean()), couverture=float(np.mean((lo <= Cv) & (Cv <= hi))),
                eta_rms_pct=float(100 * np.sqrt(np.mean(e_eta ** 2))), eta_rms_apres_1tiers_pct=float(100 * np.sqrt(np.mean(e_eta[sl] ** 2))),
                eta_1000=float(np.exp(m_eta[min(999, T - 1)])), eta_2000=float(np.exp(m_eta[min(1999, T - 1)])), eta_final=float(np.exp(m_eta[-1])),
                eta_final_sd_pct=float(100 * np.sqrt(v_eta[-1])))


if __name__ == '__main__':
    os.makedirs(os.path.join(ICI, 'sorties'), exist_ok=True)
    sortie_json = os.path.join(ICI, 'sorties', f'resultats_tout_{SCEN}_{T}.json')
    res = json.load(open(sortie_json)) if os.path.exists(sortie_json) else []
    log = lambda s: print(s, flush=True)
    lc = np.log(g.c)
    for v in VARIANTES:
        t0 = time.time()
        if v in ('tout', 'tout_sans_delta'):
            Kcs = [(1 - P_U) * transition(lc, s, f['s_saut'], f['w']) + P_U / NC for s in S_C]
            so = moteur_tout(Kcs, True, True, log, v == 'tout')
            nom = 'tout : 1-2-3 + [1][2][3][7][8]' if v == 'tout' else 'tout sans [2] : 1-2-3 + [1][3][7][8]'
            r = resultat(nom, so['pc'], so['m_eta'], so['v_eta'], list(range(M)))
            r.update(poids_s_c_finaux=so['w'][-1].tolist(), n_regrilles=so['n_regrilles'], delta_final=so['nuis']['delta'],
                     rho_vrai_final=float(rho_vrai[-1]), delta_vrai=delta_vrai,
                     fond_bleu_estime=so['nuis']['bB'].tolist(), fond_bleu_vrai=fond_vrai[0].tolist(),
                     fond_jaune_estime=so['nuis']['bY'].tolist(), fond_jaune_vrai=fond_vrai[1].tolist(), historique=so['hist'])
        elif v.startswith('abl'):                           # 1-2-3 + une seule modification : abl1, abl3, abl7, abl8
            x = v[3:]; pu = P_U if '8' in x else 0.0
            Kcs = [(1 - pu) * transition(lc, s, f['s_saut'], f['w']) + pu / NC for s in S_C]
            so = moteur_tout(Kcs, '7' in x, '3' in x, log, False, ('B',) if '1' in x else ('B', 'Y'))
            r = resultat('1-2-3 + [' + ']['.join(x) + '] seul', so['pc'], so['m_eta'], so['v_eta'], list(range(M)))
            r.update(poids_s_c_finaux=so['w'][-1].tolist(), n_regrilles=so['n_regrilles'])
        elif v == '123':                                    # même moteur, sans aucune modification : doit égaler 1-2-3
            so = moteur_tout([transition(lc, s, f['s_saut'], f['w']) for s in S_C], False, False, canaux=('B', 'Y'))
            r = resultat('1-2-3 (moteur partagé)', so['pc'], so['m_eta'], so['v_eta'], list(range(M)))
            r.update(poids_s_c_finaux=so['w'][-1].tolist())
        elif v in ('123moteur', 'ref'):
            formes = np.empty((M, T, NC)); NB = H[:, :, 0].sum(-1).astype(float); NY = H[:, :, 1].sum(-1).astype(float)
            for m in range(M):
                for j in range(T):
                    formes[m, j] = g.formes(H[m, j, 0].astype(float), H[m, j, 1].astype(float))
            ratio = Ratio(g)
            if v == '123moteur':
                L0, DD = np.meshgrid(np.log(np.geomspace(0.7, 1.5, 81)), np.linspace(-0.15, 0.15, 31), indexing='ij')
                L0, DD = L0.ravel(), DD.ravel(); param = lambda k_: L0 + DD * k_ / (T - 1)
                s = [filtrer(formes, NB, NY, ratio, transition(lc, sc, f['s_saut'], f['w']), param) for sc in S_C]
                pc, me, ve, w = moyenne_modeles(s); r = resultat('1-2-3 (moteur.py)', pc, me, ve, list(range(M)))
            else:
                le_g = np.log(g.eta); Kc = transition(lc, 0.04, f['s_saut'], f['w']); rs = []
                for m in range(M):
                    s = [filtrer(formes[m:m + 1], NB[m:m + 1], NY[m:m + 1], ratio, Kc, lambda k_: le_g, transition(le_g, se))
                         for se in f['s_eta_candidats']]
                    pc, me, ve, _ = moyenne_modeles(s); rs.append(resultat('', pc, me, ve, [m]))
                r = {k_: float(np.mean([x[k_] for x in rs])) for k_ in rs[0] if k_ != 'variante'}
                r = dict(variante='référence : η marche aléatoire, 1 ROI', **r)
        r['duree_s'] = time.time() - t0
        log({k_: (round(x_, 4) if isinstance(x_, float) else x_) for k_, x_ in r.items() if k_ != 'historique'})
        res = json.load(open(sortie_json)) if os.path.exists(sortie_json) else []     # relu : plusieurs processus écrivent
        res = [x for x in res if x['variante'] != r['variante']] + [r]
        json.dump(res, open(sortie_json, 'w'), indent=1, default=float)
