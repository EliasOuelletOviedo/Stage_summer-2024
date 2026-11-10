#!/usr/bin/env python3
"""
run_tests.py — les vérifications du modèle et du filtre (après step0_calibration.py).

Lancer depuis le dossier bayesian_filter :  python3 test/run_tests.py
"""
import os, sys, traceback
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.modele import charger_params, charger_calibration, Senseur, isotherme, x_vers_c
from src.filtre import (Grille, grille_c, grille_eta, portes_log, transition, filtre_lisseur, FiltreTempsReel,
                        moyenne_de_modeles, resume_c, resume_eta, grille_theta, noyaux_c, logs_ratio_eta, lisseur_derive)
from src.simulation import trame_modele, Mesures

P = charger_params(); CAL = charger_calibration(); F = P['filtre']
S = Senseur(CAL)
G = Grille(S, grille_c(F), grille_eta(F))
RESULTATS = []


def test(f):
    try:
        f(); RESULTATS.append((f.__name__, True, ''))
    except Exception:
        RESULTATS.append((f.__name__, False, traceback.format_exc(limit=2)))
    return f


@test
def calibration_comme_le_document():
    """Annexe (section 17) : durées de vie, λ, G, g₂, κ, isotherme de x."""
    k, x = CAL['constantes'], CAL['isothermes']['x']
    for cle, doc, tol in (('t1', 0.0868, 0.001), ('t2', 0.855, 0.01), ('tD', 3.755, 0.005), ('tA', 3.483, 0.005),
                          ('lam', 0.764, 0.01), ('G', 2.34, 0.02), ('G2r', 0.758, 0.01), ('k1', 0.658, 0.01)):
        assert abs(k[cle] - doc) < tol, (cle, k[cle], doc)
    assert abs(x['fmax'] - 0.548) < 0.002 and abs(x['fmin'] - 0.017) < 0.002 and abs(x['K'] - 62.5) < 0.5


@test
def isotherme_inversible():
    iso = CAL['isothermes']['x']
    for c in (0.5, 10.0, 131.0):
        assert abs(x_vers_c(isotherme(c, iso['fmax'], iso['fmin'], iso['K']), iso) - c) < 1e-9


@test
def amplitude_et_eta_hors_des_formes():
    """Les formes ne dépendent ni de A ni de η ; η ne touche que le jaune, A les deux canaux (section 11)."""
    b1, y1 = S.mu(10.0, 1.0, 1.0); b2, y2 = S.mu(10.0, 3.0, 1.4)
    assert np.allclose(b1 / b1.sum(), b2 / b2.sum()) and np.allclose(y1 / y1.sum(), y2 / y2.sum())
    assert np.isclose(y2.sum() / b2.sum(), 1.4 * y1.sum() / b1.sum())


@test
def decoupage_exact():
    """ln L de Poisson complète, A profilé, = morceau ratio + morceau FLIM à une constante près (section 11)."""
    rng = np.random.default_rng(0); B, Y = trame_modele(S, 12.0, 8000, rng, 1.1)
    cs, es = G.c[::37], G.eta[::15]; L = G.loglik(B, Y)[::37, ::15]; ref = np.zeros_like(L)
    for i, c in enumerate(cs):
        mB, mY = S.mu(c)
        for j, e in enumerate(es):
            A = (B.sum() + Y.sum()) / (mB.sum() + e * mY.sum())
            muB, muY = A * mB, A * e * mY
            ref[i, j] = (B * np.log(muB) - muB).sum() + (Y * np.log(muY) - muY).sum()
    d = L - ref
    assert np.ptp(d) < 1e-6 * np.abs(ref).max(), np.ptp(d)


@test
def portes_conservent_les_photons():
    Gp = portes_log(S.n_fen, 16)
    assert Gp.shape[0] == 16 and np.all(Gp.sum(0) == 1)
    g = Grille(S, G.c, G.eta, Gp)
    rng = np.random.default_rng(1); B, Y = trame_modele(S, 10.0, 5000, rng)
    assert np.isclose(g.regrouper(B).sum(), B.sum())


@test
def portes_gardent_l_information_flim():
    """Section 14 : 16 portes gardent environ 93 % de l'information FLIM sur c (formes des deux canaux)."""
    Gp = portes_log(S.n_fen, 16)
    def info(Gm):
        h = 0.01; out = 0.0
        for k in (0, 1):
            m = [np.asarray(S.mu(c)[k]) @ Gm.T for c in (10 - h, 10.0, 10 + h)]
            p = [x / x.sum() for x in m]
            out += ((p[2] - p[0]) / (2 * h)) ** 2 @ (1 / p[1])
        return out
    r = info(Gp) / info(np.eye(S.n_fen))
    assert 0.88 < r <= 1.0, r


@test
def transition_normalisee():
    K = transition(np.log(G.c), F['s'], F['s_saut'], F['w'])
    assert np.allclose(K.sum(0), 1.0)


@test
def ratio_avec_eta_connu_sans_biais():
    """Lecture du ratio sur l'histogramme attendu : c retrouvé exactement si η est connu."""
    for c in (2.0, 10.0, 40.0):
        b, y = S.mu(c, 5000.0, 1.3)
        assert abs(G.lire_ratio(b, y, eta=1.3) - c) < 0.02 * c


@test
def fond_connu_retrouve_c():
    """Avec un fond plat connu, le maximum de vraisemblance revient sur la vraie concentration."""
    b, y = S.mu(10.0, 1e5, 1.2); fb = 0.03 * b.sum() / len(b)
    gf = Grille(S, np.geomspace(5, 20, 1501), np.geomspace(1.0, 1.4, 401))
    L = gf.loglik(b + fb, y, (fb, 0.0)); i, j = np.unravel_index(np.argmax(L), L.shape)
    assert abs(gf.c[i] - 10.0) < 0.1 and abs(gf.eta[j] - 1.2) < 0.003


@test
def filtre_temps_reel_retrouve_eta():
    """Filtre par défaut (1-2-3), 10 ROI à 3 000 photons, c statique, η constant : η à mieux que 3 sd après 200 trames."""
    rng = np.random.default_rng(3); M = 10
    fl = FiltreTempsReel(G, M, F)
    cs = rng.uniform(6, 20, M)
    for _ in range(200):
        H = [trame_modele(S, c, 3000, rng, 1.25) for c in cs]
        x = fl.mettre_a_jour(np.array([h[0] for h in H]), np.array([h[1] for h in H]))
    assert abs(np.log(x['eta'] / 1.25)) < 3 * x['eta_sd'], x


@test
def eta_jamais_fige():
    """5 ROI, 10 000 photons, 150 trames : θ reste mis à jour (sa distribution bouge encore à la dernière trame)."""
    rng = np.random.default_rng(4)
    fl = FiltreTempsReel(G, 5, F)
    for _ in range(150):
        avant = fl.pe().copy()
        H = [trame_modele(S, 10.0, 10000, rng, 1.25) for _ in range(5)]
        x = fl.mettre_a_jour(np.array([h[0] for h in H]), np.array([h[1] for h in H]))
    assert not hasattr(fl, 'eta_fige') and np.abs(fl.pe() - avant).max() > 0 and x['eta_sd'] < 0.01


def _petit_probleme(n_c, T, M, graine):
    """Une grille de c de n_c valeurs, 4-5 valeurs de η0, 3 dérives, 2 vitesses de c : de quoi tout calculer à la main."""
    c = np.geomspace(3, 40, n_c); g = Grille(S, c, np.array([1.0]))
    f = dict(F, s_c_candidats=[0.05, 0.2], eta0_min=0.9, eta0_max=1.3, n_eta0=5, derive_max=20.0, n_derive=3,
             adaptatif=False, elaguer_noyaux=False)                 # numba reste : il doit être exact
    rng = np.random.default_rng(graine)
    H = [[trame_modele(S, (8.0 + 12 * m) * (1 + 0.1 * k), 2000, rng, 1.1 + 0.01 * k) for m in range(M)] for k in range(T)]
    return c, g, f, H


@test
def filtre_123_egal_a_la_force_brute():
    """2 ROI : la factorisation p(θ) Π p(c_m | θ) et la moyenne sur les vitesses de c redonnent la loi conjointe
    (c_1, c_2, θ) calculée directement, trame par trame (moyennes de c, η, poids des modèles)."""
    c, g, f, H = _petit_probleme(12, 6, 2, 1); M, T = 2, 6
    fl = FiltreTempsReel(g, M, f); out = [fl.mettre_a_jour(np.array([h[0] for h in Hk]), np.array([h[1] for h in Hk])) for Hk in H]
    le0, beta = grille_theta(f); nt, nc = len(le0), len(c); cum, pcs, etas = [], [], []
    for K in noyaux_c(c, f):
        p = np.full((nc, nc, nt), 1.0 / (nc * nc * nt)); lev = 0.0; pk, ek, ck = [], [], []
        for k in range(T):
            le = le0 + beta * k / 1000; lB, lY = logs_ratio_eta(g, le)
            L = [g.formes(*H[k][m])[:, None] + H[k][m][0].sum() * lB + H[k][m][1].sum() * lY for m in range(M)]
            if k > 0:
                p = np.einsum('ia,jb,abt->ijt', K, K, p)
            q = p * np.exp(L[0][:, None, :] + L[1][None, :, :] - L[0].max() - L[1].max()); z = q.sum()
            lev += np.log(z) + L[0].max() + L[1].max(); p = q / z
            pk.append(np.stack([p.sum((1, 2)), p.sum((0, 2))])); ek.append((p.sum((0, 1)) * le).sum()); ck.append(lev)
        cum.append(ck); pcs.append(pk); etas.append(ek)
    cum = np.array(cum); w = np.exp(cum - cum.max(0)); w /= w.sum(0)
    for k in range(T):
        c_bf = resume_c(sum(w[i, k] * pcs[i][k] for i in range(2)), c)[0]
        assert np.allclose(out[k]['c'], c_bf, rtol=1e-5), (k, out[k]['c'], c_bf)
        assert abs(np.log(out[k]['eta']) - sum(w[i, k] * etas[i][k] for i in range(2))) < 1e-6
        assert np.allclose(out[k]['poids'], w[:, k], atol=1e-6)


@test
def lisseur_derive_egal_a_la_force_brute():
    """1 ROI : le lissage par θ (lisseur_derive) redonne l'avant-arrière sur la grille conjointe (c, θ)."""
    c, g, f, H = _petit_probleme(15, 8, 1, 2); T = 8; H = [Hk[0] for Hk in H]
    fl = FiltreTempsReel(g, 1, f); [fl.mettre_a_jour(B[None], Y[None]) for B, Y in H]
    formes = np.array([g.formes(B, Y) for B, Y in H]); NB = np.array([B.sum() for B, _ in H]); NY = np.array([Y.sum() for _, Y in H])
    pc, m, sd, garde = lisseur_derive(fl, formes, NB, NY)
    le0, beta = grille_theta(f); nt, nc = len(le0), len(c); res, ev = [], []
    for K in noyaux_c(c, f):
        Ls = [formes[k][:, None] + NB[k] * logs_ratio_eta(g, le0 + beta * k / 1000)[0]
              + NY[k] * logs_ratio_eta(g, le0 + beta * k / 1000)[1] for k in range(T)]
        a = []; p = np.full((nc, nt), 1.0 / (nc * nt)); lev = 0.0
        for k in range(T):
            q = (p if k == 0 else K @ p) * np.exp(Ls[k] - Ls[k].max()); z = q.sum(); lev += np.log(z) + Ls[k].max(); p = q / z; a.append(p)
        post = [None] * T; post[-1] = a[-1]
        for k in range(T - 2, -1, -1):
            post[k] = a[k] * (K.T @ (post[k + 1] / (K @ a[k]))); post[k] /= post[k].sum()
        res.append(post); ev.append(lev)
    w = np.exp(np.array(ev) - max(ev)); w /= w.sum()
    pc_bf = np.array([sum(w[i] * res[i][k].sum(1) for i in range(2)) for k in range(T)])
    assert np.abs(pc - pc_bf).max() < 1e-6 and min(garde) > 1 - 1e-4, np.abs(pc - pc_bf).max()


@test
def filtre_123_suit_une_derive():
    """Filtre par défaut, 3 ROI à c constant (10 mM), 20 000 photons, η de 1 à 1,02 sur 600 trames : η final à moins de
    1 % et de 3 sd, et la dérive est préférée à η constant."""
    rng = np.random.default_rng(8); M, T = 3, 600; eta = np.linspace(1.0, 1.02, T)
    fl = FiltreTempsReel(G, M, F)
    for k in range(T):
        H = [trame_modele(S, 10.0, 20000, rng, eta[k]) for _ in range(M)]
        x = fl.mettre_a_jour(np.array([h[0] for h in H]), np.array([h[1] for h in H]))
    e = np.log(x['eta'] / eta[-1]); b = (fl.pe() * fl.beta).sum()
    assert abs(e) < 0.01 and abs(e) < 3 * x['eta_sd'] and b > 0.01, (x['eta'], x['eta_sd'], b)


@test
def accelerations_proches_du_filtre_exact():
    """Grille adaptative + noyaux élagués + numba (le défaut) contre le filtre sans accélération : 5 ROI, 3 000 photons,
    200 trames, η de 1 à 1,01 — c à 0,1 mM près, η à 0,2 % près, à chaque trame."""
    rng = np.random.default_rng(9); M, T = 5, 200; eta = np.linspace(1.0, 1.01, T); cs = rng.uniform(6, 25, M)
    H = [[trame_modele(S, cs[m] * (1.5 if 80 <= k < 140 else 1.0), 3000, rng, eta[k]) for m in range(M)] for k in range(T)]
    vite = FiltreTempsReel(G, M, F); exact = FiltreTempsReel(G, M, F, adaptatif=False, elaguer_noyaux=False, numba=False)
    for Hk in H:
        B, Y = np.array([h[0] for h in Hk]), np.array([h[1] for h in Hk])
        a, b = vite.mettre_a_jour(B, Y), exact.mettre_a_jour(B, Y)
        assert np.abs(a['c'] - b['c']).max() < 0.1 and abs(np.log(a['eta'] / b['eta'])) < 0.002, (a['c'], b['c'], a['eta'], b['eta'])
    assert vite.n_regrilles > 0


@test
def lisseur_couvre_la_verite():
    """Une ROI, η inconnu, 60 trames à 10 000 photons : l'intervalle lissé à 95 % contient c presque partout."""
    rng = np.random.default_rng(5); c = np.r_[np.full(30, 8.0), np.full(30, 16.0)]
    Ls = np.array([G.loglik(*trame_modele(S, x, 10000, rng, 1.2)) for x in c])
    K = transition(np.log(G.c), F['s'], F['s_saut'], F['w']); Ke = transition(np.log(G.eta), 0.0)
    pc = filtre_lisseur(lambda k: Ls[k], len(c), K, Ke)['lisse_c']; cdf = np.cumsum(pc, 1)
    lo, hi = G.c[np.argmax(cdf >= 0.025, 1)], G.c[np.argmax(cdf >= 0.975, 1)]
    assert np.mean((lo <= c) & (c <= hi)) >= 0.85


@test
def marche_sur_eta_garde_la_variance():
    """Noyau de η plus fin que la grille : après n trames, la variance ajoutée vaut n·s² (et non 0)."""
    le = np.log(np.geomspace(0.6, 1.8, 241)); K = transition(le, 5e-4)
    x = np.zeros(len(le)); x[120] = 1
    for _ in range(2000):
        x = K @ x
    assert np.isclose(x.sum(), 1) and abs(np.sqrt((x * (le - le[120]) ** 2).sum()) / (5e-4 * np.sqrt(2000)) - 1) < 0.01


@test
def lisseur_suit_une_derive_de_eta():
    """c constant, η de 1 à 1,15 sur 400 trames à 30 000 photons : la moyenne de modèles suit η — erreur rms < 1,5 %,
    partout sous 3 sd a posteriori (les bords, où le lisseur n'a de données que d'un côté, compris) — et écarte η constant."""
    from src.exemple import grille_lisseur
    rng = np.random.default_rng(7); T = 400; eta = np.linspace(1.0, 1.15, T)
    gl = grille_lisseur(S, G, F)
    H = [trame_modele(S, 10.0, 30000, rng, e) for e in eta]
    Ls = [gl.loglik(B, Y) for B, Y in H]
    K = transition(np.log(gl.c), F['s'], F['s_saut'], F['w'])
    courses = [filtre_lisseur(lambda k: Ls[k], T, K, transition(np.log(gl.eta), s)) for s in F['s_eta_candidats']]
    _, pe, w = moyenne_de_modeles(courses, causal=False)
    m, sd = resume_eta(pe, gl.eta); e = np.log(m / eta)
    assert np.sqrt(np.mean(e ** 2)) < 0.015 and np.max(np.abs(e) / sd) < 3 and w[0, 0] < 0.01, (e, w[0])


@test
def amincissement_proportionnel_aux_mesures():
    mes = Mesures(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), P['donnees']['dossier']),
                  P['donnees']['titration'], CAL)
    rng = np.random.default_rng(6); B, Y = mes.trame(10.0, 1e7, rng, 1.25); i = mes.indice(10.0)
    assert abs(Y.sum() / B.sum() / (1.25 * mes.Y[i].sum() / mes.B[i].sum()) - 1) < 0.01


if __name__ == '__main__':
    for nom, ok, msg in RESULTATS:
        print(f"{'ok   ' if ok else 'ÉCHEC'} {nom}" + (f"\n{msg}" if msg else ''))
    n_ok = sum(ok for _, ok, _ in RESULTATS)
    print(f"{n_ok}/{len(RESULTATS)} vérifications réussies")
    sys.exit(0 if n_ok == len(RESULTATS) else 1)
