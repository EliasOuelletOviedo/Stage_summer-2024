"""Les 4 variantes (et la référence) sur la même simulation : 10 ROI, 3 000 trames, 1 000 photons/ROI/trame, η 1 → 1,1."""
import sys, os, time, json, numpy as np
sys.dont_write_bytecode = True
ICI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ICI)
from moteur import *
from src.modele import charger_params, charger_calibration, Senseur
from src.filtre import Grille, grille_c, grille_eta, resume_c
from src.simulation import trame_modele
P_ = charger_params(); f = P_['filtre']; S = Senseur(charger_calibration())
g = Grille(S, grille_c(f), grille_eta(dict(f, n_eta=f['n_eta_lisseur']))); ratio = Ratio(g)
T = int(sys.argv[1]) if len(sys.argv) > 1 else 3000; M = 10; N = 1000; scen = sys.argv[2] if len(sys.argv) > 2 else 'derive'
eta_vrai = np.linspace(1.0, 1.1, T) if scen == 'derive' else np.full(T, 1.1)
rng = np.random.default_rng(20261023); k = np.arange(T, dtype=float); C = np.empty((M, T))
for m in range(M):                                       # la forme de l'exemple étiré, propre à chaque ROI
    base, amp = rng.uniform(6, 12), rng.uniform(8, 16); on = int(rng.uniform(0.2, 0.4) * T); off = on + int(rng.uniform(0.2, 0.33) * T)
    c = np.full(T, base); c[on:] = base + amp * (1 - np.exp(-(k[on:] - on) / (T / 75)))
    c[off:] = base + (c[off - 1] - base) * np.exp(-(k[off:] - off) / (T / 12)); C[m] = c
A = 1 - 0.3 * k / (T - 1)
formes = np.empty((M, T, len(g.c))); NB = np.empty((M, T)); NY = np.empty((M, T))
for m in range(M):
    for j in range(T):
        Bh, Yh = trame_modele(S, C[m, j], N * A[j], rng, eta_vrai[j])
        formes[m, j] = g.formes(Bh, Yh); NB[m, j] = Bh.sum(); NY[m, j] = Yh.sum()
lc = np.log(g.c)
Kc = {s: transition(lc, s, f['s_saut'], f['w']) for s in (0.004, 0.013, 0.04)}
# η marche aléatoire : grille de 241, 5 vitesses candidates
le_g = np.log(g.eta); Ke = {s: transition(le_g, s) for s in f['s_eta_candidats']}
rw = lambda kk: le_g
# η paramétrique : η0 81 valeurs de 0,7 à 1,5 (pas de 0,95 %) ; dérive totale D de −15 % à +15 % (ln), pas de 1 %
le0 = np.log(np.geomspace(0.7, 1.5, 81)); Dg = np.linspace(-0.15, 0.15, 31)
L0, DD = np.meshgrid(le0, Dg, indexing='ij'); L0, DD = L0.ravel(), DD.ravel()
param = lambda kk: L0 + DD * kk / (T - 1)
def resultat(nom, pc, m_eta, v_eta, ROIs):
    moy, lo, hi = resume_c(pc.reshape(-1, len(g.c)), g.c)
    moy, lo, hi = (x.reshape(len(ROIs), T) for x in (moy, lo, hi)); Cv = C[ROIs]
    err = moy - Cv; e_eta = m_eta - np.log(eta_vrai)
    sl = slice(T // 3, None)
    r = dict(variante=nom, rmse_c=float(np.sqrt(np.mean(err ** 2))), rmse_c_apres_1tiers=float(np.sqrt(np.mean(err[:, sl] ** 2))),
             biais_c=float(err.mean()), couverture=float(np.mean((lo <= Cv) & (Cv <= hi))),
             eta_rms_pct=float(100 * np.sqrt(np.mean(e_eta ** 2))), eta_rms_apres_1tiers_pct=float(100 * np.sqrt(np.mean(e_eta[sl] ** 2))),
             eta_1000=float(np.exp(m_eta[min(999, T - 1)])), eta_2000=float(np.exp(m_eta[min(1999, T - 1)])), eta_final=float(np.exp(m_eta[-1])),
             eta_final_sd_pct=float(100 * np.sqrt(v_eta[-1])))
    return r
def par_roi(nom, fn):
    """Une ROI à la fois (η propre à chaque ROI) : résultats moyennés sur les 10 ROI."""
    rs = []
    for m in range(M):
        pc, me, ve = fn(m); rs.append(resultat(nom, pc, me, ve, [m]))
    out = {k: (np.mean([r[k] for r in rs]) if k != 'variante' else nom) for k in rs[0]}
    return out
res = []; t0 = time.time()
def ref(m):
    s = [filtrer(formes[m:m + 1], NB[m:m + 1], NY[m:m + 1], ratio, Kc[0.04], rw, Ke[se]) for se in f['s_eta_candidats']]
    pc, me, ve, _ = moyenne_modeles(s); return pc, me, ve
res.append(par_roi('référence : η marche aléatoire, 1 ROI', ref)); print(res[-1], f'{time.time() - t0:.0f} s', flush=True)
def v1(m):
    pc, me, ve, _ = filtrer(formes[m:m + 1], NB[m:m + 1], NY[m:m + 1], ratio, Kc[0.04], param); return pc, me, ve
res.append(par_roi('1 : dérive paramétrique, 1 ROI', v1)); print(res[-1], f'{time.time() - t0:.0f} s', flush=True)
s = [filtrer(formes, NB, NY, ratio, Kc[0.04], rw, Ke[se]) for se in f['s_eta_candidats']]
pc, me, ve, w = moyenne_modeles(s); res.append(resultat('2 : η marche aléatoire partagé, 10 ROI', pc, me, ve, list(range(M))))
print(res[-1], 'poids finaux s_η', np.round(w[-1], 2), f'{time.time() - t0:.0f} s', flush=True)
pc, me, ve, _ = filtrer(formes, NB, NY, ratio, Kc[0.04], param); res.append(resultat('1-2 : dérive paramétrique partagée, 10 ROI', pc, me, ve, list(range(M))))
print(res[-1], f'{time.time() - t0:.0f} s', flush=True)
s = [filtrer(formes, NB, NY, ratio, Kc[sc], param) for sc in (0.004, 0.013, 0.04)]
pc, me, ve, w = moyenne_modeles(s); res.append(resultat('1-2-3 : + dynamique de c apprise', pc, me, ve, list(range(M))))
print(res[-1], 'poids finaux s_c (0,4 / 1,3 / 4 %)', np.round(w[-1], 2), f'{time.time() - t0:.0f} s', flush=True)
os.makedirs(os.path.join(ICI, 'sorties'), exist_ok=True)
json.dump(res, open(os.path.join(ICI, 'sorties', f'resultats_essai_{T}.json'), 'w'), indent=1)
