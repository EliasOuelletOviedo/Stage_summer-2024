"""
Les 5 variantes du test de dérive (essai.py : référence, 1, 2, 1-2, 1-2-3), relancées sur les mêmes données (scénario A
de tout.py = essai.py, même graine) en gardant les traces trame par trame → traces/<variante>.npz.
Usage : python3 traces.py ref 1 2 123   (123 donne 1-2 et 1-2-3 : 1-2 est le modèle s_c = 4 % de la moyenne 1-2-3)
"""
import sys, os, time, json, numpy as np
VARS = sys.argv[1:]
T_ESSAI = os.environ.get('T_ESSAI', '3000'); sys.argv = [sys.argv[0], T_ESSAI, 'A']
ICI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ICI)
import tout as t
from moteur import Ratio, filtrer, moyenne_modeles
from src.filtre import resume_c, transition

M, T, NC, g, f = t.M, t.T, t.NC, t.g, t.f
DOS = os.path.join(ICI, 'sorties', 'traces' + ('' if T_ESSAI == '3000' else f'_{T_ESSAI}')); os.makedirs(DOS, exist_ok=True)
lc = np.log(g.c); le_g = np.log(g.eta); ratio = Ratio(g)
NB = t.H[:, :, 0].sum(-1).astype(float); NY = t.H[:, :, 1].sum(-1).astype(float)

commun = os.path.join(DOS, 'commun.npz')
if not os.path.exists(commun):
    c_ratio = np.array([[g.lire_ratio(t.H[m, j, 0].astype(float), t.H[m, j, 1].astype(float)) for j in range(T)] for m in range(M)])
    c_flim = np.array([[g.lire_flim(t.H[m, j, 0].astype(float), t.H[m, j, 1].astype(float)) for j in range(T)] for m in range(M)])
    np.savez_compressed(commun, C=t.C, eta_vrai=t.eta_vrai, c_ratio=c_ratio, c_flim=c_flim)
formes = None


def les_formes():
    global formes
    if formes is None:
        formes = np.empty((M, T, NC))
        for m in range(M):
            for j in range(T):
                formes[m, j] = g.formes(t.H[m, j, 0].astype(float), t.H[m, j, 1].astype(float))
    return formes


def sauver(nom, titre, pcs, mes, ves, poids=None, etiquette_poids=''):
    """pcs : liste de (1 ou M, T, nc) ; mes, ves : (n_η, T) — une trace de η par ROI ou une seule partagée."""
    pc = np.concatenate(pcs, 0)
    moy, lo, hi = (x.reshape(M, T) for x in resume_c(pc.reshape(-1, NC), g.c))
    np.savez_compressed(os.path.join(DOS, f'{nom}.npz'), c=moy, bas=lo, haut=hi, eta=np.exp(np.asarray(mes)),
                        eta_sd=np.sqrt(np.asarray(ves)), titre=titre, poids=np.asarray(poids if poids is not None else []),
                        etiquette_poids=etiquette_poids)
    err = moy - t.C
    print(f"{titre}: rmse {np.sqrt(np.mean(err ** 2)):.4f}, biais {err.mean():+.4f}, couverture "
          f"{np.mean((lo <= t.C) & (t.C <= hi)):.4f}, η final {np.exp(np.asarray(mes)[:, -1]).mean():.4f}", flush=True)


t0 = time.time()
if 'ref' in VARS:
    Kc = transition(lc, 0.04, f['s_saut'], f['w']); pcs, mes, ves, ws = [], [], [], []
    for m in range(M):
        s = [filtrer(les_formes()[m:m + 1], NB[m:m + 1], NY[m:m + 1], ratio, Kc, lambda k: le_g, transition(le_g, se))
             for se in f['s_eta_candidats']]
        pc, me, ve, w = moyenne_modeles(s); pcs.append(pc); mes.append(me); ves.append(ve); ws.append(w[-1])
    sauver('ref', 'référence : η marche aléatoire, chaque ROI seule', pcs, mes, ves, np.mean(ws, 0), 's_η')
    print(f'{time.time() - t0:.0f} s', flush=True)
if '1' in VARS:
    L0, DD = np.meshgrid(np.log(np.geomspace(0.7, 1.5, 81)), np.linspace(-0.15, 0.15, 31), indexing='ij')
    L0, DD = L0.ravel(), DD.ravel(); param = lambda k: L0 + DD * k / (T - 1)
    Kc = transition(lc, 0.04, f['s_saut'], f['w']); pcs, mes, ves = [], [], []
    for m in range(M):
        pc, me, ve, _ = filtrer(les_formes()[m:m + 1], NB[m:m + 1], NY[m:m + 1], ratio, Kc, param)
        pcs.append(pc); mes.append(me); ves.append(ve)
    sauver('1', '1 : dérive paramétrique de η, chaque ROI seule', pcs, mes, ves)
    print(f'{time.time() - t0:.0f} s', flush=True)
if '2' in VARS:
    Kc = transition(lc, 0.04, f['s_saut'], f['w'])
    s = [filtrer(les_formes(), NB, NY, ratio, Kc, lambda k: le_g, transition(le_g, se)) for se in f['s_eta_candidats']]
    pc, me, ve, w = moyenne_modeles(s)
    sauver('2', '2 : η marche aléatoire, partagé par les 10 ROI', [pc], [me], [ve], w[-1], 's_η')
    print(f'{time.time() - t0:.0f} s', flush=True)
if '123' in VARS:                       # le moteur partagé de tout.py, sans aucune modification (identique à moteur.py)
    Kcs = [transition(lc, s, f['s_saut'], f['w']) for s in t.S_C]
    tab0 = t.TablesFiltre.__new__(t.TablesFiltre); tab0.canaux = ('B', 'Y')
    tab0.MB, tab0.MY, tab0.pB, tab0.pY, tab0.lpB, tab0.lpY = g.MB, g.MY, g.pB, g.pY, g.lpB, g.lpY
    tab0.lB, tab0.lY = ratio.lB.astype(np.float32), ratio.lY.astype(np.float32); tab0.bB = tab0.bY = np.zeros(M)
    et = t.Etat(Kcs, tab0); chat = np.zeros((M, T, len(t.QUANTILES_E)), int)
    so = dict(pc=np.empty((M, T, NC), np.float32), m_eta=np.empty(T), v_eta=np.empty(T), w=np.empty((T, 3)))
    pc12 = np.empty((M, T, NC), np.float32); me12 = np.empty(T); ve12 = np.empty(T)
    for k in range(T):
        et.trame(k, chat, so, adaptatif=False)
        mo = et.modeles[2]; lek = et.gr.ln_eta(k)                # s_c = 4 % : le modèle de 1-2
        pc12[:, k] = np.stack([mo.P[m] @ mo.pe.astype(np.float32) for m in range(M)])
        me12[k] = mo.pe @ lek; ve12[k] = mo.pe @ (lek - me12[k]) ** 2
    sauver('12', '1-2 : dérive paramétrique de η, partagée par les 10 ROI', [pc12], [me12], [ve12])
    sauver('123', '1-2-3 : dérive paramétrique partagée + dynamique de c apprise', [so['pc']], [so['m_eta']], [so['v_eta']],
           so['w'][-1], 's_c')
    print(f'{time.time() - t0:.0f} s', flush=True)
