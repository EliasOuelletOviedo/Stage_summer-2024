"""Étape 12 — banc d'essai : quelle lecture de la [Cl⁻] convient le mieux à ces données en tranche ?

Usage :  python3 step12_banc_essai.py <dossier> [--processus 6]      (~15 min)
         python3 step12_banc_essai.py <dossier> --sans-calcul         refait notes et figures depuis le cache
         python3 step12_banc_essai.py <dossier> --essai               10 somas par acquisition, pour vérifier le code

Regard neuf sur les déclins de l'étape 4 : une quarantaine de lectures (banc.py : rapports, formes du déclin sans
ajustement, ajustements, gabarits mesurés, apprentissage, modèle de l'étape 9) sont calculées sur les mêmes photons,
à plusieurs niveaux : chaque soma (cycles gardés réunis), ses deux moitiés de cycles (pairs / impairs), chaque
cycle de 10 s, la somme amincie (tirage binomial) à un nombre de photons donné, les cycles des time-lapses, et les
photons réunis de chaque champ.

Critères, fixés avant de regarder les résultats. « mM-éq. » : écart de la lecture divisé par sa pente moyenne entre
les médianes des somas à 0 et à 20 mM (la gamme utile).
  A. séparation soma par soma : aire sous la courbe ROC, 0 contre 20 mM et 20 contre 131 mM ;
  B. écart entre champs (écart-type des médianes de champ, 0 et 20 mM), dispersion entre somas d'un champ, et
     dispersion « vraie » : covariance des deux moitiés de cycles, donc sans bruit de photons ;
  C. répétabilité : bruit d'un cycle de 10 s (différences successives), et son rapport au bruit de photons attendu
     (somme amincie au nombre de photons d'un cycle) ;
  D. précision par photon : écart-type et biais à 300, 1 000, 3 000 et 10 000 photons (CFP + YFP) ;
  E. exactitude en validation croisée (isotherme par les médianes d'entraînement), plis de l'étape 10 : un champ de
     chaque condition retiré, ou somas tirés au hasard ; erreur quadratique des médianes lues ; [Cl⁻]i basale ;
  F. changement dans une même cellule : time-lapse 20→131 (4 derniers cycles contre 4 premiers), z = Δ / erreur
     type ; aussi basal→0 (fin contre 4–9 min) et 0→20 ;
  G. dérive pendant la minute d'illumination d'un champ (pente par minute) ;
  K. lecture d'un champ entier (photons réunis), validation par champ.
Classement : moyenne des rangs sur huit critères — erreur par champ à 0 et à 20 mM (E), écart entre champs (B),
AUC 0/20 (A), écart-type à 1 000 photons (D), bruit d'un cycle (C), z de 20→131 (F), |dérive| (G).

Les lectures apprises (lda, ridge, arbres) sont toujours lues hors de leurs champs d'entraînement ; les gabarits et
l'axe du phaseur, références de toutes les données sauf pour E et K. Les [Cl⁻] négatives sont gardées.

Compléments, depuis le cache : les mêmes notes sur la moitié la plus brillante des somas de chaque champ ; écarts
des champs lus par le rapport et par la forme CFP (commun ou discordant) ; time-lapses cellule par cellule. Le
classement « rang_6 » laisse de côté les deux critères qui supposent que chaque champ est à sa [Cl⁻] nominale
(erreur par champ à 20 mM, écart entre champs) : la comparaison rapport / forme montre que les champs à 20 mM
diffèrent réellement.

Sorties : results/banc_lectures.csv (une ligne par lecture), banc_somas.csv, banc_precision.csv, banc_champs.csv,
banc_correlations.csv, banc_lectures_brillants.csv, banc_champs_rapport_forme.csv, banc_time_lapses.csv,
banc_classement.png, banc_precision.png, banc_deux_etats.png, banc_champs.png
"""
import argparse
import json
import pickle
import re
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.stats import rankdata, spearmanr

import banc as BA
import config as C
import deux_canaux as DC
import figures as G
import flim
import recalage
from figures import plt
from mesures import comptes_cycles, selection
from step7_calibration import ajuster, inverser

W, FEN, CONDITIONS = BA.W, BA.FEN, BA.CONDITIONS
CONC = [0.0, 20.0, 131.0]
BUDGETS = (300, 1000, 3000, 10000)
N_AMINCIS, N_REP = 400, 6
GRAINE = 20261110                 # celle de l'étape 10 : mêmes plis aléatoires
ALPHA_NP, EROSION_PX = 0.7, 3
ANNEXES = ("R_tm", "R_np", "R_erode", "R_brillant")
CLASSEMENT = [("rmse_champ_0mM", 1), ("rmse_champ_20mM", 1), ("champs_mM", 1), ("auc_0_20", -1),
              ("sigma_1000", 1), ("cycle_mM", 1), ("z_20_to_131", -1), ("derive_abs", 1)]
FAMILLE = {"combien": ("R", "R_tard", "R_eta") + ANNEXES,
           "quand, CFP": ("tm_B", "tm_B_sf", "tm_queue_B", "porte_B", "tau_rld_B", "tphi_B", "tmod_B", "tphi2_B", "tau_mono",
                          "tau_mono_s", "tau_bi", "tau_bi_amp", "tau_3pop", "f_long", "phas_axe_B", "gab_B"),
           "quand, YFP": ("tm_Y", "tphi_Y", "tau_mono_Y", "gab_Y"),
           "deux canaux": ("gab_formes", "gab_formes_fond", "gab_joint", "lda_portes", "lda_3", "ridge_portes",
                           "gbr_portes", "pc1", "mod2_flim", "mod1_flim", "mod2_hyb", "mod1_hyb")}


def num_fov(nom):
    m = re.search(r"fov(\d+)$", nom)
    return int(m.group(1)) if m else 0


def dit(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ données
def eroder(masks, r):
    """Chaque ROI rogné de r pixels (un pixel reste si tout son voisinage porte la même étiquette)."""
    t = 2 * r + 1
    return np.where((masks > 0) & (ndi.minimum_filter(masks, size=t) == masks)
                    & (ndi.maximum_filter(masks, size=t) == masks), masks, 0)


def brillants(masks, som):
    """La moitié la plus brillante des pixels de chaque ROI (somme des deux canaux, cycles gardés)."""
    out = np.zeros_like(masks)
    for lab, sl in enumerate(ndi.find_objects(masks), 1):
        if sl is None:
            continue
        m = masks[sl] == lab
        out[sl][m & (som[sl] >= np.median(som[sl][m]))] = lab
    return out


def comptes_images(ims, lab, n, decalages):
    """Photons par ROI, par cycle et par canal (cycles, n + 1, 2), masques déplacés comme à l'étape 4."""
    out = np.zeros((len(ims), n + 1, 2))
    for c in range(len(ims)):
        lk = recalage.deplacer(lab, decalages[c]) if tuple(decalages[c]) != (0, 0) else lab
        lk = lk.ravel().astype(np.int64)
        for ch in range(2):
            out[c, :, ch] = np.bincount(lk, weights=ims[c, ch].ravel().astype(float), minlength=n + 1)
    return out


def annexes(D, a, m, n, g, k):
    """Rapports qui demandent les images : numérateur et dénominateur (cycles gardés, somas) de chaque variante."""
    z = np.load(D["projections"] / f"{a['nom']}.npz")
    lab, dk = m["etiquettes"], m["decalages"]
    cc = comptes_cycles(z, lab, n + 1, dk)[g]                      # (cy, n + 2, canal, [n, corrigé, n²/E])
    out = {"R_tm": (cc[:, k, 1, 1], cc[:, k, 0, 1])}
    aire = np.bincount(lab.ravel(), minlength=n + 2).astype(float)
    np_px = cc[:, n + 1, :, 0] / aire[n + 1]                        # neuropile : photons par pixel (cy, canal)
    out["R_np"] = tuple(cc[:, k, ch, 0] - ALPHA_NP * aire[k][None, :] * np_px[:, ch][:, None] for ch in (1, 0))
    ims = z["images"]
    masks = np.where((lab >= 1) & (lab <= n), lab, 0)
    for nom, l in (("R_erode", eroder(masks, EROSION_PX)), ("R_brillant", brillants(masks, z["somme"].sum(0)))):
        c2 = comptes_images(ims, l, n, dk)[g]
        out[nom] = (c2[:, k, 1], c2[:, k, 0])
    return out


def assembler(root, D, essai):
    N = {q: dict(B=[], Y=[], meta=[], ann=[]) for q in ("soma", "demi", "cycle", "trans", "champ")}

    def ajouter(niv, d, meta, ann=None):
        N[niv]["B"].append(d[:, 0])
        N[niv]["Y"].append(d[:, 1])
        N[niv]["meta"].append(meta)
        if ann is not None:
            N[niv]["ann"].append(ann)

    i_soma = 0
    for a in C.acquisitions(root):
        f = D["mesures"] / f"{a['nom']}.npz"
        if not f.exists():
            continue
        m = np.load(f)
        n, g = int(m["n_rois"]), np.nonzero(m["garde"])[0]
        k = selection(D, a["nom"], n)
        if essai:
            k = k[:10]
        dec = m["declins"][g][:, k].astype(np.float32)             # (cycles gardés, somas, canal, 256)
        t_s = m["debut"][g] - m["debut"][0] + m["duree_s"][g] / 2
        ncy, ns = dec.shape[:2]
        base = dict(acq=a["nom"], cond=a["condition"], fov=num_fov(a["nom"]))
        ajouter("champ", dec.sum((0, 1))[None], pd.DataFrame([dict(base, type=a["type"], somas=ns)]))
        if a["type"] == "fov":
            ann = annexes(D, a, m, n, g, k)
            with np.errstate(divide="ignore", invalid="ignore"):
                rap = lambda sel: {q: v[0][sel].sum(0) / v[1][sel].sum(0) for q, v in ann.items()}
                soma = np.arange(i_soma, i_soma + ns)
                i_soma += ns
                ajouter("soma", dec.sum(0), pd.DataFrame(dict(base, roi=k, soma=soma, n_cycles=ncy)), rap(slice(None)))
                for h in (0, 1):
                    ajouter("demi", dec[h::2].sum(0), pd.DataFrame(dict(base, roi=k, soma=soma, moitie=h)),
                            rap(slice(h, None, 2)))
                for c in range(ncy):
                    ajouter("cycle", dec[c], pd.DataFrame(dict(base, roi=k, soma=soma, cycle=c, t_s=t_s[c])),
                            rap(slice(c, c + 1)))
        else:
            t_min = t_s / 60
            if a["nom"].startswith("basal"):           # fin contre le plateau du début (4 à 9 min)
                cols = [(c, "debut", 0) for c in np.nonzero((t_min >= 4) & (t_min <= 9))[0]] \
                    + [(c, "fin", 1) for c in range(ncy - 6, ncy)]
            else:
                cols = [(c, "debut" if c < 4 else "fin" if c >= ncy - 4 else "milieu", 0) for c in range(ncy)]
            for j, (c, bloc, serie) in enumerate(cols):
                ajouter("trans", dec[c], pd.DataFrame(dict(base, roi=k, roi_i=np.arange(ns), col=j, bloc=bloc,
                                                           serie=serie, t_min=t_min[c])))
    for q, v in N.items():
        v["B"], v["Y"] = np.concatenate(v["B"]), np.concatenate(v["Y"])
        v["meta"] = pd.concat(v["meta"], ignore_index=True)
        v["ann"] = {x: np.concatenate([d[x] for d in v["ann"]]) for x in ANNEXES} if v["ann"] else {}
    return N


def amincir(N, rng):
    """Somme des cycles de somas à 0 et à 20 mM, amincie à un nombre de photons donné (et à celui d'un cycle)."""
    S, Bs, Ys = N["soma"]["meta"], N["soma"]["B"], N["soma"]["Y"]
    tot = Bs[:, W].sum(1) + Ys[:, W].sum(1)
    pris = np.concatenate([rng.choice(v, min(len(v), N_AMINCIS // 2), replace=False) for v in
                           (np.nonzero((S.cond == c).values & (tot >= 2 * max(BUDGETS)))[0] for c in ("0 mM", "20 mM"))])
    B, Y, rows = [], [], []
    for j, i in enumerate(pris):
        for b, nb in enumerate(list(BUDGETS) + [None]):
            p = 1 / S.n_cycles[i] if nb is None else nb / tot[i]
            B.append(rng.binomial(Bs[i].astype(np.int64), p, size=(N_REP, 256)))
            Y.append(rng.binomial(Ys[i].astype(np.int64), p, size=(N_REP, 256)))
            rows += [dict(acq=S.acq[i], cond=S.cond[i], fov=S.fov[i], soma=i, j=j, b=b, rep=r, p=p) for r in range(N_REP)]
    N["aminci"] = dict(B=np.concatenate(B).astype(np.float32), Y=np.concatenate(Y).astype(np.float32),
                       meta=pd.DataFrame(rows), ann={})


def instrument(root, D, N, essai):
    f = D["qc"] / ("banc_instrument_essai.json" if essai else "banc_instrument.json")
    if f.exists():
        d = json.loads(f.read_text())
    else:
        MB = flim.Modele(C.charger_irf("cfp", root), fenetre=FEN)
        MY = flim.Modele(C.charger_irf("yfp", root), fenetre=FEN)
        sB = [MB.ajuster_multi(b.astype(float), 3)["s"] for b in N["champ"]["B"]]
        sY = [MY.ajuster_multi(y.astype(float), 3)["s"] for y in N["champ"]["Y"]]
        d = dict(sB=float(np.median(sB)), sY=float(np.median(sY)), sB_ecart_type=float(np.std(sB)),
                 sY_ecart_type=float(np.std(sY)), pic=int(np.argmax(N["champ"]["B"].sum(0)[W])),
                 note="décalage de l'IRF en canaux (3 exponentielles sur les photons réunis de chaque acquisition, "
                      "médiane) ; pic du déclin CFP réuni, en canaux de la fenêtre")
        f.write_text(json.dumps(d, indent=1, ensure_ascii=False))
    taus3 = json.loads((D["qc"] / "irf.json").read_text())["taus_cfp"]
    dit(f"  IRF : décalage CFP {d['sB'] * C.DT_NS * 1e3:.1f} ± {d['sB_ecart_type'] * C.DT_NS * 1e3:.1f} ps, "
        f"YFP {d['sY'] * C.DT_NS * 1e3:.1f} ± {d['sY_ecart_type'] * C.DT_NS * 1e3:.1f} ps ; pic au canal {d['pic']}")
    return BA.Instrument(root, d["sB"], d["sY"], d["pic"], taus3)


def gabarits_de(S, Bs, Ys, masque):
    """Gabarits tirés des somas du masque : photons réunis de chaque champ à 0 et à 131 mM."""
    args = []
    for c in ("0 mM", "131 mM"):
        mc = masque & (S.cond == c).values
        ch = [mc & (S.acq == a).values for a in pd.unique(S.acq[mc])]
        args += [[Bs[m].sum(0, dtype=float) for m in ch], [Ys[m].sum(0, dtype=float) for m in ch]]
    return BA.Gabarits(*args)


def appris_de(S, Bs, Ys, masque, ins):
    return BA.Appris(Bs[masque], Ys[masque], S.cond.values[masque], S.acq.values[masque], ins)


def plis_aleatoires(S):
    """Les plis de l'étape 10 : dans chaque condition imposée, somas mélangés (même graine) puis répartis en 5."""
    rng = np.random.default_rng(GRAINE)
    pli = np.zeros(len(S), int)
    for c in CONDITIONS:
        idx = np.nonzero((S.cond == c).values)[0]
        for i, j in enumerate(rng.permutation(len(idx))):
            pli[idx[j]] = i % 5 + 1
    return pli


# ------------------------------------------------------------------ calcul des lectures
def calculer(B, Y, ins, gab, app, pool):
    o = BA.simples(B, Y, ins)
    o.update(BA.trois_pop(B, ins))
    o.update(gab.lire(B, Y))
    o.update(app.lire(B, Y, o))
    lents = {}
    for r in pool.imap(BA.lent, [(B[i:i + 250], Y[i:i + 250]) for i in range(0, len(B), 250)]):
        for q, v in r.items():
            lents.setdefault(q, []).append(v)
    o.update({q: np.concatenate(v) for q, v in lents.items()})
    return o


def eta_des_champs(gab, vch, Ch):
    """η de chaque acquisition : R des photons réunis / R prédit par leur forme (gabarits)."""
    f = vch["gab_formes"]
    return pd.Series(vch["R"] / (gab.R0 * (1 - f) + gab.R1 * f), index=Ch.acq.values)


def tout_calculer(root, D, args):
    t0 = time.time()
    N = assembler(root, D, args.essai)
    amincir(N, np.random.default_rng(GRAINE + 1))
    S, Bs, Ys = N["soma"]["meta"], N["soma"]["B"], N["soma"]["Y"]
    Ch = N["champ"]["meta"]
    dit("  niveaux : " + ", ".join(f"{q} {len(v['B'])}" for q, v in N.items()) + f" ({time.time() - t0:.0f} s)")
    ins = instrument(root, D, N, args.essai)
    imposes = S.cond.isin(CONDITIONS).values
    gab, app = gabarits_de(S, Bs, Ys, imposes), appris_de(S, Bs, Ys, imposes, ins)
    dit(f"  gabarits : R0 {gab.R0:.3f}, R1 {gab.R1:.3f}, f de {gab.f[0]:.2f} à {gab.f[-1]:.2f}")
    Bc = [Bs[(S.cond == c).values].sum(0, dtype=float) for c in CONDITIONS]
    Yc = [Ys[(S.cond == c).values].sum(0, dtype=float) for c in CONDITIONS]
    irfp = DC.ajuster_irf(DC.donnees(root, Bc, Yc, CONC, FEN))
    lecteurs = {f"mod{nt}_flim": BA.lecteur_modele(root, Bc, Yc, CONC, irfp, nt) for nt in (2, 1)}
    dit(f"  modèles de l'étape 9 calibrés ({time.time() - t0:.0f} s)")

    V = {}
    with Pool(args.processus, initializer=BA.init_lent, initargs=(ins, gab, lecteurs)) as pool:
        for q in N:
            V[q] = calculer(N[q]["B"], N[q]["Y"], ins, gab, app, pool)
            V[q].update(N[q]["ann"])
            dit(f"  lectures du niveau {q} ({time.time() - t0:.0f} s)")
    eta = eta_des_champs(gab, V["champ"], Ch)
    for q in N:
        V[q]["R_eta"] = V[q]["R"] / eta.reindex(N[q]["meta"].acq.values).values

    # lectures apprises : chaque champ imposé est lu par un apprentissage qui ne l'a pas vu
    plis = {"fov": [imposes & (S.fov.values == k) for k in range(1, 6)],
            "alea": [plis_aleatoires(S) == k for k in range(1, 6)]}
    P = {"fov": [], "alea": []}
    for mode in ("fov", "alea"):
        for k, te in enumerate(plis[mode], 1):
            tr = imposes & ~te
            gk, ak = gabarits_de(S, Bs, Ys, tr), appris_de(S, Bs, Ys, tr, ins)
            vs, vc = {}, {}
            for v, niv in ((vs, "soma"), (vc, "champ")):
                B, Y = N[niv]["B"], N[niv]["Y"]
                v.update(gk.lire(B, Y))
                v["gab_formes_fond"] = gk.lire_fond(B, Y)
                v.update(ak.lire(B, Y, V[niv]))
            ek = eta_des_champs(gk, dict(vc, R=V["champ"]["R"]), Ch)
            vs["R_eta"] = V["soma"]["R"] / ek.reindex(S.acq.values).values
            vc["R_eta"] = V["champ"]["R"] / ek.values
            P[mode].append(dict(soma={q: x for q, x in vs.items() if not q.startswith("_")},
                                champ={q: x for q, x in vc.items() if not q.startswith("_")}))
            if mode == "fov":
                for q in N:
                    lig = (N[q]["meta"].cond.isin(CONDITIONS) & (N[q]["meta"].fov == k)).values
                    if lig.any():
                        r = ak.lire(N[q]["B"][lig], N[q]["Y"][lig], {x: V[q][x][lig] for x in ("R", "tm_B", "tm_Y", "_gB", "_sB")})
                        for x in BA.APPRISES:
                            V[q][x][lig] = r[x]
        dit(f"  plis « {mode} » ({time.time() - t0:.0f} s)")
    cache = dict(meta={q: v["meta"] for q, v in N.items()}, V=V, P=P, plis=plis,
                 gab=dict(R0=gab.R0, R1=gab.R1, rho=gab.rho, p0B=gab.p0B, p1B=gab.p1B, p0Y=gab.p0Y, p1Y=gab.p1Y),
                 p20=[np.mean([N["champ"]["B" if ch == 0 else "Y"][i][W] / (N["champ"]["B"][i][W].sum() + N["champ"]["Y"][i][W].sum())
                               for i in np.nonzero((Ch.cond == "20 mM").values & (Ch.type == "fov").values)[0]], 0)
                      for ch in (0, 1)])
    with open(D["mesures"] / ("banc_cache_essai.pkl" if args.essai else "banc_cache.pkl"), "wb") as f:
        pickle.dump(cache, f)
    return cache


# ------------------------------------------------------------------ notes
def lecture_mM(x, med):
    """[Cl⁻] par l'isotherme qui passe par les trois médianes ; à défaut (courbe inexistante), interpolation
    linéaire par morceaux, prolongée. Rend (c, isotherme utilisée ?)."""
    x = np.asarray(x, float)
    med = np.asarray(med, float)
    f = ajuster(CONC, med) if np.all(np.isfinite(med)) else None
    if f is not None:
        c = np.asarray(inverser(x, *f), float)
        c[np.isnan(x)] = np.nan
        return c, True
    if not np.all(np.isfinite(med)) or med[1] == med[0]:
        return np.full(len(x), np.nan), False
    if not (np.all(np.diff(med) > 0) or np.all(np.diff(med) < 0)):
        # médianes non monotones : droite par 0 et 20 mM seulement (la gamme utile) ; rien de fiable au-delà
        return np.where(np.isnan(x), np.nan, CONC[0] + (x - med[0]) * (CONC[1] - CONC[0]) / (med[1] - med[0])), False
    s = 1.0 if med[1] > med[0] else -1.0
    xm, cm, u = s * med, np.asarray(CONC), s * x
    c = np.interp(u, xm, cm)
    c = np.where(u < xm[0], cm[0] + (u - xm[0]) * (cm[1] - cm[0]) / (xm[1] - xm[0]), c)
    c = np.where(u > xm[-1], cm[-1] + (u - xm[-1]) * (cm[-1] - cm[-2]) / (xm[-1] - xm[-2]), c)
    return np.where(np.isnan(x), np.nan, c), False


def valider(xs, tests, cond):
    """Validation croisée : xs = valeur de chaque soma sous les références de chaque pli, tests = somas lus."""
    imp = np.isin(cond, CONDITIONS)
    err, basal, iso = {c: [] for c in CONDITIONS}, [], 0
    for x, te in zip(xs, tests):
        tr = imp & ~te
        med = [np.nanmedian(x[tr & (cond == c)]) for c in CONDITIONS]
        cm, ok = lecture_mM(x, med)
        iso += ok
        for c, cl in zip(CONDITIONS, CONC):
            v = cm[te & (cond == c)]
            v = v[~np.isnan(v)]
            err[c].append(np.median(v) - cl if len(v) else np.nan)
        v = cm[cond == "basal"]
        v = v[~np.isnan(v)]
        basal.append(np.median(v) if len(v) else np.nan)
    out = {}
    for c in CONDITIONS:
        e = np.array(err[c], float)
        fin = np.all(np.isfinite(e))
        out[c.replace(" ", "")] = (float(np.sqrt(np.mean(e ** 2))) if fin else np.inf, float(np.mean(e)) if fin else np.inf)
    return out, (float(np.nanmedian(basal)), float(np.nanmin(basal)), float(np.nanmax(basal))), iso


def valider_champs(xs, Ch):
    """Comme valider, sur les photons réunis de chaque champ : le champ k de chaque condition lu hors calibration."""
    cond, fov = Ch.cond.values, Ch.fov.values
    fo = (Ch.type == "fov").values & np.isin(cond, CONDITIONS)
    out = {}
    err = {c: [] for c in CONDITIONS}
    for k, x in enumerate(xs, 1):
        med = [np.nanmedian(x[fo & (fov != k) & (cond == c)]) for c in CONDITIONS]
        cm, _ = lecture_mM(x, med)
        for c, cl in zip(CONDITIONS, CONC):
            v = cm[fo & (fov == k) & (cond == c)]
            err[c].append(v[0] - cl if len(v) else np.nan)
    for c in CONDITIONS:
        e = np.array(err[c], float)
        out[c.replace(" ", "")] = float(np.sqrt(np.mean(e ** 2))) if np.all(np.isfinite(e)) else np.inf
    return out


def auc(a, b, sens):
    """Probabilité qu'un soma de b soit lu « plus de Cl⁻ » qu'un soma de a."""
    a, b = sens * a[~np.isnan(a)], sens * b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    r = rankdata(np.r_[a, b])
    return float((r[len(a):].sum() - len(b) * (len(b) + 1) / 2) / (len(a) * len(b)))


def fini(x):
    return np.where(np.isfinite(x), x, np.nan)


def mad(v):
    v = v[np.isfinite(v)]
    return 1.4826 * np.median(np.abs(v - np.median(v))) if len(v) else np.nan


def commun(meta):
    """Indices partagés par toutes les lectures : masques, matrices soma × cycle, soma × moitié, amincis, time-lapses."""
    S, Cy, Dm, Th, Tr = meta["soma"], meta["cycle"], meta["demi"], meta["aminci"], meta["trans"]
    K = dict(cond=S.cond.values, acq=S.acq.values, fov=S.fov.values)
    K["m"] = {c: (S.cond == c).values for c in ("basal",) + CONDITIONS}
    K["bas"] = K["m"]["0 mM"] | K["m"]["20 mM"]
    K["champs_bas"] = [(S.acq == a).values for a in pd.unique(S.acq[K["bas"]])]
    K["champs"] = {a: (S.acq == a).values for a in pd.unique(S.acq)}
    ic = -np.ones((len(S), int(Cy.cycle.max()) + 1), int)
    ic[Cy.soma.values, Cy.cycle.values] = np.arange(len(Cy))
    tc = np.full(ic.shape, np.nan)
    tc[Cy.soma.values, Cy.cycle.values] = Cy.t_s.values
    idm = -np.ones((len(S), 2), int)
    idm[Dm.soma.values, Dm.moitie.values] = np.arange(len(Dm))
    ith = -np.ones((int(Th.j.max()) + 1, len(BUDGETS) + 1, N_REP), int)
    ith[Th.j.values, Th.b.values, Th.rep.values] = np.arange(len(Th))
    pth = np.zeros(ith.shape[:2])
    pth[Th.j.values, Th.b.values] = Th.p.values
    som = np.zeros(ith.shape[0], int)
    som[Th.j.values] = Th.soma.values
    K.update(ic=ic, tc=tc, idm=idm, ith=ith, pth=pth, som_th=som)
    K["trans"] = {}
    for nom, g in Tr.groupby("acq"):
        it = -np.ones((int(g.roi_i.max()) + 1, int(g.col.max()) + 1), int)
        it[g.roi_i.values, g.col.values] = g.index.values
        cols = g.drop_duplicates("col").sort_values("col")
        K["trans"][nom] = dict(idx=it, bloc=cols.bloc.values, serie=cols.serie.values)
    return K


def mat(x, idx):
    return np.where(idx >= 0, x[np.maximum(idx, 0)], np.nan)


def cov_moities(X1, X2, K):
    """Covariance entre deux lectures d'un soma, sans bruit de photons : l'une sur les cycles pairs, l'autre sur les
    cycles impairs (et l'inverse), centrées par champ, somas à plus de 4 σ robustes de leur champ écartés ; champs à 0
    et à 20 mM réunis. X = (valeurs des somas, des moitiés)."""
    (s1, d1), (s2, d2) = X1, X2
    a1, b1, a2, b2 = (fini(mat(d, K["idm"][:, h])) for d in (d1, d2) for h in (0, 1))
    num = n = 0.0
    for m in K["champs_bas"]:
        ok = m & np.isfinite(a1) & np.isfinite(b1) & np.isfinite(a2) & np.isfinite(b2) & np.isfinite(s1) & np.isfinite(s2)
        for s in (s1, s2):
            if ok.sum() >= 8:
                ok = ok & (np.abs(s - np.median(s[ok])) <= 4 * max(mad(s[ok]), 1e-300))
        if ok.sum() < 8:
            continue
        c = lambda u, v: np.mean((u[ok] - u[ok].mean()) * (v[ok] - v[ok].mean()))
        num += 0.5 * (c(a1, b2) + c(b1, a2)) * ok.sum()
        n += ok.sum()
    return num / n if n else np.nan


def noter(X, K):
    """Toutes les notes d'une lecture. X : {niveau: valeurs ou None}."""
    r = {}
    xs = X["soma"]
    med = {c: float(np.nanmedian(xs[m])) for c, m in K["m"].items()}
    sens = 1.0 if med["131 mM"] >= med["0 mM"] else -1.0
    pente = (med["20 mM"] - med["0 mM"]) / 20.0
    ap = abs(pente) if pente else np.nan
    r.update(med_basal=med["basal"], med_0=med["0 mM"], med_20=med["20 mM"], med_131=med["131 mM"],
             monotone=bool((med["20 mM"] - med["0 mM"]) * sens > 0 and (med["131 mM"] - med["20 mM"]) * sens > 0),
             echec_pct=100 * float(np.mean(~np.isfinite(xs))))
    r["auc_0_20"] = auc(xs[K["m"]["0 mM"]], xs[K["m"]["20 mM"]], sens)
    r["auc_20_131"] = auc(xs[K["m"]["20 mM"]], xs[K["m"]["131 mM"]], sens)
    # B. entre champs, entre somas
    ec = [np.nanmedian(fini(xs[m])) for m in K["champs_bas"]]
    cb = [K["cond"][m][0] for m in K["champs_bas"]]
    res = np.concatenate([np.array([e for e, c in zip(ec, cb) if c == cc]) - np.mean([e for e, c in zip(ec, cb) if c == cc])
                          for cc in ("0 mM", "20 mM")])
    r["champs_mM"] = float(np.sqrt(np.sum(res ** 2) / (len(res) - 2)) / ap)
    r["somas_mM"] = float(np.nanmedian([mad(xs[m]) for m in K["champs_bas"]]) / ap)
    rho = []
    for m in K["champs_bas"]:
        ok = m & np.isfinite(xs)
        if ok.sum() >= 8:
            rho.append(spearmanr(sens * xs[ok], K["lnph"][ok])[0])
    r["rho_photons"] = float(np.median(rho)) if rho else np.nan
    if X.get("demi") is not None:
        r["vrai_somas_mM"] = float(np.sqrt(max(cov_moities((xs, X["demi"]), (xs, X["demi"]), K), 0)) / ap)
    # C. répétabilité, G. dérive
    if X.get("cycle") is not None:
        xc = fini(mat(X["cycle"], K["ic"]))
        with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sig = np.sqrt(np.nanmean(np.diff(xc, axis=1) ** 2, axis=1) / 2)
            r["cycle_mM"] = float(np.nanmedian(sig[K["bas"]]) / ap)
            t = K["tc"] - np.nanmean(np.where(np.isfinite(xc), K["tc"], np.nan), 1, keepdims=True)
            xm = xc - np.nanmean(xc, 1, keepdims=True)
            slope = np.nansum(t * xm, 1) / np.nansum(np.where(np.isfinite(xc), t, np.nan) ** 2, 1)
            slope[np.isfinite(xc).sum(1) < 4] = np.nan
            r["derive_mM_min"] = float(np.nanmedian(slope[K["bas"]]) * 60 / pente) if pente else np.nan
            r["derive_abs"] = abs(r["derive_mM_min"])
    # D. précision par photon
    if X.get("aminci") is not None:
        xt = fini(mat(X["aminci"], K["ith"]))                       # (somas tirés, budgets + cycle, répétitions)
        with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sd = np.nanstd(xt, axis=2, ddof=1) / np.sqrt(1 - K["pth"])
            for b, nb in enumerate(BUDGETS):
                r[f"sigma_{nb}"] = float(np.nanmedian(sd[:, b]) / ap)
                r[f"biais_{nb}"] = float(np.nanmedian(np.nanmean(xt[:, b], 1) - xs[K["som_th"]]) / pente) if pente else np.nan
                r[f"echec_{nb}"] = 100 * float(np.mean(~np.isfinite(mat(X["aminci"], K["ith"][:, b]))))
            if X.get("cycle") is not None:
                r["exces_bruit"] = float(np.nanmedian(sig[K["som_th"]]) / np.nanmedian(sd[:, -1]))
    # F. changement dans une même cellule
    if X.get("trans") is not None:
        for nom, tr in K["trans"].items():
            xt = fini(mat(X["trans"], tr["idx"]))
            with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                deb, fin = xt[:, tr["bloc"] == "debut"], xt[:, tr["bloc"] == "fin"]
                d2 = np.concatenate([np.diff(xt[:, tr["serie"] == s], axis=1) ** 2 for s in np.unique(tr["serie"])], 1)
                se = np.sqrt(np.nanmean(d2, 1) / 2 * (1 / np.isfinite(deb).sum(1) + 1 / np.isfinite(fin).sum(1)))
                z = sens * (np.nanmean(fin, 1) - np.nanmean(deb, 1)) / se
                r[f"z_{nom}"] = float(np.nanmedian(z))
                r[f"delta_{nom}_mM"] = float(np.nanmedian(np.nanmean(fin, 1) - np.nanmean(deb, 1)) / pente) if pente else np.nan
    # basal par l'isotherme de toutes les données
    cm, iso = lecture_mM(xs[K["m"]["basal"]], [med["0 mM"], med["20 mM"], med["131 mM"]])
    r["basal_mM"], r["isotherme"] = float(np.nanmedian(cm)), bool(iso)
    return r, sens, pente


def tout_noter(D, cache):
    meta, V, P, plis = cache["meta"], cache["V"], cache["P"], cache["plis"]
    S, Ch = meta["soma"], meta["champ"]
    K = commun(meta)
    K["lnph"] = np.log(V["soma"]["_NB"] + V["soma"]["_NY"])
    cond = K["cond"]
    noms = [q for q in V["soma"] if not q.startswith("_")]
    lignes, sens, pentes = [], {}, {}
    for q in noms:
        X = {niv: (V[niv][q] if q in V[niv] else None) for niv in V}
        r, sens[q], pentes[q] = noter(X, K)
        # E. validation croisée ; K. champs entiers
        for mode in ("fov", "alea"):
            xs = [P[mode][k]["soma"].get(q, V["soma"][q]) for k in range(5)]
            e, b, iso = valider(xs, plis[mode], cond)
            lab = "champ" if mode == "fov" else "alea"
            for c, (rm, bi) in e.items():
                r[f"rmse_{lab}_{c}"], r[f"biais_{lab}_{c}"] = rm, bi
            r[f"basal_{lab}"], r[f"basal_{lab}_min"], r[f"basal_{lab}_max"] = b
        if q in V["champ"]:
            for c, v in valider_champs([P["fov"][k]["champ"].get(q, V["champ"][q]) for k in range(5)], Ch).items():
                r[f"rmse_reuni_{c}"] = v
        lignes.append(dict(lecture=q, **r))
    T = pd.DataFrame(lignes)

    # modèle de l'étape 9 : exactitude reprise de l'étape 10 (mêmes plis), hybrides au niveau du soma
    def de_etape10(fichier, prefixe, table):
        f = D["results"] / fichier
        if not f.exists():
            return
        e = pd.read_csv(f)
        col = "rmse" if "rmse" in e else "rmse_fov"
        for lec, cible in (("flim_2p", "mod2_flim"), ("flim_1p", "mod1_flim"), ("hybride_2p", "mod2_hyb"), ("hybride_1p", "mod1_hyb")):
            for _, x in e[e.lecture == lec].iterrows():
                table.setdefault(cible, {})[f"rmse_{prefixe}_{x.condition.replace(' ', '')}"] = x[col]
                table[cible][f"biais_{prefixe}_{x.condition.replace(' ', '')}"] = x.biais
        b = pd.read_csv(D["results"] / fichier.replace("resume", "basal")).set_index("lecture")
        for lec, cible in (("flim_2p", "mod2_flim"), ("flim_1p", "mod1_flim"), ("hybride_2p", "mod2_hyb"), ("hybride_1p", "mod1_hyb")):
            if lec in b.index:
                table.setdefault(cible, {}).update({f"basal_{prefixe}": b.loc[lec, "median"],
                                                    f"basal_{prefixe}_min": b.loc[lec, "min"], f"basal_{prefixe}_max": b.loc[lec, "max"]})
    e10 = {}
    de_etape10("validation_croisee_resume.csv", "champ", e10)
    de_etape10("validation_aleatoire_resume.csv", "alea", e10)
    for suf, nom in (("", "mod2_hyb"), ("_1t", "mod1_hyb")):
        f = D["results"] / f"deux_canaux_somas{suf}.csv"
        if f.exists():
            h = pd.read_csv(f)[["acquisition", "roi", "c_hybride"]].rename(columns={"acquisition": "acq"})
            x = S.merge(h, on=["acq", "roi"], how="left").c_hybride.values
            r, sens[nom], pentes[nom] = noter({"soma": x}, K)
            T = pd.concat([T, pd.DataFrame([dict(lecture=nom, **r)])], ignore_index=True)
            V["soma"][nom] = x
    T = T.set_index("lecture")
    for cible, d in e10.items():
        if cible in T.index:
            for c, v in d.items():
                T.loc[cible, c] = v
    T["famille"] = [next((f for f, v in FAMILLE.items() if q in v), "") for q in T.index]

    # classement : moyenne des rangs sur les huit critères, parmi les lectures qui les ont tous
    complet = T[[c for c, _ in CLASSEMENT]].notna().all(1)
    rg = pd.DataFrame({c: pd.Series(rankdata(s * T.loc[complet, c]), index=T.index[complet]) for c, s in CLASSEMENT})
    T["rang_moyen"] = rg.mean(1)
    T["rang_6"] = rg[[c for c in rg if c not in ("rmse_champ_20mM", "champs_mM")]].mean(1)
    for c in rg:
        T["rang_" + c] = rg[c]
    T = T.sort_values("rang_moyen", na_position="last")
    T.to_csv(D["results"] / "banc_lectures.csv")

    # valeurs par soma, précision par photon, champs, corrélations vraies
    pd.concat([S[["acq", "cond", "fov", "roi", "n_cycles"]], pd.DataFrame({q: V["soma"][q] for q in T.index if q in V["soma"]})],
              axis=1).to_csv(D["results"] / "banc_somas.csv", index=False)
    pd.DataFrame([dict(lecture=q, photons=nb, sigma_mM=T.loc[q, f"sigma_{nb}"], biais_mM=T.loc[q, f"biais_{nb}"],
                       echec_pct=T.loc[q, f"echec_{nb}"]) for q in T.index for nb in BUDGETS
                  if f"sigma_{nb}" in T and np.isfinite(T.loc[q, f"sigma_{nb}"])]).to_csv(D["results"] / "banc_precision.csv", index=False)
    g = cache["gab"]
    ch = Ch.copy()
    for q in ("R", "tm_B", "tm_Y", "tau_mono", "tphi_B", "gab_B", "gab_Y", "gab_formes", "gab_joint", "R_eta", "_gB", "_sB"):
        ch[q.lstrip("_")] = V["champ"][q]
    ch["f_R"] = (g["R0"] - ch.R) / (g["R0"] - g["R1"])
    ch["eta"] = ch.R / (g["R0"] * (1 - ch.gab_formes) + g["R1"] * ch.gab_formes)
    ch["photons_cfp"], ch["photons_yfp"] = V["champ"]["_NB"], V["champ"]["_NY"]
    ch.to_csv(D["results"] / "banc_champs.csv", index=False)
    cles = [q for q in ("R", "tm_B", "tau_mono", "tphi_B", "gab_B", "tm_Y", "gab_Y", "gab_formes", "lda_3") if q in V["demi"]]
    cv = {(a, b): cov_moities((V["soma"][a], V["demi"][a]), (V["soma"][b], V["demi"][b]), K) for a in cles for b in cles}
    with np.errstate(invalid="ignore"):
        cor = pd.DataFrame([[sens[a] * sens[b] * cv[(a, b)] / np.sqrt(cv[(a, a)] * cv[(b, b)]) for b in cles] for a in cles],
                           index=cles, columns=cles)
    cor.to_csv(D["results"] / "banc_correlations.csv")
    return T, cor, ch, K


def complements(D, cache, T, ch, K):
    """Trois vérifications qui ne demandent que le cache : (1) les mêmes notes sur la moitié la plus brillante des
    somas de chaque champ ; (2) écarts des champs à la moyenne de leur condition, lus par le rapport et par la forme
    CFP (part commune : l'état du capteur diffère vraiment ; différence : l'une des deux lectures se trompe) ;
    (3) time-lapses : changement de chaque cellule lu par le rapport et par la forme, bruit du rapport entre cycles."""
    meta, V, P, plis = cache["meta"], cache["V"], cache["P"], cache["plis"]
    S, cond = meta["soma"], K["cond"]
    # (1) somas brillants
    brillant = np.zeros(len(S), bool)
    for m in K["champs"].values():
        brillant[m] = K["lnph"][m] >= np.median(K["lnph"][m])
    garde = {"soma": brillant, "demi": brillant[meta["demi"].soma.values], "cycle": brillant[meta["cycle"].soma.values]}
    lig = []
    for q in T.index:
        if q not in V["cycle"]:
            continue
        X = {niv: np.where(g, V[niv][q], np.nan) for niv, g in garde.items()}
        r, _, _ = noter(X, K)
        e, b, _ = valider([np.where(brillant, P["fov"][k]["soma"].get(q, V["soma"][q]), np.nan) for k in range(5)], plis["fov"], cond)
        lig.append(dict(lecture=q, **{c: r[c] for c in ("auc_0_20", "champs_mM", "somas_mM", "vrai_somas_mM", "rho_photons",
                                                         "cycle_mM", "derive_mM_min", "basal_mM")},
                        rmse_champ_0mM=e["0mM"][0], rmse_champ_20mM=e["20mM"][0], basal_champ=b[0]))
    Tb = pd.DataFrame(lig).set_index("lecture")
    Tb.to_csv(D["results"] / "banc_lectures_brillants.csv")
    # (2) champs : rapport contre forme CFP
    fo = ch[ch.type == "fov"]
    lig = []
    for nom, s_ in list(fo.groupby("cond")) + [("0 et 20 mM", fo[fo.cond.isin(["0 mM", "20 mM"])])]:
        dR = s_.f_R - s_.groupby("cond").f_R.transform("mean")
        dB = s_.gab_B - s_.groupby("cond").gab_B.transform("mean")
        w, v = np.linalg.eigh(np.cov(dR, dB))
        lig.append(dict(champs=nom, n=len(s_), sd_rapport=dR.std(), sd_forme=dB.std(), r=np.corrcoef(dR, dB)[0, 1],
                        pente_forme_sur_rapport=v[1, 1] / v[0, 1], sd_commun=((dR + dB) / 2).std(),
                        sd_difference=((dB - dR) / 2).std(), forme_moins_rapport=(s_.gab_B - s_.f_R).mean()))
    Lh = pd.DataFrame(lig)
    Lh.to_csv(D["results"] / "banc_champs_rapport_forme.csv", index=False)
    # (3) time-lapses
    Vt, lig = V["trans"], []
    for nom, tr in K["trans"].items():
        M = lambda q: fini(mat(Vt[q], tr["idx"]))
        deb, fin = tr["bloc"] == "debut", tr["bloc"] == "fin"
        d = lambda q: np.nanmean(M(q)[:, fin], 1) - np.nanmean(M(q)[:, deb], 1)
        R = M("R")
        d2 = np.concatenate([np.diff(R[:, tr["serie"] == x], axis=1) ** 2 for x in np.unique(tr["serie"])], 1)
        poisson = np.nanmean(R * np.sqrt(1 / M("_NB") + 1 / M("_NY")), 1)
        dR = d("R")
        r = dict(time_lapse=nom, somas=len(dR), delta_R_pct=100 * np.nanmedian(dR / np.nanmean(R[:, deb], 1)),
                 bruit_R_sur_poisson=np.nanmedian(np.sqrt(np.nanmean(d2, 1) / 2)) / np.nanmedian(poisson))
        for q in ("gab_B", "tm_B", "tau_mono"):
            dq = d(q)
            ok = np.isfinite(dR) & np.isfinite(dq)
            r[f"delta_{q}"], r[f"rho_dR_d{q}"] = np.nanmedian(dq), spearmanr(dR[ok], dq[ok])[0]
        lig.append(r)
    Tl = pd.DataFrame(lig)
    Tl.to_csv(D["results"] / "banc_time_lapses.csv", index=False)
    return Tb, Lh, Tl


# ------------------------------------------------------------------ figures
COL_FAM = {"combien": G.PALETTE[1], "quand, CFP": G.PALETTE[0], "quand, YFP": "#c99a06", "deux canaux": G.PALETTE[2], "": G.INK_2}


def fig_classement(D, T):
    cols = [("rmse_champ_0mM", "erreur par\nchamp, 0 mM", 1), ("rmse_champ_20mM", "erreur par\nchamp, 20 mM", 1),
            ("champs_mM", "écart entre\nchamps", 1), ("auc_0_20", "AUC\n0 / 20 mM", -1), ("sigma_1000", "écart-type à\n1 000 photons", 1),
            ("cycle_mM", "bruit d'un\ncycle", 1), ("z_20_to_131", "z dans la\ncellule, 20→131", -1), ("derive_abs", "|dérive|\npar min", 1),
            ("vrai_somas_mM", "dispersion vraie\nentre somas", 1), ("rmse_alea_20mM", "erreur, tirage\naléatoire, 20 mM", 1),
            ("exces_bruit", "bruit / bruit\nde photons", 1), ("basal_champ", "basal\n(mM)", 0)]
    cols = [c for c in cols if c[0] in T]
    t = T[T.rang_moyen.notna()]
    M = np.full((len(t), len(cols)), np.nan)
    for j, (c, _, s) in enumerate(cols):
        v = t[c]
        if s:
            M[:, j] = (rankdata(s * np.nan_to_num(v.values.astype(float), nan=np.inf, posinf=1e300)) - 1) / max(len(t) - 1, 1)
    fig, ax = plt.subplots(figsize=(1.25 * len(cols) + 3.2, 0.36 * len(t) + 1.9), facecolor=G.SURFACE)
    ax.imshow(M, cmap="RdYlGn_r", vmin=0, vmax=1, aspect="auto", alpha=0.7)
    for i, q in enumerate(t.index):
        for j, (c, _, _) in enumerate(cols):
            v = t.loc[q, c]
            txt = "—" if pd.isna(v) else "∞" if not np.isfinite(v) else f"{v:.2f}" if c in ("auc_0_20", "exces_bruit") else f"{v:.1f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8.5, color=G.INK)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], fontsize=8.5, color=G.INK)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([f"{q}  ({t.loc[q, 'rang_moyen']:.1f})" for q in t.index], fontsize=9)
    for lab, q in zip(ax.get_yticklabels(), t.index):
        lab.set_color(COL_FAM[t.loc[q, "famille"]])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.axvline(7.5, color=G.INK, lw=1)
    fig.text(0.01, 0.005, "Lectures triées par rang moyen (entre parenthèses) sur les huit premières colonnes ; vert = meilleur de la "
             "colonne. Valeurs en mM-équivalents sauf AUC, z et rapport de bruit. Couleur du nom : famille (orange « combien », "
             "bleu forme CFP, jaune forme YFP, vert deux canaux).", fontsize=8, color=G.INK_2, wrap=True)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    G.sauver(fig, D["results"] / "banc_classement.png")


def fig_precision(D, T):
    choix = [q for q in ("R", "gab_joint", "lda_3", "gab_formes", "gab_B", "phas_axe_B", "tm_B", "tphi_B", "tau_mono",
                         "tau_mono_s", "tau_bi", "tau_3pop", "mod1_flim", "mod2_flim", "tm_Y") if q in T.index]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6))
    for ax in axes:
        G.style(ax)
    cm = plt.get_cmap("tab20")
    for i, q in enumerate(choix):
        s = [T.loc[q, f"sigma_{nb}"] for nb in BUDGETS]
        b = [T.loc[q, f"biais_{nb}"] for nb in BUDGETS]
        ls = "-" if T.loc[q, "famille"] != "combien" else "--"
        axes[0].plot(BUDGETS, s, "o-", color=cm(i), lw=1.4, ms=4, ls=ls, label=q)
        axes[1].plot(BUDGETS, b, "o-", color=cm(i), lw=1.4, ms=4, ls=ls, label=q)
    x = np.array(BUDGETS, float)
    if "R" in T.index:
        axes[0].plot(x, T.loc["R", "sigma_1000"] * np.sqrt(1000 / x), color=G.TRACE, lw=3, zorder=0, label="∝ 1/√N")
    axes[0].set_xscale("log"); axes[0].set_yscale("log")
    axes[1].set_xscale("log"); axes[1].set_yscale("symlog", linthresh=1)
    axes[1].axhline(0, color=G.TRACE)
    for ax in axes:
        ax.set_xlabel("photons du soma (CFP + YFP)", color=G.INK)
    axes[0].set_ylabel("écart-type entre tirages (mM-éq., pente 0–20 mM)", color=G.INK)
    axes[1].set_ylabel("biais par rapport à la somme complète (mM-éq.)", color=G.INK)
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    G.titre(axes[0], "Précision selon le nombre de photons", "somme des cycles de 400 somas (0 et 20 mM), amincie par tirage binomial")
    G.titre(axes[1], "Biais à bas photons", "une lecture qui se déplace quand les photons manquent fausse les somas pâles")
    fig.tight_layout()
    G.sauver(fig, D["results"] / "banc_precision.png")


def fig_deux_etats(D, cache, ch):
    g = cache["gab"]
    t = BA.T
    fig, axes = plt.subplots(2, 3, figsize=(17, 9.6))
    for ax in axes.flat:
        G.style(ax)
    # motif à 20 mM comme mélange des motifs à 0 et à 131 mM
    P20B, P20Y = cache["p20"]
    for ax, p0, p1, p20, nom in ((axes[0, 0], g["p0B"], g["p1B"], P20B / P20B.sum(), "CFP"),
                                 (axes[0, 1], g["p0Y"], g["p1Y"], P20Y / P20Y.sum(), "YFP")):
        d = p1 - p0
        f = float(np.sum((p20 - p0) * d / p20) / np.sum(d * d / p20))
        ax.plot(t, (p20 - p0) / p0 * 100, color=C.COULEUR_CONDITION["20 mM"], lw=1.2, label="20 mM")
        ax.plot(t, (p1 - p0) / p0 * 100, color=C.COULEUR_CONDITION["131 mM"], lw=1.2, label="131 mM")
        ax.plot(t, f * d / p0 * 100, color=G.INK, lw=1, ls="--", label=f"mélange : f = {f:.2f}")
        ax.axhline(0, color=G.TRACE)
        ax.set_xlabel("temps (ns)", color=G.INK)
        ax.set_ylabel("écart au motif 0 mM (%)", color=G.INK)
        ax.legend(frameon=False, fontsize=8.5)
        G.titre(ax, f"Forme du déclin {nom} : 20 mM est-il un mélange de 0 et 131 mM ?", "motifs moyens des champs, normalisés", fontsize=10.5)
    fo = ch[ch.type == "fov"]
    for ax, col, nom in ((axes[0, 2], "gab_B", "forme CFP"), (axes[1, 0], "gab_Y", "forme YFP")):
        for cond in C.CONDITIONS_FIXES:
            s = fo[fo.cond == cond]
            ax.plot(s.f_R, s[col], "o", ms=8, color=C.COULEUR_CONDITION[cond], label=cond)
            for _, r in s.iterrows():
                ax.text(r.f_R, r[col], f" {r.fov}", fontsize=8, color=C.COULEUR_CONDITION[cond], va="center")
        lim = [min(fo.f_R.min(), fo[col].min()) - 0.1, max(fo.f_R.max(), fo[col].max()) + 0.1]
        ax.plot(lim, lim, color=G.TRACE, lw=1)
        ax.set_xlabel("f d'après le rapport : (R0 − R) / (R0 − R1)", color=G.INK)
        ax.set_ylabel(f"f d'après la {nom}", color=G.INK)
        ax.legend(frameon=False, fontsize=8.5)
        G.titre(ax, f"Chaque champ : rapport contre {nom}", "sur la diagonale : un seul état, lu de deux façons", fontsize=10.5)
    ax = axes[1, 1]
    for cond in C.CONDITIONS_FIXES:
        s = fo[fo.cond == cond]
        ax.plot(s.gB, s.sB, "o", ms=8, color=C.COULEUR_CONDITION[cond], label=cond)
        for _, r in s.iterrows():
            ax.text(r.gB, r.sB, f" {r.fov}", fontsize=8, color=C.COULEUR_CONDITION[cond], va="center")
    th = np.linspace(0, np.pi, 200)
    ax.plot(0.5 + 0.5 * np.cos(th), 0.5 * np.sin(th), color=G.TRACE, lw=1)
    ax.set_xlim(fo.gB.min() - 0.03, fo.gB.max() + 0.03); ax.set_ylim(fo.sB.min() - 0.02, fo.sB.max() + 0.02)
    ax.set_xlabel("g", color=G.INK); ax.set_ylabel("s", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Phaseur CFP de chaque champ (photons réunis)", "rapporté à l'IRF ; arc : mono-exponentielles", fontsize=10.5)
    ax = axes[1, 2]
    for cond in C.CONDITIONS_FIXES:
        s = fo[fo.cond == cond]
        ax.plot(s.photons_cfp / s.somas / 1e3, s.eta, "o", ms=8, color=C.COULEUR_CONDITION[cond], label=cond)
        for _, r in s.iterrows():
            ax.text(r.photons_cfp / r.somas / 1e3, r.eta, f" {r.fov}", fontsize=8, color=C.COULEUR_CONDITION[cond], va="center")
    ax.axhline(1, color=G.TRACE)
    ax.set_xscale("log")
    ax.set_xlabel("photons CFP par soma (milliers)", color=G.INK)
    ax.set_ylabel("η du champ : R mesuré / R prédit par les formes", color=G.INK)
    G.titre(ax, "Écart hors de la courbe de chaque champ", "η = 1 : le rapport et les formes disent la même chose", fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "banc_deux_etats.png")


def fig_champs(D, T, cache, K):
    V = cache["V"]["soma"]
    choix = [q for q in list(T.index[T.rang_moyen.notna()][:5]) + ["R", "tau_mono", "tm_B", "gab_B", "mod1_flim"] if q in V]
    choix = list(dict.fromkeys(choix))[:8]
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.4))
    for ax, q in zip(axes.flat, choix):
        G.style(ax)
        x = V[q]
        cm, _ = lecture_mM(x, [T.loc[q, "med_0"], T.loc[q, "med_20"], T.loc[q, "med_131"]])
        for j, cond in enumerate(C.CONDITIONS_FIXES):
            for a, m in K["champs"].items():
                if K["cond"][m][0] != cond:
                    continue
                v = cm[m]
                v = v[~np.isnan(v)]
                me = np.median(v)
                q1, q3 = np.percentile(v, [25, 75])
                xx = j + (K["fov"][m][0] - 3) * 0.14
                ax.plot([xx, xx], [q1, min(q3, 400)], color=C.COULEUR_CONDITION[cond], lw=2, alpha=0.5)
                ax.plot(xx, me if np.isfinite(me) else 400, "o", ms=6, color=C.COULEUR_CONDITION[cond])
            if cond != "basal":
                ax.plot([j - 0.4, j + 0.4], [C.CL_MM[cond]] * 2, color=G.INK, ls="--", lw=1)
        ax.set_yscale("symlog", linthresh=20)
        ax.set_ylim(-20, 450)
        ax.set_xticks(range(4)); ax.set_xticklabels(C.CONDITIONS_FIXES)
        G.titre(ax, q, f"erreur par champ : {T.loc[q, 'rmse_champ_0mM']:.1f} mM à 0, {T.loc[q, 'rmse_champ_20mM']:.1f} à 20", fontsize=11)
    for ax in axes[:, 0]:
        ax.set_ylabel("[Cl⁻] lue (mM) : médiane et quartiles de chaque champ", color=G.INK)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "banc_champs.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--processus", type=int, default=6)
    ap.add_argument("--essai", action="store_true")
    ap.add_argument("--sans-calcul", action="store_true")
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    f = D["mesures"] / ("banc_cache_essai.pkl" if args.essai else "banc_cache.pkl")
    if args.sans_calcul:
        if not f.exists():
            sys.exit("pas de cache : lance d'abord sans --sans-calcul")
        cache = pickle.load(open(f, "rb"))
    else:
        cache = tout_calculer(root, D, args)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    T, cor, ch, K = tout_noter(D, cache)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_rows", 80)
    vue = ["rang_moyen", "rmse_champ_0mM", "rmse_champ_20mM", "champs_mM", "auc_0_20", "sigma_1000", "cycle_mM",
           "z_20_to_131", "derive_mM_min", "vrai_somas_mM", "somas_mM", "rho_photons", "exces_bruit", "basal_champ"]
    dit(T[[c for c in vue if c in T]].round(2).to_string())
    vue = ["auc_20_131", "rmse_champ_131mM", "rmse_alea_0mM", "rmse_alea_20mM", "rmse_reuni_0mM", "rmse_reuni_20mM", "sigma_300",
           "sigma_10000", "biais_300", "biais_1000", "echec_300", "delta_20_to_131_mM", "z_0_to_20", "delta_0_to_20_mM",
           "z_basal_to_0", "delta_basal_to_0_mM", "basal_mM", "monotone", "isotherme"]
    dit(T[[c for c in vue if c in T]].round(2).to_string())
    dit("\nCorrélations « vraies » entre lectures, d'un soma à l'autre dans un même champ (0 et 20 mM) :")
    dit(cor.round(2).to_string())
    Tb, Lh, Tl = complements(D, cache, T, ch, K)
    dit("\nMoitié la plus brillante des somas de chaque champ :")
    dit(Tb.loc[[q for q in T.index if q in Tb.index]].round(2).to_string())
    dit("\nÉcarts des champs à la moyenne de leur condition, rapport (f_R) contre forme CFP (gab_B), en f :")
    dit(Lh.round(3).to_string(index=False))
    dit("\nTime-lapses, cellule par cellule :")
    dit(Tl.round(3).to_string(index=False))
    dit("\nClassement sans les deux critères qui supposent la [Cl⁻] nominale de chaque champ (rang_6) :")
    dit(T.sort_values("rang_6")[["rang_6", "rang_moyen"]].dropna().round(2).head(15).T.to_string())
    fig_classement(D, T)
    fig_precision(D, T)
    fig_deux_etats(D, cache, ch)
    fig_champs(D, T, cache, K)


if __name__ == "__main__":
    main()
