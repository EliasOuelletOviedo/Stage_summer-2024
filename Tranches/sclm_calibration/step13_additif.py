"""Étape 13 — ce que le tissu ajoute dans chaque pixel : le mesurer sans tranche témoin, le retirer, voir l'effet.

Usage :  python3 step13_additif.py <dossier>        (~1 min ; après l'étape 12, dont elle reprend le cache)

Constat de l'étape 12 : dans un même champ (même [Cl⁻] imposée), les somas pâles ont plus de photons très précoces
dans le canal CFP et un temps moyen plus court ; rapport et forme dépendent de la brillance, en sens opposés.
Hypothèse : chaque pixel reçoit, en plus du capteur, un signal du tissu (autofluorescence, diffusion) qui ne dépend
pas de la quantité de capteur. Alors, pour chaque canal de temps t de chaque canal de couleur, entre les somas i
d'un champ :

    photons_i(t) = σ(t) · N_i + ι(t) · A_i        N_i : photons du soma (CFP + YFP), A_i : pixels × cycles

σ(t) est le motif du capteur seul (somme 1 sur les deux canaux) et ι(t) ce qui s'ajoute par pixel et par cycle —
à un multiple près du motif du capteur (Σ ι ≈ 0 : « plus de tissu qui ressemble au capteur » ne se distingue pas de
« plus de capteur »). Régression pondérée (bruit de Poisson + 10 % de dispersion entre cellules), itérée. Elle
suppose que les somas du champ sont dans le même état et que cet état ne dépend pas de leur brillance.

Vérification hors échantillon : les somas de chaque champ sont séparés en deux moitiés au hasard ; le ι mesuré sur
une moitié est retiré des somas de l'autre. Les lectures fermées de l'étape 12 (et le gabarit CFP, ici par moindres
carrés : un déclin corrigé peut avoir des canaux négatifs) sont notées avant et après : dépendance à la brillance,
dispersion vraie entre somas, corrélation vraie entre rapport et forme, écart entre champs, validation par champ.

Seconde partie : ce qui s'ajoute est surtout une composante prompte (elle arrive avec l'impulsion et dure moins de
0,5 ns). Sa forme étant connue, on peut la laisser libre dans chaque cellule, sans régression entre cellules :
gabs_B (mélange des motifs du capteur seul à 0 et 131 mM), gabp_B (+ composante prompte), gabp_B_fond (+ fond plat),
R_prompt (rapport après retrait des photons prompts de chaque canal) — notés comme à l'étape 12.

Sorties : results/additif_champs.csv (un champ par ligne : ce qui s'ajoute, rapport et forme du capteur seul),
results/additif_lectures.csv (notes avant / après), results/additif_correlations.csv, results/additif_prompt.csv,
results/additif.png, mesures/additif.npz (σ et ι de chaque acquisition)
"""
import pickle
import sys

import numpy as np
import pandas as pd

import banc as BA
import config as C
import figures as G
import step12_banc_essai as B12
from figures import plt
from mesures import selection

W, T, CONDITIONS = BA.W, BA.T, BA.CONDITIONS
NW = BA.NW
DISPERSION = 0.10
GRAINE = 20261111
LECTURES = ("R", "R_tard", "tm_B", "porte_B", "tphi_B", "gab_B_mc", "tm_queue_B", "tm_Y")


def regression(Z, N, A, iters=6):
    """Colonnes de Z (somas × canaux de temps) = σ·N + ι·A, poids 1 / (μ + (DISPERSION·μ)²), itéré."""
    mu = np.maximum(np.outer(N, Z.sum(0) / N.sum()), 0.5)
    for _ in range(iters):
        w = 1 / (mu + (DISPERSION * mu) ** 2)
        a, b, c = (w * (N * N)[:, None]).sum(0), (w * (N * A)[:, None]).sum(0), (w * (A * A)[:, None]).sum(0)
        u, v = (w * N[:, None] * Z).sum(0), (w * A[:, None] * Z).sum(0)
        det = a * c - b * b
        s, i = (c * u - b * v) / det, (a * v - b * u) / det
        mu = np.maximum(np.outer(N, s) + np.outer(A, i), 0.5)
    return s, i


def charger(D, nom):
    m = np.load(D["mesures"] / f"{nom}.npz")
    n, g = int(m["n_rois"]), np.nonzero(m["garde"])[0]
    k = selection(D, nom, n)
    dec = m["declins"][g][:, k].astype(float)[..., W]                  # (cycles, somas, canal, fenêtre)
    aire = np.bincount(m["etiquettes"].ravel(), minlength=n + 2)[k].astype(float)
    return dec, aire, k


def gab_mc(Zw, p0, p1):
    """Part du motif « 131 mM » dans la forme d'un déclin, par moindres carrés pondérés (tolère des canaux négatifs)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        b = Zw / Zw.sum(1, keepdims=True)
    d, w = p1 - p0, 2 / (p0 + p1)
    return ((b - p0) * d * w).sum(1) / (d * d * w).sum()


def operateur(P):
    """Mélange de motifs par moindres carrés pondérés (poids 1 / motif moyen du capteur) : poids = déclin @ M."""
    w = 1 / np.maximum(P[:2].mean(0), 1e-6)
    return (P * w).T @ np.linalg.inv((P * w) @ P.T)


def motifs_capteur(motifs, fo, pk, sauf=0):
    """Par canal : motif du capteur seul à 0 et à 131 mM (pentes des régressions, moyenne des champs), motif de la
    composante prompte (ce qui s'ajoute de −0,3 à +0,6 ns autour du pic, moyenne des champs), fond plat. `sauf` :
    numéro du champ de chaque condition imposée laissé de côté (validation par champ)."""
    garde = [(a, c) for a, c, k in zip(fo.acq, fo.cond, fo.fov) if not (c in CONDITIONS and k == sauf)]
    sl = slice(pk - 6, pk + 13)
    out = {}
    for ch in ("B", "Y"):
        moy = lambda c: np.mean([motifs[a]["s" + ch] / motifs[a]["s" + ch].sum() for a, cc in garde if cc == c], 0)
        v = [np.clip(motifs[a]["i" + ch][sl], 0, None) for a, _ in garde]
        pr = np.zeros(NW)
        pr[sl] = np.mean([x / x.sum() for x in v if x.sum() > 0], 0)
        out[ch] = np.vstack([moy("0 mM"), moy("131 mM"), pr / pr.sum(), np.full(NW, 1 / NW)])
    return out


def lire_prompt(B, Y, mot):
    """Lectures cellule par cellule qui laissent libre la composante prompte (et, pour _fond, un fond plat)."""
    Bw, Yw = B[:, W].astype(float), Y[:, W].astype(float)
    a2, a3, a4 = Bw @ operateur(mot["B"][:2]), Bw @ operateur(mot["B"][:3]), Bw @ operateur(mot["B"])
    y3 = Yw @ operateur(mot["Y"][:3])
    with np.errstate(divide="ignore", invalid="ignore"):
        return {"gabs_B": a2[:, 1] / a2.sum(1), "gabp_B": a3[:, 1] / a3[:, :2].sum(1), "gabp_B_fond": a4[:, 1] / a4[:, :2].sum(1),
                "R_prompt": (Yw.sum(1) - y3[:, 2]) / (Bw.sum(1) - a3[:, 2]), "prompt_pct": 100 * a3[:, 2] / Bw.sum(1)}


def partie_prompt(root, D, cache, motifs, fo, pk, K):
    """La composante prompte prise en compte cellule par cellule : mêmes notes qu'à l'étape 12."""
    N = B12.assembler(root, D, False)
    B12.amincir(N, np.random.default_rng(B12.GRAINE + 1))
    S = cache["meta"]["soma"]
    assert (N["soma"]["meta"].roi.values == S.roi.values).all() and len(N["aminci"]["B"]) == len(cache["meta"]["aminci"])
    cond = K["cond"]
    tous = motifs_capteur(motifs, fo, pk)
    V = {niv: lire_prompt(N[niv]["B"], N[niv]["Y"], tous) for niv in ("soma", "demi", "cycle", "aminci", "trans")}
    plis = [lire_prompt(N["soma"]["B"], N["soma"]["Y"], motifs_capteur(motifs, fo, pk, k)) for k in range(1, 6)]
    lignes, sens = [], {}
    for q in ("gabs_B", "gabp_B", "gabp_B_fond", "R_prompt"):
        r, sens[q], _ = B12.noter({niv: V[niv][q] for niv in V}, K)
        e, b, _ = B12.valider([p[q] for p in plis], cache["plis"]["fov"], cond)
        lignes.append(dict(lecture=q, **{c: r.get(c, np.nan) for c in (
            "auc_0_20", "auc_20_131", "champs_mM", "somas_mM", "vrai_somas_mM", "rho_photons", "cycle_mM", "exces_bruit",
            "sigma_300", "sigma_1000", "sigma_10000", "biais_300", "biais_1000", "z_20_to_131", "z_basal_to_0",
            "derive_mM_min", "basal_mM", "med_0", "med_20", "med_131", "med_basal")},
                           rmse_champ_0mM=e["0mM"][0], rmse_champ_20mM=e["20mM"][0], basal_champ=b[0]))
    L = pd.DataFrame(lignes).set_index("lecture")
    # corrélation vraie avec le rapport, d'un soma à l'autre
    Vc = cache["V"]
    for q in ("gabs_B", "gabp_B", "gabp_B_fond"):
        cv = lambda xa, xb: B12.cov_moities(xa, xb, K)
        R_, X_ = (Vc["soma"]["R"], Vc["demi"]["R"]), (V["soma"][q], V["demi"][q])
        L.loc[q, "corr_vraie_avec_R"] = -sens[q] * cv(R_, X_) / np.sqrt(cv(R_, R_) * cv(X_, X_))
    L.to_csv(D["results"] / "additif_prompt.csv")
    pp = pd.Series(V["soma"]["prompt_pct"]).groupby(cond).median()
    print("\nComposante prompte laissée libre dans chaque cellule (moindres carrés, motifs du capteur seul) :")
    print(L.round(2).T.to_string())
    print("  part des photons CFP attribuée à la composante prompte, médiane des somas (%) :", pp.round(2).to_dict())
    return L


def plein(Zw):
    out = np.zeros((len(Zw), C.N_CANAUX))
    out[:, W] = Zw
    return out


def main():
    root = C.racine()
    D = C.dirs(root)
    f = D["mesures"] / "banc_cache.pkl"
    if not f.exists():
        sys.exit("mesures/banc_cache.pkl absent : lance d'abord step12_banc_essai.py")
    cache = pickle.load(open(f, "rb"))
    meta, plis = cache["meta"], cache["plis"]
    S, Ch = meta["soma"], meta["champ"]
    ins = B12.instrument(root, D, None, False)
    pk = ins.pic
    pic = slice(pk - 5, pk + 6)                                       # ± 0,25 ns autour du pic
    rng = np.random.default_rng(GRAINE)

    niveaux = {v: dict(soma=[[], []], demi=[[], []]) for v in ("brut", "corrige")}       # [B, Y] par niveau
    champs, motifs = [], {}
    for a in [x for x in C.acquisitions(root) if (D["mesures"] / f"{x['nom']}.npz").exists()]:
        dec, aire, k = charger(D, a["nom"])
        ncy, ns = dec.shape[:2]
        ZB, ZY = dec.sum(0)[:, 0], dec.sum(0)[:, 1]
        N, A = (ZB + ZY).sum(1), aire * ncy
        ok = N > 0
        sB, iB = regression(ZB[ok], N[ok], A[ok])
        sY, iY = regression(ZY[ok], N[ok], A[ok])
        motifs[a["nom"]] = dict(sB=sB, sY=sY, iB=iB, iY=iY)
        med = lambda v: float(np.median(v[ok]))
        champs.append(dict(acq=a["nom"], cond=a["condition"], type=a["type"], fov=B12.num_fov(a["nom"]), somas=int(ok.sum()),
                           photons_px=med(N / A), cfp_px=med(ZB.sum(1) / A), R_reuni=ZY[ok].sum() / ZB[ok].sum(),
                           R_capteur=sY.sum() / sB.sum(),
                           ajout_cfp_pic=iB[pic].sum(), ajout_yfp_pic=iY[pic].sum(), ajout_cfp=iB.sum(), ajout_yfp=iY.sum(),
                           part_pic_cfp_pct=100 * iB[pic].sum() / med(ZB[:, pic].sum(1) / A),
                           part_pic_yfp_pct=100 * iY[pic].sum() / med(ZY[:, pic].sum(1) / A),
                           part_cfp_pct=100 * iB.sum() / med(ZB.sum(1) / A)))
        if a["type"] != "fov":
            continue
        # ι d'une moitié des somas, retiré de l'autre
        moitie = rng.permutation(ns) % 2
        iBx, iYx = np.zeros((ns, NW)), np.zeros((ns, NW))
        for h in (0, 1):
            autre = ok & (moitie != h)
            iBx[moitie == h] = regression(ZB[autre], N[autre], A[autre])[1]
            iYx[moitie == h] = regression(ZY[autre], N[autre], A[autre])[1]
        for niv, sel, cy in (("soma", slice(None), ncy), ("demi", slice(0, None, 2), len(range(0, ncy, 2))),
                             ("demi", slice(1, None, 2), len(range(1, ncy, 2)))):
            b, y = dec[sel].sum(0)[:, 0], dec[sel].sum(0)[:, 1]
            niveaux["brut"][niv][0].append(b)
            niveaux["brut"][niv][1].append(y)
            niveaux["corrige"][niv][0].append(b - iBx * (aire * cy)[:, None])
            niveaux["corrige"][niv][1].append(y - iYx * (aire * cy)[:, None])
    Chp = pd.DataFrame(champs)
    np.savez_compressed(D["mesures"] / "additif.npz", **{f"{n}|{q}": v for n, d in motifs.items() for q, v in d.items()})

    # --- lectures avant / après, notées comme à l'étape 12
    K = B12.commun(meta)
    K["lnph"] = np.log(cache["V"]["soma"]["_NB"] + cache["V"]["soma"]["_NY"])
    cond, acq = K["cond"], K["acq"]
    fo = Chp[Chp.type == "fov"].reset_index(drop=True)
    X, refs = {}, {}
    for v in ("brut", "corrige"):
        Bs, Ys = np.concatenate(niveaux[v]["soma"][0]), np.concatenate(niveaux[v]["soma"][1])
        Bd, Yd = np.concatenate(niveaux[v]["demi"][0]), np.concatenate(niveaux[v]["demi"][1])
        assert len(Bs) == len(S) and len(Bd) == len(meta["demi"]), "ordre des somas différent de celui de l'étape 12"
        # motif CFP de chaque champ : photons réunis (brut) ou motif du capteur seul (corrigé)
        if v == "brut":
            pB = {a: Bs[acq == a].sum(0) / Bs[acq == a].sum() for a in fo.acq}
        else:
            pB = {a: motifs[a]["sB"] / motifs[a]["sB"].sum() for a in fo.acq}
        refs[v] = pB
        ref = lambda c, sauf=0: np.mean([pB[a] for a, cc, k in zip(fo.acq, fo.cond, fo.fov) if cc == c and k != sauf], 0)
        os_, od = BA.simples(plein(Bs), plein(Ys), ins), BA.simples(plein(Bd), plein(Yd), ins)
        os_["gab_B_mc"], od["gab_B_mc"] = gab_mc(Bs, ref("0 mM"), ref("131 mM")), gab_mc(Bd, ref("0 mM"), ref("131 mM"))
        for q in LECTURES:
            X[(q, v)] = dict(soma=os_[q], demi=od[q],
                             plis=[gab_mc(Bs, ref("0 mM", k), ref("131 mM", k)) for k in range(1, 6)] if q == "gab_B_mc" else None)
        fo[f"gab_B_{v}"] = [float(gab_mc(pB[a][None], ref("0 mM"), ref("131 mM"))[0]) for a in fo.acq]
    lignes, sens = [], {}
    for (q, v), x in X.items():
        r, sens[(q, v)], _ = B12.noter(dict(soma=x["soma"], demi=x["demi"]), K)
        e, b, _ = B12.valider(x["plis"] or [x["soma"]] * 5, plis["fov"], cond)
        lignes.append(dict(lecture=q, version=v, auc_0_20=r["auc_0_20"], champs_mM=r["champs_mM"], somas_mM=r["somas_mM"],
                           vrai_somas_mM=r["vrai_somas_mM"], rho_photons=r["rho_photons"], rmse_champ_0mM=e["0mM"][0],
                           rmse_champ_20mM=e["20mM"][0], basal_mM=r["basal_mM"], basal_champ=b[0],
                           med_0=r["med_0"], med_20=r["med_20"], med_131=r["med_131"], med_basal=r["med_basal"]))
    L = pd.DataFrame(lignes)
    L.to_csv(D["results"] / "additif_lectures.csv", index=False)
    cor = []
    for v in ("brut", "corrige"):
        for q in ("tm_B", "porte_B", "tphi_B", "gab_B_mc", "tm_queue_B"):
            cv = lambda a, b: B12.cov_moities((X[(a, v)]["soma"], X[(a, v)]["demi"]), (X[(b, v)]["soma"], X[(b, v)]["demi"]), K)
            cor.append(dict(version=v, forme=q, correlation_vraie_avec_R=sens[("R", v)] * sens[(q, v)] * cv("R", q)
                            / np.sqrt(cv("R", "R") * cv(q, q))))
    Cor = pd.DataFrame(cor)
    Cor.to_csv(D["results"] / "additif_correlations.csv", index=False)

    # --- champs : rapport et forme, photons réunis contre capteur seul
    R0 = {v: fo[fo.cond == "0 mM"][c].mean() for v, c in (("brut", "R_reuni"), ("corrige", "R_capteur"))}
    R1 = {v: fo[fo.cond == "131 mM"][c].mean() for v, c in (("brut", "R_reuni"), ("corrige", "R_capteur"))}
    fo["f_R_brut"] = (R0["brut"] - fo.R_reuni) / (R0["brut"] - R1["brut"])
    fo["f_R_corrige"] = (R0["corrige"] - fo.R_capteur) / (R0["corrige"] - R1["corrige"])
    Chp = Chp.merge(fo[["acq", "f_R_brut", "f_R_corrige", "gab_B_brut", "gab_B_corrige"]], on="acq", how="left")
    Chp.to_csv(D["results"] / "additif_champs.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(Chp[["acq", "somas", "photons_px", "R_reuni", "R_capteur", "ajout_cfp_pic", "ajout_yfp_pic", "part_pic_cfp_pct",
               "part_pic_yfp_pct", "part_cfp_pct", "f_R_brut", "gab_B_brut", "f_R_corrige", "gab_B_corrige"]].round(3).to_string(index=False))
    print("\nÉcarts des champs à la moyenne de leur condition : rapport contre forme CFP (en f)")
    for v in ("brut", "corrige"):
        for nom, s in list(fo.groupby("cond")) + [("0 et 20 mM", fo[fo.cond.isin(["0 mM", "20 mM"])])]:
            dR = s[f"f_R_{v}"] - s.groupby("cond")[f"f_R_{v}"].transform("mean")
            dB = s[f"gab_B_{v}"] - s.groupby("cond")[f"gab_B_{v}"].transform("mean")
            print(f"  {v:8s} {nom:>10} : sd rapport {dR.std():.3f}, sd forme {dB.std():.3f}, r {np.corrcoef(dR, dB)[0, 1]:+.2f}, "
                  f"sd de la demi-différence {((dB - dR) / 2).std():.3f} ; forme − rapport {(s[f'gab_B_{v}'] - s[f'f_R_{v}']).mean():+.3f}")
    print("\nLectures avant / après retrait de ce qui s'ajoute (ι d'une moitié des somas appliqué à l'autre) :")
    print(L.pivot(index="lecture", columns="version", values=["rho_photons", "vrai_somas_mM", "champs_mM", "auc_0_20",
                                                              "rmse_champ_0mM", "rmse_champ_20mM", "basal_mM"]).round(2).to_string())
    print("\nCorrélation vraie entre le rapport et la forme CFP, d'un soma à l'autre dans un champ (0 et 20 mM) :")
    print(Cor.pivot(index="forme", columns="version", values="correlation_vraie_avec_R").round(2).to_string())
    figure(D, Chp, motifs, L, fo, K, X, pk)
    partie_prompt(root, D, cache, motifs, fo, pk, K)


def figure(D, Chp, motifs, L, fo, K, X, pk):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    for ax in axes.flat:
        G.style(ax)
    # ce qui s'ajoute, par canal de temps
    for ax, cle, nom in ((axes[0, 0], "iB", "CFP"), (axes[0, 1], "iY", "YFP")):
        for cond in C.CONDITIONS_FIXES:
            v = np.mean([motifs[a][cle] / c for a, cc, c in zip(fo.acq, fo.cond, fo.cfp_px) if cc == cond], 0)
            ax.plot(T - T[pk], 100 * v, color=C.COULEUR_CONDITION[cond], lw=1.2, label=cond)
        ax.axhline(0, color=G.TRACE)
        ax.set_xlim(-0.6, 4)
        ax.set_xlabel("temps après le pic (ns)", color=G.INK)
        ax.set_ylabel("ajout par pixel, en % du CFP d'un soma médian (par canal de 49 ps)", color=G.INK, fontsize=8.5)
        ax.legend(frameon=False, fontsize=8.5)
        G.titre(ax, f"Ce qui s'ajoute dans chaque pixel, canal {nom}", "moyenne des champs de chaque condition ; Σ ≈ 0 sur les deux canaux", fontsize=10.5)
    # brillance : part du pic et rapport, avant / après, un champ exemple par condition
    ax = axes[0, 2]
    lnph = K["lnph"]
    for v, ls in (("brut", "--"), ("corrige", "-")):
        for cond in ("0 mM", "20 mM"):
            m = K["m"][cond]
            x = X[("tm_B", v)]["soma"]
            # tendance : temps moyen centré par champ contre brillance centrée par champ, par quintile
            dx, dl = np.full(len(x), np.nan), np.full(len(x), np.nan)
            for a in pd.unique(K["acq"][m]):
                mm = K["acq"] == a
                dx[mm], dl[mm] = x[mm] - np.nanmedian(x[mm]), lnph[mm] - np.median(lnph[mm])
            ok = m & np.isfinite(dx)
            q = np.quantile(dl[ok], np.linspace(0, 1, 7))
            xs = [np.median(dl[ok & (dl >= a) & (dl <= b)]) for a, b in zip(q[:-1], q[1:])]
            ys = [np.median(dx[ok & (dl >= a) & (dl <= b)]) for a, b in zip(q[:-1], q[1:])]
            ax.plot(np.exp(xs), ys, "o" + ls, color=C.COULEUR_CONDITION[cond], ms=4, lw=1.3,
                    label=f"{cond}, {'avant' if v == 'brut' else 'après'}")
    ax.axhline(0, color=G.TRACE)
    ax.set_xscale("log")
    ax.set_xlabel("photons du soma / médiane de son champ", color=G.INK)
    ax.set_ylabel("temps moyen CFP − médiane du champ (ns)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Temps moyen CFP selon la brillance, dans un champ", "ι mesuré sur une moitié des somas, retiré de l'autre", fontsize=10.5)
    # notes avant / après
    for ax, col, nom in ((axes[1, 0], "rho_photons", "dépendance à la brillance (ρ de Spearman dans un champ)"),
                         (axes[1, 1], "vrai_somas_mM", "dispersion vraie entre somas d'un champ (mM-éq.)")):
        p = L.pivot(index="lecture", columns="version", values=col).loc[list(LECTURES)]
        y = np.arange(len(p))
        ax.barh(y - 0.18, p["brut"], 0.34, color=G.TRACE, label="avant")
        ax.barh(y + 0.18, p["corrige"], 0.34, color=G.PALETTE[0], label="après")
        ax.set_yticks(y)
        ax.set_yticklabels(p.index)
        ax.invert_yaxis()
        ax.axvline(0, color=G.INK, lw=0.8)
        ax.legend(frameon=False, fontsize=8.5)
        G.titre(ax, nom[0].upper() + nom[1:], None, fontsize=10.5)
    ax = axes[1, 2]
    for v, mk, al in (("brut", "o", 0.35), ("corrige", "D", 1.0)):
        for cond in C.CONDITIONS_FIXES:
            s = fo[fo.cond == cond]
            ax.plot(s[f"f_R_{v}"], s[f"gab_B_{v}"], mk, ms=7, color=C.COULEUR_CONDITION[cond], alpha=al,
                    label=f"{cond}, {'photons réunis' if v == 'brut' else 'capteur seul'}" if cond in ("basal", "20 mM") else None)
    for _, r in fo.iterrows():
        ax.plot([r.f_R_brut, r.f_R_corrige], [r.gab_B_brut, r.gab_B_corrige], color=C.COULEUR_CONDITION[r.cond], lw=0.7, alpha=0.6)
    lim = [-0.3, 1.4]
    ax.plot(lim, lim, color=G.TRACE)
    ax.set_xlabel("f d'après le rapport", color=G.INK)
    ax.set_ylabel("f d'après la forme CFP", color=G.INK)
    ax.legend(frameon=False, fontsize=8)
    G.titre(ax, "Chaque champ : rapport contre forme", "pâle : photons réunis ; plein : motif du capteur seul (pentes)", fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "additif.png")


if __name__ == "__main__":
    main()
