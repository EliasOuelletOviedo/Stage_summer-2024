"""Étape 9 — FLIM à deux canaux : le modèle de Simulateurs/bayesian_filter, par maximum de vraisemblance,
soma par soma (sans filtre temporel), comparé au τ mono-exponentiel.

Usage :  python3 step9_flim_deux_canaux.py <dossier>

1. **IRF** : IRF paramétrique de sclm_hybride (gaussienne ⊗ traînes exponentielles + fond) ajustée sur l'IRF
   mesurée du jeu (`irf*.sdt`, ligne de chaque canal : config.IRF_LIGNE), ramenée sur la période.
2. **Calibration** avec les données complètes : déclins bleus (CFP) et jaunes (YFP) réunis de tous les somas
   gardés de chaque [Cl⁻] imposée (0, 20, 131 mM), ajustement de Poisson conjoint du modèle (deux_canaux.py) :
   τ1, τ2, τD, τA, λ, G, g2, κ, k0 partagés ; A, f1, f2, fonds et décalages d'horloge par condition ; η = 1 par
   définition dans ces conditions. Puis les isothermes de f1 et f2 (exactes : 3 concentrations).
3. **Lecture de chaque soma** (FOV des 4 conditions), maximum de vraisemblance sur une grille de c qui descend
   sous 0 jusqu'à la limite du modèle (f1 + f2 = 1, ~−15 mM ici) et monte à 300 mM : **les [Cl⁻] négatives sont
   gardées partout** (les tronquer biaiserait les lectures basses vers le haut) ; en haut, la grille va à 10⁵ mM et
   une lecture dont le maximum reste au bord est « hors courbe » (inf), comme le τ mono au-delà de son asymptote ;
   fond libre par canal (choix de l'utilisateur, 2026-11-10) :
     - temps de vie unique : une exponentielle sur le bleu → τ̂, puis c par la courbe τ(c) des déclins de calibration ;
     - FLIM du système : formes des deux canaux, η libre par soma (le ratio n'apporte rien) ;
     - ratio : la quantité seule, à η = 1 ;
     - **hybride** (formes + quantité) : η **partagé par le FOV** (comme θ dans le filtre), η̂ par vraisemblance
       profilée sur tous les somas du champ, puis c de chaque soma à η̂ ; intervalle à ΔlnL = 1,92.
4. **Comparaison** avec le τ mono-exponentiel de l'étape 6 (relu en mM par l'isotherme de l'étape 7).

Les somas des conditions imposées sont aussi ceux de la calibration : leur relecture dit la dispersion d'une
cellule à l'autre et la précision d'une lecture, pas l'exactitude d'une calibration indépendante.

Sorties : results/deux_canaux_calibration.json, results/deux_canaux_somas.csv, results/deux_canaux_resume.csv,
results/deux_canaux_eta.csv, results/deux_canaux_calibration.png, results/deux_canaux_lectures.png,
results/deux_canaux_tau_mono.png
"""
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

import config as C
import deux_canaux as DC
import figures as G
from figures import plt
from mesures import selection

SUFFIXE = ""
FENETRE = (17, 254)              # canaux pleins (comme sclm_flim/step8) : le dernier canal est rogné par le TDC
# Grille de c : partie négative jusqu'à la limite du modèle (c_min), positive jusqu'à 10⁵ mM. Une lecture dont le
# maximum reste au bord haut n'a pas de maximum fini (forme au-delà de l'asymptote) : elle est « hors courbe » (inf).
C_POSITIFS = np.r_[0.0, np.geomspace(0.05, 1e5, 520)]
TAU_GRILLE = np.linspace(0.3, 5.0, 471)
LN_ETA = np.log(np.geomspace(0.3, 3.0, 241))
CONDITIONS = ("0 mM", "20 mM", "131 mM")
LECTURES = [("c_tau_mono", "τ mono-exp. (étape 6) → isotherme de l'étape 7"),
            ("c_tau_unique", "temps de vie unique (IRF paramétrique) → τ(c)"),
            ("c_ratio", "ratio, η = 1"),
            ("c_flim", "FLIM du système (η libre par soma)"),
            ("c_hybride", "hybride, η partagé par le FOV")]
LECTURES_FIG = [l for l in LECTURES if l[0] != "c_tau_unique"]      # figure des lectures : sans le τ unique


def somas(D, A):
    """{acq: (labels, B (n, 256), Y (n, 256))} des somas gardés, cycles gardés réunis."""
    out = {}
    for a in A:
        m = np.load(D["mesures"] / f"{a['nom']}.npz")
        n = int(m["n_rois"])
        dec = m["declins"][m["garde"]].sum(0).astype(float)          # (n + 2, 2, 256)
        k = selection(D, a["nom"], n)
        out[a["nom"]] = (k, dec[k, 0], dec[k, 1])
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--transfert", type=int, choices=(1, 2), default=2,
                    help="populations en FRET du donneur : 2 (modèle du simulateur) ou 1 (un seul temps de transfert)")
    args = ap.parse_args()
    global SUFFIXE
    SUFFIXE = "" if args.transfert == 2 else "_1t"
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    t0 = time.time()
    A = [a for a in C.acquisitions(root) if a["type"] == "fov" and (D["mesures"] / f"{a['nom']}.npz").exists()]
    so = somas(D, A)
    cond_de = {a["nom"]: a["condition"] for a in A}

    # --- 1-2. calibration sur les déclins réunis
    conc = [C.CL_MM[c] for c in CONDITIONS]
    B = [sum(so[n][1].sum(0) for n in so if cond_de[n] == c) for c in CONDITIONS]
    Y = [sum(so[n][2].sum(0) for n in so if cond_de[n] == c) for c in CONDITIONS]
    d = DC.donnees(root, B, Y, conc, FENETRE)
    irfp = DC.ajuster_irf(d)
    d["irfp"] = irfp
    for k, v in irfp.items():
        print(f"  IRF voie {k} : μ {v['mu']:.3f} ns, σ {1000 * v['sig']:.0f} ps, traînes {np.round(v['taus'], 3)} ns, "
              f"poids {np.round(v['weights'], 3)}, déviance/ddl {v['dev_per_dof']:.2f}")
    R = DC.calibrer(d, irfp, args.transfert)
    print(f"  modèle : {args.transfert} population(s) en FRET")
    iso = DC.isothermes(R, conc)
    p = R["p"]
    se = dict(zip(R["P"].names, np.sqrt(np.maximum(np.diag(R["cov"]), 0))))
    print(f"  calibration : déviance/ddl {R['dev']:.0f}/{R['dof']} = {R['dev'] / R['dof']:.2f}")
    print("  τ1 {t1:.4f}  τ2 {t2:.4f}  τD {tD:.4f}  τA {tA:.4f} ns | λ {lam:.3f}  G {G:.3f}  g2 {G2r:.3f}  "
          "κ {k1:.3f}  k0 {k0:.3f}".format(**p))
    for n_, v in iso.items():
        print(f"  {n_:2s} : valeurs {np.round(v['values'], 3)} → max {v['fmax']:.3f}, min {v['fmin']:.3f}, K {v['K']:.1f} mM")
    w = d["win"][1]
    mB, mY = R["model_fn"](R["x"])
    cal = dict(n_transfert=args.transfert, irf=C.irf_sdt(root).name, fenetre=list(FENETRE), periode_ns=C.PERIODE_NS,
               deviance=R["dev"], ddl=R["dof"], constantes={k: float(p[k]) for k in DC.CLES},
               erreurs_types={k: float(se.get(k, np.nan)) for k in DC.CLES},
               isothermes={k: dict(fmax=float(v["fmax"]), fmin=float(v["fmin"]), K=float(v["K"]),
                                   valeurs=[float(x) for x in v["values"]]) for k, v in iso.items()},
               irf_parametrique={str(k): dict(mu=float(v["mu"]), sig=float(v["sig"]), poids=list(map(float, v["weights"])),
                                              taus=list(map(float, v["taus"])), dev_par_ddl=float(v["dev_per_dof"]))
                                 for k, v in irfp.items()},
               conditions=[dict(c_mM=c, A=float(p[f"A{i}"]), f1=float(p[f"f1_{i}"]), f2=float(p[f"f2_{i}"]),
                                decalage_bleu_ps=1e3 * float(p[f"s1_{i}"]), decalage_jaune_ps=1e3 * float(p[f"s2_{i}"]),
                                part_fond_bleu=float(p[f"bgB{i}"] * w.sum() / mB[i][w].sum()),
                                part_fond_jaune=float(p[f"bgY{i}"] * w.sum() / mY[i][w].sum()),
                                photons_bleus=float(d["B"][i][w].sum()), photons_jaunes=float(d["Y"][i][w].sum()))
                           for i, c in enumerate(conc)])

    # --- 3. lecteur, courbe τ(c) du temps de vie unique
    S = DC.senseur(d, R, iso)
    c_min = DC.c_min_modele(iso)
    c_grille = np.r_[np.linspace(c_min, 0, 200, endpoint=False), C_POSITIFS]
    print(f"  grille de c : {c_min:.2f} à {c_grille[-1]:g} mM (sous {c_min:.2f} mM, f1 + f2 > 1 : le modèle n'a plus de sens)")
    cal["c_min_modele"] = c_min
    L = DC.LecteurSoma(S, c_grille, TAU_GRILLE)
    tau_cal = np.array([L.tau_unique(b[w]) for b in d["B"]])
    print(f"  temps de vie unique des déclins de calibration : {np.round(tau_cal, 4)} ns")
    cal["tau_unique_calibration"] = tau_cal.tolist()
    cal["R1"] = {f"{c:g}": float(np.exp(np.interp(c, L.c, L.lnR1))) for c in (0, 5, 20, 131)}
    monotone = np.all(np.diff(tau_cal) > 0)
    pchip = PchipInterpolator(tau_cal, conc, extrapolate=False) if monotone else None

    def c_depuis_tau(t):
        if pchip is None:
            return np.nan
        if t < tau_cal[0]:
            c = conc[0] + (t - tau_cal[0]) * (conc[1] - conc[0]) / (tau_cal[1] - tau_cal[0])
        elif t > tau_cal[-1]:
            c = conc[-1] + (t - tau_cal[-1]) * (conc[-1] - conc[-2]) / (tau_cal[-1] - tau_cal[-2])
        else:
            c = float(pchip(t))
        return float(c)                               # négatif gardé : pas de plancher

    # lecture de tous les somas
    rows, profils = [], {}
    for a in A:
        labs, Bs, Ys = so[a["nom"]]
        for lab, b, y in zip(labs, Bs, Ys):
            b, y = b[w], y[w]
            if b.sum() < C.MIN_PHOTONS_FIT or y.sum() < C.MIN_PHOTONS_FIT:
                continue
            pr = L.lire(b, y)
            profils[(a["nom"], int(lab))] = pr
            i_f = int(np.argmax(pr["formes"]))
            quant = -0.5 * (pr["ln_q"] - L.lnR1) ** 2 / pr["v"]
            tu = L.tau_unique(b)
            rows.append(dict(acquisition=a["nom"], condition=a["condition"], cl_mM=a["cl_mM"], roi=int(lab),
                             photons_cfp=float(b.sum()), photons_yfp=float(y.sum()),
                             c_flim=DC.argmax_affine(pr["formes"], L.c),
                             eta_soma=float(np.exp(pr["ln_q"][i_f] - L.lnR1[i_f])),
                             c_ratio=DC.argmax_affine(quant, L.c),
                             tau_unique=tu, c_tau_unique=c_depuis_tau(tu),
                             fond_bleu=float(pr["bB"][i_f] / b.sum()), fond_jaune=float(pr["bY"][i_f] / y.sum())))
    T = pd.DataFrame(rows)
    haut = L.c[-1] * (1 - 1e-9)
    for col in ("c_flim", "c_ratio"):
        T.loc[T[col] >= haut, col] = np.inf

    # η partagé par le FOV : vraisemblance profilée
    eta_rows = []
    for nom, s in T.groupby("acquisition"):
        prof = [profils[(nom, int(r))] for r in s.roi]
        prof_fov = lambda grille: sum(np.array([L.hybride(pr, x).max() for x in grille]) for pr in prof)
        le0 = LN_ETA[int(np.argmax(prof_fov(LN_ETA)))]          # grille grossière, puis fine (pas 0,05 %)
        fine = le0 + np.linspace(-0.05, 0.05, 201)
        tot = prof_fov(fine)
        le = DC.argmax_affine(tot, fine)
        lo, hi = DC.intervalle(tot, fine)
        eta_rows.append(dict(acquisition=nom, condition=cond_de[nom], somas=len(s), eta=np.exp(le),
                             eta_bas=np.exp(lo), eta_haut=np.exp(hi)))
        for idx, pr in zip(s.index, prof):
            Lh = L.hybride(pr, le)
            ch = DC.argmax_affine(Lh, L.c)
            T.loc[idx, "c_hybride"] = np.inf if ch >= L.c[-1] * (1 - 1e-9) else ch
            T.loc[idx, ["c_hybride_bas", "c_hybride_haut"]] = DC.intervalle(Lh, L.c)
            T.loc[idx, "eta_fov"] = np.exp(le)
    E = pd.DataFrame(eta_rows)

    # --- 4. τ mono-exponentiel (étape 6) et sa relecture (étape 7)
    s6 = pd.read_csv(D["results"] / "somas.csv")
    s6 = s6[s6.type == "soma"][["acquisition", "roi", "tau_mono", "R", "tau_m"]]
    T = T.merge(s6, on=["acquisition", "roi"], how="left")
    c7 = pd.read_csv(D["results"] / "calibration.csv").set_index("lecture")
    X0, Xinf, K = (float(c7.loc["tau_mono", k]) for k in ("X0", "Xinf", "K_mM"))
    with np.errstate(divide="ignore", invalid="ignore"):
        cm = K * (X0 - T.tau_mono) / (T.tau_mono - Xinf)
    T["c_tau_mono"] = np.where((T.tau_mono - Xinf) * (X0 - Xinf) > 0, cm, np.inf)

    # résumé
    res = []
    for cond in C.CONDITIONS_FIXES:
        s = T[T.condition == cond]
        for col, nom in LECTURES:
            v = s[col]                     # hors courbe = +inf, gardé : l'exclure biaiserait la médiane vers le bas
            bord = 100 * float((s[col] <= L.c[0] + 1e-9).mean()) if col in ("c_flim", "c_ratio", "c_hybride") else np.nan
            res.append(dict(condition=cond, lecture=col, nom=nom, somas=len(s), mediane=v.median(), au_bord_bas_pct=bord,
                            q1=v.quantile(0.25), q3=v.quantile(0.75), iqr=v.quantile(0.75) - v.quantile(0.25),
                            hors_courbe_pct=100 * float(np.isinf(s[col]).mean() + s[col].isna().mean()),
                            mediane_fov_ecart=s.groupby("acquisition")[col].median().replace(np.inf, np.nan).std()))
        if len(s):
            res.append(dict(condition=cond, lecture="largeur_ic_hybride", nom="largeur de l'IC 95 % (hybride)",
                            somas=len(s), mediane=(s.c_hybride_haut - s.c_hybride_bas).median()))
    RS = pd.DataFrame(res)
    T.to_csv(D["results"] / f"deux_canaux_somas{SUFFIXE}.csv", index=False)
    RS.to_csv(D["results"] / f"deux_canaux_resume{SUFFIXE}.csv", index=False)
    E.to_csv(D["results"] / f"deux_canaux_eta{SUFFIXE}.csv", index=False)
    (D["results"] / f"deux_canaux_calibration{SUFFIXE}.json").write_text(json.dumps(cal, indent=1, ensure_ascii=False))
    print(RS[RS.lecture != "largeur_ic_hybride"].pivot(index="lecture", columns="condition", values="mediane")
          .round(2).to_string())
    print(RS[RS.lecture != "largeur_ic_hybride"].pivot(index="lecture", columns="condition", values="iqr")
          .round(2).to_string())
    print(E.round(3).to_string(index=False))
    print(f"  -> results/deux_canaux_* ({time.time() - t0:.0f} s)")
    figures(D, d, R, iso, T, E, L, tau_cal, conc)


def figures(D, d, R, iso, T, E, L, tau_cal, conc):
    t = C.temps_canaux()
    w = d["win"][1]
    mB, mY, comp = R["model_fn"](R["x"], parts=True)
    fig, axes = plt.subplots(3, 4, figsize=(19, 11), gridspec_kw=dict(height_ratios=[3, 1.1, 3]))
    for i, c in enumerate(conc):
        ax, axr = axes[0, i], axes[1, i]
        for x in (ax, axr):
            G.style(x)
        for y, m, col, nom in ((d["B"][i], mB[i], G.COULEUR_CANAL["cfp"], "CFP"), (d["Y"][i], mY[i], G.COULEUR_CANAL["yfp"], "YFP")):
            ax.semilogy(t[w], y[w], ".", ms=2, color=col, label=f"{nom}")
            ax.semilogy(t[w], m[w], color=G.INK, lw=0.8)
            axr.plot(t[w], SH_resid(y[w], m[w]), color=col, lw=0.6)
        cp = comp[i]
        for k, ls in (("fuite", ":"), ("sensibilisee", "--"), ("directe", "-.")):
            ax.semilogy(t[w], np.maximum(cp[k][w], 1e-3), ls, color=G.INK_2, lw=0.8, label=f"jaune : {k}")
        ax.set_ylim(max(min(d["B"][i][w].min(), d["Y"][i][w].min()), 1) * 0.5, max(d["B"][i].max(), d["Y"][i].max()) * 2)
        G.titre(ax, f"{c:g} mM — somas réunis", f"f1 {R['p'][f'f1_{i}']:.3f}, f2 {R['p'][f'f2_{i}']:.3f}", fontsize=10.5)
        axr.axhline(0, color=G.TRACE)
        axr.set_xlabel("temps (ns)", color=G.INK)
        if i == 0:
            ax.legend(frameon=False, fontsize=7.5)
            axr.set_ylabel("résidus de déviance", color=G.INK)
    axes[0, 3].axis("off"); axes[1, 3].axis("off")
    k = R["p"]
    axes[0, 3].text(0, 1, "Constantes ajustées\n\n" + "\n".join([
        f"τ1 = {k['t1']:.3f} ns", f"τ2 = {k['t2']:.3f} ns", f"τD = {k['tD']:.3f} ns", f"τA = {k['tA']:.3f} ns",
        f"λ = {k['lam']:.3f}", f"G = {k['G']:.3f}, g2 = {k['G2r']:.3f}", f"κ = {k['k1']:.3f}, k0 = {k['k0']:.3f}",
        f"déviance/ddl = {R['dev'] / R['dof']:.2f}"]), va="top", fontsize=10.5, color=G.INK, family="monospace")
    cc = np.r_[np.linspace(L.c[0], 0, 50, endpoint=False), np.geomspace(0.05, 300, 300)]
    ax = axes[2, 0]; G.style(ax)
    for n_, col in (("f1", G.PALETTE[0]), ("f2", G.PALETTE[1]), ("x", G.INK)):
        v = iso[n_]
        ax.plot(conc, v["values"], "o", color=col)
        ax.plot(cc, DC.SH.isotherm(cc, v["fmax"], v["fmin"], v["K"]), color=col, lw=1.2, label=f"{n_} : K {v['K']:.1f} mM")
    ax.set_xscale("symlog", linthresh=1); ax.set_xlabel("[Cl⁻] (mM)", color=G.INK); ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Isothermes des populations", "exactes : 3 concentrations", fontsize=10.5)
    ax = axes[2, 1]; G.style(ax)
    ax.plot(L.c, np.exp(L.lnR1), color=G.INK); ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("[Cl⁻] (mM)", color=G.INK); ax.set_ylabel("R₁(c) = M_Y / M_B", color=G.INK)
    G.titre(ax, "Ratio du modèle à η = 1", None, fontsize=10.5)
    ax = axes[2, 2]; G.style(ax)
    ax.plot(conc, tau_cal, "o-", color=G.COULEUR_CANAL["cfp"])
    ax.set_xlabel("[Cl⁻] (mM)", color=G.INK); ax.set_ylabel("τ unique (ns)", color=G.INK)
    G.titre(ax, "Temps de vie unique des déclins de calibration", None, fontsize=10.5)
    ax = axes[2, 3]; G.style(ax)
    for cond in C.CONDITIONS_FIXES:
        e = E[E.condition == cond]
        ax.errorbar(np.full(len(e), C.CONDITIONS_FIXES.index(cond)) + np.linspace(-0.2, 0.2, len(e)), e.eta,
                    yerr=[np.maximum(e.eta - e.eta_bas, 0), np.maximum(e.eta_haut - e.eta, 0)], fmt="o", color=C.COULEUR_CONDITION[cond], capsize=2)
    ax.axhline(1, color=G.TRACE)
    ax.set_xticks(range(4)); ax.set_xticklabels(C.CONDITIONS_FIXES)
    ax.set_ylabel("η̂ du FOV", color=G.INK)
    G.titre(ax, "η partagé par chaque FOV", "1 = conditions de la calibration réunies", fontsize=10.5)
    fig.suptitle("Calibration FLIM à deux canaux (modèle de bayesian_filter) sur les somas réunis — "
                 + ("une population en FRET" if R.get("n_transfert", 2) == 1 else "deux populations en FRET"), color=G.INK, fontsize=13)
    fig.tight_layout()
    G.sauver(fig, D["results"] / f"deux_canaux_calibration{SUFFIXE}.png")

    # lectures par condition : chaque boîte sans les lectures hors courbe ni celles à plus de 5 σ de sa moyenne
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.4), sharey=True)
    bornes = []                                   # moyenne ± 5 σ de chaque boîte (après retrait)
    for ax, cond in zip(axes, C.CONDITIONS_FIXES):
        G.style(ax)
        s = T[T.condition == cond]
        data, retires = [], []
        for col, _ in LECTURES_FIG:
            v = s[col].to_numpy(float)
            v = v[np.isfinite(v)]
            garde = np.ones(len(v), bool)
            while True:                                   # retrait itéré : moyenne et σ recalculés sans les retirés
                m, sd = v[garde].mean(), v[garde].std()
                nouveau = garde & (np.abs(v - m) <= 5 * sd)
                if nouveau.sum() == garde.sum():
                    break
                garde = nouveau
            bornes.append((m - 5 * sd, m + 5 * sd))
            data.append(v[garde])
            retires.append(len(s) - garde.sum())
        bp = ax.boxplot(data, vert=True, widths=0.6, showfliers=False, patch_artist=True)
        for b in bp["boxes"]:
            b.set(facecolor=C.COULEUR_CONDITION[cond], alpha=0.35, edgecolor=G.INK_2)
        for j, v in enumerate(data):
            ax.plot(np.full(len(v), j + 1) + np.random.default_rng(j).normal(0, 0.06, len(v)), v, ".", ms=2,
                    color=C.COULEUR_CONDITION[cond], alpha=0.4)
        ax.axhline(0, color=G.TRACE, lw=0.8)
        if cond != "basal":
            ax.axhline(C.CL_MM[cond], color=G.INK, lw=1, ls="--")
        ax.set_xticks(range(1, len(LECTURES_FIG) + 1))
        ax.set_xticklabels([f"{n}\n(−{r})" for n, r in zip(["τ mono", "ratio", "FLIM\nsystème", "hybride"], retires)],
                           fontsize=8.5)
        G.titre(ax, f"{cond} — {len(s)} somas", ("tireté : [Cl⁻] imposée ; " if cond != "basal" else "") +
                "médianes : " + ", ".join(f"{np.median(v):.1f}" for v in data) + " mM", fontsize=10.5)
    axes[0].set_ylabel("[Cl⁻] lue (mM)", color=G.INK)
    axes[0].set_ylim(-50, 150)                    # axe commun fixé (choix de l'utilisateur, 2026-11-10)
    fig.text(0.01, 0.005, "Axe y commun de −50 à 150 mM (des points peuvent sortir du cadre ; les boîtes sont calculées sur toutes "
             "les valeurs gardées). Retirés de chaque boîte (entre parenthèses) : lectures hors courbe, puis lectures à plus de 5 σ de la "
             "moyenne, retrait répété jusqu'à ce que plus rien ne sorte. Les valeurs négatives sont gardées.",
             color=G.INK_2, fontsize=8.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    G.sauver(fig, D["results"] / f"deux_canaux_lectures{SUFFIXE}.png")

    # comparaison au τ mono-exponentiel
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    for ax in axes:
        G.style(ax)
    for cond in C.CONDITIONS_FIXES:
        s = T[T.condition == cond]
        col = C.COULEUR_CONDITION[cond]
        axes[0].plot(s.tau_mono, s.c_hybride, "o", ms=3, alpha=0.5, color=col, mec="none", label=cond)
        fin = np.isfinite(s.c_tau_mono)
        axes[1].plot(s.c_tau_mono[fin], s.c_hybride[fin], "o", ms=3, alpha=0.5, color=col, mec="none", label=cond)
        axes[2].plot(s.tau_mono, s.tau_unique, "o", ms=3, alpha=0.5, color=col, mec="none", label=cond)
    axes[0].set_yscale("symlog", linthresh=1); axes[0].set_xlabel("τ CFP mono-exp. (étape 6, ns)", color=G.INK)
    axes[0].set_ylabel("[Cl⁻] hybride (mM)", color=G.INK)
    G.titre(axes[0], "Hybride contre τ mono-exponentiel", "chaque point : un soma")
    axes[1].set_xscale("symlog", linthresh=1); axes[1].set_yscale("symlog", linthresh=1)
    axes[1].plot([-30, 3000], [-30, 3000], color=G.TRACE)
    axes[1].set_xlabel("[Cl⁻] par τ mono (isotherme de l'étape 7, mM)", color=G.INK)
    axes[1].set_ylabel("[Cl⁻] hybride (mM)", color=G.INK)
    G.titre(axes[1], "Les deux lectures en mM", "diagonale : accord ; échelle linéaire de −1 à 1, log au-delà (négatifs gardés)")
    axes[2].set_xlabel("τ mono-exp. (IRF mesurée, étape 6, ns)", color=G.INK)
    axes[2].set_ylabel("τ unique (IRF paramétrique, ns)", color=G.INK)
    G.titre(axes[2], "Deux mono-exponentielles, deux IRF", "fond libre dans les deux")
    axes[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    G.sauver(fig, D["results"] / f"deux_canaux_tau_mono{SUFFIXE}.png")


def SH_resid(y, m):
    return DC.SH.dev_resid(y, m)


if __name__ == "__main__":
    main()
