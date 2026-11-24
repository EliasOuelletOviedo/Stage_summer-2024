"""Étape 16 — biais sur [Cl⁻] dû au signal du tissu, dans les tranches : la figure A.2 refaite avec les déclins réels.

Usage :  python3 step16_autofluorescence.py <dossier> [--tirages 600] [--N 1e4]     (après l'étape 13)

La figure A.2 du simulateur (Simulateurs/FLIM-FRET_2-channels, étape 20) mesure le biais qu'une autofluorescence
de la littérature (NADH dans le bleu, flavines dans le jaune, 5 % des photons) cause à chaque lecture. Ici, tout
vient des tranches :

  - le capteur : motifs du capteur seul (pentes des régressions de l'étape 13, σ, les deux canaux ensemble) à 0 et
    131 mM, moyennés sur les champs ; entre les deux, mélange en parts de photons h(c) = c(131 + K)/(131(c + K)),
    K tiré du motif à 20 mM ;
  - le tissu : ce que l'étape 13 voit s'ajouter dans chaque pixel (ι), réduit à sa partie prompte et positive
    (−0,3 à +0,6 ns autour du pic, la seule qui ne dépende pas de la convention Σι ≈ 0) ; niveau d'un soma médian,
    médiane des 20 champs (basal, 0, 20, 131 mM), indépendant du Cl⁻.

Lectures, chacune calibrée sur le capteur seul puis appliquée au capteur + tissu (IRF, décalages, connus) :
  rapport YFP/CFP brut ; τ du donneur (mono-exponentielle, CFP seul, fond libre) ; 2 canaux, motifs du capteur
  (part h et échelle commune) sans terme pour le tissu ; le même avec un terme par canal, de forme mesurée (plein)
  ou de forme NADH / flavines de la littérature (pointillé, comme le « mal caractérisé » de la figure A.2).

Contrôle sur les cellules réelles : le biais prédit au basal est comparé à ce que le retrait du signal du tissu
change dans les somas (étape 13 : rapport 4,4 → 3,5 mM, formes CFP 2,4–5,4 → 5,9–9,1 mM).

Sorties : results/autofluorescence_tranches.png, autofluorescence_tranches_biais.csv,
autofluorescence_tranches_precision.csv, autofluorescence_tranches_monde.json
"""
import argparse
import json

import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

import banc as BA
import config as C
import figures as G
import flim
from figures import plt

W, NW, T = BA.W, BA.NW, BA.T
C.MIN_PHOTONS_FIT = 1
NADH = dict(amps=[0.75, 0.25], taus=[0.4, 2.0])     # littérature (deux photons), comme la figure A.2
FAD = dict(amps=[0.75, 0.25], taus=[0.3, 2.6])
COUL = {"ratio": "#2a78d6", "mono": "#eb6834", "deux": "#008300", "deux_auto": "#4a3aa7"}
NOMS = {"ratio": "rapport YFP/CFP (brut)", "mono": "τ du donneur, mono-exponentielle (CFP seul)",
        "deux": "2 canaux, motifs du capteur, sans terme pour le tissu",
        "deux_auto": "2 canaux + terme pour le tissu (forme mesurée)",
        "deux_lit": "2 canaux + terme NADH / flavines de la littérature"}
VIRG = FuncFormatter(lambda x, _: f"{x:g}".replace(".", ",").replace("-", "−"))


def fr(x, nd, signe=False):
    """Nombre à la française : virgule décimale, vrai signe moins."""
    s = f"{x:+.{nd}f}" if signe else f"{x:.{nd}f}"
    return s.replace(".", ",").replace("-", "−")


def irls(y, X, it=60):
    """Poisson, lien identité : μ = X·θ (θ de signe libre), maximum de vraisemblance par scores de Fisher."""
    th = np.linalg.lstsq(X, y, rcond=None)[0]
    mu = X @ th
    if np.any(mu <= 0):
        th = np.linalg.lstsq(X, np.maximum(y, 0.5), rcond=None)[0]
    for _ in range(it):
        mu = np.maximum(X @ th, 1e-9)
        wX = X / mu[:, None]
        step = np.linalg.solve(X.T @ wX, wX.T @ (y - mu))
        a = 1.0
        while np.any(X @ (th + a * step) <= 0) and a > 1e-4:
            a /= 2
        th = th + a * step
        if np.max(np.abs(a * step)) < 1e-10 * max(1.0, np.max(np.abs(th))):
            break
    return th


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dossier", nargs="?")
    ap.add_argument("--tirages", type=int, default=600)
    ap.add_argument("--N", type=float, default=1e4)
    args = ap.parse_args()
    root = C.racine([args.dossier] if args.dossier else [])
    D = C.dirs(root)
    z = np.load(D["mesures"] / "additif.npz")
    fo = pd.read_csv(D["results"] / "additif_champs.csv")
    fo = fo[fo.type == "fov"]
    ins = json.loads((D["qc"] / "banc_instrument.json").read_text())
    sB, sY, pk = ins["sB"], ins["sY"], ins["pic"]
    MB = flim.Modele(C.charger_irf("cfp", root), fenetre=BA.FEN)
    MY = flim.Modele(C.charger_irf("yfp", root), fenetre=BA.FEN)

    # --- le capteur : motifs du capteur seul, les deux canaux ensemble (somme 1)
    motif = {c: np.mean([np.r_[z[f"{a}|sB"], z[f"{a}|sY"]] for a in fo[fo.cond == c].acq], 0)
             for c in ("basal", "0 mM", "20 mM", "131 mM")}
    u0, u1 = motif["0 mM"], motif["131 mM"]
    w8, d8 = 1 / ((u0 + u1) / 2), u1 - u0
    proj = lambda m: float(((m - u0) * d8 * w8).sum() / (d8 * d8 * w8).sum())
    h20, hb = proj(motif["20 mM"]), proj(motif["basal"])
    K = 20 * 131 * (1 - h20) / (131 * h20 - 20)
    h_de = lambda c: c * (131 + K) / (131 * (c + K))
    c_de = lambda h: 131 * K * h / (131 + K - 131 * h)
    c_basal = c_de(hb)
    print(f"  capteur : part 131 mM du motif à 20 mM {h20:.3f} → K = {K:.1f} mM ; basal h = {hb:.3f} → {c_basal:.1f} mM")

    # --- le tissu : partie prompte et positive de ι, niveau d'un soma médian (part des photons du soma, B + Y)
    sl = slice(pk - 6, pk + 13)
    formes, niv = {"B": [], "Y": []}, {"B": [], "Y": []}
    for _, r in fo.iterrows():
        for ch in ("B", "Y"):
            p = np.zeros(NW)
            p[sl] = np.clip(z[f"{r.acq}|i{ch}"][sl], 0, None)
            if p.sum() > 0:
                formes[ch].append(p / p.sum())
            niv[ch].append(p.sum() / r.photons_px)
    tissu = {ch: np.mean(formes[ch], 0) for ch in ("B", "Y")}
    for ch in ("B", "Y"):
        tissu[ch] /= tissu[ch].sum()
    a_med = {ch: float(np.median(niv[ch])) for ch in ("B", "Y")}
    a_rng = {ch: (float(np.min(niv[ch])), float(np.max(niv[ch]))) for ch in ("B", "Y")}
    lit = {"B": MB.forme(NADH["taus"], NADH["amps"], sB)[W], "Y": MY.forme(FAD["taus"], FAD["amps"], sY)[W]}
    lit = {ch: v / v.sum() for ch, v in lit.items()}
    monde = lambda c: (1 - h_de(c)) * u0 + h_de(c) * u1
    partB20 = monde(20.0)[:NW].sum()
    print(f"  tissu, soma médian : CFP {100 * a_med['B']:.2f} % des photons du soma ({100 * a_med['B'] / partB20:.1f} % du CFP à 20 mM), "
          f"YFP {100 * a_med['Y']:.2f} % ({100 * a_med['Y'] / (1 - partB20):.1f} % du YFP) ; champs : CFP "
          f"{100 * a_rng['B'][0]:.2f}–{100 * a_rng['B'][1]:.2f} %, YFP {100 * a_rng['Y'][0]:.2f}–{100 * a_rng['Y'][1]:.2f} %")

    def ajoute(v, aB, aY):
        out = v.copy()
        out[:NW] += aB * tissu["B"]
        out[NW:] += aY * tissu["Y"]
        return out

    # --- lectures
    plein = lambda b: np.r_[np.zeros(BA.FEN[0]), b, np.zeros(C.N_CANAUX - BA.FEN[1])]
    zB, zY = np.zeros(NW), np.zeros(NW)
    X2 = np.c_[u0, u1]
    X4m = np.c_[u0, u1, np.r_[tissu["B"], zY], np.r_[zB, tissu["Y"]]]
    X4l = np.c_[u0, u1, np.r_[lit["B"], zY], np.r_[zB, lit["Y"]]]

    def lire(v):
        out = dict(ratio=v[NW:].sum() / v[:NW].sum(), mono=MB.ajuster_mono(plein(v[:NW]), s=sB, tau0=1.6)[0])
        for cle, X in (("deux", X2), ("deux_auto", X4m), ("deux_lit", X4l)):
            th = irls(v, X)
            out[cle] = th[1] / (th[0] + th[1])
        return out

    ECH = 1e6                                          # photons du soma pour les lectures sans bruit
    c_bas = c_de(-0.25)
    grille = np.linspace(c_bas, 400.0, 500)
    cal = pd.DataFrame([lire(ECH * monde(c)) for c in grille])
    for k in cal:
        dd = np.diff(cal[k].values)
        if not (np.all(dd > 0) or np.all(dd < 0)):
            print(f"  attention : calibration non monotone pour {k}")

    def inverser(k, x):
        xs, cs = cal[k].values, grille
        if xs[-1] < xs[0]:
            xs, cs = xs[::-1], cs[::-1]
        return np.interp(x, xs, cs, left=np.nan, right=np.nan)

    controle = max(abs(inverser(k, lire(ECH * monde(c))[k]) - c) for k in cal for c in (0.0, 5.0, 20.0, 80.0))
    print(f"  contrôle sans tissu : écart max {controle:.1e} mM")

    # --- biais : tissu du CFP seul, du YFP seul, des deux ; niveau du soma médian
    lignes = []
    cs = np.arange(0, 80.01, 2.5)
    for sc, (fB, fY) in (("CFP", (1, 0)), ("YFP", (0, 1)), ("les_deux", (1, 1))):
        for c in cs:
            lu = lire(ECH * ajoute(monde(c), fB * a_med["B"], fY * a_med["Y"]))
            for k, x in lu.items():
                lignes.append(dict(scenario=sc, cl_vraie=c, lecture=k, cl_lue=inverser(k, x), biais=inverser(k, x) - c))
    Bi = pd.DataFrame(lignes)
    Bi.to_csv(D["results"] / "autofluorescence_tranches_biais.csv", index=False)

    # contrôle sur les cellules réelles : basal, tissu des deux canaux, au niveau médian des 20 champs et à celui
    # des champs basaux (imagés à plus faible puissance, plus de signal du tissu par photon du capteur)
    est_basal = (fo.cond == "basal").values
    a_bas = {ch: float(np.median(np.asarray(niv[ch])[est_basal])) for ch in ("B", "Y")}
    pred, pred_bas = {}, {}
    for d_, a_ in ((pred, a_med), (pred_bas, a_bas)):
        lu = lire(ECH * ajoute(monde(c_basal), a_["B"], a_["Y"]))
        d_.update({k: inverser(k, x) - c_basal for k, x in lu.items()})
    print(f"  basal ({c_basal:.1f} mM), tissu des deux canaux, biais prédit (niveau médian des champs) : " +
          ", ".join(f"{k} {v:+.2f} mM" for k, v in pred.items()))
    print(f"  idem, niveau des champs basaux (CFP {100 * a_bas['B']:.2f} %, YFP {100 * a_bas['Y']:.2f} %) : " +
          ", ".join(f"{k} {v:+.2f} mM" for k, v in pred_bas.items()))

    # --- précision à N photons CFP (à 20 mM), tissu des deux canaux
    rng = np.random.default_rng(20261123)
    k_N = args.N / partB20
    prec = []
    for c in (0.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0):
        mu = k_N * ajoute(monde(c), a_med["B"], a_med["Y"])
        lus = {k: [] for k in ("ratio", "mono", "deux", "deux_auto")}
        for _ in range(args.tirages):
            y = rng.poisson(mu).astype(float)
            r = {"ratio": y[NW:].sum() / y[:NW].sum(), "mono": MB.ajuster_mono(plein(y[:NW]), s=sB, tau0=1.6)[0]}
            for cle, X in (("deux", X2), ("deux_auto", X4m)):
                th = irls(y, X)
                r[cle] = th[1] / (th[0] + th[1])
            for k in lus:
                lus[k].append(inverser(k, r[k]))
        for k, v in lus.items():
            v = np.asarray(v)
            ok = v[np.isfinite(v)]
            prec.append(dict(cl_vraie=c, lecture=k, N_cfp=args.N, ecart_type_robuste=np.subtract(*np.percentile(ok, [75, 25])) / 1.349,
                             hors_echelle=int((~np.isfinite(v)).sum())))
        print(f"  précision {c:5.1f} mM : " + " · ".join(f"{p['lecture']} {p['ecart_type_robuste']:.1f}" for p in prec[-4:]))
    P = pd.DataFrame(prec)
    P.to_csv(D["results"] / "autofluorescence_tranches_precision.csv", index=False)
    json.dump(dict(K_mM=K, h20=h20, h_basal=hb, cl_basal_mM=c_basal, tissu_soma_median=a_med, tissu_champs_min_max=a_rng,
                   tissu_champs_basaux=a_bas, part_cfp_20mM=float(partB20), biais_predit_basal=pred,
                   biais_predit_basal_niveau_basal=pred_bas, NADH=NADH, FAD=FAD),
              open(D["results"] / "autofluorescence_tranches_monde.json", "w"), indent=2, ensure_ascii=False)

    # --- figure : même mise en page que la figure A.2
    fig = plt.figure(figsize=(15.5, 12.4))
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.98, top=0.9, bottom=0.22, hspace=0.4, wspace=0.18)
    tt = T - T[pk]
    sub = gs[0, 0].subgridspec(2, 1, hspace=0.55)
    yt_log = ([1e-3, 1e-2, 1e-1, 1], ["0,001", "0,01", "0,1", "1"])
    for j, (ch, titre_c, s0, s1, nom_lit) in enumerate((
            ("B", "(a) Formes des déclins · canal CFP", u0[:NW], u1[:NW], "NADH (littérature)"),
            ("Y", "canal YFP", u0[NW:], u1[NW:], "flavines (littérature)"))):
        ax = fig.add_subplot(sub[j])
        G.style(ax)
        for v, lab, col, ls in ((s0, "capteur, 0 mM", G.INK, "-"), (s1, "capteur, 131 mM", G.INK_2, "--"),
                                (tissu[ch], "tissu, mesuré", COUL["deux_auto"], "-"), (lit[ch], nom_lit, COUL["deux_auto"], ":")):
            ok = v > 1e-4 * v.max()
            ax.semilogy(tt[ok], v[ok] / v.max(), color=col, ls=ls, lw=1.6, label=lab)
        ax.set_ylim(1e-3, 1.5)
        ax.set_yticks(*yt_log)
        ax.minorticks_off()
        ax.set_xlim(-0.8, 10)
        ax.set_ylabel("normalisé (max = 1)", color=G.INK_2)
        if j == 1:
            ax.set_xlabel("temps après le pic (ns)", color=G.INK_2)
        ax.legend(frameon=False, fontsize=8.5, loc="lower right", ncol=2)
        G.titre(ax, titre_c, ("tissu : ce qui s'ajoute dans chaque pixel (étape 13), partie prompte (−0,3 à +0,6 ns)"
                              if j == 0 else None), fontsize=12)

    axp = fig.add_subplot(gs[0, 1])
    G.style(axp)
    for k in ("ratio", "mono", "deux", "deux_auto"):
        p = P[P.lecture == k]
        axp.semilogy(p.cl_vraie, p.ecart_type_robuste, marker="o", color=COUL[k], ls="--" if k == "deux" else "-",
                     lw=3 if k == "deux_auto" else 2, ms=4.5)
    yt = [v for v in (0.2, 0.5, 1, 2, 5, 10, 20, 50) if P.ecart_type_robuste.min() / 1.3 <= v <= P.ecart_type_robuste.max() * 1.3]
    axp.set_yticks(yt, [f"{v:g}".replace(".", ",") for v in yt])
    axp.minorticks_off()
    axp.set_xlabel("[Cl⁻] vraie (mM)", color=G.INK_2)
    axp.set_ylabel("écart-type de la [Cl⁻] lue (mM)", color=G.INK_2)
    G.titre(axp, f"(b) Le prix : précision à {args.N:.0e} photons CFP".replace("1e+04", "10⁴"),
            f"tissu du soma médian dans les deux canaux ; écart-type robuste de {args.tirages} tirages de Poisson", fontsize=12)

    def panneau(pos, sc, titre_p, part):
        ax = fig.add_subplot(pos)
        G.style(ax)
        ax.axhline(0, color=G.TRACE, lw=1)
        for k in ("deux_auto", "ratio", "deux", "deux_lit", "mono"):
            b = Bi[(Bi.scenario == sc) & (Bi.lecture == k)]
            kw = dict(color=COUL.get(k, COUL["deux_auto"]), lw=2)
            if k == "deux":
                kw.update(ls="--", lw=1.8)
            elif k == "deux_auto":
                kw.update(lw=3.2)
            elif k == "deux_lit":
                kw.update(ls=":", lw=2.4)
            elif k == "mono":
                kw.update(ls="--", lw=2.6, zorder=5)
            ax.plot(b.cl_vraie, b.biais, label=NOMS[k], **kw)
        ax.set_xlim(0, 80)
        ax.set_xlabel("[Cl⁻] vraie (mM)", color=G.INK_2)
        ax.set_ylabel("biais : [Cl⁻] lue − vraie (mM)", color=G.INK_2)
        G.titre(ax, titre_p, part, fontsize=12)
        return ax

    axc = panneau(gs[1, 0], "CFP", "(c) Signal du tissu dans le canal CFP",
                  f"soma médian : {fr(100 * a_med['B'], 2)} % des photons du soma ({fr(100 * a_med['B'] / partB20, 1)} % du CFP "
                  f"à 20 mM), indépendant du Cl⁻")
    panneau(gs[1, 1], "YFP", "(d) Signal du tissu dans le canal YFP",
            f"soma médian : {fr(100 * a_med['Y'], 2)} % des photons du soma ({fr(100 * a_med['Y'] / (1 - partB20), 1)} % du YFP "
            f"à 20 mM), indépendant du Cl⁻")
    h_, l_ = axc.get_legend_handles_labels()
    ordre = [l_.index(NOMS[k]) for k in ("ratio", "mono", "deux", "deux_auto", "deux_lit")]
    fig.legend([h_[i] for i in ordre], [l_[i] for i in ordre], loc="upper center", bbox_to_anchor=(0.5, 0.165), ncol=3,
               frameon=False, fontsize=10)
    import textwrap
    phrases = [
        f"Monde des tranches du 2026-11-09 : motifs du capteur seul (étape 13) à 0 et 131 mM, mélange h(c), K = {fr(K, 0)} mM "
        f"(motif à 20 mM) ; basal lu {fr(c_basal, 1)} mM.",
        f"Tissu : partie prompte de ce qui s'ajoute par pixel, soma médian des 20 champs (selon le champ : CFP "
        f"{fr(100 * a_rng['B'][0], 2)}–{fr(100 * a_rng['B'][1], 2)} %, YFP {fr(100 * a_rng['Y'][0], 2)}–"
        f"{fr(100 * a_rng['Y'][1], 2)} % des photons du soma) ; biais proportionnels à ce niveau.",
        "Lectures calibrées sur le capteur seul ; IRF et décalages connus. Pointillé : terme de forme NADH (CFP) / flavines (YFP) "
        "de la littérature au lieu de la forme mesurée.",
        f"Contrôle sur les somas réels (basal, tissu des deux canaux) : biais prédit du rapport {fr(pred['ratio'], 1, signe=True)} mM "
        f"({fr(pred_bas['ratio'], 1, signe=True)} au niveau des champs basaux), du τ du donneur {fr(pred['mono'], 1, signe=True)} mM "
        f"({fr(pred_bas['mono'], 1, signe=True)}) ; à l'étape 13, retirer le signal du tissu déplace le basal de −0,9 mM (rapport) "
        f"et de +3,5 à +3,7 mM (formes CFP)."]
    note = "\n".join(textwrap.fill(p_, 190) for p_ in phrases)
    fig.text(0.07, 0.012, note, fontsize=9, color=G.INK_2, va="bottom")
    for ax in fig.axes:
        ax.xaxis.set_major_formatter(VIRG)
        if ax.get_yscale() == "linear":
            ax.yaxis.set_major_formatter(VIRG)
    fig.suptitle("Biais sur [Cl⁻] dû au signal du tissu · tranches SCLM du 2026-11-09 · 850 nm", x=0.06, y=0.965,
                 ha="left", fontsize=14, color=G.INK)
    G.sauver(fig, D["results"] / "autofluorescence_tranches.png")


if __name__ == "__main__":
    main()
