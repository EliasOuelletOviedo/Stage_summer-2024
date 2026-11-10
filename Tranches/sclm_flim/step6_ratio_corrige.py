"""Étape 6 — ratio ch2 / ch1 image par image, corrigé de sa descente graduelle.

Usage :  python3 step6_ratio_corrige.py <dossier>   (après les étapes 1 à 4)

Deux figures au format du panneau c de `qc/par_image.png` (une valeur par image de
60 ms, moyenne glissante sur 1 s) : le champ entier, puis les somas Cellpose retenus à
l'étape 3 (photons de tous leurs pixels réunis).

Correction. Entre les stims, le ratio descend lentement ; à chaque stim, il chute puis
remonte en partie pendant une vingtaine de secondes. La descente graduelle est estimée
là où aucune stim n'agit : la base avant la 1re stim et, après chaque stim, de
EXCLURE_APRES_STIM_S secondes après sa fin jusqu'au début de la suivante. Modèle :

    ln(ratio) = a_k + β · t        une pente β commune, un niveau a_k par intervalle

Les niveaux libres laissent les marches des stims hors de la pente : β ne mesure que la
descente à l'intérieur des intervalles. Le ratio corrigé vaut ratio · exp(−β · t), le
temps compté depuis la première image (il garde donc l'échelle du ratio brut au départ).

Un ajustement lisse sur toute la session (exponentielle passant par toutes les bases)
aurait absorbé les marches : il aurait fabriqué, après chaque stim, un faux retour vers
la base. C'est ce que ce modèle évite.

Chaque figure a sa propre pente : les somas ne dérivent pas comme le champ entier.

Sorties : results/ratio_corrige_champ.png, results/ratio_corrige_cellpose.png,
results/derive_ratio.csv, results/reponse_ratio_corrige.csv,
results/ratio_corrige_par_image.csv
"""
import json

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

import config
from config import EXCLURE_APRES_STIM_S, PRE_WINDOW_S, RESPONSE_WINDOW_S
from figures import CLAIR, LISSAGE, bandes_stim, plt, sauver, style, ylim_robuste


# ------------------------------------------------------------------ dérive
def intervalles(tm, fen, exclure):
    """Fenêtres sans effet de stim : base, puis de fin + exclure à la stim suivante."""
    fin = tm[-1] + 1.0
    if not fen:
        return [(tm[0] - 1.0, fin)]
    iv = [(tm[0] - 1.0, fen[0][0])]
    for i, (a, b) in enumerate(fen):
        suivant = fen[i + 1][0] if i + 1 < len(fen) else fin
        iv.append((b + exclure, suivant))
    return iv


def ajuster_derive(tm, ratio, fen, exclure=EXCLURE_APRES_STIM_S):
    """β commun et niveaux a_k de ln(ratio) = a_k + β t, par moindres carrés intra-intervalle."""
    lr = np.log(ratio)
    groupes = []
    for a, b in intervalles(tm, fen, exclure):
        m = (tm >= a) & (tm < b) & np.isfinite(lr)
        if m.sum() >= 20:
            groupes.append((a, b, m))
    num = den = 0.0
    for a, b, m in groupes:
        x, y = tm[m] - tm[m].mean(), lr[m] - lr[m].mean()
        num += (x * y).sum()
        den += (x * x).sum()
    beta = num / den
    niveaux = [(a, b, float(lr[m].mean() - beta * tm[m].mean()),
                float(np.polyfit(tm[m], lr[m], 1)[0]), int(m.sum())) for a, b, m in groupes]
    return {"beta": float(beta), "niveaux": niveaux, "exclure": exclure}


def reponses(tm, rc, fen):
    """Ratio corrigé autour de chaque stim : 10 s avant, 10 s après, et 10 s avant la
    stim suivante (la fin de la session pour la dernière), en % du niveau d'avant."""
    rows = []
    for i, (a, b) in enumerate(fen):
        suivant = fen[i + 1][0] if i + 1 < len(fen) else tm[-1] + 1e-9
        av = rc[(tm >= a - PRE_WINDOW_S) & (tm < a)].mean()
        ap = rc[(tm > b) & (tm <= b + RESPONSE_WINDOW_S)].mean()
        fi = rc[(tm >= suivant - PRE_WINDOW_S) & (tm < suivant)].mean()
        rows.append({"stim": i + 1, "avant": av, "apres_10s": ap, "avant_suivante": fi,
                     "chute_pct": (ap / av - 1) * 100, "reste_pct": (fi / av - 1) * 100,
                     "delai_s": round(suivant - b, 1)})
    return rows


# ------------------------------------------------------------------ figure
def figure(tm, c1, c2, fen, fit, etiquette, out, T=CLAIR):
    lis = lambda v: ndi.uniform_filter1d(v.astype(float), LISSAGE)
    r = c2 / np.maximum(c1, 1)
    rl = lis(c2) / lis(c1)
    corr = np.exp(-fit["beta"] * (tm - tm[0]))
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True, facecolor=T["surface"])
    for ax in axes:
        style(ax, T)
        bandes_stim(ax, fen, T)
        for a, b in fen:
            ax.axvspan(b, b + fit["exclure"], color=T["neutre"], alpha=T["exclu_alpha"], lw=0, zorder=0)

    axes[0].plot(tm, r, color=T["ink_2"], lw=0.3, alpha=0.35)
    axes[0].plot(tm, rl, color=T["ink"], lw=1.4)
    for k, (a, b, niv, _, _) in enumerate(fit["niveaux"]):
        x = np.array([max(a, tm[0]), min(b, tm[-1])])
        axes[0].plot(x, np.exp(niv + fit["beta"] * x), color=T["palette"][1], lw=2.0, ls="--",
                     label="dérive ajustée" if k == 0 else None, zorder=4)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left", labelcolor=T["ink"])
    ylim_robuste(axes[0], r)
    axes[0].set_ylabel(f"ratio ch2 / ch1\n(brut, {etiquette})", color=T["ink"], fontsize=10)

    rc, rcl = r * corr, rl * corr
    axes[1].plot(tm, rc, color=T["ink_2"], lw=0.3, alpha=0.35)
    axes[1].plot(tm, rcl, color=T["ink"], lw=1.4)
    if fen:
        base = np.median(rc[tm < fen[0][0]])
        axes[1].axhline(base, color=T["ink_2"], lw=0.9, ls="--", zorder=2)
    ylim_robuste(axes[1], rc)
    axes[1].set_ylabel(f"ratio ch2 / ch1\n(corrigé, {etiquette})", color=T["ink"], fontsize=10)
    axes[1].set_xlabel("temps (s, horloge des cartes SPC)", color=T["ink"], fontsize=10.5)

    pentes = [p * 1e4 for *_, p, _ in fit["niveaux"]]
    axes[0].set_title(f"{etiquette[0].upper() + etiquette[1:]}, image par image (60 ms) — ratio corrigé de sa "
                      "descente graduelle", color=T["ink"], fontsize=12, loc="left", pad=36)
    axes[0].text(0, 1.02,
                 f"dérive : exp(β·t), β = {fit['beta'] * 1e4:+.2f} % par 100 s (de {min(pentes):+.2f} à "
                 f"{max(pentes):+.2f} selon l'intervalle), estimée hors stims en écartant les "
                 f"{fit['exclure']:.0f} s qui suivent chacune (gris)\n"
                 "un niveau par intervalle, pour ne pas absorber les marches ; trait épais : moyenne "
                 "glissante sur 1 s ; tirets du bas : base avant la 1re stim",
                 transform=axes[0].transAxes, color=T["ink_2"], fontsize=9, va="bottom", linespacing=1.5)
    fig.tight_layout()
    sauver(fig, out, dpi=180, T=T)
    return rc


# ------------------------------------------------------------------ main
def donnees(root):
    """(dossiers, temps des images, fenêtres de stim, familles) — chaque famille est
    (clé, étiquette, photons ch1, photons ch2) image par image : champ entier, puis les
    somas retenus réunis. La dernière image, coupée par l'arrêt, est écartée."""
    D = config.dirs(root)
    info = json.loads((D["projections"] / "session.json").read_text())
    f = np.load(D["projections"] / "flim.npz")
    d = np.load(D["mesures"] / "declins.npz")
    if "par_image_cellpose" not in d.files:
        raise SystemExit("mesures/declins.npz sans photons image par image : relance step4_declins.py")
    kept = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
    t = f["trames_t"]
    n = len(t) - 1
    tm = (t + np.median(np.diff(t)) / 2)[:n]
    fen = [(s["t0"], s["t1"]) for s in info["stims"]]
    pc = d["par_image_cellpose"][kept].sum(0)          # [n_tr, 2], somas retenus réunis
    familles = [("champ", "champ entier", f["par_trame_ch1"][:n].astype(float),
                 f["par_trame_ch2"][:n].astype(float)),
                ("cellpose", f"{len(kept)} somas Cellpose réunis", pc[:n, 0], pc[:n, 1])]
    return D, tm, fen, familles


def main():
    root = config.racine()
    D, tm, fen, familles = donnees(root)
    if fen and fen[0][0] < PRE_WINDOW_S:
        print(f"{root.name} — pas de ligne de base (1re stim à {fen[0][0]:.1f} s) : "
              "la descente graduelle ne peut pas être séparée des réponses ; étape 6 sautée.")
        return

    rows_fit, rows_rep = [], []
    serie = pd.DataFrame({"t_s": np.round(tm, 4),
                          "stim": np.any([(tm >= a) & (tm < b) for a, b in fen], axis=0) if fen else False})
    print(f"{root.name} — {len(fen)} stims, pente estimée hors des {EXCLURE_APRES_STIM_S:.0f} s "
          "qui suivent chaque stim")
    for cle, etiquette, c1, c2 in familles:
        c1, c2 = c1.astype(float), c2.astype(float)
        fit = ajuster_derive(tm, c2 / np.maximum(c1, 1), fen)
        rc = figure(tm, c1, c2, fen, fit, etiquette, D["results"] / f"ratio_corrige_{cle}.png")
        serie[f"ratio_{cle}"] = c2 / np.maximum(c1, 1)
        serie[f"ratio_{cle}_corrige"] = rc
        for k, (a, b, niv, p, nn) in enumerate(fit["niveaux"]):
            rows_fit.append({"famille": cle, "intervalle": k, "t0_s": round(max(a, tm[0]), 2),
                             "t1_s": round(min(b, tm[-1]), 2), "n_images": nn,
                             "pente_propre_pct_100s": p * 1e4,
                             "pente_commune_pct_100s": fit["beta"] * 1e4,
                             "exclure_apres_stim_s": fit["exclure"]})
        rep = reponses(tm, rc, fen)
        for x in rep:
            rows_rep.append({"famille": cle, **x})
        print(f"\n  {etiquette} : β = {fit['beta'] * 1e4:+.2f} %/100 s, soit "
              f"{(np.exp(fit['beta'] * (tm[-1] - tm[0])) - 1) * 100:+.1f} % sur la session")
        for x in rep:
            print(f"    stim {x['stim']} : 10 s après {x['chute_pct']:+.2f} %, "
                  f"{x['delai_s']:.0f} s après la fin {x['reste_pct']:+.2f} % (ratio corrigé, "
                  "par rapport aux 10 s d'avant)")
    pd.DataFrame(rows_fit).to_csv(D["results"] / "derive_ratio.csv", index=False)
    pd.DataFrame(rows_rep).to_csv(D["results"] / "reponse_ratio_corrige.csv", index=False)
    serie.to_csv(D["results"] / "ratio_corrige_par_image.csv", index=False)
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
