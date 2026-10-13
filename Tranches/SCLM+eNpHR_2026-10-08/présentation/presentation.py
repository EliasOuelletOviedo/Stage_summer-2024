"""Figures de présentation, sur fond noir — écrites dans ce dossier, à côté de ce script.

Usage :  python3 presentation.py            toutes les figures
         python3 presentation.py retour     seulement celles dont le nom contient « retour »
                                            (plusieurs mots possibles : overlay ratio …)

Pour retoucher une figure, tout est dans ce fichier : le bloc RÉGLAGES ci-dessous
(couleurs, polices, tailles, et tous les textes : titres, titres d'axes, légendes,
notes), puis une fonction par figure. Les calculs, eux, viennent du pipeline
(Tranches/sclm_flim/ : étape 6 pour la dérive, étape 7 pour les retours) : on ne change
ici que l'apparence.

Le jeu de données est le dossier parent de celui-ci. Pour un autre jeu, copier ce fichier
dans son propre dossier `présentation/`.
"""
import json
import sys
from pathlib import Path

ICI = Path(__file__).resolve().parent              # <jeu>/présentation/
ROOT = ICI.parent                                  # <jeu>/
sys.path.insert(0, str(ROOT.parent / "sclm_flim"))  # Tranches/sclm_flim/

import warnings  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

import config  # noqa: E402
import slide_resultats as SR  # noqa: E402  (même dossier : données corrigées des slides)
import step6_ratio_corrige as step6  # noqa: E402
import step7_retour as step7  # noqa: E402


# =====================================================================================
# RÉGLAGES
# =====================================================================================

COULEURS = dict(
    fond="#000000",
    texte="#ffffff",             # titres, titres d'axes, traits principaux
    texte_2="#c3c2b7",           # sous-titres, notes, graduations
    grille="#2f2f2c",
    stims=["#3987e5", "#d95926", "#199e70", "#9085e9"],   # une couleur par stim (1 à 4)
    cfp="#3987e5",
    yfp="#c98500",
    ratio="#ffffff",             # moyenne glissante du ratio
    ratio_images="#c3c2b7",      # valeurs image par image, en fond
    bande_stim="#c9a227", alpha_stim=0.30,
    bande_exclue="#8d8d86", alpha_exclue=0.22,   # secondes écartées pour estimer la dérive
    derive="#d95926",
    ajustement="#ffffff",        # courbes ajustées des panneaux du bas (retour)
    soma="#3987e5",              # somas, figure de distribution
    bruit="#8d8d86",             # témoin « bruit seul »
    reperes="#4a4a46",           # quartiles, figure de distribution
    contour_cellpose="#ff3b30", numero_cellpose="#ffd60a",
    contour_flimapp="#00e5ff",
    alpha_points=0.8,
    duree_de_vie="#199e70",      # durée de vie de la CFP
    redescente="#d95926",        # redescentes ajustées de la durée de vie
    temoin="#8d8d86",            # témoins (fausses stims), figure de corrélation
    droite="#ffffff",            # droite de tendance, figure de corrélation
)
# Les 4 couleurs des stims ont été vérifiées sur fond noir : distinctes pour les
# daltoniens et contraste ≥ 3:1. En changer : garder des teintes bien séparées.

POLICE = dict(titre=13, sous_titre=10, axe=11, graduations=10, annotation=9, numeros=7)

FORMAT = dict(            # tailles des figures, en pouces (largeur, hauteur)
    dpi=200,
    image=(9, 9.6),         # images avec titre (IMAGES_NUES = False)
    image_nue=(8, 8),       # images nues : l'image remplit tout le PNG (8 po x 200 dpi = 1600 px)
    par_image=(11, 7.6),
    ratio_corrige=(11, 6.6),
    retour=(13, 8.4),
    distribution=(16, 6.4),
    redescente=(11, 4.6),
    correlation=(12.5, 5.4),
)

# Images de fond des figures « overlay ». `source` : "flim" (photons CFP + YFP de la
# session) ou une image Bliq déjà recalée dans le repère FLIM (projections/fovs.npz :
# 800_nm_avant_C2, 950_nm_avant_C2, 800_nm_après_C2, 950_nm_après_C2…).
# `titre` : titre sans ROI ; `suffixe` : ajouté au titre quand il y a des ROIs ;
# `note` : sous-titre quand il y a des ROIs ("" : aucun).
FONDS = {
    "flim": dict(source="flim", titre="Champ imagé — photons CFP + YFP sur toute la session",
                 suffixe="", note=""),
    "950_avant": dict(source="950_nm_avant_C2", titre="950 nm, avant la session (Bliq, canal C2)",
                      suffixe=" — sur le 950 nm d'avant",
                      note="image prise à 15:17, à un autre plan focal que la session :\n"
                           "les contours ne tombent pas exactement sur les mêmes cellules"),
}

IMAGES_NUES = True        # images « overlay » : seulement l'image, les contours et la barre —
                          # pas de titre, pas de marge, fond transparent. False : titre et marge noire.
BARRE_UM = 50             # barre d'échelle des images, en µm (0 : pas de barre)
NUMEROS_ROIS = True       # numéros sur les contours des images
LISSAGE_S = 1.0           # moyenne glissante des courbes image par image, en s
TAU_AXE = (2, 300)        # bornes de l'axe des τ (s, échelle log)

# Textes. Les {…} sont remplis par le script (valeurs disponibles indiquées à côté) ;
# pour une accolade littérale, l'écrire deux fois : {{ }}. « \n » = retour à la ligne.
TEXTES = {
    "image": dict(        # le titre sans ROI est celui du fond (FONDS ci-dessus)
        titre_cellpose="{n} somas segmentés (Cellpose)",                       # n
        titre_flimapp="{n_f} ROIs dessinées pendant l'expérience (FLIMApp)",  # n_f
        titre_les_deux="{n} somas (Cellpose, rouge) et {n_f} ROIs FLIMApp (bleu)",  # n, n_f
        barre="{um} µm",          # étiquette de la barre ("" : la barre seule)  # um
    ),
    "par_image": dict(
        titre="{etiquette}, image par image (60 ms)",                         # etiquette, n
        sous_titre="trait épais : moyenne glissante sur 1 s ; bandes jaunes : stimulations",
        y_cfp="CFP\nphotons / image",
        y_yfp="YFP\nphotons / image",
        y_ratio="ratio YFP / CFP",
        x="temps (s)",
    ),
    "ratio_corrige": dict(
        titre="{etiquette} — ratio corrigé de sa descente graduelle",
        sous_titre="dérive : β = {beta:+.2f} % par 100 s, estimée hors stims "
                   "(gris : {exclure:.0f} s écartées après chacune) ; tirets du bas : base avant la 1re stim",
        y_brut="ratio YFP / CFP\n(brut)",
        y_corrige="ratio YFP / CFP\n(corrigé)",
        legende_derive="dérive ajustée",
        x="temps (s)",
    ),
    "retour": dict(
        titre="{etiquette} — retour à l'équilibre après chaque stim : C + A·exp(−t/τ)",
        y_haut="ratio YFP / CFP\n(corrigé)",
        x_haut="temps (s)",
        tau_courbe="τ = {tau:.1f} s",                                          # tau
        panneau_stim="stim {i}",                                               # i
        panneau_moyenne="moyenne des 4",
        titre_panneau="{nom}\nτ = {tau:.1f} ± {se:.1f} s",                     # nom, tau, se
        note_panneau="A = {A:+.2f} %\néquilibre {C:+.2f} %",                  # A, C
        y_bas="ΔR / R avant la stim (%)",
        x_bas="s après la stim",
        pied="points : tranches de 1 s ; tirets : niveau des 10 s avant la stim ; "
             "pointillés : équilibre ajusté ; ± : erreur type",
    ),
    "distribution": dict(
        titre="Temps de remontée du ratio — somas Cellpose ; trait vertical (A) : "
              "somas réunis, {tau:.1f} ± {se:.1f} s",                          # tau, se
        a_titre="A · Chaque soma, moyenne de ses 4 retours",
        a_x="τ de remontée (s)",
        a_y="ROI",
        a_note="{n_mes} mesurables sur {n}\n({n_non} sans baisse nette)\ntrait : IC 95 %",
        b_titre="B · Plus étalé que le bruit de mesure ?",
        b_x="τ de remontée (s)",
        b_y="fraction des somas (cumulée)",
        b_bruit="bruit seul : tous à {tau:.1f} s",                             # tau
        b_obs="observé ({n} somas)",                                           # n
        b_note="médiane [Q1–Q3]\nobservé : {o1:.1f} [{o0:.1f}–{o2:.1f}] s\n"
               "bruit seul : {s1:.1f} [{s0:.1f}–{s2:.1f}] s\nKolmogorov-Smirnov p = {ks:.2f}",
        c_titre="C · Stim par stim",
        c_y="τ de remontée (s)",
        c_graduation="stim {i}\n{n}/{n_tot}",                                  # i, n, n_tot
        c_note="points et boîte : chaque soma (retour isolé, très bruité)\n"
               "◆ : somas réunis ± IC 95 %\nsous l'axe : retours mesurables",
    ),
    "redescente": dict(           # la figure seule, sans titre
        x="temps (s)",
        y="Δτ CFP (ps)",
    ),
    "correlation": dict(
        a_titre="À la ligne de base, d'une cellule à l'autre",
        a_x="ratio YFP / CFP, R₀",
        a_y="durée de vie de la CFP, τ₀ (ns)",
        a_note="ρ = {rho:+.2f}, p = {p:.2g}",                                  # rho, p
        b_titre="En réponse aux stimulations, par cellule",
        b_x="ΔR / R₀ (points, moyenne des stims)",
        b_y="Δτ CFP (ps, moyenne des stims)",
        b_stims="stims : ρ = {rho:+.2f}, p = {p:.2g}",                         # rho, p
        b_temoins="témoins (fausses stims) : ρ = {rho:+.2f}, p = {p:.2g}",
    ),
}


# =====================================================================================
# Outils communs
# =====================================================================================
C = COULEURS


def style(ax):
    ax.set_facecolor(C["fond"])
    ax.grid(True, color=C["grille"], linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C["grille"])
    ax.tick_params(colors=C["texte_2"], length=0, labelsize=POLICE["graduations"])


def axe_titre(ax, x=None, y=None):
    if x is not None:
        ax.set_xlabel(x, color=C["texte"], fontsize=POLICE["axe"])
    if y is not None:
        ax.set_ylabel(y, color=C["texte"], fontsize=POLICE["axe"])


def titre(ax, texte, sous_titre=None):
    lignes = sous_titre.count("\n") + 1 if sous_titre else 0
    ax.set_title(texte, color=C["texte"], fontsize=POLICE["titre"], loc="left",
                 pad=8 + lignes * 1.6 * POLICE["sous_titre"])
    if sous_titre:
        ax.text(0, 1.02, sous_titre, transform=ax.transAxes, color=C["texte_2"],
                fontsize=POLICE["sous_titre"], va="bottom")


def bandes_stim(ax, fen):
    for a, b in fen:
        ax.axvspan(a, b, color=C["bande_stim"], alpha=C["alpha_stim"], lw=0, zorder=0)


def ylim_robuste(ax, v, marge=0.12):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size > 10:
        lo, hi = np.percentile(v, (0.5, 99.5))
        ax.set_ylim(lo - marge * (hi - lo), hi + marge * (hi - lo))


def axe_log(ax, axe="x"):
    """Échelle log lisible : graduations 2, 5, 10, 20… en clair plutôt qu'en 10^n."""
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000) if TAU_AXE[0] <= t <= TAU_AXE[1]]
    a = ax.xaxis if axe == "x" else ax.yaxis
    (ax.set_xscale if axe == "x" else ax.set_yscale)("log")
    (ax.set_xlim if axe == "x" else ax.set_ylim)(*TAU_AXE)
    a.set_major_locator(FixedLocator(ticks))
    a.set_major_formatter(FixedFormatter([str(t) for t in ticks]))
    a.set_minor_locator(NullLocator())


def lisse(v, tm):
    k = max(int(round(LISSAGE_S / np.median(np.diff(tm)))), 1)
    return ndi.uniform_filter1d(np.asarray(v, float), k)


def sauver(fig, nom, nue=False):
    out = ICI / f"{nom}.png"
    if nue:       # l'axe occupe toute la figure : ni marge, ni fond
        fig.savefig(out, dpi=FORMAT["dpi"], transparent=True, pad_inches=0)
    else:
        fig.savefig(out, dpi=FORMAT["dpi"], facecolor=C["fond"], bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {out.name}")


# =====================================================================================
# Données (calculées par le pipeline)
# =====================================================================================
def charger():
    D, tm, fen, familles = step6.donnees(ROOT)
    _, etiquette, c1, c2 = next(x for x in familles if x[0] == "cellpose")
    info = json.loads((D["projections"] / "session.json").read_text())
    kept = json.loads((D["selection"] / "cellpose_selection.json").read_text())["kept"]
    fovs = list(info.get("fovs", {}).values())
    # Pixel FLIM en x : pixel Bliq (OME) divisé par l'échelle x mesurée au recalage.
    um_par_px = config.PIXEL_SIZE_UM / (fovs[0]["echelle_x"] if fovs else 1.0)
    return dict(D=D, tm=tm, fen=fen, etiquette=etiquette, c1=c1, c2=c2, kept=kept,
                n=len(kept), um_par_px=um_par_px,
                somme=np.load(D["projections"] / "flim.npz")["somme"].astype(float),
                fovs=dict(np.load(D["projections"] / "fovs.npz")) if (D["projections"] / "fovs.npz").exists() else {},
                masks=np.load(D["masks"] / "cellpose_masks.npy"),
                rois_flimapp=config.session(ROOT)["rois_flimapp"],
                sr=SR.charger(),                                     # blocs de 3 s corrigés (slides)
                par_roi=pd.read_csv(D["results"] / "resume_par_roi.csv"))


# =====================================================================================
# Figures
# =====================================================================================
def image(d, contours, nom, fond="flim"):
    """Image de fond (FONDS) ; contours = "aucun", "cellpose", "flimapp" ou "les_deux"."""
    T, F = TEXTES["image"], FONDS[fond]
    img = d["somme"] if F["source"] == "flim" else d["fovs"][F["source"]].astype(float)
    lo, hi = np.nanpercentile(img, (1, 99.7))
    if IMAGES_NUES:
        fig = plt.figure(figsize=FORMAT["image_nue"])
        ax = fig.add_axes([0, 0, 1, 1])
    else:
        fig, ax = plt.subplots(figsize=FORMAT["image"], facecolor=C["fond"])
    cmap = plt.get_cmap("gray").copy()
    cmap.set_bad(C["fond"], alpha=0 if IMAGES_NUES else 1)
    ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi)
    ax.set_axis_off()
    if contours in ("cellpose", "les_deux"):
        for lab in d["kept"]:
            m = d["masks"] == lab
            ax.contour(m, levels=[0.5], colors=C["contour_cellpose"], linewidths=1.2)
            if NUMEROS_ROIS:
                ys, xs = np.nonzero(m)
                ax.text(xs.mean(), ys.mean(), str(lab), color=C["numero_cellpose"],
                        fontsize=POLICE["numeros"], ha="center", va="center", fontweight="bold")
    if contours in ("flimapp", "les_deux"):
        for r in d["rois_flimapp"]:
            ax.plot(r["xs"], r["ys"], color=C["contour_flimapp"], lw=1.5, ls="--")
            if NUMEROS_ROIS:
                if contours == "les_deux":      # au-dessus du polygone : le centre porte déjà un numéro de soma
                    k = int(np.argmin(r["ys"]))
                    x, y, va = r["xs"][k], r["ys"][k] - 6, "bottom"
                else:
                    x, y, va = np.mean(r["xs"][:-1]), np.mean(r["ys"][:-1]), "center"
                ax.text(x, y, f"F{r['index']}", color=C["contour_flimapp"],
                        fontsize=POLICE["numeros"] + 1, ha="center", va=va, fontweight="bold")
    if contours == "aucun":
        texte, note = F["titre"], ""
    else:
        texte = T[f"titre_{contours}"].format(n=d["n"], n_f=len(d["rois_flimapp"])) + F["suffixe"]
        note = F["note"]
    if BARRE_UM:
        long = BARRE_UM / d["um_par_px"]
        x0, y0 = 0.04 * img.shape[1], 0.955 * img.shape[0]
        ax.plot([x0, x0 + long], [y0, y0], color=C["texte"], lw=4, solid_capstyle="butt")
        if T["barre"]:
            ax.text(x0 + long / 2, y0 - 0.015 * img.shape[0], T["barre"].format(um=BARRE_UM),
                    color=C["texte"], fontsize=POLICE["annotation"], ha="center", va="bottom")
    if not IMAGES_NUES:
        titre(ax, texte, note or None)
    sauver(fig, nom, nue=IMAGES_NUES)


def par_image(d, nom):
    """Photons CFP, photons YFP et ratio, image par image, somas réunis."""
    T = TEXTES["par_image"]
    tm, c1, c2, fen = d["tm"], d["c1"], d["c2"], d["fen"]
    fig, axes = plt.subplots(3, 1, figsize=FORMAT["par_image"], sharex=True, facecolor=C["fond"])
    for ax in axes:
        style(ax)
        bandes_stim(ax, fen)
    for ax, v, coul, lab in ((axes[0], c1, C["cfp"], T["y_cfp"]), (axes[1], c2, C["yfp"], T["y_yfp"])):
        ax.plot(tm, v, color=coul, lw=0.3, alpha=0.35)
        ax.plot(tm, lisse(v, tm), color=coul, lw=1.4)
        axe_titre(ax, y=lab)
        ylim_robuste(ax, v)
    r = c2 / np.maximum(c1, 1)
    axes[2].plot(tm, r, color=C["ratio_images"], lw=0.3, alpha=0.35)
    axes[2].plot(tm, lisse(c2, tm) / lisse(c1, tm), color=C["ratio"], lw=1.4)
    ylim_robuste(axes[2], r)
    axe_titre(axes[2], x=T["x"], y=T["y_ratio"])
    titre(axes[0], T["titre"].format(etiquette=d["etiquette"], n=d["n"]), T["sous_titre"])
    fig.tight_layout()
    sauver(fig, nom)


def ratio_corrige(d, nom):
    """Ratio brut avec la dérive ajustée, puis corrigé (méthode de l'étape 6)."""
    T = TEXTES["ratio_corrige"]
    tm, c1, c2, fen = d["tm"], d["c1"], d["c2"], d["fen"]
    r = c2 / np.maximum(c1, 1)
    fit = step6.ajuster_derive(tm, r, fen)
    rl = lisse(c2, tm) / lisse(c1, tm)
    corr = np.exp(-fit["beta"] * (tm - tm[0]))
    fig, axes = plt.subplots(2, 1, figsize=FORMAT["ratio_corrige"], sharex=True, facecolor=C["fond"])
    for ax in axes:
        style(ax)
        bandes_stim(ax, fen)
        for a, b in fen:
            ax.axvspan(b, b + fit["exclure"], color=C["bande_exclue"], alpha=C["alpha_exclue"], lw=0, zorder=0)
    axes[0].plot(tm, r, color=C["ratio_images"], lw=0.3, alpha=0.35)
    axes[0].plot(tm, rl, color=C["ratio"], lw=1.4)
    for k, (a, b, niv, _, _) in enumerate(fit["niveaux"]):
        x = np.array([max(a, tm[0]), min(b, tm[-1])])
        axes[0].plot(x, np.exp(niv + fit["beta"] * x), color=C["derive"], lw=2.0, ls="--",
                     label=T["legende_derive"] if k == 0 else None, zorder=4)
    axes[0].legend(frameon=False, fontsize=POLICE["annotation"], loc="lower left", labelcolor=C["texte"])
    ylim_robuste(axes[0], r)
    axe_titre(axes[0], y=T["y_brut"])
    rc = r * corr
    axes[1].plot(tm, rc, color=C["ratio_images"], lw=0.3, alpha=0.35)
    axes[1].plot(tm, rl * corr, color=C["ratio"], lw=1.4)
    if fen:
        axes[1].axhline(np.median(rc[tm < fen[0][0]]), color=C["texte_2"], lw=0.9, ls="--", zorder=2)
    ylim_robuste(axes[1], rc)
    axe_titre(axes[1], x=T["x"], y=T["y_corrige"])
    titre(axes[0], T["titre"].format(etiquette=d["etiquette"], n=d["n"]),
          T["sous_titre"].format(beta=fit["beta"] * 1e4, exclure=fit["exclure"]))
    fig.tight_layout()
    sauver(fig, nom)


def retour(R, nom):
    """Retour à l'équilibre des somas réunis : les 4 ajustements, puis chacun en détail."""
    T = TEXTES["retour"]
    tm, c1, c2, fen, beta = R["tm"], R["pc"][:, 0], R["pc"][:, 1], R["fen"], R["beta"]
    pal = C["stims"]
    fig = plt.figure(figsize=FORMAT["retour"], facecolor=C["fond"])
    gs = fig.add_gridspec(2, len(fen) + 1, height_ratios=[1.25, 1], hspace=0.45, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    style(ax)
    bandes_stim(ax, fen)
    corr = np.exp(-beta * (tm - tm[0]))
    rc = c2 / np.maximum(c1, 1) * corr
    ax.plot(tm, rc, color=C["ratio_images"], lw=0.3, alpha=0.3)
    ax.plot(tm, lisse(c2, tm) / lisse(c1, tm) * corr, color=C["ratio"], lw=1.2)
    for i, ((a, b), (xc, y, r_av), p) in enumerate(zip(fen, R["rets"], R["fits"])):
        xx = np.linspace(step7.GARDE_S, xc[-1], 300)
        coul = pal[i % len(pal)]
        ax.plot(b + xx, r_av * (1 + p["C"] + p["A"] * np.exp(-xx / p["tau"])), color=coul, lw=2.6, zorder=5)
        ax.text(b + 6, r_av * (1 + p["A"] * 1.35), T["tau_courbe"].format(tau=p["tau"]), color=coul,
                fontsize=POLICE["annotation"] + 1, fontweight="bold", va="top")
    ylim_robuste(ax, rc, marge=0.1)
    axe_titre(ax, x=T["x_haut"], y=T["y_haut"])
    titre(ax, T["titre"].format(etiquette=R["etiquette"], n=len(R["kept"])))

    panneaux = [(T["panneau_stim"].format(i=i + 1), r[0], r[1], p, pal[i % len(pal)])
                for i, (r, p) in enumerate(zip(R["rets"], R["fits"]))]
    panneaux.append((T["panneau_moyenne"], R["rets"][0][0],
                     np.nanmean([r[1] for r in R["rets"]], axis=0), R["moyen"], C["texte"]))
    lims = np.nanpercentile(np.concatenate([p[2] for p in panneaux]) * 100, (1, 99))
    for j, (nom_p, xc, y, p, coul) in enumerate(panneaux):
        a2 = fig.add_subplot(gs[1, j])
        style(a2)
        a2.scatter(xc, y * 100, s=7, color=coul, alpha=C["alpha_points"], lw=0)
        xx = np.linspace(step7.GARDE_S, xc[-1], 300)
        a2.plot(xx, (p["C"] + p["A"] * np.exp(-xx / p["tau"])) * 100, color=C["ajustement"], lw=1.8)
        a2.axhline(0, color=C["texte_2"], lw=0.8, ls="--")
        a2.axhline(p["C"] * 100, color=coul, lw=0.8, ls=":")
        a2.set_ylim(lims[0] - 0.4, lims[1] + 0.4)
        a2.set_title(T["titre_panneau"].format(nom=nom_p, tau=p["tau"], se=p["se_tau"]),
                     color=C["texte"], fontsize=POLICE["axe"])
        a2.text(0.97, 0.05, T["note_panneau"].format(A=p["A"] * 100, C=p["C"] * 100),
                transform=a2.transAxes, ha="right", va="bottom", fontsize=POLICE["annotation"] - 1,
                color=C["texte_2"])
        axe_titre(a2, x=T["x_bas"], y=T["y_bas"] if j == 0 else None)
    fig.text(0.01, 0.005, T["pied"], color=C["texte_2"], fontsize=POLICE["annotation"])
    sauver(fig, nom)


def distribution(R, nom):
    """Distribution des temps de remontée : par soma, contre le bruit seul, stim par stim."""
    T = TEXTES["distribution"]
    df, sims, ks = R["df"], R["sims"], R["ks"]
    tau_c, se_c = R["moyen"]["tau"], R["moyen"]["se_tau"]
    roi = df[(df["niveau"] == "roi") & (df["stim"] == "moyenne")]
    mes = roi[roi["mesurable"]].sort_values("tau_s")
    som = df[df["niveau"] == "somas réunis"]
    fig, axes = plt.subplots(1, 3, figsize=FORMAT["distribution"], facecolor=C["fond"],
                             gridspec_kw=dict(width_ratios=[1.15, 1, 1]))

    # A — chaque soma, moyenne de ses retours, avec son intervalle à 95 %
    ax = axes[0]
    style(ax)
    lo_c, hi_c = step7.ic95(tau_c, se_c)
    ax.axvspan(lo_c, hi_c, color=C["texte"], alpha=0.10, lw=0)
    ax.axvline(tau_c, color=C["texte"], lw=1.4)
    lo, hi = step7.ic95(mes["tau_s"].to_numpy(), mes["tau_se_s"].to_numpy())
    yy = np.arange(len(mes))
    ax.hlines(yy, np.clip(lo, *TAU_AXE), np.clip(hi, *TAU_AXE), color=C["soma"], lw=1.6, alpha=0.6)
    ax.scatter(mes["tau_s"], yy, color=C["soma"], s=26, zorder=3)
    axe_log(ax, "x")
    ax.set_yticks(yy)
    ax.set_yticklabels(mes["roi"], fontsize=POLICE["graduations"] - 2.5, color=C["texte_2"])
    axe_titre(ax, x=T["a_x"], y=T["a_y"])
    ax.set_title(T["a_titre"], color=C["texte"], fontsize=POLICE["axe"], loc="left")
    ax.text(0.98, 0.02, T["a_note"].format(n_mes=len(mes), n=len(roi), n_non=len(roi) - len(mes)),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=POLICE["annotation"], color=C["texte_2"])

    # B — observé contre bruit seul
    ax = axes[1]
    style(ax)
    obs, sim = np.sort(mes["tau_s"].to_numpy()), np.sort(sims)
    ax.step(sim, np.arange(1, len(sim) + 1) / len(sim), where="post", color=C["bruit"], lw=2.2,
            label=T["b_bruit"].format(tau=tau_c))
    ax.step(obs, np.arange(1, len(obs) + 1) / len(obs), where="post", color=C["soma"], lw=2.2,
            label=T["b_obs"].format(n=len(obs)))
    for q, ls in ((0.25, ":"), (0.5, "--"), (0.75, ":")):
        ax.axhline(q, color=C["reperes"], lw=0.8, ls=ls)
    axe_log(ax, "x")
    ax.set_ylim(0, 1.02)
    axe_titre(ax, x=T["b_x"], y=T["b_y"])
    ax.set_title(T["b_titre"], color=C["texte"], fontsize=POLICE["axe"], loc="left")
    qo, qs = np.percentile(obs, [25, 50, 75]), np.percentile(sim, [25, 50, 75])
    ax.text(0.98, 0.02, T["b_note"].format(o0=qo[0], o1=qo[1], o2=qo[2], s0=qs[0], s1=qs[1], s2=qs[2], ks=ks),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=POLICE["annotation"], color=C["texte_2"])
    ax.legend(frameon=False, fontsize=POLICE["annotation"], loc="upper left", labelcolor=C["texte"])

    # C — stim par stim
    ax = axes[2]
    style(ax)
    rs = df[(df["niveau"] == "roi") & (df["stim"] != "moyenne") & df["mesurable"]]
    rng = np.random.default_rng(1)
    stims = sorted(rs["stim"].unique(), key=int)
    graduations = []
    for i in map(int, stims):
        v = rs.loc[rs["stim"] == str(i), "tau_s"].to_numpy()
        coul = C["stims"][(i - 1) % len(C["stims"])]
        if len(v) >= 3:
            ax.boxplot([v], positions=[i - 0.08], widths=0.32, showfliers=False, whis=0,
                       showcaps=False, patch_artist=True, manage_ticks=False,
                       boxprops=dict(facecolor=coul, alpha=0.22, edgecolor=coul),
                       medianprops=dict(color=coul, lw=2.2), whiskerprops=dict(lw=0))
        ax.scatter(np.full(len(v), i - 0.08) + rng.uniform(-0.13, 0.13, len(v)), v, s=14,
                   color=coul, alpha=C["alpha_points"], lw=0, zorder=3)
        p = som[som["stim"] == str(i)].iloc[0]
        lo, hi = step7.ic95(p["tau_s"], p["tau_se_s"])
        ax.errorbar(i + 0.3, p["tau_s"], yerr=[[p["tau_s"] - lo], [hi - p["tau_s"]]], fmt="D",
                    color=C["texte"], ms=6, capsize=3, zorder=5)
        n_tot = ((df["niveau"] == "roi") & (df["stim"] == str(i))).sum()
        graduations.append(T["c_graduation"].format(i=i, n=len(v), n_tot=n_tot))
    axe_log(ax, "y")
    ax.set_xlim(0.5, len(stims) + 0.6)
    ax.set_xticks([int(i) for i in stims])
    ax.set_xticklabels(graduations, color=C["texte"], fontsize=POLICE["graduations"] - 1)
    axe_titre(ax, y=T["c_y"])
    ax.set_title(T["c_titre"], color=C["texte"], fontsize=POLICE["axe"], loc="left")
    ax.text(0.98, 0.98, T["c_note"], transform=ax.transAxes, ha="right", va="top",
            fontsize=POLICE["annotation"], color=C["texte_2"])
    fig.suptitle(T["titre"].format(tau=tau_c, se=se_c), color=C["texte"], fontsize=POLICE["titre"],
                 x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    sauver(fig, nom)


def ajuster_redescentes(t, med, fen):
    """Après chaque stimulation, C_k + A_k·exp(−(t − fin_k)/τ) jusqu'à la stimulation
    suivante : un τ commun, une amplitude et un équilibre par stim, moindres carrés.
    Une à une, les redescentes de la durée de vie sont trop bruitées pour être ajustées
    (τ de 29 ± 122 s à la stim 1, à la borne de 300 s à la stim 4) ; ensemble, si.
    Rend τ, son erreur type, et [(A_k, C_k, fin_k, premier, dernier temps ajusté)]."""
    segs = []
    for i, (_, b) in enumerate(fen):
        nxt = fen[i + 1][0] if i + 1 < len(fen) else t[-1] + 1
        m = (t > b) & (t < nxt) & np.isfinite(med)
        segs.append((b, t[m] - b, med[m]))

    def residus(p):
        return np.concatenate([p[2 + 2 * k] + p[1 + 2 * k] * np.exp(-x / p[0]) - y
                               for k, (_, x, y) in enumerate(segs)])

    p0 = [15.0] + [v for _, x, y in segs for v in (y[:3].mean() - y[-10:].mean(), y[-10:].mean())]
    n = len(segs)
    r = least_squares(residus, p0, bounds=([0.3] + [-200, -100] * n, [300] + [200, 100] * n))
    s2 = (r.fun ** 2).sum() / max(len(r.fun) - len(r.x), 1)
    se = np.sqrt(np.diag(np.linalg.inv(r.jac.T @ r.jac) * s2))
    return r.x[0], se[0], [(r.x[1 + 2 * k], r.x[2 + 2 * k], b, x.min(), x.max())
                           for k, (b, x, _) in enumerate(segs)]


def redescente(d, nom):
    """Durée de vie de la CFP sur la session, somas (médiane ± écart interquartile, corrigée
    de la descente graduelle comme dans les slides), blocs de stim compris ; et les
    redescentes après chaque stim, ajustées hors stim sur la médiane non lissée."""
    T, s = TEXTES["redescente"], d["sr"]
    t, fen = s["t"], s["fen"]
    fig, ax = plt.subplots(figsize=FORMAT["redescente"], facecolor=C["fond"])
    style(ax)
    bandes_stim(ax, fen)
    M = SR.lisser(s["dT_tout"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q1, med, q3 = (np.nanpercentile(M, q, axis=0) for q in (25, 50, 75))
        med_fit = np.nanmedian(s["dT"], axis=0)
    ax.fill_between(t, q1, q3, color=C["duree_de_vie"], alpha=0.25, lw=0, zorder=2)
    ax.plot(t, med, color=C["duree_de_vie"], lw=2.2, zorder=3)
    ax.axhline(0, color=C["texte_2"], lw=0.9, ls="--", zorder=1)
    tau, se, courbes = ajuster_redescentes(t, med_fit, fen)
    for A_, C_, b, x0, x1 in courbes:
        xx = np.linspace(x0, x1, 300)
        ax.plot(b + xx, C_ + A_ * np.exp(-xx / tau), color=C["redescente"], lw=2.2, zorder=4)
    ylim_robuste(ax, np.r_[q1, q3], marge=0.08)
    axe_titre(ax, x=T["x"], y=T["y"])
    sauver(fig, nom)
    print(f"     redescentes : τ commun = {tau:.1f} ± {se:.1f} s ; "
          f"A par stim : {', '.join(f'{a:+.1f}' for a, *_ in courbes)} ps")


def correlation(d, nom):
    """Ratio et durée de vie de la CFP, deux lectures du même FRET : à la ligne de base,
    d'une cellule à l'autre ; puis en réponse aux stims, par cellule (moyenne des stims),
    avec les témoins des mêmes cellules (étape 5)."""
    T = TEXTES["correlation"]
    c = d["par_roi"]
    c = c[c["famille"] == "cellpose"]
    fig, axes = plt.subplots(1, 2, figsize=FORMAT["correlation"], facecolor=C["fond"])

    def droite(ax, x, y):
        res = stats.theilslopes(y, x)
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 50)
        ax.plot(xx, res.intercept + res.slope * xx, color=C["droite"], lw=1.6, zorder=2)

    ax = axes[0]
    style(ax)
    x, y = c["R0"].to_numpy(), c["tau0_ch1_ns"].to_numpy()
    ax.scatter(x, y, s=38, color=C["soma"], alpha=C["alpha_points"], lw=0, zorder=3)
    droite(ax, x, y)
    rho, p = stats.spearmanr(x, y)
    ax.text(0.97, 0.95, T["a_note"].format(rho=rho, p=p), transform=ax.transAxes, ha="right",
            va="top", color=C["texte_2"], fontsize=POLICE["annotation"] + 1)
    axe_titre(ax, x=T["a_x"], y=T["a_y"])
    ax.set_title(T["a_titre"], color=C["texte"], fontsize=POLICE["axe"], loc="left")

    ax = axes[1]
    style(ax)
    ax.axhline(0, color=C["texte_2"], lw=0.8, ls="--", zorder=1)
    ax.axvline(0, color=C["texte_2"], lw=0.8, ls="--", zorder=1)
    xt, yt = c["dR_R0_pct_ecart_temoins"].to_numpy(), c["dtau_ch1_ps_ecart_temoins"].to_numpy()
    xs, ys = c["dR_R0_pct_ecart_stims"].to_numpy(), c["dtau_ch1_ps_ecart_stims"].to_numpy()
    rho, p = stats.spearmanr(xs, ys)
    rho_t, p_t = stats.spearmanr(xt, yt)
    ax.scatter(xt, yt, s=28, color=C["temoin"], alpha=0.7, lw=0, zorder=2,
               label=T["b_temoins"].format(rho=rho_t, p=p_t))
    ax.scatter(xs, ys, s=38, color=C["soma"], alpha=C["alpha_points"], lw=0, zorder=3,
               label=T["b_stims"].format(rho=rho, p=p))
    droite(ax, xs, ys)
    ax.legend(frameon=False, fontsize=POLICE["annotation"] + 1, loc="lower left", labelcolor=C["texte"])
    axe_titre(ax, x=T["b_x"], y=T["b_y"])
    ax.set_title(T["b_titre"], color=C["texte"], fontsize=POLICE["axe"], loc="left")
    fig.tight_layout()
    sauver(fig, nom)
    print(f"     base : ρ = {stats.spearmanr(x, y)[0]:+.2f} ; réponses : stims ρ = {rho:+.2f} "
          f"(p = {p:.2g}), témoins ρ = {rho_t:+.2f} (p = {p_t:.2g})")


# =====================================================================================
# Les figures produites : nom du fichier (sans .png) -> (fonction, données)
# Commenter une ligne pour ne plus la produire.
# =====================================================================================
FIGURES = {
    "cellpose_overlay_sans_roi":      (lambda d, R, n: image(d, "aucun", n), "d"),
    "cellpose_overlay_rois_cellpose": (lambda d, R, n: image(d, "cellpose", n), "d"),
    "cellpose_overlay_rois_flimapp":  (lambda d, R, n: image(d, "flimapp", n), "d"),
    "cellpose_overlay_rois_les_deux": (lambda d, R, n: image(d, "les_deux", n), "d"),
    "cellpose_overlay_950nm_avant_sans_roi":      (lambda d, R, n: image(d, "aucun", n, "950_avant"), "d"),
    "cellpose_overlay_950nm_avant_rois_cellpose": (lambda d, R, n: image(d, "cellpose", n, "950_avant"), "d"),
    "cellpose_overlay_950nm_avant_rois_flimapp":  (lambda d, R, n: image(d, "flimapp", n, "950_avant"), "d"),
    "cellpose_overlay_950nm_avant_rois_les_deux": (lambda d, R, n: image(d, "les_deux", n, "950_avant"), "d"),
    "par_image_cellpose":             (lambda d, R, n: par_image(d, n), "d"),
    "ratio_corrige_cellpose":         (lambda d, R, n: ratio_corrige(d, n), "d"),
    "retour_equilibre_cellpose":      (lambda d, R, n: retour(R, n), "R"),
    "temps_de_remontee":              (lambda d, R, n: distribution(R, n), "R"),
    "duree_de_vie_redescente":        (lambda d, R, n: redescente(d, n), "d"),
    "correlation_ratio_duree_de_vie": (lambda d, R, n: correlation(d, n), "d"),
}


def main():
    filtres = [a for a in sys.argv[1:] if not a.startswith("-")]
    choix = {k: v for k, v in FIGURES.items() if not filtres or any(f in k for f in filtres)}
    if not choix:
        sys.exit(f"Aucune figure ne correspond à {filtres}. Noms : {', '.join(FIGURES)}")
    print(f"{ROOT.name} — figures de présentation ({len(choix)})")
    d = charger()
    # Les retours (étape 7) ne sont recalculés que si une figure en a besoin.
    R = step7.calculer(ROOT) if any(v[1] == "R" for v in choix.values()) else None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for nom, (f, _) in choix.items():
            f(d, R, nom)
    print(f"\nFigures -> {ICI}")


if __name__ == "__main__":
    main()
