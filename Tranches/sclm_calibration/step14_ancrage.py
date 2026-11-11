"""Étape 14 — le rapport seul contre le rapport ancré par la forme du déclin, en validation croisée par champ.

Usage :  python3 step14_ancrage.py <dossier>        (~1 min ; après l'étape 13)

Le rapport YFP/CFP est précis mais son échelle dépend de η (transmission du jaune relative au bleu), qui peut changer
d'une tranche, d'un réglage ou d'une séance à l'autre ; la forme du déclin CFP n'en dépend pas. Ancrage : pour chaque
champ, la forme des photons réunis donne la part f du motif « 131 mM », d'où le rapport attendu R(f) ; η̂ = rapport
mesuré / R(f), et chaque soma du champ est lu à R / η̂.

Calibration, sur les champs d'entraînement seulement (le k-ième champ de chaque condition imposée est lu, pli k) :
  - ancrage, sur les deux conditions au clamp fiable (0 et 131 mM) : gabarits CFP = motif du capteur seul (pente de la
    régression entre somas de l'étape 13, composante prompte du tissu retirée), moyenne des champs ; droite
    R(f) = R₀·(1 − f) + R₁·f (capteur à deux états) entre leurs f et leurs rapports du capteur seul ;
  - conversion en mM, la même pour les deux lectures : isotherme sur les médianes des somas d'entraînement à 0, 20 et
    131 mM (rapports bruts pour le rapport seul ; divisés par le η̂ de leur champ pour le rapport ancré).

Deux situations pour le champ lu (et pour les champs du basal, lus par chaque pli) :
  - données telles quelles : ce que coûte l'ancrage quand η ne change pas ;
  - η changé (0,7, 0,85, 1,15, 1,3) : canal jaune aminci (η < 1) ou canal bleu aminci (η > 1) par tirage binomial,
    photon par photon, puis régression de l'étape 13 refaite sur le champ modifié — une autre tranche ou un autre
    réglage, simulé ; la forme des déclins ne change pas, seule l'échelle du rapport.

Critères : erreur de la médiane de chaque champ lu par rapport à la [Cl⁻] imposée (biais et erreur quadratique sur
les 5 plis ; à 20 mM, la référence est elle-même incertaine — étape 12) ; [Cl⁻]i basale et son déplacement quand η
change ; erreur de η̂.

Sorties : results/ancrage_champs.csv, results/ancrage_resume.csv, results/ancrage.png
"""
import sys

import numpy as np
import pandas as pd

import banc as BA
import config as C
import figures as G
import step12_banc_essai as B12
from figures import plt
from step7_calibration import ajuster, inverser
from step13_additif import charger, regression

W, CONDITIONS = BA.W, BA.CONDITIONS
CONC = [0.0, 20.0, 131.0]
ETAS = (1.0, 0.7, 0.85, 1.15, 1.3)
GRAINE = 20261112


def forme_f(sB, p0, p1):
    """Part du motif « 131 mM » dans un motif CFP (moindres carrés pondérés : les motifs du capteur seul peuvent avoir
    des canaux légèrement négatifs)."""
    b = sB / sB.sum()
    d, w = p1 - p0, 2 / np.maximum(p0 + p1, 1e-9)
    return float(((b - p0) * d * w).sum() / (d * d * w).sum())


def champ(dec, aire, eta, rng):
    """Photons de chaque soma (cycles gardés réunis), η imposé par amincissement binomial ; motif du capteur seul."""
    d = dec.sum(0)                                           # (somas, canal, fenêtre)
    if eta != 1.0:
        d = d.astype(np.int64).copy()
        c, p = (1, eta) if eta < 1 else (0, 1 / eta)
        d[:, c] = rng.binomial(d[:, c], p)
    B, Y = d[:, 0].astype(float), d[:, 1].astype(float)
    N, A = (B + Y).sum(1), aire * len(dec)
    ok = N > 0
    sB, _ = regression(B[ok], N[ok], A[ok])
    sY, _ = regression(Y[ok], N[ok], A[ok])
    return dict(R=np.where(B.sum(1) > 0, Y.sum(1) / np.maximum(B.sum(1), 1), np.nan), sB=sB, R_capteur=sY.sum() / sB.sum())


def main():
    root = C.racine()
    D = C.dirs(root)
    A = [a for a in C.acquisitions(root) if a["type"] == "fov" and (D["mesures"] / f"{a['nom']}.npz").exists()]
    rng = np.random.default_rng(GRAINE)
    donnees = {a["nom"]: (a["condition"], B12.num_fov(a["nom"])) + charger(D, a["nom"])[:2] for a in A}
    # champs tels quels, une fois
    base = {n: champ(dec, aire, 1.0, rng) for n, (_, _, dec, aire) in donnees.items()}
    lignes = []
    for k in range(1, 6):
        train = [n for n, (c, f, _, _) in donnees.items() if c in CONDITIONS and f != k]
        lus = [n for n, (c, f, _, _) in donnees.items() if (c in CONDITIONS and f == k) or c == "basal"]
        # --- ancrage : gabarits et droite R(f) sur les champs d'entraînement à 0 et 131 mM
        ext = {c: [n for n in train if donnees[n][0] == c] for c in ("0 mM", "131 mM")}
        p0 = np.mean([base[n]["sB"] / base[n]["sB"].sum() for n in ext["0 mM"]], 0)
        p1 = np.mean([base[n]["sB"] / base[n]["sB"].sum() for n in ext["131 mM"]], 0)
        fx = np.array([forme_f(base[n]["sB"], p0, p1) for n in ext["0 mM"] + ext["131 mM"]])
        Rx = np.array([base[n]["R_capteur"] for n in ext["0 mM"] + ext["131 mM"]])
        b1, b0 = np.polyfit(fx, Rx, 1)
        R_f = lambda f: b0 + b1 * f
        eta_hat = lambda ch: ch["R_capteur"] / R_f(forme_f(ch["sB"], p0, p1))
        # --- conversion en mM : isotherme sur les médianes des somas d'entraînement
        med = lambda ancre: [np.nanmedian(np.concatenate(
            [base[n]["R"] / (eta_hat(base[n]) if ancre else 1.0) for n in train if donnees[n][0] == c])) for c in CONDITIONS]
        iso = {False: ajuster(CONC, med(False)), True: ajuster(CONC, med(True))}
        for n in lus:
            cond, fov, dec, aire = donnees[n]
            for eta in ETAS:
                ch = base[n] if eta == 1.0 else champ(dec, aire, eta, rng)
                eh = eta_hat(ch)
                for ancre in (False, True):
                    f = iso[ancre]
                    x = ch["R"] / (eh if ancre else 1.0)
                    c = inverser(x[np.isfinite(x)], *f) if f else np.full(np.isfinite(x).sum(), np.nan)
                    lignes.append(dict(pli=k, acq=n, cond=cond, fov=fov, eta=eta, ancre=ancre, eta_hat=eh,
                                       f_forme=forme_f(ch["sB"], p0, p1), mediane_mM=float(np.median(c)),
                                       K_mM=f[2] if f else np.nan))
        print(f"  pli {k} : R(f) = {b0:.3f} + {b1:.3f}·f ; K {iso[False][2]:.1f} mM (seul), {iso[True][2]:.1f} mM (ancré)")
    L = pd.DataFrame(lignes)
    L["cl_mM"] = L.cond.map(C.CL_MM)
    L["erreur"] = L.mediane_mM - L.cl_mM
    L.to_csv(D["results"] / "ancrage_champs.csv", index=False)

    res = []
    for (cond, eta, ancre), g in L.groupby(["cond", "eta", "ancre"]):
        r = dict(cond=cond, eta=eta, lecture="rapport ancré" if ancre else "rapport seul",
                 eta_hat_moyen=g.eta_hat.mean(), eta_hat_relatif_sd=(g.eta_hat / eta).std())
        if cond == "basal":
            ref = L[(L.cond == "basal") & (L.eta == 1.0) & (L.ancre == ancre)].set_index(["pli", "acq"]).mediane_mM
            d = g.set_index(["pli", "acq"]).mediane_mM - ref.reindex(g.set_index(["pli", "acq"]).index)
            r.update(basal_mM=g.groupby("pli").mediane_mM.median().median(), deplacement_basal_mM=d.groupby(level=0).median().mean())
        else:
            e = g.erreur
            r.update(biais_mM=e.mean(), rmse_mM=float(np.sqrt((e ** 2).mean())))
        res.append(r)
    RS = pd.DataFrame(res)
    RS.to_csv(D["results"] / "ancrage_resume.csv", index=False)
    pd.set_option("display.width", 220)
    for q in ("rmse_mM", "biais_mM"):
        print(f"\n{q} (médiane du champ lu − [Cl⁻] imposée, 5 plis)")
        print(RS[RS.cond != "basal"].pivot_table(index=["cond", "lecture"], columns="eta", values=q).round(2).to_string())
    print("\n[Cl⁻]i basale (médiane des plis) et son déplacement quand η change (mM)")
    print(RS[RS.cond == "basal"].pivot_table(index="lecture", columns="eta", values=["basal_mM", "deplacement_basal_mM"]).round(2).to_string())
    print("\nη̂ / η vrai : écart-type entre champs, par condition (rapport ancré)")
    print(RS[RS.lecture == "rapport ancré"].pivot_table(index="cond", columns="eta", values="eta_hat_relatif_sd").round(3).to_string())
    figure(D, RS, L)


def figure(D, RS, L):
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.8))
    for ax in axes:
        G.style(ax)
    for ax, cond in zip(axes[:3], CONDITIONS):
        for ancre, col, nom in ((False, G.PALETTE[1], "rapport seul"), (True, G.PALETTE[0], "rapport ancré")):
            g = RS[(RS.cond == cond) & (RS.lecture == nom)].sort_values("eta")
            ax.plot(g.eta, g.rmse_mM, "o-", color=col, label=nom)
            ax.plot(g.eta, g.biais_mM.abs(), "o:", color=col, alpha=0.6, ms=3)
        ax.axvline(1, color=G.TRACE)
        ax.set_yscale("symlog", linthresh=5)
        ax.set_xlabel("η imposé au champ lu", color=G.INK)
        ax.set_ylabel("erreur de la médiane du champ (mM)", color=G.INK)
        ax.legend(frameon=False, fontsize=8.5)
        G.titre(ax, f"{cond} — champ lu hors calibration",
                "plein : erreur quadratique sur 5 plis ; pointillé : |biais|", fontsize=10.5)
    ax = axes[3]
    for nom, col in (("rapport seul", G.PALETTE[1]), ("rapport ancré", G.PALETTE[0])):
        g = RS[(RS.cond == "basal") & (RS.lecture == nom)].sort_values("eta")
        ax.plot(g.eta, g.basal_mM, "o-", color=col, label=nom)
    ax.axvline(1, color=G.TRACE)
    ax.set_xlabel("η imposé aux champs du basal", color=G.INK)
    ax.set_ylabel("[Cl⁻]i basale lue (mM)", color=G.INK)
    ax.legend(frameon=False, fontsize=8.5)
    G.titre(ax, "Basal", "médiane des somas, médiane des plis", fontsize=10.5)
    fig.tight_layout()
    G.sauver(fig, D["results"] / "ancrage.png")


if __name__ == "__main__":
    main()
