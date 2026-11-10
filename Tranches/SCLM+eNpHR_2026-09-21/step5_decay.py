"""Étape 5 — les meilleures réponses, toutes données confondues, sur toute l'acquisition.

On retient les `N_BEST_TRACES` réponses (ROI x stim) les plus fortes de tout le jeu,
jugées sur `dR_R0_apres_immediat`, et on les trace avec exactement la disposition des
figures `trace_<FOV>_<acq>.png` : ratio C2/C3 brut en haut, dR/R0 en bas, sur toute la
durée de l'acquisition. L'ajustement `A·exp(-t/tau) + C` de la médiane après la stim —
le retour attendu vers la base — est calculé et imprimé, mais pas tracé.

Un même neurone peut apparaître plusieurs fois s'il répond fort à plusieurs stims.
Les acquisitions laissant moins de `MIN_POST_S` après la stim sont écartées
(FOV1/stim-002, dont la durée et le calendrier diffèrent).

**Biais de sélection.** Les réponses sont choisies sur leur amplitude juste après la
stim puis tracées sur cette même fenêtre : cette amplitude est surestimée (régression
vers la moyenne), et une partie des réponses retenues sont des ROIs bruités plutôt que
de vrais répondeurs. La forme du retour, et donc tau, en souffre beaucoup moins.

Usage :  python3 step5_decay.py   (après step4_ratio.py)
"""
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from config import MIN_POST_S, N_BEST_TRACES, dirs
from step4_ratio import trace_figure

COLOR = "#2a78d6"       # médiane, comme pour la première acquisition des figures trace
MIN_COVER = 0.8         # part des traits requise pour qu'un point entre dans l'ajustement
TAU_BOUNDS = (0.5, 400.0)


# ------------------------------------------------------------------ sélection
def post_duration(df):
    """Durée disponible après la stim, par acquisition."""
    out = {}
    for (f, a), g in df.groupby(["fov", "acquisition"]):
        fin = g.loc[g["stim"], "t_s"].max()
        out[(f, a)] = float(g.loc[g["phase"] == "apres", "t_s"].max() - fin)
    return out


def pick(D, df):
    """Les `N_BEST_TRACES` réponses les plus négatives, parmi les acquisitions candidates."""
    per_roi = pd.read_csv(D["results"] / "resume_par_roi.csv")
    dur = post_duration(df)
    per_roi["apres_s"] = [dur[(f, a)] for f, a in zip(per_roi["fov"], per_roi["acquisition"])]
    for (f, a), d in per_roi.groupby(["fov", "acquisition"])["apres_s"].first().items():
        if d < MIN_POST_S:
            print(f"  ! {f}/{a} écartée : {d:.0f} s après la stim "
                  f"(< {MIN_POST_S:.0f} s requises)")
    cand = per_roi[per_roi["apres_s"] >= MIN_POST_S]
    return cand.nsmallest(N_BEST_TRACES, "dR_R0_apres_immediat").reset_index(drop=True)


def gather(df, sel):
    """Décours complets des réponses retenues ; `roi` devient un identifiant de trait."""
    key = df.set_index(["fov", "acquisition", "roi"]).index
    want = pd.MultiIndex.from_frame(sel[["fov", "acquisition", "roi"]])
    sub = df[key.isin(want)].copy()
    sub["roi_origine"] = sub["roi"]
    sub["roi"] = sub["fov"] + "/" + sub["acquisition"] + "#" + sub["roi"].astype(str)
    # Les fenêtres de stim diffèrent d'une image d'une acquisition à l'autre (29–89,
    # 30–89, 30–90…). Sans masque commun, les images de bord ne seraient masquées que
    # pour une partie des traits et la médiane y ferait un saut. On masque donc partout
    # la réunion des fenêtres.
    union = sub.loc[sub["stim"], "image"].unique()
    sub["stim"] = sub["image"].isin(union)
    return sub


# ------------------------------------------------------------------ ajustement
def _model(t, A, tau, C):
    return A * np.exp(-t / tau) + C


def fit_decay(t, y):
    """A·exp(-t/tau) + C par moindres carrés. (A, tau, C, R2) ou None."""
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = np.asarray(t)[ok], np.asarray(y)[ok]
    if len(t) < 8:
        return None
    try:
        p, _ = curve_fit(_model, t, y, p0=(y[0] - y[-1], 20.0, y[-1]), maxfev=20000,
                         bounds=([-200, TAU_BOUNDS[0], -100], [200, TAU_BOUNDS[1], 100]))
    except (RuntimeError, ValueError):
        return None
    ss = np.sum((y - y.mean()) ** 2)
    r2 = 1 - np.sum((y - _model(t, *p)) ** 2) / ss if ss > 0 else np.nan
    return (*p, r2)


def median_after(sub):
    """Médiane des traits image par image après la stim, là où assez de traits restent.

    Les acquisitions durent 150 ou 180 images : au-delà de la plus courte, la médiane ne
    porte plus que sur une partie des traits et changerait de composition en cours de
    route. L'ajustement s'arrête donc là où au moins MIN_COVER des traits contribuent.
    """
    n = sub["roi"].nunique()
    fin = sub.loc[sub["stim"], "t_s"].max()
    a = sub[(sub["phase"] == "apres") & ~sub["stim"]]
    g = a.groupby("image").agg(t=("t_s", "median"), y=("dR_R0_pct", "median"),
                               n=("dR_R0_pct", "count"))
    g = g[g["n"] >= MIN_COVER * n]
    return g["t"].to_numpy() - fin, g["y"].to_numpy(), fin


# ------------------------------------------------------------------ main
def main():
    D = dirs()
    path = D["results"] / "ratio_temps.csv"
    if not path.exists():
        sys.exit(f"{path} absent — lance d'abord step4_ratio.py")
    df = pd.read_csv(path)

    sel = pick(D, df)
    sub = gather(df, sel)
    n_neur = sel[["fov", "roi"]].drop_duplicates().shape[0]
    print(f"\n{len(sel)} meilleures réponses ({sel['dR_R0_apres_immediat'].min():+.1f} à "
          f"{sel['dR_R0_apres_immediat'].max():+.1f} %), {n_neur} neurones distincts :")
    print(sel.groupby(["fov", "acquisition"]).size().to_string())

    t, y, fin = median_after(sub)
    fit = fit_decay(t, y)
    trace_figure(sub, f"Les {len(sel)} meilleures réponses, toutes stims confondues",
                 f"{len(sel)} réponses (ROI × stim), {n_neur} neurones, "
                 f"{sel['fov'].nunique()} FOVs", COLOR,
                 D["results"] / "retour_base_meilleurs_neurones.png")

    fits = []
    for trace, g in sub.groupby("roi", sort=False):
        f_ = g.loc[g["stim"], "t_s"].max()
        a = g[(g["phase"] == "apres") & ~g["stim"]]
        r = fit_decay(a["t_s"].to_numpy() - f_, a["dR_R0_pct"].to_numpy())
        fits.append(dict(fov=g["fov"].iloc[0], acquisition=g["acquisition"].iloc[0],
                         roi=int(g["roi_origine"].iloc[0]),
                         A=r[0] if r else np.nan, tau=r[1] if r else np.nan,
                         C=r[2] if r else np.nan, R2=r[3] if r else np.nan))
    fits = pd.DataFrame(fits)
    fits.to_csv(D["results"] / "retour_base_ajustements.csv", index=False)
    sel.drop(columns="apres_s").to_csv(D["results"] / "retour_base_neurones.csv", index=False)

    if fit:
        A, tau, C, r2 = fit
        print(f"\najustement de la médiane après la stim (jusqu'à {t.max():.0f} s) : "
              f"A={A:+.2f}  tau={tau:.1f} s  C={C:+.2f}  R2={r2:.2f}  "
              f"-> {abs(C / (A + C)):.0%} de la réponse subsiste")
    good = fits[(fits["R2"] > 0.2)]
    print(f"ajustements par trait : {len(good)}/{len(fits)} exploitables (R2 > 0,2), "
          f"tau médian {good['tau'].median():.0f} s")
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
