"""FOV11 — train de stimulations, analysé à part du pipeline principal.

FOV11 contient des trains de stimulations (`train_stim-002`, `train_stim-003`) : 270
images (275 s) et **quatre stimulations** d'une dizaine de secondes, espacées d'une
minute. Tous les décours de FOV11 sont traités, chacun avec ses propres figures. Les étapes 1 à 6 ne
savent traiter qu'une stimulation par acquisition ; ce script reprend leurs fonctions
en ajoutant ce qu'il faut pour un train, et écrit tout dans `FOV11_train/` — les
résultats des autres FOVs ne sont pas touchés.

Étapes enchaînées :
  1. fenêtres des stimulations (détection multi-pulses), recalage, projections ;
  2. segmentation Cellpose, lancée dans l'environnement dédié — seulement si les masques
     n'existent pas encore, pour ne pas renuméroter les ROIs d'un tri déjà fait ;
  3. sélection par défaut (tous les ROIs), seulement si aucune sélection n'existe ;
  4. ratio C2/C3 dans le temps, figure trace ;
  5. ΔR après chaque stimulation en fonction du R basal, comme deltaR_vs_Rbasal.png,
     une figure par train, et la reproductibilité d'un train à l'autre ;
  6. pour chaque train, les N_BEST neurones à la réponse typique la plus forte ;
  7. avant / après chaque stimulation, comme avant_apres_stim.png, un panneau par train ;
     et, pour les trains de PAR_PULSE, un panneau par stimulation.

Le tri manuel des ROIs se fait avec l'outil habituel, redirigé vers ce dossier :
    SCLM_OUT=FOV11_train SCLM_FOVS=FOV11 python3 step3_select.py FOV11
puis relancer ce script.

Usage :  python3 train_FOV11.py  [--resegmenter]
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOV = "FOV11"
# Avant tout import de config : toutes les sorties vont dans FOV11_train/.
os.environ["SCLM_OUT"] = str(HERE / "FOV11_train")
os.environ["SCLM_FOVS"] = FOV

import json            # noqa: E402
import subprocess      # noqa: E402

import numpy as np     # noqa: E402
import matplotlib      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import pandas as pd    # noqa: E402
from scipy.ndimage import binary_dilation, label   # noqa: E402
from scipy.stats import spearmanr, theilslopes      # noqa: E402

import config                                        # noqa: E402
from config import (MIN_VALID_FRAC, STIM_EDGE_W, STIM_MIN_RISE, STIM_PAD,  # noqa: E402
                    STIM_STEP_K, dirs, slug, survey)
import step1_align                                   # noqa: E402
import step4_ratio                                   # noqa: E402
import step6_deltaR_vs_Rbasal as step6               # noqa: E402

CELLPOSE_PY = Path.home() / ".venvs" / "cellpose" / "bin" / "python"
MIN_PULSE = STIM_EDGE_W + 1   # images entre bords : une image aberrante isolée, étalée
                              # par la fenêtre, donne deux bords distants d'une image
WIN = step6.WIN        # 10 s avant / après chaque stimulation
N_BEST = 5             # neurones montrés par train à l'étape 6
PAR_PULSE = ("train_stim-003",)   # trains aux pulses de durées différentes : une figure
                                  # avant / après avec un panneau par stimulation
BEST_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#9b51e0", "#d64550", "#00867d",
               "#c9a227", "#6b6b66"]


# ------------------------------------------------------------------ 1. stimulations
def train_window(c2):
    """Images contaminées par la lumière de stim, pour un train de pulses.

    Même principe que la détection du pipeline principal (step1_align.stim_window) :
    la fuite dans C2 fait des plateaux nets sur la moyenne pleine image, repérés par
    leurs bords sur le contraste fenêtré `edge_response`. Mais au lieu du seul couple
    (bord montant le plus fort, bord descendant le plus fort), on garde **tous** les
    bords au-dessus du seuil et on apparie chaque bord montant au bord descendant qui
    le suit. Chaque pulse est ensuite délimité par un seuil à mi-chemin entre son
    plateau et la ligne de base qui l'entoure, puis élargi de STIM_PAD image.
    """
    m = c2.reshape(len(c2), -1).mean(axis=1)
    n = len(m)
    none = np.zeros(n, bool)
    diag = {"bord_sur_bruit": 0.0, "niveau": None, "hausse": None, "pulses": []}
    g = step1_align.edge_response(m)
    ok = np.isfinite(g)
    if not ok.any():
        return none, m, diag
    scale = 1.4826 * np.median(np.abs(g[ok] - np.median(g[ok]))) or 1e-12
    gz = np.where(ok, g, 0.0)

    def edges(sign):
        """Position du maximum de chaque région où sign·g dépasse le seuil."""
        lab, k = label(sign * gz > STIM_STEP_K * scale)
        return [int(np.flatnonzero(lab == j)[np.argmax(sign * gz[lab == j])])
                for j in range(1, k + 1)]

    ups, downs = edges(+1), edges(-1)
    hit, ratios, rises, W = none.copy(), [], [], STIM_EDGE_W
    for r in ups:
        # r = première image éclairée. Un bord montant déjà couvert par un pulse (fuite
        # qui monte par paliers) n'en ouvre pas un nouveau.
        if hit[r]:
            continue
        lo = float(np.median(m[max(r - 2 * W, 0):max(r - 1, 1)]))
        h = gz[r]
        # Fin du pulse : premier bord descendant après lequel C2 est revenu près de la
        # base. Le plateau n'est pas toujours plat (la fuite peut décroître pendant la
        # stim) : un simple « bord descendant suivant » s'arrêterait sur une marche
        # interne.
        f = next((d for d in downs if d > r and np.median(m[d:d + W]) < lo + 0.5 * h), None)
        if f is None or f - r < MIN_PULSE:
            continue
        hi = float(np.median(m[r + 1:f - 1]))
        base = float(np.median(np.r_[m[max(r - 2 * W, 0):max(r - 1, 1)], m[f + 1:f + 2 * W]]))
        # Le plateau doit réellement être haut : une image aberrante isolée, étalée par
        # la fenêtre, a deux bords mais un intérieur à la hauteur de la base.
        if hi < lo + 0.5 * h or base <= 0 or hi / base - 1 < STIM_MIN_RISE:
            continue
        span = np.zeros(n, bool)
        span[max(r - W, 0):min(f + W, n)] = True
        lab, _ = label(span & (m > (hi + base) / 2))
        mid = (r + f - 1) // 2
        # Au moins STIM_EDGE_W images franchement hautes : deux images aberrantes
        # voisines passent les autres tests mais ne font qu'un « pulse » de 1-2 images.
        if not lab[mid] or (lab == lab[mid]).sum() < STIM_EDGE_W:
            continue
        hit |= binary_dilation(lab == lab[mid], np.ones(2 * STIM_PAD + 1, bool))
        ratios.append(min(gz[r], -gz[f]) / scale)
        rises.append(hi / base - 1)
    lab, k = label(hit)
    diag["pulses"] = [[int(np.flatnonzero(lab == j)[0]) + 1, int(np.flatnonzero(lab == j)[-1]) + 1]
                      for j in range(1, k + 1)]
    if ratios:
        diag.update(bord_sur_bruit=float(min(ratios)), hausse=float(np.median(rises)))
    return hit, m, diag


def align(D):
    # process_fov appelle stim_window : on lui substitue la détection multi-pulses.
    step1_align.stim_window = train_window
    recs = survey()[FOV]
    log = ["acquisition\tetape\tdy / detail\tdx\tNCC\tstatut"]
    print(f"  {FOV} ({len(recs)} acquisitions)")
    step1_align.process_fov(FOV, recs, D, log)
    (D["qc"] / "alignement.tsv").write_text("\n".join(log) + "\n")
    info = json.loads((D["projections"] / f"{FOV}_align.json").read_text())["info"]
    for rec, i in info.items():
        pulses = (i.get("stim_detection") or {}).get("pulses") or []
        if pulses:
            t = i["times"]
            print(f"    {rec} : {len(pulses)} stimulations — " + ", ".join(
                f"{a}–{b} ({t[a-1]:.0f}–{t[b-1]:.0f} s)" for a, b in pulses))


# ------------------------------------------------------------------ 2-3. ROIs
def segment(D, force):
    if (D["masks"] / f"{FOV}_masks.npy").exists() and not force:
        print(f"  masques existants conservés (--resegmenter pour relancer Cellpose)")
        return
    if not CELLPOSE_PY.exists():
        sys.exit(f"environnement Cellpose introuvable : {CELLPOSE_PY}")
    subprocess.run([str(CELLPOSE_PY), str(HERE / "step2_cellpose.py"), FOV],
                   check=True, env=os.environ.copy(), cwd=HERE)


def default_selection(D):
    path = D["selection"] / f"{FOV}_selection.json"
    if path.exists():
        kept = json.loads(path.read_text())["kept"]
        print(f"  sélection existante conservée ({len(kept)} ROIs)")
        return
    n = int(np.load(D["masks"] / f"{FOV}_masks.npy").max())
    path.write_text(json.dumps({"fov": FOV, "n_rois": n, "kept": list(range(1, n + 1)),
                                "keep": {str(l): True for l in range(1, n + 1)}}, indent=1))
    print(f"  {n} ROIs gardés (sélection par défaut — tri manuel : voir l'en-tête du script)")


# ------------------------------------------------------------------ 4. ratio dans le temps
def measure(D):
    kept = json.loads((D["selection"] / f"{FOV}_selection.json").read_text())["kept"]
    align = json.loads((D["projections"] / f"{FOV}_align.json").read_text())
    masks = np.load(D["masks"] / f"{FOV}_masks.npy")
    n_lab = int(masks.max())
    labels = np.where(np.isin(masks, kept), masks, 0).astype(np.int32)
    rows = []
    for rec, n_frames, date in survey()[FOV]:
        info = align["info"][rec]
        if n_frames < config.MIN_FRAMES_TIMELAPSE or not any(info["stim"]):
            continue
        mean, _, frac, npx = step4_ratio.measure_recording(FOV, rec, info, labels, n_lab)
        stim = np.asarray(info["stim"], bool)
        # Numéro du pulse que chaque image suit : 0 avant le premier, k après le k-ième.
        lab, _ = label(stim)
        starts = [np.flatnonzero(lab == j)[0] for j in range(1, lab.max() + 1)]
        pulse = np.searchsorted(starts, np.arange(n_frames), side="right")
        for lab_roi in kept:
            c2, c3 = mean["C2"][:, lab_roi], mean["C3"][:, lab_roi]
            ok = frac[:, lab_roi] >= MIN_VALID_FRAC
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio = np.where(ok & (c3 > 0), c2 / c3, np.nan)
            for i in range(n_frames):
                rows.append(dict(fov=FOV, acquisition=rec, date=date, roi=int(lab_roi),
                                 image=i + 1, t_s=round(info["times"][i], 4),
                                 mean_C2=c2[i], mean_C3=c3[i], ratio=ratio[i],
                                 stim=bool(stim[i]), pulse=int(pulse[i]),
                                 frac_valide=round(float(frac[i, lab_roi]), 4)))
    df = pd.DataFrame(rows)
    # R0 = médiane avant le premier pulse (add_normalised prend la première image de stim).
    df["note"] = ""
    df = step4_ratio.add_normalised(df)
    # add_normalised ne connaît qu'une stim : les images entre deux pulses sont renommées.
    df["phase"] = np.select([df["stim"], df["pulse"] == 0],
                            ["pendant", "avant"], default="apres")
    df.to_csv(D["results"] / "ratio_temps.csv", index=False)
    return df


# ------------------------------------------------------------------ 5. ΔR contre R basal
def delta_vs_basal(df):
    """ΔR autour de chaque pulse ; R basal pris avant le premier, commun aux quatre."""
    rows = []
    for rec, g in df.groupby("acquisition"):
        t_first = g.loc[g["stim"], "t_s"].min()
        bas = g[(g["pulse"] == 0) & ~g["stim"] & (g["t_s"] < t_first - WIN)
                & (g["image"] > step4_ratio.BASE_SKIP)].groupby("roi")["ratio"].mean()
        for k, s in g[g["stim"]].groupby("pulse"):
            t0, t1 = s["t_s"].min(), s["t_s"].max()
            av = g[~g["stim"] & (g["t_s"] < t0) & (g["t_s"] >= t0 - WIN)]
            ap = g[~g["stim"] & (g["t_s"] > t1) & (g["t_s"] <= t1 + WIN)]
            w = pd.DataFrame({"R_basal": bas, "R_avant": av.groupby("roi")["ratio"].mean(),
                              "R_apres": ap.groupby("roi")["ratio"].mean()})
            w["dR"] = w["R_apres"] - w["R_avant"]
            w["dR_R0_pct"] = 100 * w["dR"] / w["R_avant"]
            rows.append(w.reset_index().assign(fov=FOV, acquisition=rec, pulse=int(k),
                                               debut_s=round(t0, 1)))
    x = pd.concat(rows, ignore_index=True).dropna(subset=["R_basal", "dR"])
    return x[["fov", "acquisition", "pulse", "debut_s", "roi", "R_basal", "R_avant",
              "R_apres", "dR", "dR_R0_pct"]]


# ------------------------------------------------------------------ 6. meilleurs neurones
def best_neurons(xr, n=N_BEST):
    """Les n neurones à la réponse typique la plus forte d'un train.

    Réponse typique = médiane, sur les pulses du train, du ΔR/R0 (10 s après - 10 s
    avant). La médiane écarte un pulse isolé exceptionnel : un neurone n'est retenu que
    s'il répond fort à la plupart des pulses.
    """
    w = xr.pivot(index="roi", columns="pulse", values="dR_R0_pct")
    out = pd.DataFrame({"typique": w.median(axis=1),
                        "pulses_en_baisse": (w < 0).sum(axis=1)})
    out = out.join(w.add_prefix("pulse_"))
    return out.nsmallest(n, "typique")


def best_figure(g, best, rec, D, colors):
    """Disposition des figures trace (ratio en haut, ΔR/R0 en bas, toute la durée),
    un neurone par couleur : valeurs brutes en trait fin, lissées en trait épais."""
    SMOOTH = step4_ratio.SMOOTH
    fig, axes = plt.subplots(2, 1, figsize=(9.6, 7.3), sharex=True,
                             facecolor=step4_ratio.SURFACE)
    for ax in axes:
        step4_ratio._style(ax)
        step4_ratio._stim_band(ax, g)
    axes[1].axhline(0, color=step4_ratio.INK_2, linewidth=0.9, linestyle="--", zorder=1)
    lows = []
    for j, (roi, row) in enumerate(best.iterrows()):
        s = g[g["roi"] == roi].sort_values("image")
        t, color = s["t_s"].to_numpy(), colors[roi]
        for ax, col in ((axes[0], "ratio"), (axes[1], "dR_R0_pct")):
            y = np.where(s["stim"], np.nan, s[col].to_numpy())
            ax.plot(t, y, color=color, linewidth=0.6, alpha=1, zorder=2)
            # ax.plot(t, step4_ratio._smooth(y[None, :])[0], color=color, linewidth=1.9,
                    # zorder=3, label=f"ROI {roi}  ({row['typique']:+.1f} %)" if col == "ratio" else None)
            if col == "dR_R0_pct":
                lows.append(y)
    step4_ratio._robust_ylim(axes[1], np.concatenate(lows))
    axes[0].set_ylabel("Ratio C2 / C3", color=step4_ratio.INK, fontsize=11)
    axes[1].set_ylabel("$\\Delta$R / R$_0$  (%)", color=step4_ratio.INK, fontsize=11)
    axes[1].set_xlabel("temps (s)", color=step4_ratio.INK, fontsize=11)
    axes[0].set_title(f"{FOV} — {rec} : les {len(best)} neurones aux réponses typiques "
                      f"les plus fortes", color=step4_ratio.INK, fontsize=13, pad=40, loc="left")
    axes[0].text(0, 1.015, f"réponse typique = médiane des {g['pulse'].max()} pulses du "
                 f"ΔR/R$_0$ (10 s après − 10 s avant), entre parenthèses\ntraits épais lissés "
                 f"sur {SMOOTH} images   ·   bandes jaunes = stimulations, images masquées",
                 transform=axes[0].transAxes, color=step4_ratio.INK_2, fontsize=9.5,
                 va="bottom", linespacing=1.5)
    fig.legend(loc="lower center", ncol=len(best), frameon=False, fontsize=10,
               labelcolor=step4_ratio.INK_2, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = D["results"] / f"meilleurs_neurones_{slug(rec)}.png"
    fig.savefig(out, dpi=200, facecolor=step4_ratio.SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


# ------------------------------------------------------------------ 7. avant / après
def pre_post(df):
    """dR/R0 moyen par ROI, 10 s avant le début et 10 s après la fin de chaque pulse.

    R0 est la ligne de base avant le premier pulse : quand les réponses s'additionnent
    (second train), le « avant » des derniers pulses est déjà sous zéro.
    """
    rows = []
    for rec, g in df.groupby("acquisition"):
        for k, s in g[g["stim"]].groupby("pulse"):
            t0, t1 = s["t_s"].min(), s["t_s"].max()
            av = g[~g["stim"] & (g["t_s"] < t0) & (g["t_s"] >= t0 - WIN)]
            ap = g[~g["stim"] & (g["t_s"] > t1) & (g["t_s"] <= t1 + WIN)]
            w = pd.DataFrame({"dR_R0_avant_pct": av.groupby("roi")["dR_R0_pct"].mean(),
                              "dR_R0_apres_pct": ap.groupby("roi")["dR_R0_pct"].mean()})
            rows.append(w.reset_index().assign(fov=FOV, acquisition=rec, pulse=int(k)))
    pp = pd.concat(rows, ignore_index=True)
    return pp[["fov", "acquisition", "pulse", "roi", "dR_R0_avant_pct", "dR_R0_apres_pct"]]


def pre_post_figure(pp, D):
    """Même présentation que avant_apres_stim.png du pipeline principal."""
    S = step4_ratio
    recs = list(dict.fromkeys(pp["acquisition"]))
    cols = ["dR_R0_avant_pct", "dR_R0_apres_pct"]
    labels = [f"{WIN:.0f} s\navant", f"{WIN:.0f} s\naprès"]
    fig, axes = plt.subplots(1, len(recs), figsize=(3.4 * len(recs) + 1.6, 5.4),
                             sharey=True, facecolor=S.SURFACE, squeeze=False)
    for ax, rec in zip(axes[0], recs):
        sub = pp[pp["acquisition"] == rec]
        w = sub[cols].dropna().reset_index(drop=True)
        S._paired(ax, w, labels)
        ax.axhline(0, color=S.INK_2, linewidth=0.9, linestyle="--", zorder=2)
        med = w.median()
        ax.annotate(f"{med.iloc[1]:+.1f} %", (1, med.iloc[1]), textcoords="offset points",
                    xytext=(10, 0), va="center", color=S.INK, fontsize=10.5, zorder=5)
        ax.set_title(f"{rec}\n{sub['pulse'].nunique()} stimulations, n = {len(w)}",
                     color=S.INK, fontsize=11, pad=6)
    v = pp[cols].to_numpy(); v = v[np.isfinite(v)]
    pad = 0.05 * (v.max() - v.min())
    axes[0][0].set_ylim(v.min() - pad, v.max() + pad)      # tous les points visibles
    axes[0][0].set_ylabel("$\\Delta$R / R$_0$ moyen  (%)", color=S.INK, fontsize=11)
    fig.suptitle(f"{FOV} — changement du ratio C2/C3, {WIN:.0f} s avant chaque stimulation "
                 f"et {WIN:.0f} s après\nun trait par ROI et par stimulation ; trait noir = "
                 f"médiane", color=S.INK, fontsize=12, x=0.01, ha="left", y=0.985,
                 linespacing=1.4)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    out = D["results"] / "avant_apres_stim.png"
    fig.savefig(out, dpi=200, facecolor=S.SURFACE); plt.close(fig)
    print(f"  -> {out.name}")


def pulse_durations(g):
    """Durée de chaque pulse d'un train, en s : images de stim moins les images de garde."""
    dt = float(np.median(np.diff(np.sort(g["t_s"].unique()))))
    n = g[g["stim"]].groupby("pulse")["image"].nunique()
    return ((n - 2 * STIM_PAD) * dt).round().astype(int)


def pre_post_pulses_figure(pp, df, rec, D):
    """Avant / après pour un seul train, un panneau par stimulation, axe y commun."""
    S = step4_ratio
    sub = pp[pp["acquisition"] == rec]
    dur = pulse_durations(df[df["acquisition"] == rec])
    pulses = sorted(sub["pulse"].unique())
    cols = ["dR_R0_avant_pct", "dR_R0_apres_pct"]
    labels = [f"{WIN:.0f} s\navant", f"{WIN:.0f} s\naprès"]
    fig, axes = plt.subplots(1, len(pulses), figsize=(3.0 * len(pulses) + 1.6, 5.4),
                             sharey=True, facecolor=S.SURFACE, squeeze=False)
    for ax, k in zip(axes[0], pulses):
        w = sub.loc[sub["pulse"] == k, cols].dropna().reset_index(drop=True)
        S._paired(ax, w, labels)
        ax.axhline(0, color=S.INK_2, linewidth=0.9, linestyle="--", zorder=2)
        med = w.median()
        ax.annotate(f"{med.iloc[1]:+.1f} %", (1, med.iloc[1]), textcoords="offset points",
                    xytext=(10, 0), va="center", color=S.INK, fontsize=10.5, zorder=5)
        ax.set_title(f"stimulation {k}  ·  ≈ {dur[k]} s\nn = {len(w)}",
                     color=S.INK, fontsize=11, pad=6)
    v = sub[cols].to_numpy(); v = v[np.isfinite(v)]
    pad = 0.05 * (v.max() - v.min())
    axes[0][0].set_ylim(v.min() - pad, v.max() + pad)
    axes[0][0].set_ylabel("$\\Delta$R / R$_0$ moyen  (%)", color=S.INK, fontsize=11)
    fig.suptitle(f"{FOV} — {rec} : changement du ratio C2/C3, {WIN:.0f} s avant et "
                 f"{WIN:.0f} s après chaque stimulation\nun trait par ROI ; trait noir = "
                 f"médiane ; R$_0$ = base avant la première stimulation", color=S.INK,
                 fontsize=12, x=0.01, ha="left", y=0.985, linespacing=1.4)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    out = D["results"] / f"avant_apres_stim_{slug(rec)}_par_stimulation.png"
    fig.savefig(out, dpi=200, facecolor=S.SURFACE); plt.close(fig)
    print(f"  -> {out.name}")
    return dur


# ------------------------------------------------------------------ main
def main():
    D = dirs()
    print(f"Sorties -> {config.OUT_ROOT}\n\n1. stimulations et recalage")
    align(D)
    print("\n2. segmentation")
    segment(D, "--resegmenter" in sys.argv[1:])
    print("\n3. sélection")
    default_selection(D)

    print("\n4. ratio dans le temps")
    df = measure(D)
    for rec, g in df.groupby("acquisition"):
        step4_ratio.trace_figure(g, f"{FOV} — {rec} : ratio FRET dans le temps",
                                 f"{g['roi'].nunique()} ROIs", step4_ratio.PALETTE[0],
                                 D["results"] / f"trace_{FOV}_{slug(rec)}.png")

    print("\n5. ΔR en fonction du R basal")
    x = delta_vs_basal(df)
    x.to_csv(D["results"] / "deltaR_vs_Rbasal.csv", index=False)
    per_cell = {}
    for rec, xr in x.groupby("acquisition"):
        print(f"  {rec}")
        st = step6.within(xr, "dR", group="pulse")
        step6.figure(xr, st, D, ylabel="ΔR après stim  (ratio C2 / C3)",
                     out=D["results"] / f"deltaR_vs_Rbasal_{slug(rec)}.png")
        for k, g in xr.groupby("pulse"):
            print(f"    stimulation {k} ({g['debut_s'].iloc[0]:.0f} s) : ΔR médian "
                  f"{g['dR'].median():+.3f}, ΔR/R0 médian {g['dR_R0_pct'].median():+.2f} %  "
                  f"({int((g['dR'] < 0).sum())}/{len(g)} ROIs en baisse)")
        # Les 4 points d'une cellule partagent le même R basal : le test porte sur la
        # moyenne des pulses par cellule, pas sur les points.
        c = xr.groupby("roi").agg(R=("R_basal", "first"), dR=("dR", "mean"),
                                  dRR=("dR_R0_pct", "mean"))
        per_cell[rec] = c
        for col, unit in (("dR", "ratio"), ("dRR", "%")):
            r, p = spearmanr(c["R"], c[col])
            sl, _, lo, hi = theilslopes(c[col], c["R"])
            print(f"    par cellule (n = {len(c)}), {'ΔR absolu' if col == 'dR' else 'ΔR/R0   '} : "
                  f"rho = {r:+.2f} ({step6._pfmt(p)}), pente {sl:+.3g} {unit}/unité de R "
                  f"[IC95 {lo:+.3g}, {hi:+.3g}]")
    if len(per_cell) > 1:
        # Mêmes cellules d'un train à l'autre : la réponse et le R basal se reproduisent-ils ?
        recs = list(per_cell)
        for a, b in zip(recs[:-1], recs[1:]):
            m = per_cell[a].join(per_cell[b], lsuffix="_a", rsuffix="_b", how="inner").dropna()
            print(f"  {a} -> {b} (n = {len(m)} cellules) : R basal rho = "
                  f"{spearmanr(m['R_a'], m['R_b'])[0]:+.2f} (médian {m['R_a'].median():.2f} -> "
                  f"{m['R_b'].median():.2f}) ; ΔR/R0 rho = "
                  f"{spearmanr(m['dRR_a'], m['dRR_b'])[0]:+.2f} "
                  f"({step6._pfmt(spearmanr(m['dRR_a'], m['dRR_b'])[1])})")
    print(f"\n6. les {N_BEST} neurones aux réponses typiques les plus fortes")
    bests = {rec: best_neurons(xr) for rec, xr in x.groupby("acquisition")}
    # Une couleur par ROI, la même d'un train à l'autre.
    union = list(dict.fromkeys(r for b in bests.values() for r in b.index))
    colors = {r: BEST_COLORS[i % len(BEST_COLORS)] for i, r in enumerate(union)}
    rows = []
    for rec, xr in x.groupby("acquisition"):
        best = bests[rec]
        best_figure(df[df["acquisition"] == rec], best, rec, D, colors)
        print(f"  {rec} : " + ", ".join(
            f"ROI {r} ({v['typique']:+.1f} %, {int(v['pulses_en_baisse'])}/{xr['pulse'].nunique()} pulses en baisse)"
            for r, v in best.iterrows()))
        rows.append(best.reset_index().assign(acquisition=rec))
    pd.concat(rows, ignore_index=True).to_csv(D["results"] / "meilleurs_neurones.csv", index=False)

    print(f"\n7. avant / après chaque stimulation")
    pp = pre_post(df)
    pp.to_csv(D["results"] / "avant_apres_stim.csv", index=False)
    pre_post_figure(pp, D)
    for rec, g in pp.dropna().groupby("acquisition"):
        dd = g["dR_R0_apres_pct"] - g["dR_R0_avant_pct"]
        print(f"  {rec} : {g['dR_R0_avant_pct'].median():+.2f} % -> "
              f"{g['dR_R0_apres_pct'].median():+.2f} %  (écart médian {dd.median():+.2f} pts, "
              f"{int((dd < 0).sum())}/{len(dd)} en baisse)")
        for k, gk in g.groupby("pulse"):
            dk = gk["dR_R0_apres_pct"] - gk["dR_R0_avant_pct"]
            print(f"      pulse {k} : {gk['dR_R0_avant_pct'].median():+.2f} -> "
                  f"{gk['dR_R0_apres_pct'].median():+.2f} %  (écart médian {dk.median():+.2f} pts, "
                  f"{int((dk < 0).sum())}/{len(dk)} en baisse)")
    for rec in PAR_PULSE:
        if (pp["acquisition"] == rec).any():
            dur = pre_post_pulses_figure(pp, df, rec, D)
            print(f"  {rec} : durées des stimulations " + ", ".join(f"{d} s" for d in dur))
    print(f"\nRésultats -> {D['results']}")


if __name__ == "__main__":
    main()
