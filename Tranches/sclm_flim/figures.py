"""Style commun des figures — celui des tranches 2026-09-21 (step4_ratio.py)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

from config import NOMS_CANAUX  # noqa: E402

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
TRACE = "#b4b4ac"
STIM_BAND = "#f4d35e"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#9b51e0", "#d64550", "#00867d"]
COULEUR_CANAL = {1: "#2a78d6", 2: "#c99a06"}      # CFP bleu, YFP jaune foncé
LISSAGE = 17          # images, ~1 s à 16,7 Hz : moyenne glissante des figures image par image

# Thème : les fonctions qui l'acceptent (style, bandes_stim, sauver, figure_par_image et
# les figures des étapes 6 et 7) dessinent avec CLAIR par défaut ; un autre dict aux mêmes
# clés change l'apparence. Les figures de présentation ont leurs propres réglages, dans
# <jeu>/présentation/presentation.py.
CLAIR = dict(surface=SURFACE, ink=INK, ink_2=INK_2, grid=GRID, stim=STIM_BAND, stim_alpha=0.35,
             palette=PALETTE, neutre="#8d8d86", repere="#d9d9d4", alpha_points=0.55,
             canaux=COULEUR_CANAL, exclu_alpha=0.13)


def style(ax, T=None):
    T = T or CLAIR
    ax.set_facecolor(T["surface"])
    ax.grid(True, color=T["grid"], linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(T["grid"])
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=T["ink_2"], length=0)


def bandes_stim(ax, fenetres, T=None, **kw):
    T = T or CLAIR
    for a, b in fenetres:
        ax.axvspan(a, b, color=T["stim"], alpha=kw.get("alpha", T["stim_alpha"]), lw=0, zorder=0)


def ylim_robuste(ax, M, pad=0.12):
    v = np.asarray(M, float)
    v = v[np.isfinite(v)]
    if v.size < 10:
        return
    lo, hi = np.percentile(v, (0.5, 99.5))
    m = (hi - lo) * pad or 1.0
    ax.set_ylim(lo - m, hi + m)


def contraste(im, lo=1, hi=99.7):
    v = np.asarray(im, float)
    v = v[np.isfinite(v)]
    a, b = np.percentile(v, (lo, hi)) if v.size else (0, 1)
    return dict(vmin=a, vmax=max(b, a + 1e-9))


def sauver(fig, path, dpi=180, T=None, **kw):
    fig.savefig(path, dpi=dpi, facecolor=(T or CLAIR)["surface"], **kw)
    plt.close(fig)
    print(f"  -> {path.name}")


def figure_par_image(tm, c1, c2, fen, etiquette, out, T=None):
    """Photons ch1, photons ch2 et ratio ch2 / ch1, une valeur par image (60 ms), moyenne
    glissante sur 1 s ; `etiquette` dit sur quels pixels (« champ entier », « 34 somas…»)."""
    T = T or CLAIR
    fig, axes = plt.subplots(3, 1, figsize=(11, 7.6), sharex=True, facecolor=T["surface"])
    for ax in axes:
        style(ax, T)
        bandes_stim(ax, fen, T)
    lis = lambda v: ndi.uniform_filter1d(np.asarray(v, float), LISSAGE)
    for ax, k, v in zip(axes[:2], (1, 2), (c1, c2)):
        ax.plot(tm, v, color=T["canaux"][k], lw=0.3, alpha=0.35)
        ax.plot(tm, lis(v), color=T["canaux"][k], lw=1.4)
        ax.set_ylabel(f"ch{k} ({NOMS_CANAUX[k]})\nphotons / image", color=T["ink"], fontsize=10)
        ylim_robuste(ax, v)
    r = np.asarray(c2, float) / np.maximum(c1, 1)
    axes[2].plot(tm, r, color=T["ink_2"], lw=0.3, alpha=0.35)
    axes[2].plot(tm, lis(c2) / lis(c1), color=T["ink"], lw=1.4)
    ylim_robuste(axes[2], r)
    axes[2].set_ylabel(f"ratio ch2 / ch1\n({etiquette})", color=T["ink"], fontsize=10)
    axes[2].set_xlabel("temps (s, horloge des cartes SPC)", color=T["ink"], fontsize=10.5)
    axes[0].set_title(f"{etiquette[0].upper() + etiquette[1:]}, image par image (60 ms)",
                      color=T["ink"], fontsize=12, loc="left", pad=22)
    axes[0].text(0, 1.02, "trait épais : moyenne glissante sur 1 s ; bandes jaunes : commande de "
                 "stimulation", transform=axes[0].transAxes, color=T["ink_2"], fontsize=9.5,
                 va="bottom")
    fig.tight_layout()
    sauver(fig, out, T=T)
