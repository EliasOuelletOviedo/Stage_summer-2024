"""Style commun des figures — celui des tranches (sclm_flim/figures.py)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e6e2"
TRACE = "#b4b4ac"
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#9b51e0", "#d64550", "#00867d"]
COULEUR_CANAL = {"cfp": "#2a78d6", "yfp": "#c99a06"}      # CFP bleu, YFP jaune foncé


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, length=0)


def titre(ax, texte, sous_titre=None, fontsize=12):
    ax.set_title(texte, color=INK, fontsize=fontsize, loc="left", pad=22 if sous_titre else 8)
    if sous_titre:
        ax.text(0, 1.02, sous_titre, transform=ax.transAxes, color=INK_2, fontsize=9.2, va="bottom")


def contraste(im, lo=1, hi=99.7):
    v = np.asarray(im, float)
    v = v[np.isfinite(v)]
    a, b = np.percentile(v, (lo, hi)) if v.size else (0, 1)
    return dict(vmin=a, vmax=max(b, a + 1e-9))


def sauver(fig, path, dpi=160, **kw):
    fig.savefig(path, dpi=dpi, facecolor=SURFACE, **kw)
    plt.close(fig)
    print(f"  -> {path.parent.name}/{path.name}")
