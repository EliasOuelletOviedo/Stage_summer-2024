"""Étape 3 — sélection manuelle des ROIs, un FOV à la fois.

Usage :  python3 step3_select.py  [Slice1_FOV3 ...]

Commandes
  clic sur un ROI  garder / rejeter
  a / n            tout garder / tout rejeter
  b                change le fond (somme, Basal, High_K, Return)
  l                affiche ou masque les numéros de ROI
  flèche droite    enregistre et passe au FOV suivant
  flèche gauche    enregistre et revient au FOV précédent
  q                enregistre et quitte
"""
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("macosx")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from config import CONDITIONS, SLICES, N_FOV, fov_key, dirs

D = dirs()
BACKGROUNDS = ["summed"] + [f"{c}_C2" for c in CONDITIONS]
BG_LABEL = {"summed": "somme (C2+C3 x 3 conditions)", "Basal_C2": "Basal C2",
            "High_K_C2": "High_K C2", "Return_C2": "Return C2"}


class Selector:
    def __init__(self, keys):
        self.keys = keys
        self.i = 0
        self.bg = 0
        self.labels_on = True
        self.fig, self.ax = plt.subplots(figsize=(11, 11))
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.load()
        plt.show()

    # ---------- données ----------
    def load(self):
        key = self.keys[self.i]
        self.key = key
        self.npz = np.load(D["registered"] / f"{key}.npz")
        self.masks = np.load(D["masks"] / f"{key}_masks.npy")
        self.n = int(self.masks.max())
        self.centroids = {}
        for lab in range(1, self.n + 1):
            ys, xs = np.nonzero(self.masks == lab)
            if len(ys):
                self.centroids[lab] = (xs.mean(), ys.mean())

        path = D["selection"] / f"{key}_selection.json"
        if path.exists():
            saved = json.loads(path.read_text())
            self.keep = {int(k): v for k, v in saved["keep"].items()}
            for lab in range(1, self.n + 1):
                self.keep.setdefault(lab, True)
        else:
            self.keep = {lab: True for lab in range(1, self.n + 1)}   # tout gardé par défaut
        self.draw()

    def save(self):
        path = D["selection"] / f"{self.key}_selection.json"
        path.write_text(json.dumps(
            {"fov": self.key, "n_rois": self.n,
             "kept": sorted(l for l, v in self.keep.items() if v),
             "keep": {str(k): bool(v) for k, v in self.keep.items()}}, indent=1))

    # ---------- affichage ----------
    def draw(self):
        self.ax.clear()
        name = BACKGROUNDS[self.bg]
        img = self.npz[name] if name in self.npz.files else self.npz["summed"]
        img = np.nan_to_num(img, nan=0.0)
        lo, hi = np.percentile(img, (1, 99.7))
        self.ax.imshow(img, cmap="gray", vmin=lo, vmax=max(hi, lo + 1))

        for lab in range(1, self.n + 1):
            keep = self.keep.get(lab, True)
            col = "#2ecc71" if keep else "#e74c3c"
            self.ax.contour(self.masks == lab, levels=[0.5], colors=[col],
                            linewidths=1.4 if keep else 0.7,
                            linestyles="solid" if keep else "dotted")
            if self.labels_on and lab in self.centroids:
                x, y = self.centroids[lab]
                self.ax.text(x, y, str(lab), color=col, fontsize=7,
                             ha="center", va="center", fontweight="bold",
                             bbox=dict(boxstyle="round,pad=0.1", fc=to_rgba("k", .45), ec="none"))

        nk = sum(self.keep.values())
        miss = [c for c in CONDITIONS if f"{c}_C2" not in self.npz.files]
        warn = f"   [manque : {', '.join(miss)}]" if miss else ""
        self.ax.set_title(
            f"{self.key}   ({self.i+1}/{len(self.keys)})   fond : {BG_LABEL.get(name, name)}\n"
            f"{nk}/{self.n} ROIs gardés{warn}", fontsize=12)
        self.ax.set_xlabel("clic = garder/rejeter   |   a/n = tout/rien   |   b = fond   |   "
                           "l = numéros   |   ← → = FOV   |   q = quitter", fontsize=9)
        self.ax.set_xticks([]); self.ax.set_yticks([])
        self.fig.canvas.draw_idle()

    # ---------- interactions ----------
    def on_click(self, ev):
        if ev.inaxes is not self.ax or ev.xdata is None:
            return
        lab = int(self.masks[int(round(ev.ydata)), int(round(ev.xdata))])
        if lab > 0:
            self.keep[lab] = not self.keep.get(lab, True)
            self.draw()

    def go(self, step):
        self.save()
        j = self.i + step
        if 0 <= j < len(self.keys):
            self.i = j
            self.load()
        else:
            print("Fin de la liste.")

    def on_key(self, ev):
        if ev.key == "a":
            self.keep = {l: True for l in self.keep}; self.draw()
        elif ev.key == "n":
            self.keep = {l: False for l in self.keep}; self.draw()
        elif ev.key == "b":
            self.bg = (self.bg + 1) % len(BACKGROUNDS); self.draw()
        elif ev.key == "l":
            self.labels_on = not self.labels_on; self.draw()
        elif ev.key in ("right", "enter"):
            self.go(+1)
        elif ev.key == "left":
            self.go(-1)
        elif ev.key == "q":
            self.save(); print(f"Sélections enregistrées dans {D['selection']}"); plt.close(self.fig)


def main():
    keys = sys.argv[1:] or [fov_key(s, i) for s in SLICES for i in range(N_FOV)]
    keys = [k for k in keys if (D["masks"] / f"{k}_masks.npy").exists()]
    if not keys:
        sys.exit("Aucun masque trouvé — lance d'abord step2_cellpose.py")
    print(__doc__)
    Selector(keys)


if __name__ == "__main__":
    main()
