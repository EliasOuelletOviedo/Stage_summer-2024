"""Étape 3 — sélection manuelle des ROIs, un FOV à la fois.

Usage :  python3 step3_select.py  [FOV1 ...]
         python3 step3_select.py --tout   (sans fenêtre : garde tous les ROIs)

Commandes
  clic sur un ROI  garder / rejeter
  a / n            tout garder / tout rejeter
  b                change le fond (somme, puis chaque acquisition)
  l                affiche ou masque les numéros de ROI
  flèche droite    enregistre et passe au FOV suivant
  flèche gauche    enregistre et revient au FOV précédent
  q                enregistre et quitte

Tous les ROIs sont gardés par défaut : il suffit de cliquer ceux à rejeter.
La sélection vaut pour **toutes** les acquisitions du FOV, puisque les ROIs y sont
définis dans un repère commun (étape 1).
"""
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("macosx")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

from config import dirs, slug, survey

D = dirs()


class Selector:
    def __init__(self, keys):
        self.keys = keys
        self.data = survey()
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
        self.npz = np.load(D["projections"] / f"{key}.npz")
        self.masks = np.load(D["masks"] / f"{key}_masks.npy")
        self.n = int(self.masks.max())
        self.backgrounds = [("somme de toutes les acquisitions", "summed")]
        for rec, n, _ in self.data.get(key, []):
            name = f"proj_{slug(rec)}"
            if name in self.npz.files:
                self.backgrounds.append((f"{rec} ({n} images)", name))
        self.bg = min(self.bg, len(self.backgrounds) - 1)

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
        label, name = self.backgrounds[self.bg]
        img = np.nan_to_num(self.npz[name], nan=0.0)
        lo, hi = np.percentile(img, (1, 99.7))
        self.ax.imshow(img, cmap="gray", vmin=lo, vmax=max(hi, lo + 1e-6))

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
        self.ax.set_title(
            f"{self.key}   ({self.i+1}/{len(self.keys)})   fond : {label}\n"
            f"{nk}/{self.n} ROIs gardés", fontsize=12)
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
            self.bg = (self.bg + 1) % len(self.backgrounds); self.draw()
        elif ev.key == "l":
            self.labels_on = not self.labels_on; self.draw()
        elif ev.key in ("right", "enter"):
            self.go(+1)
        elif ev.key == "left":
            self.go(-1)
        elif ev.key == "q":
            self.save(); print(f"Sélections enregistrées dans {D['selection']}"); plt.close(self.fig)


def keep_all(keys):
    """Sélection par défaut sans passer par la fenêtre : tous les ROIs gardés.

    Permet de faire tourner l'étape 4 tout de suite ; le tri manuel peut être fait
    ensuite (le fichier est relu et complété), puis l'étape 4 relancée.
    """
    for key in keys:
        masks = np.load(D["masks"] / f"{key}_masks.npy")
        n = int(masks.max())
        (D["selection"] / f"{key}_selection.json").write_text(json.dumps(
            {"fov": key, "n_rois": n, "kept": list(range(1, n + 1)),
             "keep": {str(l): True for l in range(1, n + 1)}}, indent=1))
        print(f"  {key:10s} {n} ROIs gardés (sélection par défaut)")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    keys = args or list(survey())
    keys = [k for k in keys if (D["masks"] / f"{k}_masks.npy").exists()]
    if not keys:
        sys.exit("Aucun masque trouvé — lance d'abord step2_cellpose.py")
    if "--tout" in sys.argv[1:]:
        keep_all(keys)
        return
    print(__doc__)
    Selector(keys)


if __name__ == "__main__":
    main()
