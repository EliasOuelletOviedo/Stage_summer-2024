"""Étape 3 — sélection manuelle des ROIs Cellpose.

Usage :  python3 step3_select.py <dossier>
         python3 step3_select.py <dossier> --tout   (sans fenêtre : garde tous les ROIs)

Commandes
  clic sur un ROI  garder / rejeter
  a / n            tout garder / tout rejeter
  b                change le fond (FLIM ch1 + ch2, ch1, ch2, ratio, images Bliq)
  l                affiche ou masque les numéros de ROI
  q                enregistre et quitte

Tous les ROIs sont gardés par défaut : il suffit de cliquer ceux à rejeter.
Les ROIs FLIMApp (cyan) ne se trient pas : ils sont toujours mesurés.
"""
import json
import sys

import numpy as np
import matplotlib

import config


def charger(root):
    D = config.dirs(root)
    masks = np.load(D["masks"] / "cellpose_masks.npy")
    return D, masks, int(masks.max())


def ecrire(D, n, keep):
    (D["selection"] / "cellpose_selection.json").write_text(json.dumps(
        {"n_rois": n, "kept": sorted(l for l, v in keep.items() if v),
         "keep": {str(k): bool(v) for k, v in keep.items()}}, indent=1))


def keep_all(root):
    """Sélection par défaut sans passer par la fenêtre : tous les ROIs gardés."""
    D, masks, n = charger(root)
    ecrire(D, n, {l: True for l in range(1, n + 1)})
    print(f"  {n} ROIs gardés (sélection par défaut)")


class Selector:
    def __init__(self, root):
        matplotlib.use("macosx")
        import matplotlib.pyplot as plt
        from matplotlib.colors import to_rgba
        self.plt, self.to_rgba = plt, to_rgba
        self.D, self.masks, self.n = charger(root)
        self.polys = config.session(root)["rois_flimapp"]
        f = np.load(self.D["projections"] / "flim.npz")
        red = lambda im: np.kron(im.reshape(256, 4, 256, 4).sum((1, 3)), np.ones((4, 4)))
        self.backgrounds = [("FLIM ch1 + ch2", f["somme"]), ("FLIM ch1 (CFP)", f["intensite_ch1"]),
                            ("FLIM ch2 (YFP)", f["intensite_ch2"]),
                            ("ratio ch2 / ch1", red(f["intensite_ch2"]) / np.maximum(red(f["intensite_ch1"]), 1))]
        p = self.D["projections"] / "fovs.npz"
        if p.exists():
            for k, v in np.load(p).items():
                if k.endswith("C2"):
                    self.backgrounds.append((f"Bliq {k.replace('_', ' ')}", v))
        self.bg, self.labels_on = 0, True
        self.centroids = {}
        for lab in range(1, self.n + 1):
            ys, xs = np.nonzero(self.masks == lab)
            if len(ys):
                self.centroids[lab] = (xs.mean(), ys.mean())
        path = self.D["selection"] / "cellpose_selection.json"
        self.keep = {lab: True for lab in range(1, self.n + 1)}
        if path.exists():
            saved = json.loads(path.read_text())
            self.keep.update({int(k): v for k, v in saved["keep"].items()})
        self.fig, self.ax = plt.subplots(figsize=(11, 11))
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.draw()
        plt.show()

    def draw(self):
        self.ax.clear()
        label, img = self.backgrounds[self.bg]
        img = np.nan_to_num(img, nan=0.0)
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
                self.ax.text(x, y, str(lab), color=col, fontsize=7, ha="center", va="center",
                             fontweight="bold",
                             bbox=dict(boxstyle="round,pad=0.1", fc=self.to_rgba("k", .45), ec="none"))
        for r in self.polys:
            self.ax.plot(r["xs"], r["ys"], color="#00e5ff", lw=1.0, ls="--")
        nk = sum(self.keep.values())
        self.ax.set_title(f"fond : {label}\n{nk}/{self.n} ROIs gardés   (cyan : ROIs FLIMApp)",
                          fontsize=12)
        self.ax.set_xlabel("clic = garder/rejeter   |   a/n = tout/rien   |   b = fond   |   "
                           "l = numéros   |   q = quitter", fontsize=9)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.fig.canvas.draw_idle()

    def on_click(self, ev):
        if ev.inaxes is not self.ax or ev.xdata is None:
            return
        lab = int(self.masks[int(round(ev.ydata)), int(round(ev.xdata))])
        if lab > 0:
            self.keep[lab] = not self.keep.get(lab, True)
            self.draw()

    def on_key(self, ev):
        if ev.key == "a":
            self.keep = {l: True for l in self.keep}
            self.draw()
        elif ev.key == "n":
            self.keep = {l: False for l in self.keep}
            self.draw()
        elif ev.key == "b":
            self.bg = (self.bg + 1) % len(self.backgrounds)
            self.draw()
        elif ev.key == "l":
            self.labels_on = not self.labels_on
            self.draw()
        elif ev.key == "q":
            ecrire(self.D, self.n, self.keep)
            print(f"Sélection enregistrée dans {self.D['selection']}")
            self.plt.close(self.fig)


def main():
    root = config.racine()
    if not (root / "masks" / "cellpose_masks.npy").exists():
        sys.exit("Aucun masque trouvé — lance d'abord step2_cellpose.py")
    if "--tout" in sys.argv[1:]:
        keep_all(root)
        return
    print(__doc__)
    Selector(root)


if __name__ == "__main__":
    main()
