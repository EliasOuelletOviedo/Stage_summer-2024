"""Étape 3 — sélection manuelle des ROIs Cellpose, acquisition par acquisition.

Usage :  python3 step3_select.py <dossier> [acquisition]
         python3 step3_select.py <dossier> --tout   (sans fenêtre : garde tous les ROIs)
         python3 step3_select.py <dossier> --top-photons 50   (sans fenêtre : dans chaque acquisition, les
             ROIs dont le total de photons CFP + YFP est dans les 50 % les plus élevés)

Commandes
  clic sur un ROI  garder / rejeter
  a / n            tout garder / tout rejeter
  b                change le fond (√photons, CFP, YFP, rapport YFP/CFP)
  l                affiche ou masque les numéros
  ← →              acquisition précédente / suivante (enregistre au passage)
  q                enregistre et quitte

Tous les ROIs sont gardés par défaut : il suffit de cliquer ceux à rejeter. La sélection est
écrite dans selection/<acq>_selection.json ; relancer l'étape 6 suffit ensuite (l'étape 4
mesure tous les ROIs une fois pour toutes).
"""
import json
import sys

import numpy as np
import matplotlib

import config as C


def ecrire(D, nom, n, keep):
    (D["selection"] / f"{nom}_selection.json").write_text(json.dumps(
        {"n_rois": n, "kept": sorted(l for l, v in keep.items() if v),
         "keep": {str(k): bool(v) for k, v in keep.items()}}, indent=1))


def lire(D, nom, n):
    p = D["selection"] / f"{nom}_selection.json"
    keep = {l: True for l in range(1, n + 1)}
    if p.exists():
        keep.update({int(k): v for k, v in json.loads(p.read_text())["keep"].items()})
    return keep


class Selecteur:
    def __init__(self, root, A, depart=0):
        matplotlib.use("macosx")
        import matplotlib.pyplot as plt
        from matplotlib.colors import to_rgba
        self.plt, self.to_rgba = plt, to_rgba
        self.D, self.A, self.i = C.dirs(root), A, depart
        self.bg, self.labels_on = 0, True
        self.fig, self.ax = plt.subplots(figsize=(15, 8.5))
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.charger()
        plt.show()

    def charger(self):
        nom = self.A[self.i]["nom"]
        self.masks = np.load(self.D["masks"] / f"{nom}_masks.npy")
        self.n = int(self.masks.max())
        self.keep = lire(self.D, nom, self.n)
        s = np.load(self.D["projections"] / f"{nom}.npz")["somme"]
        from scipy import ndimage as ndi
        lis = [ndi.uniform_filter(s[c], 5) for c in range(2)]
        self.fonds = [("√photons CFP + YFP", np.sqrt(s.sum(0))), ("CFP", s[0]), ("YFP", s[1]),
                      ("rapport YFP / CFP (lissé 5 px)", np.where(lis[0] > 0.5, lis[1] / np.maximum(lis[0], 1e-9), np.nan))]
        self.centres = {}
        for lab in range(1, self.n + 1):
            ys, xs = np.nonzero(self.masks == lab)
            if len(ys):
                self.centres[lab] = (xs.mean(), ys.mean())
        self.draw()

    def draw(self):
        self.ax.clear()
        label, img = self.fonds[self.bg]
        v = img[np.isfinite(img)]
        lo, hi = np.percentile(v, (1, 99.7))
        self.ax.imshow(np.nan_to_num(img, nan=lo), cmap="gray" if self.bg < 3 else "viridis", vmin=lo, vmax=max(hi, lo + 1e-6))
        for lab in range(1, self.n + 1):
            k = self.keep.get(lab, True)
            col = "#2ecc71" if k else "#e74c3c"
            self.ax.contour(self.masks == lab, levels=[0.5], colors=[col], linewidths=1.3 if k else 0.7,
                            linestyles="solid" if k else "dotted")
            if self.labels_on and lab in self.centres:
                x, y = self.centres[lab]
                self.ax.text(x, y, str(lab), color=col, fontsize=7, ha="center", va="center",
                             fontweight="bold", bbox=dict(boxstyle="round,pad=0.1", fc=self.to_rgba("k", .45), ec="none"))
        nk = sum(self.keep.values())
        self.ax.set_title(f"{self.A[self.i]['nom']} ({self.i + 1}/{len(self.A)}) — fond : {label}\n"
                          f"{nk}/{self.n} ROIs gardés", fontsize=12)
        self.ax.set_xlabel("clic = garder/rejeter | a/n = tout/rien | b = fond | l = numéros | "
                           "← → = acquisition | q = quitter", fontsize=9)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.fig.canvas.draw_idle()

    def enregistrer(self):
        ecrire(self.D, self.A[self.i]["nom"], self.n, self.keep)

    def on_click(self, ev):
        if ev.inaxes is not self.ax or ev.xdata is None:
            return
        lab = int(self.masks[int(round(ev.ydata)), int(round(ev.xdata))])
        if lab > 0:
            self.keep[lab] = not self.keep.get(lab, True)
            self.draw()

    def on_key(self, ev):
        if ev.key in ("a", "n"):
            self.keep = {l: ev.key == "a" for l in self.keep}
            self.draw()
        elif ev.key == "b":
            self.bg = (self.bg + 1) % len(self.fonds)
            self.draw()
        elif ev.key == "l":
            self.labels_on = not self.labels_on
            self.draw()
        elif ev.key in ("left", "right"):
            self.enregistrer()
            self.i = (self.i + (1 if ev.key == "right" else -1)) % len(self.A)
            self.charger()
        elif ev.key == "q":
            self.enregistrer()
            print(f"Sélections enregistrées dans {self.D['selection']}")
            self.plt.close(self.fig)


def main():
    root = C.racine()
    D = C.dirs(root)
    A = [a for a in C.acquisitions(root) if (D["masks"] / f"{a['nom']}_masks.npy").exists()]
    if not A:
        sys.exit("Aucun masque trouvé — lance d'abord step2_cellpose.py")
    if "--tout" in sys.argv[1:]:
        for a in A:
            n = int(np.load(D["masks"] / f"{a['nom']}_masks.npy").max())
            ecrire(D, a["nom"], n, {l: True for l in range(1, n + 1)})
        print(f"  {len(A)} acquisitions : tous les ROIs gardés (sélection par défaut)")
        return
    if "--top-photons" in sys.argv[1:]:
        pct = float(sys.argv[sys.argv.index("--top-photons") + 1])
        for a in A:
            masks = np.load(D["masks"] / f"{a['nom']}_masks.npy")
            n = int(masks.max())
            somme = np.load(D["projections"] / f"{a['nom']}.npz")["somme"].sum(0)      # CFP + YFP, cycles gardés
            ph = np.bincount(masks.ravel(), weights=somme.ravel(), minlength=n + 1)[1:]
            seuil = np.percentile(ph, 100 - pct)
            ecrire(D, a["nom"], n, {l: bool(ph[l - 1] >= seuil) for l in range(1, n + 1)})
            print(f"  {a['nom']:>12} : {int((ph >= seuil).sum())}/{n} ROIs gardés (photons ≥ {seuil:.0f})")
        return
    noms = [a["nom"] for a in A]
    demande = [x for x in sys.argv[1:] if x in noms]
    print(__doc__)
    Selecteur(root, A, noms.index(demande[0]) if demande else 0)


if __name__ == "__main__":
    main()
