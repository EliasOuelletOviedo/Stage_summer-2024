from tkinter import filedialog
from tkinter import *
from PIL import Image as Im
from os import walk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as LSC


nb = int(input("Nombre de fichiers : "))

# Dictionnaire de correspondance entre les valeurs RGB et CMY
normal_conversion_dict = {"R": (1, 0, 0),
                   "G": (0, 1, 0),
                   "B": (0, 0, 1),
                   "C": (0, 0.5, 0.5),
                   "M": (0.5, 0, 0.5),
                   "Y": (0.5, 0.5, 0),
                   }

conversion_dict = {"R": (1, 0, 0),
                   "G": (0, 1, 0),
                   "B": (0, 0, 1),
                   "C": (0, 1, 1),
                   "M": (1, 0, 1),
                   "Y": (1, 1, 0),
                   }

# Couleurs voulues pour l'affichage
colors = []

paths = []

for n in range(nb):
    color = input(f"Couleur associée au fichier {n+1} : ")
    colors.append(color)
    
    root = Tk()
    root.withdraw()
    paths.append(filedialog.askdirectory())

fig, axes = plt.subplots(nrows=nb, ncols=nb)

for i, path in enumerate(paths):
    files = next(walk(path))[2]
    N = len(files)

    if i == 0:
        w, h = Im.open(path + "/" + files[0]).size
        stack = np.zeros((h, w, nb), float)
        RGB_stack = np.zeros((h, w, 3), float)

    single = np.zeros((h, w),float)

    for file in files:
        stack[...,i] = stack[...,i] + np.array(Im.open(path + "/" + file), dtype = float)/N

    for k, v in enumerate(conversion_dict[colors[i]]):
        RGB_stack[...,k] = RGB_stack[...,k] + v*stack[...,i]

    # Affichage d'un channel en format Matplotlib
    cmap = LSC.from_list("Custom", [(0, 0, 0), conversion_dict[colors[i]]])
    im = axes[0, i].imshow(stack[...,i], cmap = cmap)
    plt.colorbar(im, ax = axes[0, i])

stack = stack[...,0] - stack[...,1]
cmap = LSC.from_list("Custom", [conversion_dict[colors[1]], (0, 0, 0), conversion_dict[colors[0]]])
im = axes[1, 1].imshow(stack, cmap = cmap)
im.set_clim(-200,200)
plt.colorbar(im, ax = axes[1, 1])

RGB_stack_int = np.array(np.round(RGB_stack), dtype = np.uint8)
axes[1, 0].imshow(RGB_stack_int)

plt.show()