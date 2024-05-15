from tkinter import filedialog
from tkinter import *
from PIL import Image as Im
from os import walk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as LSC


nb = int(input("Nombre de fichiers : "))

# Dictionnaire de correspondance entre les valeurs RGB et CMY
conversion_dict = {"R": (1, 0, 0),
                   "G": (0, 1, 0),
                   "B": (0, 0, 1),
                   "C": (0, 0.5, 0.5),
                   "M": (0.5, 0, 0.5),
                   "Y": (0.5, 0.5, 0),
                   }

# conversion_dict = {"R": [1, 0, 0],
#                    "G": [0, 1, 0],
#                    "B": [0, 0, 1],
#                    "C": [0, 0.5, 0.5],
#                    "M": [0.5, 0, 0.5],
#                    "Y": [0.5, 0.5, 0],
#                    }

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

    if i == 0:
        w, h = Im.open(path + "/" + files[0]).size
        N = len(files)
        stack = np.zeros((h, w, nb), float)

        print(np.max(np.array(Im.open(path + "/" + files[0]), dtype = float)))
    else:
        stack[...,i] = single

    single = np.zeros((h, w),float)

    for file in files:
        single = single + np.array(Im.open(path + "/" + file), dtype = float)

    # Affichage d'un channel en format Matplotlib
    cmap = LSC.from_list("Custom", [(0, 0, 0), conversion_dict[colors[i]]])
    im = axes[0, i].imshow(single, cmap = cmap)
    plt.colorbar(im, ax = axes[0, i])

plt.show()