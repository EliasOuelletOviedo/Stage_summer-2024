from PIL import Image as Im
from os import walk
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/VMS/Pollen/"
folders = next(walk(path))[1]

# Dictionnaire de correspondance entre les valeurs RGB et CMY
conversion_dict = {"R": [1, 0, 0],
                   "G": [0, 1, 0],
                   "B": [0, 0, 1],
                   "C": [0, 0.5, 0.5],
                   "M": [0.5, 0, 0.5],
                   "Y": [0.5, 0.5, 0],
                   }

# Couleurs voulues pour l'affichage
colors = ["R", "G"]

for i, folder in enumerate(folders):
    files = next(walk(path + folder + "/"))[2]

    if i == 0:
        w, h = Im.open(path + folder + "/" + files[0]).size
        N = len(files)
        stack = np.zeros((h, w, 3), float)
    else:
        stack += single

    single = np.zeros((h, w, 3),float)

    for file in files:
        single_arr = np.array(Im.open(path + folder + "/" + file), dtype = float)

        for j, v in enumerate(conversion_dict[colors[i]]):
            single[..., j] = single[..., j] + v*single_arr/N

    single_int = np.array(np.round(single), dtype = np.uint8)

    # Affichage d'un channel en format PNG
    single_im = Im.fromarray(single_int, mode = "RGB")
    # single_im.show()

    # Affichage d'un channel en format Matplotlib
    plt.imshow(single_int)
    # plt.show()

stack += single
stack_int = np.array(np.round(stack), dtype = np.uint8)

# Affichage de la superposition des channels en format PNG
stack_im = Im.fromarray(stack_int, mode = "RGB")
# stack_im.show()

# Affichage de la superposition des channels en format Matplotlib
plt.imshow(stack_int)
plt.show()
