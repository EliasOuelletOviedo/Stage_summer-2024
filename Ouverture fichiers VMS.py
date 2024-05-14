from PIL import Image as Im
from os import walk
import numpy as np

path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Données/2024-05-13/VMS/2024-05-13_pollen/pollen-001/Bliq VMS/"
folders = next(walk(path))[1]

conversion_dict = {"R": [1, 0, 0],
                   "G": [0, 1, 0],
                   "B": [0, 0, 1],
                   "C": [0, 0.5, 0.5],
                   "M": [0.5, 0, 0.5],
                   "Y": [0.5, 0.5, 0],
                   }

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
    single_im = Im.fromarray(single_int, mode = "RGB")
    single_im.show()

stack += single
stack_int = np.array(np.round(stack), dtype = np.uint8)

stack_im = Im.fromarray(stack_int, mode = "RGB")
stack_im.show()