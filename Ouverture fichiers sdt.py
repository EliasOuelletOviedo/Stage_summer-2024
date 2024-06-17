from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Billes/billes.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_1/pollen.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_2/pollen_2.sdt"
path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-24/Neurone1/2024-05-24-neurone_1.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers/1.sdt"

file = SdtFile(path)

print(np.array(file.data).shape)

array = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256)).sum(2)[::2,]
array = array[~np.all(array == 0, axis=1)]
array = array[-513:-1:,]

plt.figure()
ax = plt.gca()
im = ax.imshow(np.matrix(array), cmap='gray')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.2)

plt.colorbar(im, cax=cax)
plt.show()
