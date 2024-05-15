from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt

# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Billes/billes.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_1/pollen.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_2/pollen_2.sdt"
path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_3/pollen_3.sdt"

file = SdtFile(path)

array = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256)).sum(2)[::2,]
array = array[~np.all(array == 0, axis=1)]
array = array[-513:-1:,]

plt.imshow(np.matrix(array), cmap='gray')
plt.colorbar()
plt.show()

