from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from read_roi import read_roi_file



# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Billes/billes.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_1/pollen.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_2/pollen_2.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-24/Neurone1/2024-05-24-neurone_1.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers/1.sdt"

path = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie2/Données/2025-10-10/0/Iris/512/10/image.sdt"
roipath = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie2/Données/2025-10-10/ROI_512.roi"

file = SdtFile(path)
roi = read_roi_file(roipath)


print(roi)

# print(np.array(file.data).shape)

array = np.sum(np.squeeze(np.array(file.data)).reshape((-1, 1024, 256))[::2,:,:], axis = (0,1))
print(array.shape)
# array = array[~np.all(array == 0, axis=1)]
# array = array[-513:-1:,]

# plt.figure(figsize=(10,6))
# plt.imshow(array, cmap='gray', aspect='auto')
# plt.axis('off')
# plt.show()