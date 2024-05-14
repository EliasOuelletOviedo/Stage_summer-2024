from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt

# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Billes/billes.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_1/pollen.sdt"
# path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_2/pollen_2.sdt"
path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-05-13/FLIM/Pollen_3/pollen_3.sdt"

file = SdtFile(path)

if len(file.data[0]) == 1048576:
    array = np.squeeze(np.array(file.data)).reshape((1024, 1024, 256))

    plt.imshow(np.matrix(array.mean(2)[::2,:]), cmap='gray')
    plt.colorbar()
    plt.show()

elif len(file.data[0]) == 2097152:
    array = np.squeeze(np.array(file.data)).reshape((2048, 1024, 256))

    plt.imshow(np.matrix(array.mean(2)[24:1047:2,:]), cmap='gray')
    plt.colorbar()
    plt.show()

else:
    print("Erreur de format d'image")
