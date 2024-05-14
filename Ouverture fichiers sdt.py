from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt

path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Données/2024-05-13/FLIM/2024-05-13_demo/pollen_2.sdt"
path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Données/2024-05-13/FLIM/2024-05-13_demo/billes.sdt"

file = SdtFile(path)

print(file)

array = np.squeeze(np.array(file.data)).reshape((1024, 1024, 256)).mean(2)[::2,:]

# plt.imshow(np.matrix(array.mean(2)), cmap='gray')
# plt.colorbar()
# plt.show()

plt.imshow(np.matrix(array), cmap='gray')
plt.colorbar()
plt.show()