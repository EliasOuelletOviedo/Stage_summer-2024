import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

path = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Stage_summer-2024/SDT.csv"

data = np.append([1], np.array(pd.read_csv(path)))

data = np.reshape(data, (1024, 1024, 256))

print(data.shape)

plt.figure()
ax = plt.gca()
im = ax.imshow(np.matrix(np.sum(data[::2,:,:], axis=2)), cmap='gray')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.2)

plt.colorbar(im, cax=cax)
plt.show()