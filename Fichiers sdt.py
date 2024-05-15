from tkinter import filedialog
from tkinter import *
from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt
import os

root = Tk()
root.withdraw()
path = filedialog.askopenfile()
file = SdtFile(os.path.abspath(path.name))

array = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256)).sum(2)[::2,]
array = array[~np.all(array == 0, axis=1)]
array = array[-513:-1:,]

plt.imshow(np.matrix(array), cmap='gray')
plt.colorbar()
plt.show()
