from tkinter import *
from tkinter import ttk, filedialog
from tkinter.filedialog import askopenfile
from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt
import os

def openfile():
    path = filedialog.askopenfile()
    file = SdtFile(os.path.abspath(path.name))

    array = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256)).sum(2)[::2,]
    array = array[~np.all(array == 0, axis=1)]
    array = array[-513:-1:,]

    plt.imshow(np.matrix(array), cmap='gray')
    plt.colorbar()
    plt.show()

root = Tk()
root.geometry("700x350")

label = Label(root, text="Cliquez pour sélectrionner un fichier", font=('Cambria 13'))
label.pack(pady=10)

ttk.Button(root, text="Browse", command=openfile).pack(pady=20)

root.mainloop()
