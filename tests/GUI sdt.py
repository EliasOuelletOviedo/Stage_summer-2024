from tkinter import *
from tkinter import ttk, filedialog
from tkinter.filedialog import askopenfile
from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
import os

def openfile():
    path = filedialog.askopenfile()
    file = SdtFile(os.path.abspath(path.name))

    array = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256)).sum(2)[::2,]
    array = array[~np.all(array == 0, axis=1)]
    array = array[-513:-1:,]

    fig = Figure(figsize = (12, 7), dpi = 100) 

    plot = fig.add_subplot(111) 
  
    im = plot.imshow(np.matrix(array), cmap='gray')
    plt.colorbar(im)

    canvas = FigureCanvasTkAgg(fig, master = root)   
    canvas.draw() 
  
    canvas.get_tk_widget().pack() 
  
    toolbar = NavigationToolbar2Tk(canvas, root) 
    toolbar.update() 
  
    canvas.get_tk_widget().pack(pady = 10)


root = Tk()
root.geometry("1440x900")

# ttk.Button(root, text="Delete", command=delete).pack(side = BOTTOM, pady = 10)
ttk.Button(root, text="Open", command=openfile).pack(side = BOTTOM, pady = 0)

root.mainloop()
