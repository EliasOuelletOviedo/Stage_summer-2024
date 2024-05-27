using Gtk4
using PyCall
using FFTW
using Statistics
using Distributed
using DelimitedFiles


lifetimes = readdlm("lifetimes.csv", ',', Float64, '\n')


if isfile("file_filepath.txt")
    filepath = open(f->read(f, String), "file_filepath.txt")
else
    filepath = open_dialog("File: Choose .sdt file to open.")
    open("file_filepath.txt", "w") do io
        write(io, filepath)
    end
end

py"""
from sdtfile import SdtFile
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

filepath = $filepath

file = SdtFile(filepath)
pyarray = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256))[::2,:,:]
pyarray = pyarray[:512]

lifetimes = $lifetimes

plt.figure()
ax = plt.gca()


# im = ax.imshow(np.matrix(lifetimes), cmap='gray')
im = ax.imshow(np.sum(pyarray, axis=2), cmap='gray')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.2)
plt.colorbar(im, cax=cax)

plt.show()
"""