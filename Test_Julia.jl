using Gtk4
using PyCall
using FFTW
using Statistics
using Distributed
using DelimitedFiles

time = readdlm("lifetimes.csv", ',', Float64, '\n')
# time = readdlm("times.csv", ',', Float64, '\n')
# approx_time = zeros(0)
# for (i, x) in time
#     append!(approx_time, x*length(time)/i)
# end

py"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

data = $time
# data_2 = $approx_time
plt.figure()
ax = plt.gca()


im = ax.imshow(np.matrix(data), cmap='gray')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.2)
plt.colorbar(im, cax=cax)


# im = ax.plot(np.matrix(data))


plt.show()


# im = ax.plot(np.matrix(data_2))


# plt.show()
"""