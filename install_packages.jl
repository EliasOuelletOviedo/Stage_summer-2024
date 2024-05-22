using Pkg

# Install necessary Julia packages
Pkg.add(["PyCall", "Conda", "FFTW", "Statistics", "Distributed", "Plots", "GR_jll", "GR"])

# Install necessary Python packages
using Conda
Conda.add("scipy")