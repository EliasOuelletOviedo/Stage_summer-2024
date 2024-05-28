using Pkg

# Install necessary Julia packages
Pkg.add(["PyCall", "Conda", "FFTW", "Statistics", "Distributed", "ProgressBars", "DPS"])

# Install necessary Python packages
using Conda
Conda.add("scipy")