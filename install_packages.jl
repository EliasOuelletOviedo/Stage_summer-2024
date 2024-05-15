using Pkg

# Install necessary Julia packages
Pkg.add(["PyCall", "Conda", "FFTW", "Statistics", "Distributed"])

# Install necessary Python packages
using Conda
Conda.add("scipy")