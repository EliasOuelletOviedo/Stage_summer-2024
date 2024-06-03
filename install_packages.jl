using Pkg

# Install necessary Julia packages
Pkg.add(["PyCall", "Conda", "FFTW", "Statistics", "Distributed", "ProgressBars", "Plots"])

# Install necessary Python packages
using Conda
Conda.add("scipy")
Conda.add("sdtfile")

Pkg.build("CodecZlib")