using Pkg

# Install necessary Julia packages
Pkg.add(["PyCall", "Conda", "FFTW", "Statistics", "Distributed", "ProgressBars", "StructTypes", "Base64", "DelimitedFiles", "ZipFile", "Gtk4", "DSP", "BenchmarkTools", "Plots", "GR"])

# Install necessary Python packages
using Conda
Conda.add("scipy")
Conda.add("sdtfile")
