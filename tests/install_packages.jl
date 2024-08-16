using Pkg

# Install necessary Julia packages
Pkg.add(["PyCall", "Conda", "FFTW", "Statistics", "Distributed", "ProgressBars", "StructTypes", "Base64", "DelimitedFiles", "ZipFile", "Gtk4", "DSP", "BenchmarkTools", "Plots", "GR", "GLMakie", "Observables", "Dates", "Colors"])
Pkg.add(["GLMakie", "PyCall", "FFTW", "Statistics", "NativeFileDialog", "ZipFile", "Observables", "Dates", "DelimitedFiles"])

# Install necessary Python packages
using Conda
Conda.add("scipy")
Conda.add("sdtfile")
