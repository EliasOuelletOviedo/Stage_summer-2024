using Pkg

Pkg.add(["Colors", "Base64", "StructTypes", "DSP", "ProgressBars", "PyCall", "FFTW", "Statistics", "Distributed", "NativeFileDialog", "ZipFile", "Observables", "Dates", "DelimitedFiles", "LibSerialPort", "DataStructures", "BenchmarkTools", "Conda"])

using Conda
Conda.add("scipy")
# Conda.add("scipy.optimize")
Conda.add("sdtfile")