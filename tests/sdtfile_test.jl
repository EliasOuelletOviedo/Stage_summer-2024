using Plots
using Images
using PyCall
using BenchmarkTools

sdt = pyimport("sdtfile")

include("SdtFile.jl")
using .SdtFile

filepath = "/Users/eliasouellet-oviedo/Documents/Stages/Stage2/Partie2/Données/2026-03-12/FOVs/fov_3.sdt"
# filepath = "/Users/eliasouellet-oviedo/Documents/Stages/Stage2/Partie2/Données/2026-03-12/cell1/60_1/single_c1.sdt"
# filepath = "/Users/eliasouellet-oviedo/Desktop/puffs_20ms_2_10X.sdt"

# @benchmark pyfile = sdt.SdtFile(filepath) samples=10 seconds=30
jlfile = SdtFile.read_sdt(filepath)

data = jlfile.data[1]

vector = sum(data, dims=2)
println(size(vector))

matrix = reshape(vector, 1024, :)
println(size(matrix))

# show matrix image
heatmap(matrix[:,begin:2:end], aspect_ratio=:equal)

# println("Python file: ", size(pyfile.data))
# println("Julia file: ", sum(jlfile.data[1]))
# println("Julia file: ", sum(jlfile.data[2]))


