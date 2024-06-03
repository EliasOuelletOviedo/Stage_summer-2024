using Plots

gr()

file = open("lifetimes.csv", "r")

data = read(file)

heatmap(data, c=cgrad([:black, :white]), yflip = true)
