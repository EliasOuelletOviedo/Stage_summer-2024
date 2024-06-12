using GLMakie
using Observables
using Threads

set_theme!(backgroundcolor = :gray90)

# fig = Figure(size = (600, 400))

# empty!(fig)

# fig

# sleep(10)

# button_IRF = buttongrid[1, 1] = Button(fig, label="Import IRF", width = 100)
# button_folder = buttongrid[1, 2] = Button(fig, label="Select folder", width = 100)

fig = Figure(size = (1440, 900))

fig[1, 1] = buttongrid = GridLayout(tellwidth = false)

button_start = buttongrid[1, 1] = Button(fig, label="START", width = 200)
button_stop = buttongrid[1, 2] = Button(fig, label="STOP", width = 200)

fig[2, 1] = axisgrid = GridLayout(tellwidth = false)

axis_count = axisgrid[1, 1] = Axis(fig, xlabel = "Time [ns]", ylabel = "y label", title = "Photon counts", aspect = 2, xgridvisible=false, ygridvisible=false)
axis_histogram = axisgrid[1, 2] = Axis(fig, xlabel = "Time [ns]", ylabel = "y label", title = "Counts histogram", aspect = 2, xgridvisible=false, ygridvisible=false)
axis_lifetime = axisgrid[2, 1] = Axis(fig, xlabel = "Time [ns]", ylabel = "y label", title = "Lifetime", aspect = 2, xgridvisible=false, ygridvisible=false)
axis_concentration = axisgrid[2, 2] = Axis(fig, xlabel = "Time [ns]", ylabel = "y label", title = "Chlorine concentration", aspect = 2, xgridvisible=false, ygridvisible=false)

fig