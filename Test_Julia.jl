using GLMakie

fig = Figure()

ax = Axis(fig[1, 1], xlabel = "x label", ylabel = "y label",
title = "Title")

fig[2, 1] = buttongrid = GridLayout(tellwidth = false)

button_IRF = buttongrid[1, 1] = Button(fig, label="Import IRF")
button_folder = buttongrid[1, 2] = Button(fig, label="Select folder")
button_play = buttongrid[2, 1:2] = Button(fig, label="Select folder")

fig