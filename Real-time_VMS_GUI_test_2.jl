using PyCall
using FFTW
using Statistics
using Distributed
using DelimitedFiles
using ProgressBars
using DSP
using StructTypes
using Base64
using ZipFile
using BenchmarkTools
using Makie
using GLMakie
using Observables
using Dates
using Colors
using FileIO
using Images
using Clipper
using Luxor

function start_function(C1_path, C2_path, axis_C1_plot, axis_C2_plot)
    running[] = true

    intensity_C1 = Observable(Float64[])
    intensity_C2 = Observable(Float64[])
    real_time_C1 = Observable(Float64[])
    real_time_C2 = Observable(Float64[])

    intensity_plot_C1 = lines!(axis_C1_plot, real_time_C1, intensity_C1, color = Colors.JULIA_LOGO_COLORS.red, linewidth=2)
    intensity_plot_C2 = lines!(axis_C2_plot, real_time_C2, intensity_C2, color = Colors.JULIA_LOGO_COLORS.red, linewidth=2)

    println("Button pressed")
    @async get_intensity(C1_path, C2_path, intensity_C1, intensity_C2, real_time_C1, real_time_C2, axis_C1_plot)
end

function get_intensity(C1_path, C2_path, intensity_C1, intensity_C2, real_time_C1, real_time_C2, axis_C1_plot)
    previous_files_C1 = readdir(C1_path)
    previous_files_C2 = readdir(C2_path)

    t_0 = now()

    i = 0

    println("Search started")

    while running[]
        current_files_C1 = readdir(C1_path)
        current_files_C2 = readdir(C2_path)

        new_files_C1 = setdiff(current_files_C1, previous_files_C1)
        new_files_C2 = setdiff(current_files_C2, previous_files_C2)

        if !isempty(new_files_C1)

            println("File detected!")

            for new_file_C1 in new_files_C1

            # new_file_C1 = new_files_C1[end]

                file_C1 = rawview(channelview(load(joinpath(C1_path, new_file_C1))))'

                if i == 0
                    t_0 = now()
                    i += 1
                end

                present_time = (now() - t_0)/Millisecond(1)/1000

                default_img[] = file_C1
                push!(intensity_C1[], sum(file_C1))
                push!(real_time_C1[], present_time)

                notify(default_img)
                notify(intensity_C1)
                notify(real_time_C1)

                println(minimum(intensity_C1[]))

                # limits!(axis_C1_plot, max(0, present_time-10), max(10, present_time), minimum(intensity_C1[]), maximum(intensity_C1[]))
                # limits!(axis_C1_plot,  max(0, present_time-10), max(10, present_time), maximum(intensity_C1[])*0.9, maximum(intensity_C1[])*1.1)

                println(maximum(intensity_C1[]))

                sleep(0.01)
            end
        end

        if !isempty(new_files_C2)

            println("File detected!")

            for new_file_C2 in new_files_C2

            # new_file_C2 = new_files_C2[end]

                file_C2 = rawview(channelview(load(joinpath(C2_path, new_file_C2))))'

                if i == 0
                    t_0 = now()
                    i += 1
                end

                present_time = (now() - t_0)/Millisecond(1)/1000

                default_img[] = file_C2
                push!(intensity_C2[], sum(file_C2))
                push!(real_time_C2[], present_time)

                notify(default_img)
                notify(intensity_C2)
                notify(real_time_C2)

                # println(i)
                sleep(0.01)
            end
        end

        previous_files_C1 = current_files_C1
        previous_files_C2 = current_files_C2

        # if save[] == true
        #     filename = tb_file_name.displayed_string[]*".csv"

        #     if filename == "Enter a file name.csv"
        #         filename = split(path, "/")[end]*".csv"
        #     end

        #     open(filename, "w") do io
        #         writedlm(io, [string(now())])
        #         writedlm(io, [""])
        #         writedlm(io, [real_time[] lifetime[] photons[] concentration[]], ",")
        #     end

        #     save[] = false
        # end

        # println("1")

        sleep(0.1)
    end
end


image_colors = [RGB(1,0.25,0.25), RGB(0.25,1,0.25), RGB(0.25,0.25,1), RGB(0.25,1,1), RGB(1,0.25,1), RGB(1,1,0.25)]
border_colors = [RGBA(1,0.25,0.25,0.9), RGBA(0.25,1,0.25,0.9), RGBA(0.25,0.25,1,0.9), RGBA(0.25,1,1,0.9), RGBA(1,0.25,1,0.9), RGBA(1,1,0.25,0.9)]
fill_colors = [RGBA(1,0.25,0.25,0.1), RGBA(0.25,1,0.25,0.1), RGBA(0.25,0.25,1,0.1), RGBA(0.25,1,1,0.1), RGBA(1,0.25,1,0.1), RGBA(1,1,0.25,0.1)]
reverse_border_colors = [RGBA(0.25,1,1,0.9), RGBA(1,0.25,1,0.9), RGBA(1,1,0.25,0.9), RGBA(1,0.25,0.25,0.9), RGBA(0.25,1,0.25,0.9), RGBA(0.25,0.25,1,0.9)]
reverse_fill_colors = [RGBA(0.25,1,1,0.1), RGBA(1,0.25,1,0.1), RGBA(1,1,0.25,0.1), RGBA(1,0.25,0.25,0.1), RGBA(0.25,1,0.25,0.1), RGBA(0.25,0.25,1,0.1)]

set_theme!(backgroundcolor = :gray95)

default_img = Observable(zeros(1024,512))
default_img[][1,1] = 1

folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-18_TestSoftware/test save online-001/Bliq VMS"

C1_path = joinpath(folderpath, "C1")
C2_path = joinpath(folderpath, "C2")

fig = Figure(size = (1440, 847))

fig[1, 1] = figgrid = GridLayout()
fig[1, 2] = plotgrid = GridLayout()

axis_C1_fig = figgrid[1, 1][1, 1] = GLMakie.Axis(fig, yreversed = true, aspect = DataAspect(), xrectzoom=false, yrectzoom=false)
axis_C2_fig = figgrid[2, 1][1, 1] = GLMakie.Axis(fig, yreversed = true, aspect = DataAspect(), xrectzoom=false, yrectzoom=false)
axis_C1_plot = plotgrid[1, 1] = GLMakie.Axis(fig, title="Real-time plot of the total photon counts", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=20, xlabelsize=18, ylabelsize=18)
axis_C2_plot = plotgrid[2, 1] = GLMakie.Axis(fig, tellheight=false)

limits!(axis_C1_plot, 0, 10, 0, 1000000)

colorbar_C1 = figgrid[1, 1][1, 2] = Colorbar(fig, halign=:left, limits = (0, 255), colormap=[:black, RGB(1,0.25,0.25)])
colorbar_C2 = figgrid[2, 1][1, 2] = Colorbar(fig, halign=:left, limits = (0, 255), colormap=[:black, RGB(0.25,1,0.25)])

C1_fig = image!(axis_C1_fig, default_img, colormap=[:black, RGB(1,0.25,0.25)])
C2_fig = image!(axis_C2_fig, default_img, colormap=[:black, RGB(0.25,1,0.25)])

region_points = Observable(Point2f[])

region_plot_C1 = poly!(axis_C1_fig, region_points, color=RGBA(0.25,1,0.25,0.1), strokecolor=RGBA(0.25,1,0.25,0.9), strokewidth=2)
region_plot_C2 = poly!(axis_C2_fig, region_points, color=RGBA(1,0.25,0.25,0.1), strokecolor=RGBA(1,0.25,0.25,0.9), strokewidth=2)

region_scatter = scatter!(axis_C1_fig, region_points, markersize = 12, overdraw=true, color=RGBA(0.25,1,0.25,0.9))

fig[2, 1] = figbuttongrid = GridLayout(tellwidth=false)
fig[2, 2] = plotbuttongrid = GridLayout(tellwidth=false)

button_start = figbuttongrid[1, 1] = Button(fig, label="START", font=:bold, height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)
button_reset = figbuttongrid[2, 1] = Button(fig, label="RESET", font=:bold, height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)

menu_color_C1 = figbuttongrid[1, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow"], default = "Red", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black, tellheight=false)
menu_color_C2 = figbuttongrid[2, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow"], default = "Green", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black, tellheight=false)

button_ROI = figbuttongrid[1, 3] = Button(fig, label="Set ROI", height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)
button_noise = figbuttongrid[2, 3] = Button(fig, label="Set noise", height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)

toggle_ROI = figbuttongrid[1, 4] = Toggle(fig, height = 20, width = 40, buttoncolor=RGBf(0.4, 0.4, 0.4), framecolor_active=RGBf(0.6, 0.6, 0.6), framecolor_inactive=RGBf(0.8, 0.8, 0.8))
toggle_noise = figbuttongrid[2, 4] = Toggle(fig, height = 20, width = 40, buttoncolor=RGBf(0.4, 0.4, 0.4), framecolor_active=RGBf(0.6, 0.6, 0.6), framecolor_inactive=RGBf(0.8, 0.8, 0.8))

label_ROI = figbuttongrid[1, 5] = Label(fig, "ROI")
label_noise = figbuttongrid[2, 5] = Label(fig, "Noise")

tb_file_name = plotbuttongrid[1, 1] = Textbox(fig, placeholder="Enter a file name", textcolor=RGBf(0.4, 0.4, 0.4), textcolor_placeholder=RGBf(0.6, 0.6, 0.6),height = 40, width = 397, borderwidth=2, bordercolor=RGBf(0.8, 0.8, 0.8), bordercolor_focused=RGBf(0.4, 0.4, 0.4), bordercolor_hover=RGBf(0.6, 0.6, 0.6), boxcolor=:white, boxcolor_focused=:white, boxcolor_hover=:white, fontsize=20, cornerradius=8)
button_save = plotbuttongrid[1, 2] = Button(fig, label="Save", height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)

dragging = false
idx = 1

ROI_matrix = Observable(zeros(512,1024))
noise_matrix = Observable(zeros(512,1024))

plots = Observable(zeros(3))

running = Observable(false)

on(button_start.clicks) do _
    start_function(C1_path, C2_path, axis_C1_plot, axis_C2_plot)
end

on(button_reset.clicks) do _
    running[] = false
    println("Button pressed")
end

on(button_ROI.clicks) do _
    if !isempty(region_points[])
        poly = deepcopy(region_points[])

        ROI_plot_C1 = poly!(axis_C1_fig, poly, color=reverse_fill_colors[menu_color_C1.i_selected[]], strokecolor=reverse_border_colors[menu_color_C1.i_selected[]], strokewidth=2)
        ROI_plot_C2 = poly!(axis_C2_fig, poly, color=RGBA(1,0.25,1,0.1), strokecolor=RGBA(1,0.25,1,0.9), strokewidth=2)

        polygon_coords = []

        for _ in 1:length(region_points[])
            x, y = pop!(region_points[])
            # println(pop!(region_points[]))
            push!(polygon_coords, (x, y))
        end

        polygon_points = [Luxor.Point(x, y) for (x, y) in polygon_coords]

        for i in CartesianIndices(ROI_matrix[])
            point = Luxor.Point(i[2], i[1])

            if isinside(point, polygon_points)
                ROI_matrix[][i] = 1
            end
        end

        println(image_colors[menu_color_C1.i_selected[]])

        # image!(axis_C1_plot, img_C1.*ROI_matrix[]', colormap=[:black, :white])

        push!(plots[], 1)
        notify(region_points)
    end
end

on(button_noise.clicks) do _
    if !isempty(region_points[])
        poly = deepcopy(region_points[])

        noise_plot_C1 = poly!(axis_C1_fig, poly, color=RGBA(1,1,1,0.1), strokecolor=RGBA(1,1,1,0.9), strokewidth=2)
        noise_plot_C2 = poly!(axis_C2_fig, poly, color=RGBA(1,1,1,0.1), strokecolor=RGBA(1,1,1,0.9), strokewidth=2)
        
        polygon_coords = []

        for _ in 1:length(region_points[])
            x, y = pop!(region_points[])
            push!(polygon_coords, (x, y))
        end

        polygon_points = [Luxor.Point(x, y) for (x, y) in polygon_coords]

        for i in CartesianIndices(noise_matrix[])
            point = Luxor.Point(i[2], i[1])

            if isinside(point, polygon_points)
                noise_matrix[][i] = 1
            end
        end

        # image!(axis_C2_plot, img_C1.*noise_matrix[]', colormap=[:black, :white])

        push!(plots[], 0)
        notify(region_points)
    end
end

on(events(fig).mousebutton, priority = 2) do event
    global dragging, idx
    if running[]
        if event.button == Mouse.left
            if event.action == Mouse.press
                plt, i = pick(fig)

                if Keyboard.d in events(fig).keyboardstate && plt == region_scatter
                    deleteat!(region_points[], i)
                    notify(region_points)

                    return Consume(true)
                elseif Keyboard.a in events(fig).keyboardstate
                    push!(region_points[], mouseposition(axis_C1_fig))
                    notify(region_points)

                    return Consume(true)
                else
                    dragging = plt == region_scatter
                    idx = i

                    return Consume(dragging)
                end
            elseif event.action == Mouse.release
                dragging = false

                return Consume(false)
            end
        end
        return Consume(false)
    end
end

on(events(fig).mouseposition, priority = 2) do mp
    if dragging
        region_points[][idx] = mouseposition(axis_C1_fig)
        notify(region_points)

        return Consume(true)
    end
    return Consume(false)
end

display(fig)