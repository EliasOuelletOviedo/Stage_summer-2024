using Makie
using GLMakie
using Observables
using Dates
using FileIO
using Images
using Luxor
using Clipper
using NativeFileDialog
using Distributed

function validator(str)
    # if tryparse(UInt16, str) !== nothing
    #     if 2 <= parse(UInt16, str) <= 999
    #         return true
    #     end
    # end
    return tryparse(UInt16, str) !== nothing
end

@everywhere function get_folder_path()
    if isfile("VMS_folder_path.txt")
        folderpath = open(f->read(f, String), "VMS_folder_path.txt")

        if isempty(folderpath)
            folderpath = pick_folder()

            open("VMS_folder_path.txt", "w") do io
                write(io, folderpath)
            end
        end
    else
        folderpath = pick_folder()

        open("VMS_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    end

    if isempty(folderpath)
        println("Folder does not exist")

        folderpath = get_folder_path()
    end

    return folderpath
end

@everywhere function get_new_folder_path()
    folderpath = pick_folder()

    if !isempty(folderpath)
        open("VMS_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    else
        folderpath = open(f->read(f, String), "VMS_folder_path.txt")
    end

    return folderpath
end

@everywhere function verify_folder(folderpath)
    if isdir(joinpath(folderpath, "C1")) || isdir(joinpath(folderpath, "C2")) || isdir(joinpath(folderpath, "C3"))
        folderpath_correct = true
    elseif isdir(joinpath(folderpath, "Bliq VMS", "C1")) || isdir(joinpath(folderpath, "Bliq VMS", "C2")) || isdir(joinpath(folderpath, "Bliq VMS", "C3"))
        folderpath = joinpath(folderpath, "Bliq VMS")

        folderpath_correct = true
    else
        println("Folder error")
        folderpath_correct = false
        running[] = false

        button_start.label = "START"
        button_reset.label = "CLEAR"
    end

    return folderpath, folderpath_correct
end

@everywhere function get_file_metadata(filepath)
    file = String(read(filepath, 3680))
    framerate = parse(Float64, match(r"<M K=\"FrameRate\">(.*?)</M>", file).captures[1])
    sizex = parse(Int64, match(r" SizeX=\"(.*?)\"", file).captures[1])
    sizey = parse(Int64, match(r" SizeY=\"(.*?)\"", file).captures[1])

    return framerate, sizex, sizey
end

@everywhere function get_file_number(filename)
    # try
    #     println(fieldnames(typeof(match(r"T(.*?).tif", filename))))
    #     println(match(r"T(.*?).tif", filename).captures[1])
    #     println(parse(Int64, match(r"T(.*?).tif", filename).captures[1]))
    # catch e
    #     println(e)
    # end

    return match(r"T(.*?).tif", filename) !== nothing ? parse(Int64, match(r"T(.*?).tif", filename).captures[1]) : nothing
end

@everywhere function find_common_files(previous_files, available_channels)
    filtered_files = [previous_files[i] for i in 1:3 if available_channels[i]]

    max_file_num = maximum(reduce(intersect, get_file_number.(filtered_file) for filtered_file in filtered_files))

    common_files = ["", "", ""]

    for i in 1:3
        if available_channels[i]
            common_files[i] = previous_files[i][findfirst(x -> get_file_number(x) == max_file_num, previous_files[i])]
        else
            common_files[i] = ""
        end
    end

    return common_files

    # return [filtered_file[findfirst(x -> get_file_number(x) == max_file_num, filtered_file)] for filtered_file in filtered_files]
end

@everywhere function start_real_time2(folderpath)
    println("Button pressed")

    existing_folders = readdir(folderpath)
    folderpath_correct = false

    println("Waiting")

    while running[] && !folderpath_correct
        new_folders = setdiff(readdir(folderpath), existing_folders)

        if !isempty(new_folders)
            println("Folder found")
            folderpath = joinpath(folderpath, new_folders[end], "Bliq VMS")
            folderpath_correct = true
        end

        sleep(0.1)
    end

    println(folderpath)

    sleep(0.1)

    if !initialized[] && folderpath_correct
        global intensity_C1 = Observable(Float64[])
        global intensity_C2 = Observable(Float64[])
        global intensity_C3 = Observable(Float64[])
        global intensity_ratio = Observable(Float64[])

        global real_time_C1 = Observable(Float64[])
        global real_time_C2 = Observable(Float64[])
        global real_time_C3 = Observable(Float64[])
        global real_time_ratio = Observable(Float64[])

        global ROI_intensity_C1 = Array([])
        global ROI_intensity_C2 = Array([])
        global ROI_intensity_C3 = Array([])
        global ROI_intensity_ratio = Array([])

        global ROI_real_time_C1 = Array([])
        global ROI_real_time_C2 = Array([])
        global ROI_real_time_C3 = Array([])
        global ROI_real_time_ratio = Array([])

        lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
        lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
        lines!(axis_C3_plot, real_time_C3, intensity_C3, color=:black, linewidth=2)
        lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

        if last(split(folderpath, folderpath[1])) != "Bliq VMS"
            println("Working")
            folderpath = joinpath(folderpath, "Bliq VMS")
        end

        global C1_path = joinpath(folderpath, "C1")
        global C2_path = joinpath(folderpath, "C2")
        global C3_path = joinpath(folderpath, "C3")

        global previous_files_C1 = readdir(C1_path)
        global previous_files_C2 = readdir(C2_path)
        global previous_files_C3 = readdir(C3_path)

        global previous_time = now()

        global framerate = 31.25
        global sizex = 1024
        global sizey = 512

        global i = 1
        global j = 1
        global N = 1
    end

    while running[] && folderpath_correct
        if new_ROI[] != 0
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_C3, Observable(Float64[]))
            push!(ROI_real_time_C3, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_real_time_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_real_time_C2)], ROI_intensity_C2[length(ROI_intensity_C2)], color=ROI_colors[mod1(length(ROI_real_time_C2), length(ROI_colors)-1)])
            lines!(axis_C3_plot, ROI_real_time_C3[length(ROI_real_time_C3)], ROI_intensity_C3[length(ROI_intensity_C3)], color=ROI_colors[mod1(length(ROI_real_time_C3), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_real_time_ratio)], ROI_intensity_ratio[length(ROI_intensity_ratio)], color=ROI_colors[mod1(length(ROI_real_time_ratio), length(ROI_colors)-1)])

            new_ROI[] -= 1
        end
        
        current_files_C1 = readdir(C1_path)
        new_files_C1 = setdiff(current_files_C1, previous_files_C1)
        
        current_files_C2 = readdir(C2_path)
        new_files_C2 = setdiff(current_files_C2, previous_files_C2)
        
        current_files_C3 = readdir(C3_path)
        new_files_C3 = setdiff(current_files_C3, previous_files_C3)

        if !isempty(new_files_C1) && !isempty(new_files_C2) && !isempty(new_files_C3)
            previous_files_C1 = current_files_C1
            previous_files_C2 = current_files_C2
            previous_files_C3 = current_files_C3

            # println(new_files_C1[end])

            # time = (now())/Millisecond(1)/1000
            # time = parse(UInt16, new_files_C1[end][end-7:end-4])

            file_C1 = load(joinpath(C1_path, new_files_C1[end]))' * 255
            file_C2 = load(joinpath(C2_path, new_files_C2[end]))' * 255
            file_C3 = load(joinpath(C3_path, new_files_C3[end]))' * 255

            if i == 1
                array_C1[] = file_C1
                array_C2[] = file_C2
                array_C3[] = file_C3
    
                framerate, sizex, sizey = get_file_metadata(joinpath(C1_path, current_files_C1[1]))
            elseif slider_averaging.value[] == 1
                array_C1[] = file_C1
                array_C2[] = file_C2
                array_C3[] = file_C3
            else
                if N < slider_averaging.value[]
                    array_C1[] = 1/(N + 1) .* (N .* array_C1[] .+ file_C1)
                    array_C2[] = 1/(N + 1) .* (N .* array_C2[] .+ file_C2)
                    array_C3[] = 1/(N + 1) .* (N .* array_C3[] .+ file_C3)
    
                    N += 1
                else
                    if N > slider_averaging.value[]
                        for k in 1:N-slider_averaging.value[]
                            array_C1[] -= 1/N * rawview(channelview(load(joinpath(folderpath, "C1", current_files_C1[i-N+k]))))'
                            array_C2[] -= 1/N * rawview(channelview(load(joinpath(folderpath, "C2", current_files_C2[i-N+k]))))'
                            array_C3[] -= 1/N * rawview(channelview(load(joinpath(folderpath, "C3", current_files_C3[i-N+k]))))'
                        end
    
                        array_C1[] *= N/slider_averaging.value[]
                        array_C2[] *= N/slider_averaging.value[]
                        array_C3[] *= N/slider_averaging.value[]
    
                        N = slider_averaging.value[]
                    end
                
                    array_C1[] += 1/N * file_C1 - 1/N * rawview(channelview(load(joinpath(folderpath, "C1", current_files_C1[i-N+1]))))'
                    array_C2[] += 1/N * file_C2 - 1/N * rawview(channelview(load(joinpath(folderpath, "C2", current_files_C2[i-N+1]))))'
                    array_C3[] += 1/N * file_C3 - 1/N * rawview(channelview(load(joinpath(folderpath, "C3", current_files_C3[i-N+1]))))'
                end
            end

            time = (parse(UInt16, new_files_C1[end][end-6:end-4]))/framerate
    
            if noise_size[] != 0
                noise_C1 = sum(file_C1.*(noise_matrix[])')/noise_size[]
                noise_C2 = sum(file_C2.*(noise_matrix[])')/noise_size[]
                noise_C3 = sum(file_C3.*(noise_matrix[])')/noise_size[]
            else
                noise_C1 = 0
                noise_C2 = 0
                noise_C3 = 0
            end
    
            for n in 1:length(ROI_intensity_C1)
                ROI_size = ROI_sizes[][n+1]
                push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/ROI_size - noise_C1)
                push!(ROI_real_time_C1[n][], time)
    
                push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/ROI_size - noise_C2)
                push!(ROI_real_time_C2[n][], time)
    
                push!(ROI_intensity_C3[n][], sum(file_C3.*(ROI_matrix[][n+1,:,:])')/ROI_size - noise_C3)
                push!(ROI_real_time_C3[n][], time)
    
                notify(ROI_intensity_C1[n])
                notify(ROI_real_time_C1[n])
    
                notify(ROI_intensity_C2[n])
                notify(ROI_real_time_C2[n])
    
                notify(ROI_intensity_C3[n])
                notify(ROI_real_time_C3[n])
            end
    
            push!(intensity_C1[], sum(file_C1)/(sizex*sizey) - noise_C1)
            push!(real_time_C1[], time)
    
            push!(intensity_C2[], sum(file_C2)/(sizex*sizey) - noise_C2)
            push!(real_time_C2[], time)
    
            push!(intensity_C3[], sum(file_C3)/(sizex*sizey) - noise_C3)
            push!(real_time_C3[], time)
    
            if menu_ratio.selection[] == "C1/C2"
                if intensity_C2[][end] != 0
                    push!(intensity_ratio[], intensity_C1[][end]/intensity_C2[][end])
                else
                    push!(intensity_ratio[], 0)
                end
    
                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C1[n][]) && !isempty(ROI_intensity_C2[n][])
                        if ROI_intensity_C2[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C1[n][][end]/ROI_intensity_C2[n][][end])
                            push!(ROI_real_time_ratio[n][], time)
    
                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
    
            elseif menu_ratio.selection[] == "C1/C3"
                if intensity_C3[][end] != 0
                    push!(intensity_ratio[], intensity_C1[][end]/intensity_C3[][end])
                else
                    push!(intensity_ratio[], 0)
                end
    
                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C1[n][]) && !isempty(ROI_intensity_C3[n][])
                        if ROI_intensity_C3[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C1[n][][end]/ROI_intensity_C3[n][][end])
                            push!(ROI_real_time_ratio[n][], time)
    
                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
    
            elseif menu_ratio.selection[] == "C2/C1"
                if intensity_C1[][end] != 0
                    push!(intensity_ratio[], intensity_C2[][end]/intensity_C1[][end])
                else
                    push!(intensity_ratio[], 0)
                end
    
                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C2[n][]) && !isempty(ROI_intensity_C1[n][])
                        if ROI_intensity_C1[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C2[n][][end]/ROI_intensity_C1[n][][end])
                            push!(ROI_real_time_ratio[n][], time)
    
                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
    
            elseif menu_ratio.selection[] == "C2/C3"
                if intensity_C3[][end] != 0
                    push!(intensity_ratio[], intensity_C2[][end]/intensity_C3[][end])
                else
                    push!(intensity_ratio[], 0)
                end
    
                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C2[n][]) && !isempty(ROI_intensity_C3[n][])
                        if ROI_intensity_C3[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C2[n][][end]/ROI_intensity_C3[n][][end])
                            push!(ROI_real_time_ratio[n][], time)
    
                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
    
            elseif menu_ratio.selection[] == "C3/C1"
                if intensity_C1[][end] != 0
                    push!(intensity_ratio[], intensity_C3[][end]/intensity_C1[][end])
                else
                    push!(intensity_ratio[], 0)
                end
    
                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C3[n][]) && !isempty(ROI_intensity_C1[n][])
                        if ROI_intensity_C1[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C3[n][][end]/ROI_intensity_C1[n][][end])
                            push!(ROI_real_time_ratio[n][], time)
    
                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
    
            elseif menu_ratio.selection[] == "C3/C2"
                if intensity_C2[][end] != 0
                    push!(intensity_ratio[], intensity_C3[][end]/intensity_C2[][end])
                else
                    push!(intensity_ratio[], 0)
                end
    
                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C3[n][]) && !isempty(ROI_intensity_C2[n][])
                        if ROI_intensity_C2[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C3[n][][end]/ROI_intensity_C2[n][][end])
                            push!(ROI_real_time_ratio[n][], time)
    
                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
            end
    
            push!(real_time_ratio[], time)
    
            notify(array_C1)
            notify(intensity_C1)
            notify(real_time_C1)
    
            notify(array_C2)
            notify(intensity_C2)
            notify(real_time_C2)
    
            notify(array_C3)
            notify(intensity_C3)
            notify(real_time_C3)
    
            notify(intensity_ratio)
            notify(real_time_ratio)
    
            autolimits!(axis_C1_plot)
            xlims!(axis_C1_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))
    
            autolimits!(axis_C2_plot)
            xlims!(axis_C2_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))
    
            autolimits!(axis_C3_plot)
            xlims!(axis_C3_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))
    
            autolimits!(axis_ratio_plot)
            xlims!(axis_ratio_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))
    
            i += 1
    
            if mod1(i, 10) == 10
                current_time = now()
                label_frequency.text = "Frequency = $(string(round(1/((current_time - previous_time)/Millisecond(1)/10000), sigdigits=3))) Hz"
                previous_time = current_time
            end
        end
    
        sleep(0.01)
    end
end

@everywhere function start_real_time(folderpath)
    folderpath, folderpath_correct = verify_folder(folderpath)

    if !initialized[] && folderpath_correct
        println("Folder correct")

        arrays = [zeros(Float64, 1024, 512), zeros(Float64, 1024, 512), zeros(Float64, 1024, 512)]
        intensity = [Observable(Float64[]), Observable(Float64[]), Observable(Float64[]), Observable(Float64[])]
        real_time = Observable(Float64[])
        # ROI_intensity = []
        ROI_intensity = [[], [], [], []]
        ROI_real_time = []

        for n in 1:4
            lines!(plotgrid.content[n].content, real_time, intensity[n], color=:black, linewidth=2)
        end 

        channels_paths = [joinpath(folderpath, "C1"), joinpath(folderpath, "C2"), joinpath(folderpath, "C3")]
        available_channels = [true, true, true]
        # previous_files = readdir.(channels_paths)
        previous_files = [[], [], []]
        current_files = [[], [], []]
        new_files = [[], [], []]
        common_files = []

        for n in 1:3
            try
                previous_files[n] = readdir(channels_paths[n])
            catch e
                available_channels[n] = false
            end
        end

        i = 0
        j = 0
        N = [1, 1, 1]
        noise = [0.0, 0.0, 0.0]
        ROIs = 0

        previous_time = now()

        framerate = 31.25
        sizex = 1024
        sizey = 512

        time = 0

        initialized[] = true
    end

    while running[]
        while new_ROI[] != 0
            ROIs += 1

            # push!(ROI_intensity, [Observable(Float64[]), Observable(Float64[]), Observable(Float64[]), Observable(Float64[])])
            push!(ROI_real_time, Observable(Float64[]))

            for n in 1:4
                push!(ROI_intensity[n], Observable(Float64[]))
                lines!(plotgrid.content[n].content, ROI_real_time[ROIs], ROI_intensity[n][ROIs], color=ROI_colors[mod1(length(ROI_intensity[1]), length(ROI_colors)-1)])
            end

            # println(length(ROI_intensity[1]))
            new_ROI[] -= 1
        end

        i += 1

        current_files = [available_channels[n] ? readdir(channels_paths[n]) : [] for n in 1:3]

        new_files = setdiff.(current_files, previous_files)

        if all(isempty.(new_files) .== .!available_channels)
            println("Files detected")
            common_files = find_common_files(new_files, available_channels)
        
            println(common_files)

            files = [zeros(UInt8, sizex, sizey), zeros(UInt8, sizex, sizey), zeros(UInt8, sizex, sizey)]

            for channel in 1:3
                if available_channels[channel]
                    file = joinpath(channels_paths[channel], common_files[channel])

                    index = findfirst(==(common_files[channel]), current_files[channel])

                    println("Index : ",index)
                    
                    files[channel] = round.(UInt8, load(file)'*255)

                    if i == 1
                        arrays[channel] = files[channel]

                        if j == 0
                            framerate, sizex, sizey = get_file_metadata(file)
                        end
                    elseif slider_averaging.value[] == 1
                        arrays[channel] = files[channel]
                    else
                        if N[channel] < slider_averaging.value[]
                            arrays[channel] = 1/(N[channel] + 1) .* (N[channel] .* arrays[channel] .+ files[channel])
            
                            N[channel] += 1
                        else
                            array = arrays[channel]

                            if N[channel] > slider_averaging.value[]
                                for k in 1:N[channel]-slider_averaging.value[]
                                    array -= 255/N[channel] * load(joinpath(channels_paths[channel], current_files[channel][index-N[channel]+k]))'
                                end
            
                                array *= N[channel]/slider_averaging.value[]
            
                                N[channel] = slider_averaging.value[]
                            end
                        
                            array += 1/N[channel] * files[channel] - 255/N[channel] * load(joinpath(channels_paths[channel], current_files[channel][index-N[channel]+1]))'

                            arrays[channel] = array
                        end
                    end

                    # time = (get_file_number(file) + j)/framerate
                    time = i
                    
                    if noise_size[] != 0
                        noise[channel] = sum(files[channel].*(noise_matrix[])')/noise_size[]
                    else
                        noise[channel] = 0
                    end

                    push!(intensity[channel][], max(0, sum(files[channel])/(sizex*sizey) - noise[channel]))

                    for ROI in eachindex(ROI_intensity[channel])
                        push!(ROI_intensity[channel][ROI][], max(0, sum(files[channel].*(ROI_matrix[][ROI+1,:,:])')/ROI_sizes[][ROI+1] - noise[channel]))
                    end
                else
                    push!(intensity[channel][], 0)

                    for ROI in eachindex(ROI_intensity[channel])
                        push!(ROI_intensity[channel][ROI][], 0)
                    end
                end
            end

            num = parse(UInt8, menu_ratio.selection[][2])
            denum = parse(UInt8, menu_ratio.selection[][5])

            if intensity[denum][][end] > 1
                push!(intensity[4][], intensity[num][][end]/intensity[denum][][end])
            else
                push!(intensity[4][], 0)
            end

            for ROI in eachindex(ROI_intensity[4])
                if !isempty(ROI_intensity[num][ROI][]) && !isempty(ROI_intensity[denum][ROI][])
                    if ROI_intensity[denum][ROI][][end] > 1
                        push!(ROI_intensity[4][ROI][], ROI_intensity[num][ROI][][end]/ROI_intensity[denum][ROI][][end])
                        push!(ROI_real_time[ROI][], time)
                    else
                        push!(ROI_intensity[4][ROI][], 0)
                        push!(ROI_real_time[ROI][], time)
                    end
                end
            end

            push!(real_time[], time)

            if mod1(i, slider_refresh.value[]) == slider_refresh.value[]
                if toggle_display.active[]
                    for n in eachindex(images)
                        images[n][] = arrays[n]
                    end
                end
    
                for ROI in eachindex(ROI_intensity[4])
                    notify(ROI_real_time[ROI])
                end
    
                notify(real_time)
    
                for n in 1:4
                    autolimits!(plotgrid.content[n].content)
                    xlims!(plotgrid.content[n].content, max(0, time - slider_time.value[]+1), max(slider_time.value[], time+1))
                end
            end
    
            if mod1(i, 10+slider_refresh.value[]) == 10+slider_refresh.value[]
                current_time = now()
                label_frequency.text = "Frequency : $(string(round(1/((current_time - previous_time)/Millisecond(1))*1000*(10+slider_refresh.value[]), sigdigits=3))) Hz"
                previous_time = current_time
    
                label_number.text = "Files : $(round(UInt64, i+j))"
            end

            previous_files = current_files
    
            i += 1
            sleep(0.001)
        end

        sleep(1)
    end
end

@everywhere function start_playback2(folderpath)
    println("Button pressed")

    if isdir(joinpath(folderpath, "C1"))
        folderpath_correct = true
    elseif isdir(joinpath(folderpath, "Bliq VMS", "C1"))
        folderpath = joinpath(folderpath, "Bliq VMS")

        folderpath_correct = true
    else
        println("Folder error")
        folderpath_correct = false
        running[] = false

        button_start.label = "START"
        button_reset.label = "CLEAR"
    end

    if !initialized[] && folderpath_correct
        global intensity_C1 = Observable(Float64[])
        global intensity_C2 = Observable(Float64[])
        global intensity_C3 = Observable(Float64[])
        global intensity_ratio = Observable(Float64[])

        global real_time_C1 = Observable(Float64[])
        global real_time_C2 = Observable(Float64[])
        global real_time_C3 = Observable(Float64[])
        global real_time_ratio = Observable(Float64[])

        global ROI_intensity_C1 = Array([])
        global ROI_intensity_C2 = Array([])
        global ROI_intensity_C3 = Array([])
        global ROI_intensity_ratio = Array([])

        global ROI_real_time_C1 = Array([])
        global ROI_real_time_C2 = Array([])
        global ROI_real_time_C3 = Array([])
        global ROI_real_time_ratio = Array([])

        lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
        lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
        lines!(axis_C3_plot, real_time_C3, intensity_C3, color=:black, linewidth=2)
        lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

        global C1_path = joinpath(folderpath, "C1")
        global C2_path = joinpath(folderpath, "C2")
        global C3_path = joinpath(folderpath, "C3")

        global C1_files = readdir(C1_path)
        global C2_files = readdir(C2_path)
        global C3_files = readdir(C3_path)

        if C1_files[1] == ".DS_Store"
            deleteat!(C1_files, 1)
        end

        if C2_files[1] == ".DS_Store"
            deleteat!(C2_files, 1)
        end

        if C3_files[1] == ".DS_Store"
            deleteat!(C3_files, 1)
        end

        global i = 1
        global j = 0
        global N = 1

        global previous_time = now()

        global framerate = 31
        global sizex = 1024
        global sizey = 512

        initialized[] = true
    end

    while running[] && folderpath_correct
        if new_ROI[] != 0
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_C3, Observable(Float64[]))
            push!(ROI_real_time_C3, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_real_time_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_real_time_C2)], ROI_intensity_C2[length(ROI_intensity_C2)], color=ROI_colors[mod1(length(ROI_real_time_C2), length(ROI_colors)-1)])
            lines!(axis_C3_plot, ROI_real_time_C3[length(ROI_real_time_C3)], ROI_intensity_C3[length(ROI_intensity_C3)], color=ROI_colors[mod1(length(ROI_real_time_C3), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_real_time_ratio)], ROI_intensity_ratio[length(ROI_intensity_ratio)], color=ROI_colors[mod1(length(ROI_real_time_ratio), length(ROI_colors)-1)])

            new_ROI[] -= 1
        end

        try
            C1_files[i]
        catch e
            println("Restarting")
            j += i - 1
            i = 1
            N = 1
        end

        file_C1 = load(joinpath(C1_path, C1_files[i]))' * 255
        file_C2 = load(joinpath(C2_path, C2_files[i]))' * 255
        file_C3 = load(joinpath(C3_path, C3_files[i]))' * 255

        if i == 1 && j == 0
            array_C1[] = file_C1
            array_C2[] = file_C2
            array_C3[] = file_C3

            framerate, sizex, sizey = get_file_metadata(joinpath(C1_path, C1_files[1]))
        elseif slider_averaging.value[] == 1
            array_C1[] = file_C1
            array_C2[] = file_C2
            array_C3[] = file_C3
        else
            if N < slider_averaging.value[]
                array_C1[] = 1/(N + 1) .* (N .* array_C1[] .+ file_C1)
                array_C2[] = 1/(N + 1) .* (N .* array_C2[] .+ file_C2)
                array_C3[] = 1/(N + 1) .* (N .* array_C3[] .+ file_C3)

                N += 1
            else
                if N > slider_averaging.value[]
                    for k in 1:N-slider_averaging.value[]
                        array_C1[] -= 255/N * load(joinpath(C1_path, C1_files[i-N+k]))'
                        array_C2[] -= 255/N * load(joinpath(C2_path, C2_files[i-N+k]))'
                        array_C3[] -= 255/N * load(joinpath(C3_path, C3_files[i-N+k]))'
                    end

                    array_C1[] *= N/slider_averaging.value[]
                    array_C2[] *= N/slider_averaging.value[]
                    array_C3[] *= N/slider_averaging.value[]

                    N = slider_averaging.value[]
                end
            
                array_C1[] += 1/N * file_C1 - 255/N * load(joinpath(C1_path, C1_files[i-N+1]))'
                array_C2[] += 1/N * file_C2 - 255/N * load(joinpath(C2_path, C2_files[i-N+1]))'
                array_C3[] += 1/N * file_C3 - 255/N * load(joinpath(C3_path, C3_files[i-N+1]))'
            end
        end

        # println(C1_files[i])
        # println(get_file_number(C1_files[i]))

        time = (get_file_number(C1_files[i]) + j)/framerate

        if noise_size[] != 0
            noise_C1 = sum(file_C1.*(noise_matrix[])')/noise_size[]
            noise_C2 = sum(file_C2.*(noise_matrix[])')/noise_size[]
            noise_C3 = sum(file_C3.*(noise_matrix[])')/noise_size[]

            println(noise_C1)
        else
            noise_C1 = 0
            noise_C2 = 0
            noise_C3 = 0
        end

        for n in 1:length(ROI_intensity_C1)
            ROI_size = ROI_sizes[][n+1]
            push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/ROI_size - noise_C1)
            push!(ROI_real_time_C1[n][], time)

            push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/ROI_size - noise_C2)
            push!(ROI_real_time_C2[n][], time)

            push!(ROI_intensity_C3[n][], sum(file_C3.*(ROI_matrix[][n+1,:,:])')/ROI_size - noise_C3)
            push!(ROI_real_time_C3[n][], time)

            notify(ROI_intensity_C1[n])
            notify(ROI_real_time_C1[n])

            notify(ROI_intensity_C2[n])
            notify(ROI_real_time_C2[n])

            notify(ROI_intensity_C3[n])
            notify(ROI_real_time_C3[n])
        end

        push!(intensity_C1[], sum(file_C1)/(sizex*sizey) - noise_C1)
        push!(real_time_C1[], time)

        push!(intensity_C2[], sum(file_C2)/(sizex*sizey) - noise_C2)
        push!(real_time_C2[], time)

        push!(intensity_C3[], sum(file_C3)/(sizex*sizey) - noise_C3)
        push!(real_time_C3[], time)

        if menu_ratio.selection[] == "C1/C2"
            if intensity_C2[][end] != 0
                push!(intensity_ratio[], intensity_C1[][end]/intensity_C2[][end])
            else
                push!(intensity_ratio[], 0)
            end

            for n in 1:length(ROI_intensity_ratio)
                if !isempty(ROI_intensity_C1[n][]) && !isempty(ROI_intensity_C2[n][])
                    if ROI_intensity_C2[n][][end] != 0
                        push!(ROI_intensity_ratio[n][], ROI_intensity_C1[n][][end]/ROI_intensity_C2[n][][end])
                        push!(ROI_real_time_ratio[n][], time)

                        notify(ROI_intensity_ratio[n])
                        notify(ROI_real_time_ratio[n])
                    end
                end
            end

        elseif menu_ratio.selection[] == "C1/C3"
            if intensity_C3[][end] != 0
                push!(intensity_ratio[], intensity_C1[][end]/intensity_C3[][end])
            else
                push!(intensity_ratio[], 0)
            end

            for n in 1:length(ROI_intensity_ratio)
                if !isempty(ROI_intensity_C1[n][]) && !isempty(ROI_intensity_C3[n][])
                    if ROI_intensity_C3[n][][end] != 0
                        push!(ROI_intensity_ratio[n][], ROI_intensity_C1[n][][end]/ROI_intensity_C3[n][][end])
                        push!(ROI_real_time_ratio[n][], time)

                        notify(ROI_intensity_ratio[n])
                        notify(ROI_real_time_ratio[n])
                    end
                end
            end

        elseif menu_ratio.selection[] == "C2/C1"
            if intensity_C1[][end] != 0
                push!(intensity_ratio[], intensity_C2[][end]/intensity_C1[][end])
            else
                push!(intensity_ratio[], 0)
            end

            for n in 1:length(ROI_intensity_ratio)
                if !isempty(ROI_intensity_C2[n][]) && !isempty(ROI_intensity_C1[n][])
                    if ROI_intensity_C1[n][][end] != 0
                        push!(ROI_intensity_ratio[n][], ROI_intensity_C2[n][][end]/ROI_intensity_C1[n][][end])
                        push!(ROI_real_time_ratio[n][], time)

                        notify(ROI_intensity_ratio[n])
                        notify(ROI_real_time_ratio[n])
                    end
                end
            end

        elseif menu_ratio.selection[] == "C2/C3"
            if intensity_C3[][end] != 0
                push!(intensity_ratio[], intensity_C2[][end]/intensity_C3[][end])
            else
                push!(intensity_ratio[], 0)
            end

            for n in 1:length(ROI_intensity_ratio)
                if !isempty(ROI_intensity_C2[n][]) && !isempty(ROI_intensity_C3[n][])
                    if ROI_intensity_C3[n][][end] != 0
                        push!(ROI_intensity_ratio[n][], ROI_intensity_C2[n][][end]/ROI_intensity_C3[n][][end])
                        push!(ROI_real_time_ratio[n][], time)

                        notify(ROI_intensity_ratio[n])
                        notify(ROI_real_time_ratio[n])
                    end
                end
            end

        elseif menu_ratio.selection[] == "C3/C1"
            if intensity_C1[][end] != 0
                push!(intensity_ratio[], intensity_C3[][end]/intensity_C1[][end])
            else
                push!(intensity_ratio[], 0)
            end

            for n in 1:length(ROI_intensity_ratio)
                if !isempty(ROI_intensity_C3[n][]) && !isempty(ROI_intensity_C1[n][])
                    if ROI_intensity_C1[n][][end] != 0
                        push!(ROI_intensity_ratio[n][], ROI_intensity_C3[n][][end]/ROI_intensity_C1[n][][end])
                        push!(ROI_real_time_ratio[n][], time)

                        notify(ROI_intensity_ratio[n])
                        notify(ROI_real_time_ratio[n])
                    end
                end
            end

        elseif menu_ratio.selection[] == "C3/C2"
            if intensity_C2[][end] != 0
                push!(intensity_ratio[], intensity_C3[][end]/intensity_C2[][end])
            else
                push!(intensity_ratio[], 0)
            end

            for n in 1:length(ROI_intensity_ratio)
                if !isempty(ROI_intensity_C3[n][]) && !isempty(ROI_intensity_C2[n][])
                    if ROI_intensity_C2[n][][end] != 0
                        push!(ROI_intensity_ratio[n][], ROI_intensity_C3[n][][end]/ROI_intensity_C2[n][][end])
                        push!(ROI_real_time_ratio[n][], time)

                        notify(ROI_intensity_ratio[n])
                        notify(ROI_real_time_ratio[n])
                    end
                end
            end
        end

        push!(real_time_ratio[], time)

        notify(array_C1)
        notify(intensity_C1)
        notify(real_time_C1)

        notify(array_C2)
        notify(intensity_C2)
        notify(real_time_C2)

        notify(array_C3)
        notify(intensity_C3)
        notify(real_time_C3)

        notify(intensity_ratio)
        notify(real_time_ratio)

        autolimits!(axis_C1_plot)
        xlims!(axis_C1_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))

        autolimits!(axis_C2_plot)
        xlims!(axis_C2_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))

        autolimits!(axis_C3_plot)
        xlims!(axis_C3_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))

        autolimits!(axis_ratio_plot)
        xlims!(axis_ratio_plot, max(0, time - slider_time.value[]), max(slider_time.value[], time))

        if mod1(i, 10) == 10
            current_time = now()
            label_frequency.text = "Frequency : $(string(round(1/((current_time - previous_time)/Millisecond(1)/10000), sigdigits=3))) Hz"
            previous_time = current_time

            label_number.text = "Files : $(length(intensity_C1[]))/$(i+j)"
        end

        i += 1
        
        sleep(0.001)
    end
end

@everywhere function start_playback(folderpath)
    println("Button pressed")

    folderpath, folderpath_correct = verify_folder(folderpath)

    if !initialized[] && folderpath_correct
        println("Folder correct")

        arrays = [zeros(Float64, 1024, 512), zeros(Float64, 1024, 512), zeros(Float64, 1024, 512)]
        intensity = [Observable(Float64[]), Observable(Float64[]), Observable(Float64[]), Observable(Float64[])]
        real_time = Observable(Float64[])
        # ROI_intensity = []
        ROI_intensity = [[], [], [], []]
        ROI_real_time = []

        for n in 1:4
            lines!(plotgrid.content[n].content, real_time, intensity[n], color=:black, linewidth=2)
        end 

        channels_paths = [joinpath(folderpath, "C1"), joinpath(folderpath, "C2"), joinpath(folderpath, "C3")]
        available_channels = [true, true, true]
        channels_files = []

        for n in 1:3
            try
                channels_file = readdir(channels_paths[n])

                if channels_file[1] == ".DS_Store"
                    deleteat!(channels_file, 1)
                end
    
                sort!(channels_file, by=get_file_number)

                push!(channels_files, channels_file)
            catch e
                push!(channels_files, [])

                available_channels[n] = false
            end
        end

        i = 1
        j = 0
        N = [1, 1, 1]
        noise = [0.0, 0.0, 0.0]
        ROIs = 0

        previous_time = now()

        framerate = 31.25
        sizex = 1024
        sizey = 512

        time = 0

        initialized[] = true
    end

    println("Starting")

    println(available_channels)

    while running[] && folderpath_correct
        while new_ROI[] != 0
            ROIs += 1

            # push!(ROI_intensity, [Observable(Float64[]), Observable(Float64[]), Observable(Float64[]), Observable(Float64[])])
            push!(ROI_real_time, Observable(Float64[]))

            for n in 1:4
                push!(ROI_intensity[n], Observable(Float64[]))
                lines!(plotgrid.content[n].content, ROI_real_time[ROIs], ROI_intensity[n][ROIs], color=ROI_colors[mod1(length(ROI_intensity[1]), length(ROI_colors)-1)])
            end

            # println(length(ROI_intensity[1]))
            new_ROI[] -= 1
        end

        try
            channels_files[1][i]
        catch _
            try
                channels_files[2][i]
            catch _
                try
                    channels_files[3][i]
                catch _
                    println("Restarting")
                    println(i)
                    j = time*framerate
                    i = 1
                    N = [1, 1, 1]
                end
            end
        end

        files = [zeros(UInt8, sizex, sizey), zeros(UInt8, sizex, sizey), zeros(UInt8, sizex, sizey)]

        for channel in 1:3
            if available_channels[channel]
                if i == 1 && j == 0
                    framerate, sizex, sizey = get_file_metadata(joinpath(channels_paths[channel], channels_files[channel][1]))
                end

                files[channel] = round.(UInt8, load(joinpath(channels_paths[channel], channels_files[channel][i]))'*255)

                if toggle_display.active[]
                    if i == 1 
                        arrays[channel] = files[channel]
                    elseif slider_averaging.value[] == 1
                        arrays[channel] = files[channel]
                    else
                        if N[channel] < slider_averaging.value[]
                            arrays[channel] = 1/(N[channel] + 1) .* (N[channel] .* arrays[channel] .+ files[channel])
            
                            N[channel] += 1
                        else
                            array = arrays[channel]

                            if N[channel] > slider_averaging.value[]
                                for k in 1:N[channel]-slider_averaging.value[]
                                    array -= 255/N[channel] * load(joinpath(channels_paths[channel], channels_files[channel][i-N[channel]+k]))'
                                end
            
                                array *= N[channel]/slider_averaging.value[]
            
                                N[channel] = slider_averaging.value[]
                            end
                        
                            array += 1/N[channel] * files[channel] - 255/N[channel] * load(joinpath(channels_paths[channel], channels_files[channel][i-N[channel]+1]))'

                            arrays[channel] = array
                        end
                    end
                end

                time = (get_file_number(channels_files[channel][i]) + j)/framerate
                
                if noise_size[] != 0
                    noise[channel] = sum(files[channel].*(noise_matrix[])')/noise_size[]
                else
                    noise[channel] = 0
                end

                push!(intensity[channel][], max(0, sum(files[channel])/(sizex*sizey) - noise[channel]))

                for ROI in eachindex(ROI_intensity[channel])
                    push!(ROI_intensity[channel][ROI][], max(0, sum(files[channel].*(ROI_matrix[][ROI+1,:,:])')/ROI_sizes[][ROI+1] - noise[channel]))
                end
            else
                push!(intensity[channel][], 0)

                for ROI in eachindex(ROI_intensity[channel])
                    push!(ROI_intensity[channel][ROI][], 0)
                end
            end
        end

        num = parse(UInt8, menu_ratio.selection[][2])
        denum = parse(UInt8, menu_ratio.selection[][5])

        if intensity[denum][][end] > 1
            push!(intensity[4][], intensity[num][][end]/intensity[denum][][end])
        else
            push!(intensity[4][], 0)
        end

        for ROI in eachindex(ROI_intensity[4])
            if !isempty(ROI_intensity[num][ROI][]) && !isempty(ROI_intensity[denum][ROI][])
                if ROI_intensity[denum][ROI][][end] > 1
                    push!(ROI_intensity[4][ROI][], ROI_intensity[num][ROI][][end]/ROI_intensity[denum][ROI][][end])
                    push!(ROI_real_time[ROI][], time)
                else
                    push!(ROI_intensity[4][ROI][], 0)
                    push!(ROI_real_time[ROI][], time)
                end
            end
        end
        
        push!(real_time[], time)

        if mod1(i, slider_refresh.value[]) == slider_refresh.value[]
            if toggle_display.active[]
                for n in eachindex(images)
                    images[n][] = arrays[n]
                end
            end

            for ROI in eachindex(ROI_intensity[4])
                notify(ROI_real_time[ROI])
            end

            notify(real_time)

            for n in 1:4
                autolimits!(plotgrid.content[n].content)
                xlims!(plotgrid.content[n].content, max(0, time - slider_time.value[]+1), max(slider_time.value[], time+1))
            end
        end

        if mod1(i, 9+slider_refresh.value[]) == 9+slider_refresh.value[]
            current_time = now()
            label_frequency.text = "Frequency : $(string(round(1/((current_time - previous_time)/Millisecond(1))*1000*(9+slider_refresh.value[]), sigdigits=3))) Hz"
            previous_time = current_time

            label_number.text = "Files : $(round(UInt64, i))"
        end

        i += 1
        sleep(0.001)
    end
end 

image_colors = [RGB(1,0.25,0.25), RGB(0.25,1,0.25), RGB(0.25,0.5,1), RGB(0.25,1,1), RGB(1,0.25,1), RGB(1,1,0.25), RGB(1, 1, 1)]
ROI_colors = [RGB(1,0,0), RGB(1,0.5,0), RGB(1,1,0), RGB(0,1,0), RGB(0,1,1), RGB(0,0,1), RGB(0.5,0,1), RGB(1,0,1), RGB(1,1,1)]
ROI_colors_dimmed = [RGBA(1,0,0, 0.1), RGBA(1,0.5,0, 0.1), RGBA(1,1,0, 0.1), RGBA(0,1,0, 0.1), RGBA(0,1,1, 0.1), RGBA(0,0,1, 0.1), RGBA(0.5,0,1, 0.1), RGBA(1,0,1, 0.1), RGBA(1,1,1, 0.1)]

set_theme!(backgroundcolor = :gray95)

display_averaging_limits = 100

folderpath = get_folder_path()

images = [Observable(zeros(Float64, 1024, 512)), Observable(zeros(Float64, 1024, 512)), Observable(zeros(Float64, 1024, 512))]

# fig = Figure(size = (1440, 847))
fig = Figure(size = (1920, 1026))

fig[1, 1] = buttongrid = GLMakie.GridLayout(tellheight=false)
fig[1, 2] = figgrid = GLMakie.GridLayout()
fig[1, 3] = plotgrid = GLMakie.GridLayout()

colsize!(fig.layout, 1, GLMakie.Fixed(300))

# progress_bar = buttongrid[0, 1:2] = Slider(fig, range = 1:0.01:100, snap=false, linewidth = 14, color_active=RGBf(0.8, 0.8, 0.8), color_active_dimmed=RGBf(0.8, 0.8, 0.8), color_inactive=:gray95)
progress_bar_1 = buttongrid[0, 1:2] = Box(fig, cornerradius = 7, color=RGBf(0.8, 0.8, 0.8), strokecolor=:transparent, height = 14, width=260, halign=:left)
progress_bar_2 = buttongrid[0, 1:2] = Box(fig, cornerradius = 7, color=RGBf(0.6, 0.6, 0.6), strokecolor=:transparent, height = 14, width=130, halign=:left)
# println(fieldnames(typeof(progress_bar.layoutobservables)))
# println(progress_bar.layoutobservables)
# progress_bar.blockscene.visible[] = false

label_frequency = buttongrid[1, 1] = Label(fig, "Frequency : 0 Hz", halign=:left, fontsize=12)
label_number = buttongrid[1, 2] = Label(fig, "Files : 0/0", halign=:left, fontsize=12)

button_start = buttongrid[2, 1] = Button(fig, label="START", font=:bold, width = 120, height = 32, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)
button_reset = buttongrid[2, 2] = Button(fig, label="CLEAR", font=:bold, width = 120, height = 32, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)

button_load = buttongrid[3, 1] = Button(fig, label="LOAD", font=:bold, width = 120, height = 32, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)
button_save = buttongrid[3, 2] = Button(fig, label="SAVE", font=:bold, width = 120, height = 32, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)

button_ROI = buttongrid[4, 1] = Button(fig, label="Set ROI", width = 120, height = 32, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)
button_noise = buttongrid[4, 2] = Button(fig, label="Set noise", width = 120, height = 32, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)

rowsize!(buttongrid, 4, 48)

label_mode = buttongrid[5, 1] = Label(fig, "Mode :", halign=:right, fontsize=16)
label_color_C1 = buttongrid[6, 1] = Label(fig, "Color C1 :", halign=:right, fontsize=16)
label_color_C2 = buttongrid[7, 1] = Label(fig, "Color C2 :", halign=:right, fontsize=16)
label_color_C3 = buttongrid[8, 1] = Label(fig, "Color C3 :", halign=:right, fontsize=16)
label_ratio = buttongrid[9, 1] = Label(fig, "Ratio : ", halign=:right, fontsize=16)

menu_mode = buttongrid[5, 2] = Menu(fig, options=["Playback", "Real-time"], default = "Playback", width = 120, height = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=16, textcolor=:black)
# menu_mode = buttongrid[5, 2] = Menu(fig, options=["Playback", "Real-time"], default = "Real-time", width = 120, height = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=16, textcolor=:black)
menu_color_C1 = buttongrid[6, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "White"], default = "Red", width = 120, height = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=16, textcolor=:black)
menu_color_C2 = buttongrid[7, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "White"], default = "Green", width = 120, height = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=16, textcolor=:black)
menu_color_C3 = buttongrid[8, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "White"], default = "Blue", width = 120, height = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=16, textcolor=:black)
menu_ratio = buttongrid[9, 2] = Menu(fig, options=["C1/C2", "C1/C3", "C2/C1", "C2/C3", "C3/C1", "C3/C2"], default = "C1/C2", width = 120, height = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=16, textcolor=:black)

rowsize!(buttongrid, 5, 48)
rowsize!(buttongrid, 9, 48)

toggle_display = buttongrid[11, 1:2][1, 1] = Toggle(fig, height = 14, active=true, rimfraction=0, buttoncolor=RGBf(0.4, 0.4, 0.4), framecolor_active=RGBf(0.6, 0.6, 0.6), framecolor_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false)

label_refresh = buttongrid[10, 1:2] = Label(fig, "\nDisplay speed : 1", fontsize=16, height=10)
slider_refresh = buttongrid[11, 1:2][1, 2] = Slider(fig, range = 1:1:16, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_averaging = buttongrid[12, 1:2] = Label(fig, "\nDisplay averaging : 1", fontsize=16, height=10)
slider_averaging = buttongrid[13, 1:2] = Slider(fig, range = 1:1:display_averaging_limits, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_limit_C1 = buttongrid[14, 1:2] = Label(fig, "\nDisplay limit C1 : 50", fontsize=16, height=10)
slider_limit_C1 = buttongrid[15, 1:2] = Slider(fig, range = 3:1:99, startvalue = 50, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_limit_C2 = buttongrid[16, 1:2] = Label(fig, "\nDisplay limit C2 : 50", fontsize=16, height=10)
slider_limit_C2 = buttongrid[17, 1:2] = Slider(fig, range = 3:1:99, startvalue = 50, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_limit_C3 = buttongrid[18, 1:2] = Label(fig, "\nDisplay limit C3 : 50", fontsize=16, height=10)
slider_limit_C3 = buttongrid[19, 1:2] = Slider(fig, range = 6:1:99, startvalue = 50, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

axis_C1_fig = figgrid[1, 1][1, 1] = GLMakie.Axis(fig, xrectzoom=false, yrectzoom=false, title="Channel 1", titlesize=16, backgroundcolor=:black, aspect=2)
axis_C2_fig = figgrid[2, 1][1, 1] = GLMakie.Axis(fig, xrectzoom=false, yrectzoom=false, title="Channel 2", titlesize=16, backgroundcolor=:black, aspect=2)
axis_C3_fig = figgrid[3, 1][1, 1] = GLMakie.Axis(fig, xrectzoom=false, yrectzoom=false, title="Channel 3", titlesize=16, backgroundcolor=:black, aspect=2)

axis_C1_plot = plotgrid[1, 1] = GLMakie.Axis(fig, title="Channel 1", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))
axis_C2_plot = plotgrid[2, 1] = GLMakie.Axis(fig, title="Channel 2", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))
axis_C3_plot = plotgrid[3, 1] = GLMakie.Axis(fig, title="Channel 3", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))
axis_ratio_plot = plotgrid[4, 1] = GLMakie.Axis(fig, title="Ratio of C1/C2", xlabel = "Time [s]", ylabel = "Ratio", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))

slider_time = plotgrid[5, 1][1, 1] = Slider(fig, range = 1:1:30, startvalue = 10, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))
textbox_time = plotgrid[5, 1][1, 2] = Textbox(fig, placeholder = " ", height = 30, width = 40, textpadding=(8,6,6,7), halign=:center, valign=:center, bordercolor=:transparent, bordercolor_focused=:transparent, bordercolor_focused_invalid=:transparent, bordercolor_hover=:transparent, textcolor_placeholder=RGBf(0.4, 0.4, 0.4), textcolor=RGBf(0.4, 0.4, 0.4), validator=UInt16)
textbox_time.displayed_string[] = "30"
textbox_time.stored_string[] = "30"

limits!(axis_C1_fig, 0, 1024, 512, 0)
limits!(axis_C2_fig, 0, 1024, 512, 0)
limits!(axis_C3_fig, 0, 1024, 512, 0)

hidedecorations!(axis_C1_fig)
hidedecorations!(axis_C2_fig)
hidedecorations!(axis_C3_fig)

colorbar_C1 = figgrid[1, 1][1, 2] = Colorbar(fig, halign=:left, limits=(0, 50), colormap=[:black, image_colors[menu_color_C1.i_selected[]]])
colorbar_C2 = figgrid[2, 1][1, 2] = Colorbar(fig, halign=:left, limits=(0, 50), colormap=[:black, image_colors[menu_color_C2.i_selected[]]])
colorbar_C3 = figgrid[3, 1][1, 2] = Colorbar(fig, halign=:left, limits=(0, 50), colormap=[:black, image_colors[menu_color_C3.i_selected[]]])

img_C1 = image!(axis_C1_fig, images[1], colorrange=[0,50], colormap=[:black, image_colors[menu_color_C1.i_selected[]]])
img_C2 = image!(axis_C2_fig, images[2], colorrange=[0,50], colormap=[:black, image_colors[menu_color_C2.i_selected[]]])
img_C3 = image!(axis_C3_fig, images[3], colorrange=[0,50], colormap=[:black, image_colors[menu_color_C3.i_selected[]]])

dragging = false
idx = 1

running = Observable(false)
ran = Observable(false)
saving = Observable(false)
saved = false
initialized = Observable(false)
current_ROI = Observable(1)
current_noise = Observable(1)
region_points = Observable(Point2f[])

new_ROI = Observable(0)

ROI_matrix = Observable(ones(1, 512,1024))
ROI_sizes = Observable([512*1024])
noise_matrix = Observable(zeros(512,1024))
noise_size = Observable(0)

region_plot_C1 = poly!(axis_C1_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
region_plot_C2 = poly!(axis_C2_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
region_plot_C3 = poly!(axis_C3_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)

region_scatter = scatter!(axis_C1_fig, region_points, markersize = 12, overdraw=true, color=strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)])

hideydecorations!(axis_C1_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
hideydecorations!(axis_C2_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
hideydecorations!(axis_C3_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
hideydecorations!(axis_ratio_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)

on(slider_refresh.value) do value
    label_refresh.text = "\nDisplay speed : $value"
end

on(slider_averaging.value) do value
    label_averaging.text = "\nDisplay averaging : $value"
end

on(slider_limit_C1.value) do value
    label_limit_C1.text = "Display limit C1 : $value"
    axis_C1_fig.scene.plots[1].colorrange = [0, value]
    colorbar_C1.limits = (0, round(Int, value))
end

on(slider_limit_C2.value) do value
    label_limit_C2.text = "Display limit C2 : $value"
    axis_C2_fig.scene.plots[1].colorrange = [0, value]
    colorbar_C2.limits = (0, round(Int, value))
end

on(slider_limit_C3.value) do value
    if value == 27
        label_limit_C3.text = "Display limit C3 : 27"
        axis_C3_fig.scene.plots[1].colorrange = [0, 26]
        colorbar_C3.limits = (0, 26)
    else
        label_limit_C3.text = "Display limit C3 : $value"
        axis_C3_fig.scene.plots[1].colorrange = [0, value]
        colorbar_C3.limits = (0, round(Int, value))
    end
end

on(button_start.clicks) do _
    if !running[]
        running[] = true
        ran[] = true
        saved = false

        button_start.label = "PAUSE"
        button_reset.label = "STOP"

        if menu_mode.selection[] == "Playback"
            @async start_playback(folderpath)
        else
            @async start_real_time(folderpath)
        end
    else
        running[] = false

        button_start.label = "START"
        button_reset.label = "CLEAR"
    end
end

on(button_reset.clicks) do _
    running[] = false
    ran[] = false
    saved = false
    initialized[] = false

    button_start.label = "START"
    button_reset.label = "CLEAR"

    current_ROI[] = 1
    current_noise[] = 1
    region_points[] = Point2f[]
    noise_matrix[] = zeros(512,1024)
    noise_size = Observable(0)

    for image in images
        image[] = zeros(UInt, 1024,512)
    end
    
    ROI_matrix[] = reshape(ROI_matrix[][1, :, :], 1, size(ROI_matrix[][1, :, :])...)
    ROI_sizes[] = [ROI_sizes[][1]]

    empty!(axis_C1_fig)
    empty!(axis_C2_fig)
    empty!(axis_C3_fig)
    empty!(axis_C1_plot)
    empty!(axis_C2_plot)
    empty!(axis_C3_plot)
    empty!(axis_ratio_plot)

    image!(axis_C1_fig, images[1], colorrange=[0,50], colormap=[:black, image_colors[menu_color_C1.i_selected[]]])
    image!(axis_C2_fig, images[2], colorrange=[0,50], colormap=[:black, image_colors[menu_color_C2.i_selected[]]])
    image!(axis_C3_fig, images[3], colorrange=[0,50], colormap=[:black, image_colors[menu_color_C3.i_selected[]]])

    poly!(axis_C1_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
    poly!(axis_C2_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
    poly!(axis_C3_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)

    scatter!(axis_C1_fig, region_points, markersize = 12, overdraw=true, color=strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)])
end

on(button_load.clicks) do _
    if !initialized[]
        global folderpath = get_new_folder_path()

        println(folderpath)
    end
end

on(button_save.clicks) do _
    if !running[] && ran[] && !saved
        saved = true
        println("Saved")
    elseif running[]
        saving[] = true
        println("Saving")
    end
end

on(button_ROI.clicks) do _
    if length(region_points[]) > 2
        poly!(axis_C1_fig, region_points[], color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
        poly!(axis_C2_fig, region_points[], color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
        poly!(axis_C3_fig, region_points[], color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)

        polygon_coords = []

        for _ in 1:length(region_points[])
            x, y = pop!(region_points[])
            push!(polygon_coords, (x, y))
        end

        polygon_points = [Luxor.Point(x, y) for (x, y) in polygon_coords]

        current_ROI_matrix = zeros(1, 512,1024)

        for i in CartesianIndices(current_ROI_matrix)
            point = Luxor.Point(i[3], i[2])

            if isinside(point, polygon_points, allowonedge=true)
                current_ROI_matrix[i] = 1
            end
        end

        ROI_matrix[] = cat(ROI_matrix[], current_ROI_matrix, dims=1)
        push!(ROI_sizes[], sum(ROI_matrix[][current_ROI[]+1, :, :]))
        # println(ROI_sizes[])

        current_ROI[] += 1

        # image!(axis_C1_plot, ROI_matrix[][current_ROI[], :, :]', colormap=[:black, :white])

        axis_C1_fig.scene.plots[2].color = ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)]
        axis_C2_fig.scene.plots[2].color = ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)]
        axis_C3_fig.scene.plots[2].color = ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)]

        axis_C1_fig.scene.plots[2].strokecolor = strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)]
        axis_C2_fig.scene.plots[2].strokecolor = strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)]
        axis_C3_fig.scene.plots[2].strokecolor = strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)]

        axis_C1_fig.scene.plots[3].color = strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)]

        notify(region_points)

        new_ROI[] += 1
    else
        println("3 points necessary")
    end
end

on(button_noise.clicks) do _
    if !isempty(region_points[])
        poly!(axis_C1_fig, region_points[], color=ROI_colors_dimmed[9], strokecolor=ROI_colors[9], strokewidth=2)
        poly!(axis_C2_fig, region_points[], color=ROI_colors_dimmed[9], strokecolor=ROI_colors[9], strokewidth=2)
        poly!(axis_C3_fig, region_points[], color=ROI_colors_dimmed[9], strokecolor=ROI_colors[9], strokewidth=2)

        polygon_coords = []

        for _ in 1:length(region_points[])
            x, y = pop!(region_points[])
            push!(polygon_coords, (x, y))
        end

        polygon_points = [Luxor.Point(x, y) for (x, y) in polygon_coords]

        for i in CartesianIndices(noise_matrix[])
            point = Luxor.Point(i[2], i[1])

            if isinside(point, polygon_points, allowonedge=true)
                noise_matrix[][i] = 1
            end
        end

        noise_size[] = sum(noise_matrix[])

        current_noise[] += 1

        # image!(axis_C2_plot, noise_matrix[]', colormap=[:black, :white])

        notify(region_points)
    else
        println("3 points necessary")
    end
end

on(menu_color_C1.i_selected) do selected_index
    colorbar_C1.colormap = [:black, image_colors[selected_index]]
    axis_C1_fig.scene.plots[1].colormap = [:black, image_colors[selected_index]]
end

on(menu_color_C2.i_selected) do selected_index
    colorbar_C2.colormap = [:black, image_colors[selected_index]]
    axis_C2_fig.scene.plots[1].colormap = [:black, image_colors[selected_index]]
end

on(menu_color_C3.i_selected) do selected_index
    colorbar_C3.colormap = [:black, image_colors[selected_index]]
    axis_C3_fig.scene.plots[1].colormap = [:black, image_colors[selected_index]]
end

on(menu_ratio.selection) do selection
    axis_ratio_plot.title = "Ratio of $selection"
end

on(textbox_time.displayed_string) do string
    if length(string) > 3
        textbox_time.displayed_string[] = string[begin:end-1]
        string = string[begin:end-1]
    end
end

on(textbox_time.stored_string) do string
    if !(parse(UInt16, string) > 1)
        textbox_time.stored_string[] = "2"
        textbox_time.displayed_string[] = "2"
        string = "2"
    end

    slider_time.range = 1:1:parse(UInt16, string)
end

on(events(fig).mousebutton, priority = 2) do event
    global dragging, idx
    if running[]
        if event.button == Mouse.left
            if event.action == Mouse.press
                plt, i = pick(fig)

                if Keyboard.d in events(fig).keyboardstate && plt == axis_C1_fig.scene.plots[3]
                    deleteat!(region_points[], i)
                    notify(region_points)

                    return Consume(true)
                elseif Keyboard.a in events(fig).keyboardstate
                    push!(region_points[], mouseposition(axis_C1_fig))
                    notify(region_points)

                    return Consume(true)
                else
                    dragging = plt == axis_C1_fig.scene.plots[3]
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

fig
