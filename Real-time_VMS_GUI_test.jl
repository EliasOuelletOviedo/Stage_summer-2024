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
using FileIO
using Images
using Clipper
using Luxor

@everywhere function start_real_time(folderpath, display_x_limits)
    println("Button pressed")

    intensity_C1 = Observable(Float64[])
    intensity_C2 = Observable(Float64[])
    intensity_C3 = Observable(Float64[])
    intensity_ratio = Observable(Float64[])

    real_time_C1 = Observable(Float64[])
    real_time_C2 = Observable(Float64[])
    real_time_C3 = Observable(Float64[])
    real_time_ratio = Observable(Float64[])

    ROI_intensity_C1 = Array([])
    ROI_intensity_C2 = Array([])
    ROI_intensity_C3 = Array([])
    ROI_intensity_ratio = Array([])

    ROI_real_time_C1 = Array([])
    ROI_real_time_C2 = Array([])
    ROI_real_time_C3 = Array([])
    ROI_real_time_ratio = Array([])

    lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
    lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
    lines!(axis_C3_plot, real_time_C3, intensity_C3, color=:black, linewidth=2)
    lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

    C1_path = joinpath(folderpath, "C1")
    C2_path = joinpath(folderpath, "C2")
    C3_path = joinpath(folderpath, "C3")

    previous_files_C1 = readdir(C1_path)
    previous_files_C2 = readdir(C2_path)
    previous_files_C3 = readdir(C3_path)

    t_0 = now()

    i = 1
    N = 1

    println("Search started")

    while running[]
        if new_ROI[] != 0
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_intensity_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_intensity_C1)], ROI_intensity_C2[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_intensity_C1)], ROI_intensity_ratio[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])

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

            if i == 1
                t_0 = now()
            end

            time = (now() - t_0)/Millisecond(1)/1000

            file_C1 = rawview(channelview(load(joinpath(C1_path, new_files_C1[end]))))'
            file_C2 = rawview(channelview(load(joinpath(C2_path, new_files_C2[end]))))'
            file_C3 = rawview(channelview(load(joinpath(C3_path, new_files_C3[end]))))'

            if i == 1 || slider_averaging.value[] == 1
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
                            array_C1[] -= 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i-N+k]))))'
                            array_C2[] -= 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i-N+k]))))'
                            array_C3[] -= 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i-N+k]))))'
                        end
    
                        array_C1[] *= N/slider_averaging.value[]
                        array_C2[] *= N/slider_averaging.value[]
                        array_C3[] *= N/slider_averaging.value[]
    
                        N = slider_averaging.value[]
                    end
                
                    array_C1[] += 1/N * file_C1 - 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i-N+1]))))'
                    array_C2[] += 1/N * file_C2 - 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i-N+1]))))'
                    array_C3[] += 1/N * file_C3 - 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i-N+1]))))'
                end
            end
    
            if sum(noise_matrix[]) != 0
                noise_C1 = sum(file_C1.*(noise_matrix[])')/sum(noise_matrix[])
                noise_C2 = sum(file_C2.*(noise_matrix[])')/sum(noise_matrix[])
                noise_C3 = sum(file_C3.*(noise_matrix[])')/sum(noise_matrix[])
            else
                noise_C1 = 0
                noise_C2 = 0
                noise_C3 = 0
            end
    
            for n in 1:length(ROI_intensity_C1)
                push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C1)
                push!(ROI_real_time_C1[n][], time)
    
                push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C2)
                push!(ROI_real_time_C2[n][], time)
    
                push!(ROI_intensity_C3[n][], sum(file_C3.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C3)
                push!(ROI_real_time_C3[n][], time)
    
                notify(ROI_intensity_C1[n])
                notify(ROI_real_time_C1[n])
    
                notify(ROI_intensity_C2[n])
                notify(ROI_real_time_C2[n])
    
                notify(ROI_intensity_C3[n])
                notify(ROI_real_time_C3[n])
            end
    
            push!(intensity_C1[], sum(file_C1)/sum(ROI_matrix[][1,:,:]) - noise_C1)
            push!(real_time_C1[], time)
    
            push!(intensity_C2[], sum(file_C2)/sum(ROI_matrix[][1,:,:]) - noise_C2)
            push!(real_time_C2[], time)
    
            push!(intensity_C3[], sum(file_C3)/sum(ROI_matrix[][1,:,:]) - noise_C3)
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
            xlims!(axis_C1_plot, max(0, time - display_x_limits), max(display_x_limits, time))
    
            autolimits!(axis_C2_plot)
            xlims!(axis_C2_plot, max(0, time - display_x_limits), max(display_x_limits, time))
    
            autolimits!(axis_C3_plot)
            xlims!(axis_C3_plot, max(0, time - display_x_limits), max(display_x_limits, time))
    
            autolimits!(axis_ratio_plot)
            xlims!(axis_ratio_plot, max(0, time - display_x_limits), max(display_x_limits, time))
    
            i += 1
    
            try
                label_frequency.text = "Frequency = $(round(1/(real_time_C1[][end] - real_time_C1[][end-1]), sigdigits=2)) Hz"
            catch e
                label_frequency.text = "Frequency = $(round(1/(real_time_C1[][end]), sigdigits=2)) Hz"
            end
        end
    
        sleep(0.001)
    end
end

@everywhere function start_real_time2(C1_path, C2_path, C3_path)
    println("Button pressed")

    intensity_C1 = Observable(Float64[])
    intensity_C2 = Observable(Float64[])
    intensity_ratio = Observable(Float64[])

    real_time_C1 = Observable(Float64[])
    real_time_C2 = Observable(Float64[])
    real_time_ratio = Observable(Float64[])

    ROI_intensity_C1 = Array([])
    ROI_intensity_C2 = Array([])
    ROI_intensity_ratio = Array([])

    ROI_real_time_C1 = Array([])
    ROI_real_time_C2 = Array([])
    ROI_real_time_ratio = Array([])

    averaged_array_C1 = zeros(10, 1024,512)
    averaged_array_C2 = zeros(10, 1024,512)
    averaged_array_C3 = zeros(10, 1024,512)

    lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
    lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
    lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

    previous_files_C1 = readdir(C1_path)
    previous_files_C2 = readdir(C2_path)
    previous_files_C3 = readdir(C3_path)

    newfile = false

    t_0 = now()

    i = 0
    j = 0
    k = 0
    l = 0

    println("Search started")

    while running[]
        if new_ROI[] != 0
            # println("New ROI")
            
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_intensity_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_intensity_C1)], ROI_intensity_C2[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_intensity_C1)], ROI_intensity_ratio[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])

            new_ROI[] -= 1

            # println("New ROI done")
        end
        
        current_files_C1 = readdir(C1_path)
        new_files_C1 = setdiff(current_files_C1, previous_files_C1)

        if !isempty(new_files_C1)

            # println("File detected!")

            newfile = true

            if i == 0 && j == 0 && k == 0 && l == 0
                t_0 = now()
            end

            # new_file_C1 = new_files_C1[end]

            for new_file_C1 in new_files_C1

                file_C1 = rawview(channelview(load(joinpath(C1_path, new_file_C1))))'

                if i == 0
                    for n in 1:10
                        averaged_array_C1[n,:,:] = file_C1
                    end

                    i += 1
                else
                    for n in 1:9
                        averaged_array_C1[n,:,:] = averaged_array_C1[n+1,:,:]
                    end

                    averaged_array_C1[end,:,:] = file_C1
                end

                present_time = (now() - t_0)/Millisecond(1)/1000

                if slider_averaging.value[] == 1
                    array_C1[] = averaged_array_C1[end,:,:]
                else
                    array_C1[] = sum(averaged_array_C1[end-(slider_averaging.value[]-1):end,:,:], dims=1)[1,:,:]/slider_averaging.value[]
                end

                if sum(noise_matrix[]) != 0
                    noise_C1 = sum(file_C1.*(noise_matrix[])')/sum(noise_matrix[])
                else
                    noise_C1 = 0
                end

                for n in 1:length(ROI_intensity_C1)
                    push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C1)
                    push!(ROI_real_time_C1[n][], present_time)

                    notify(ROI_intensity_C1[n])
                    notify(ROI_real_time_C1[n])
                end

                push!(intensity_C1[], sum(file_C1)/sum(ROI_matrix[][1,:,:]) - noise_C1)
                push!(real_time_C1[], present_time)

                notify(array_C1)
                notify(intensity_C1)
                notify(real_time_C1)

                # limits!(axis_C1_plot, max(0, present_time-10), max(10, present_time), minimum(intensity_C1[])*0.9, maximum(intensity_C1[])*1.1)
                limits!(axis_C1_plot, max(0, present_time-10), max(10, present_time), 0, 10)

                # println("C1 done")

                sleep(0.01)
            end
        end

        previous_files_C1 = current_files_C1

        current_files_C2 = readdir(C2_path)
        new_files_C2 = setdiff(current_files_C2, previous_files_C2)

        if !isempty(new_files_C2)

            # println("File detected!")

            newfile = true

            if i == 0 && j == 0 && k == 0 && l == 0
                t_0 = now()
            end

            for new_file_C2 in new_files_C2

            # new_file_C2 = new_files_C2[end]

                file_C2 = rawview(channelview(load(joinpath(C2_path, new_file_C2))))'

                if j == 0
                    for n in 1:10
                        averaged_array_C2[n,:,:] = file_C2
                    end

                    j += 1
                else
                    for n in 1:9
                        averaged_array_C2[n,:,:] = averaged_array_C2[n+1,:,:]
                    end

                    averaged_array_C2[end,:,:] = file_C2
                end

                present_time = (now() - t_0)/Millisecond(1)/1000

                if slider_averaging.value[] == 1
                    array_C2[] = averaged_array_C2[end,:,:]
                else
                    array_C2[] = sum(averaged_array_C2[end-(slider_averaging.value[]-1):end,:,:], dims=1)[1,:,:]/slider_averaging.value[]
                end

                if sum(noise_matrix[]) != 0
                    noise_C2 = sum(file_C1.*(noise_matrix[])')/sum(noise_matrix[])
                else
                    noise_C2 = 0
                end

                for n in 1:length(ROI_intensity_C2)
                    push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:]) - noise_C2)
                    push!(ROI_real_time_C2[n][], present_time)

                    notify(ROI_intensity_C2[n])
                    notify(ROI_real_time_C2[n])
                end

                push!(intensity_C2[], sum(file_C2)/sum(ROI_matrix[][1,:,:]) - noise_C2)
                push!(real_time_C2[], present_time)

                notify(array_C2)
                notify(intensity_C2)
                notify(real_time_C2)

                # limits!(axis_C2_plot, max(0, present_time-10), max(10, present_time), minimum(intensity_C2[])*0.9, maximum(intensity_C2[])*1.1)
                limits!(axis_C2_plot,  max(0, present_time-10), max(10, present_time), 0, 10)

                # println("C2 done")

                sleep(0.01)
            end
        end

        previous_files_C2 = current_files_C2

        current_files_C3 = readdir(C3_path)
        new_files_C3 = setdiff(current_files_C3, previous_files_C3)

        if !isempty(new_files_C3)

            # println("File detected!")

            # if i == 0 && j == 0 && k == 0 && l == 0
            #     t_0 = now()
            # end

            for new_file_C3 in new_files_C3

            # new_file_C3 = new_files_C3[end]

                file_C3 = rawview(channelview(load(joinpath(C3_path, new_file_C3))))'

                if k == 0
                    for n in 1:10
                        averaged_array_C3[n,:,:] = file_C3
                    end

                    k += 1
                else
                    for n in 1:9
                        averaged_array_C3[n,:,:] = averaged_array_C3[n+1,:,:]
                    end

                    averaged_array_C3[end,:,:] = file_C3
                end

                if slider_averaging.value[] == 1
                    array_C3[] = averaged_array_C3[end,:,:]
                else
                    array_C3[] = sum(averaged_array_C3[end-(slider_averaging.value[]-1):end,:,:], dims=1)[1,:,:]/slider_averaging.value[]
                end

                notify(array_C3)

                sleep(0.01)
            end
        end

        previous_files_C3 = current_files_C3

        if newfile && i == 1 && j == 1
            # println("Ratio")

            present_time = (now() - t_0)/Millisecond(1)/1000

            push!(real_time_ratio[], present_time)

            if menu_ratio.selection[] == "C1/C2"
                push!(intensity_ratio[], intensity_C1[][end]/intensity_C2[][end])

                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C1[n][]) && !isempty(ROI_intensity_C2[n][])
                        if ROI_intensity_C2[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C1[n][][end]/ROI_intensity_C2[n][][end])
                            push!(ROI_real_time_ratio[n][], present_time)

                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end

            elseif menu_ratio.selection[] == "C2/C1"
                push!(intensity_ratio[], intensity_C2[][end]/intensity_C1[][end])

                for n in 1:length(ROI_intensity_ratio)
                    if !isempty(ROI_intensity_C1[n][]) && !isempty(ROI_intensity_C2[n][])
                        if ROI_intensity_C1[n][][end] != 0
                            push!(ROI_intensity_ratio[n][], ROI_intensity_C2[n][][end]/ROI_intensity_C1[n][][end])
                            push!(ROI_real_time_ratio[n][], present_time)

                            notify(ROI_intensity_ratio[n])
                            notify(ROI_real_time_ratio[n])
                        end
                    end
                end
            end

            notify(intensity_ratio)
            notify(real_time_ratio)

            limits!(axis_ratio_plot,  max(0, present_time-10), max(10, present_time), 0, 5)

            newfile = false
        end

        sleep(0.01)
    end
end

@everywhere function start_playback(folderpath, display_x_limits)
    println("Button pressed")

    intensity_C1 = Observable(Float64[])
    intensity_C2 = Observable(Float64[])
    intensity_C3 = Observable(Float64[])
    intensity_ratio = Observable(Float64[])

    real_time_C1 = Observable(Float64[])
    real_time_C2 = Observable(Float64[])
    real_time_C3 = Observable(Float64[])
    real_time_ratio = Observable(Float64[])

    ROI_intensity_C1 = Array([])
    ROI_intensity_C2 = Array([])
    ROI_intensity_C3 = Array([])
    ROI_intensity_ratio = Array([])

    ROI_real_time_C1 = Array([])
    ROI_real_time_C2 = Array([])
    ROI_real_time_C3 = Array([])
    ROI_real_time_ratio = Array([])

    lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
    lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
    lines!(axis_C3_plot, real_time_C3, intensity_C3, color=:black, linewidth=2)
    lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

    C1_files = readdir(joinpath(folderpath, "C1"))
    C2_files = readdir(joinpath(folderpath, "C2"))
    C3_files = readdir(joinpath(folderpath, "C3"))

    if C1_files[1] == ".DS_Store"
        deleteat!(C1_files, 1)
    end

    if C2_files[1] == ".DS_Store"
        deleteat!(C2_files, 1)
    end

    if C3_files[1] == ".DS_Store"
        deleteat!(C3_files, 1)
    end

    files = true

    i = 1
    j = 0
    N = 1

    previous_time = now()

    fps = 31

    println("Starting")

    while files && running[]
        if new_ROI[] != 0
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_C3, Observable(Float64[]))
            push!(ROI_real_time_C3, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_intensity_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_intensity_C1)], ROI_intensity_C2[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C3_plot, ROI_real_time_C3[length(ROI_intensity_C1)], ROI_intensity_C3[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_intensity_C1)], ROI_intensity_ratio[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])

            new_ROI[] -= 1
        end

        time = (i - 1 + j)/fps

        try
            C1_files[i]
        catch e
            println("Restarting")
            j += i - 1
            i = 1
            N = 1
        end

        file_C1 = rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i]))))'
        file_C2 = rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i]))))'
        file_C3 = rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i]))))'

        if i == 1 && j == 0 || slider_averaging.value[] == 1
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
                        array_C1[] -= 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i-N+k]))))'
                        array_C2[] -= 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i-N+k]))))'
                        array_C3[] -= 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i-N+k]))))'
                    end

                    array_C1[] *= N/slider_averaging.value[]
                    array_C2[] *= N/slider_averaging.value[]
                    array_C3[] *= N/slider_averaging.value[]

                    N = slider_averaging.value[]
                end
            
                array_C1[] += 1/N * file_C1 - 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i-N+1]))))'
                array_C2[] += 1/N * file_C2 - 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i-N+1]))))'
                array_C3[] += 1/N * file_C3 - 1/N * rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i-N+1]))))'
            end
        end

        if sum(noise_matrix[]) != 0
            noise_C1 = sum(file_C1.*(noise_matrix[])')/sum(noise_matrix[])
            noise_C2 = sum(file_C2.*(noise_matrix[])')/sum(noise_matrix[])
            noise_C3 = sum(file_C3.*(noise_matrix[])')/sum(noise_matrix[])
        else
            noise_C1 = 0
            noise_C2 = 0
            noise_C3 = 0
        end

        for n in 1:length(ROI_intensity_C1)
            push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C1)
            push!(ROI_real_time_C1[n][], time)

            push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C2)
            push!(ROI_real_time_C2[n][], time)

            push!(ROI_intensity_C3[n][], sum(file_C3.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C3)
            push!(ROI_real_time_C3[n][], time)

            notify(ROI_intensity_C1[n])
            notify(ROI_real_time_C1[n])

            notify(ROI_intensity_C2[n])
            notify(ROI_real_time_C2[n])

            notify(ROI_intensity_C3[n])
            notify(ROI_real_time_C3[n])
        end

        push!(intensity_C1[], sum(file_C1)/sum(ROI_matrix[][1,:,:]) - noise_C1)
        push!(real_time_C1[], time)

        push!(intensity_C2[], sum(file_C2)/sum(ROI_matrix[][1,:,:]) - noise_C2)
        push!(real_time_C2[], time)

        push!(intensity_C3[], sum(file_C3)/sum(ROI_matrix[][1,:,:]) - noise_C3)
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
        xlims!(axis_C1_plot, max(0, time - display_x_limits), max(display_x_limits, time))

        autolimits!(axis_C2_plot)
        xlims!(axis_C2_plot, max(0, time - display_x_limits), max(display_x_limits, time))

        autolimits!(axis_C3_plot)
        xlims!(axis_C3_plot, max(0, time - display_x_limits), max(display_x_limits, time))

        autolimits!(axis_ratio_plot)
        xlims!(axis_ratio_plot, max(0, time - display_x_limits), max(display_x_limits, time))

        i += 1

        current_time = now()
        label_frequency.text = "Frequency = $(string(round(1/((current_time - previous_time)/Millisecond(1)/1000), sigdigits=2))) Hz"
        previous_time = current_time

        sleep(0.001)
    end
end

@everywhere function start_playback1(folderpath, display_x_limits, display_averaging_limits)
    println("Button pressed")

    intensity_C1 = Observable(Float64[])
    intensity_C2 = Observable(Float64[])
    intensity_C3 = Observable(Float64[])
    intensity_ratio = Observable(Float64[])

    real_time_C1 = Observable(Float64[])
    real_time_C2 = Observable(Float64[])
    real_time_C3 = Observable(Float64[])
    real_time_ratio = Observable(Float64[])

    ROI_intensity_C1 = Array([])
    ROI_intensity_C2 = Array([])
    ROI_intensity_C3 = Array([])
    ROI_intensity_ratio = Array([])

    ROI_real_time_C1 = Array([])
    ROI_real_time_C2 = Array([])
    ROI_real_time_C3 = Array([])
    ROI_real_time_ratio = Array([])

    averaged_array_C1 = zeros(display_averaging_limits, 1024,512)
    averaged_array_C2 = zeros(display_averaging_limits, 1024,512)
    averaged_array_C3 = zeros(display_averaging_limits, 1024,512)

    lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
    lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
    lines!(axis_C3_plot, real_time_C3, intensity_C3, color=:black, linewidth=2)
    lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

    C1_files = readdir(joinpath(folderpath, "C1"))
    C2_files = readdir(joinpath(folderpath, "C2"))
    C3_files = readdir(joinpath(folderpath, "C3"))

    if C1_files[1] == ".DS_Store"
        deleteat!(C1_files, 1)
    end

    if C2_files[1] == ".DS_Store"
        deleteat!(C2_files, 1)
    end

    if C3_files[1] == ".DS_Store"
        deleteat!(C3_files, 1)
    end

    files = true

    i = 1
    j = 0
    k = 0

    previous_time = now()

    fps = 31

    println("Starting")

    while files && running[]
        if new_ROI[] != 0
            # println("New ROI")
            
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_C3, Observable(Float64[]))
            push!(ROI_real_time_C3, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_intensity_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_intensity_C1)], ROI_intensity_C2[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C3_plot, ROI_real_time_C3[length(ROI_intensity_C1)], ROI_intensity_C3[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_intensity_C1)], ROI_intensity_ratio[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])

            new_ROI[] -= 1

            # println("New ROI done")
        end

        time = (i - 1 + j)/fps

        index = mod1(k, display_averaging_limits)

        try
            C1_files[i]
        catch e
            println("Restarting")
            j += i - 1
            i = 1
        end

        file_C1 = rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i]))))'
        file_C2 = rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i]))))'
        file_C3 = rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i]))))'

        if i == 1 && j == 0
            for n in 1:display_averaging_limits
                averaged_array_C1[n,:,:] = file_C1
                averaged_array_C2[n,:,:] = file_C2
                averaged_array_C3[n,:,:] = file_C3
            end
        else
            averaged_array_C1[index,:,:] = file_C1
            averaged_array_C2[index,:,:] = file_C2
            averaged_array_C3[index,:,:] = file_C3
        end

        if slider_averaging.value[] == 1
            array_C1[] = averaged_array_C1[index,:,:]
            array_C2[] = averaged_array_C2[index,:,:]
            array_C3[] = averaged_array_C3[index,:,:]
        else
            if slider_averaging.value[] <= index
                array_C1[] = sum(averaged_array_C1[index - (slider_averaging.value[] - 1):index,:,:], dims=1)[1,:,:]/slider_averaging.value[]
                array_C2[] = sum(averaged_array_C2[index - (slider_averaging.value[] - 1):index,:,:], dims=1)[1,:,:]/slider_averaging.value[]
                array_C3[] = sum(averaged_array_C3[index - (slider_averaging.value[] - 1):index,:,:], dims=1)[1,:,:]/slider_averaging.value[]
            else
                array_C1[] = (sum(averaged_array_C1[1:index,:,:], dims=1)[1,:,:] + sum(averaged_array_C1[index - (slider_averaging.value[] - 1) + display_averaging_limits:end,:,:], dims=1)[1,:,:])/slider_averaging.value[]
                array_C2[] = (sum(averaged_array_C2[1:index,:,:], dims=1)[1,:,:] + sum(averaged_array_C2[index - (slider_averaging.value[] - 1) + display_averaging_limits:end,:,:], dims=1)[1,:,:])/slider_averaging.value[]
                array_C3[] = (sum(averaged_array_C3[1:index,:,:], dims=1)[1,:,:] + sum(averaged_array_C3[index - (slider_averaging.value[] - 1) + display_averaging_limits:end,:,:], dims=1)[1,:,:])/slider_averaging.value[]
            end
        end

        if sum(noise_matrix[]) != 0
            noise_C1 = sum(file_C1.*(noise_matrix[])')/sum(noise_matrix[])
            noise_C2 = sum(file_C2.*(noise_matrix[])')/sum(noise_matrix[])
            noise_C3 = sum(file_C3.*(noise_matrix[])')/sum(noise_matrix[])
        else
            noise_C1 = 0
            noise_C2 = 0
            noise_C3 = 0
        end

        for n in 1:length(ROI_intensity_C1)
            push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C1)
            push!(ROI_real_time_C1[n][], time)
            push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C2)
            push!(ROI_real_time_C2[n][], time)
            push!(ROI_intensity_C3[n][], sum(file_C3.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C3)
            push!(ROI_real_time_C3[n][], time)

            notify(ROI_intensity_C1[n])
            notify(ROI_real_time_C1[n])
            notify(ROI_intensity_C2[n])
            notify(ROI_real_time_C2[n])
            notify(ROI_intensity_C3[n])
            notify(ROI_real_time_C3[n])
        end

        push!(intensity_C1[], sum(file_C1)/sum(ROI_matrix[][1,:,:]) - noise_C1)
        push!(real_time_C1[], time)
        push!(intensity_C2[], sum(file_C2)/sum(ROI_matrix[][1,:,:]) - noise_C2)
        push!(real_time_C2[], time)
        push!(intensity_C3[], sum(file_C3)/sum(ROI_matrix[][1,:,:]) - noise_C3)
        push!(real_time_C3[], time)

        notify(array_C1)
        notify(intensity_C1)
        notify(real_time_C1)
        notify(array_C2)
        notify(intensity_C2)
        notify(real_time_C2)
        notify(array_C3)
        notify(intensity_C3)
        notify(real_time_C3)

        push!(real_time_ratio[], time)

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

        notify(intensity_ratio)
        notify(real_time_ratio)

        # limits!(axis_C1_plot, max(0, present_time-10), max(10, present_time), minimum(intensity_C1[])*0.9, maximum(intensity_C1[])*1.1)
        autolimits!(axis_C1_plot)
        xlims!(axis_C1_plot, max(0, time - display_x_limits), max(display_x_limits, time))
        autolimits!(axis_C2_plot)
        xlims!(axis_C2_plot, max(0, time - display_x_limits), max(display_x_limits, time))
        autolimits!(axis_C3_plot)
        xlims!(axis_C3_plot, max(0, time - display_x_limits), max(display_x_limits, time))
        autolimits!(axis_ratio_plot)
        xlims!(axis_ratio_plot, max(0, time - display_x_limits), max(display_x_limits, time))

        i += 1

        current_time = now()

        # println(round(1/((current_time - previous_time)/Millisecond(1)/1000), sigdigits=2), " Hz")

        label_frequency.text = "Frequency = $(string(round(1/((current_time - previous_time)/Millisecond(1)/1000), sigdigits=2))) Hz"

        previous_time = current_time

        k += 1

        sleep(0.001)
    end
end

@everywhere function start_playback2(folderpath, display_x_limits, display_averaging_limits)
    println("Button pressed")

    intensity_C1 = Observable(Float64[])
    intensity_C2 = Observable(Float64[])
    intensity_C3 = Observable(Float64[])
    intensity_ratio = Observable(Float64[])

    real_time_C1 = Observable(Float64[])
    real_time_C2 = Observable(Float64[])
    real_time_C3 = Observable(Float64[])
    real_time_ratio = Observable(Float64[])

    ROI_intensity_C1 = Array([])
    ROI_intensity_C2 = Array([])
    ROI_intensity_C3 = Array([])
    ROI_intensity_ratio = Array([])

    ROI_real_time_C1 = Array([])
    ROI_real_time_C2 = Array([])
    ROI_real_time_C3 = Array([])
    ROI_real_time_ratio = Array([])

    averaged_array_C1 = zeros(display_averaging_limits, 1024,512)
    averaged_array_C2 = zeros(display_averaging_limits, 1024,512)
    averaged_array_C3 = zeros(display_averaging_limits, 1024,512)

    lines!(axis_C1_plot, real_time_C1, intensity_C1, color=:black, linewidth=2)
    lines!(axis_C2_plot, real_time_C2, intensity_C2, color=:black, linewidth=2)
    lines!(axis_C3_plot, real_time_C3, intensity_C3, color=:black, linewidth=2)
    lines!(axis_ratio_plot, real_time_ratio, intensity_ratio, color=:black, linewidth=2)

    C1_files = readdir(joinpath(folderpath, "C1"))
    C2_files = readdir(joinpath(folderpath, "C2"))
    C3_files = readdir(joinpath(folderpath, "C3"))

    if C1_files[1] == ".DS_Store"
        deleteat!(C1_files, 1)
    end

    if C2_files[1] == ".DS_Store"
        deleteat!(C2_files, 1)
    end

    if C3_files[1] == ".DS_Store"
        deleteat!(C3_files, 1)
    end

    files = true

    i = 1
    j = 0

    previous_time = now()

    fps = 31

    println("Starting")

    while files && running[]
        if new_ROI[] != 0
            # println("New ROI")
            
            push!(ROI_intensity_C1, Observable(Float64[]))
            push!(ROI_real_time_C1, Observable(Float64[]))

            push!(ROI_intensity_C2, Observable(Float64[]))
            push!(ROI_real_time_C2, Observable(Float64[]))

            push!(ROI_intensity_C3, Observable(Float64[]))
            push!(ROI_real_time_C3, Observable(Float64[]))

            push!(ROI_intensity_ratio, Observable(Float64[]))
            push!(ROI_real_time_ratio, Observable(Float64[]))

            lines!(axis_C1_plot, ROI_real_time_C1[length(ROI_intensity_C1)], ROI_intensity_C1[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C2_plot, ROI_real_time_C2[length(ROI_intensity_C1)], ROI_intensity_C2[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_C3_plot, ROI_real_time_C3[length(ROI_intensity_C1)], ROI_intensity_C3[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])
            lines!(axis_ratio_plot, ROI_real_time_ratio[length(ROI_intensity_C1)], ROI_intensity_ratio[length(ROI_intensity_C1)], color=ROI_colors[mod1(length(ROI_intensity_C1), length(ROI_colors)-1)])

            new_ROI[] -= 1

            # println("New ROI done")
        end

        time = (i - 1 + j)/fps

        try
            C1_files[i]
        catch e
            println("Restarting")
            j += i - 1
            i = 1
        end

        file_C1 = rawview(channelview(load(joinpath(joinpath(folderpath, "C1"), C1_files[i]))))'
        file_C2 = rawview(channelview(load(joinpath(joinpath(folderpath, "C2"), C2_files[i]))))'
        file_C3 = rawview(channelview(load(joinpath(joinpath(folderpath, "C3"), C3_files[i]))))'

        if i == 1 && j == 0
            for n in 1:display_averaging_limits
                averaged_array_C1[n,:,:] = file_C1
                averaged_array_C2[n,:,:] = file_C2
                averaged_array_C3[n,:,:] = file_C3
            end
        else
            for n in 1:display_averaging_limits-1
                averaged_array_C1[n,:,:] = averaged_array_C1[n+1,:,:]
                averaged_array_C2[n,:,:] = averaged_array_C2[n+1,:,:]
                averaged_array_C3[n,:,:] = averaged_array_C3[n+1,:,:]
            end

            averaged_array_C1[end,:,:] = file_C1
            averaged_array_C2[end,:,:] = file_C2
            averaged_array_C3[end,:,:] = file_C3
        end

        if slider_averaging.value[] == 1
            array_C1[] = averaged_array_C1[end,:,:]
            array_C2[] = averaged_array_C2[end,:,:]
            array_C3[] = averaged_array_C3[end,:,:]
        else
            array_C1[] = sum(averaged_array_C1[end-(slider_averaging.value[]-1):end,:,:], dims=1)[1,:,:]/slider_averaging.value[]
            array_C2[] = sum(averaged_array_C2[end-(slider_averaging.value[]-1):end,:,:], dims=1)[1,:,:]/slider_averaging.value[]
            array_C3[] = sum(averaged_array_C3[end-(slider_averaging.value[]-1):end,:,:], dims=1)[1,:,:]/slider_averaging.value[]
        end

        if sum(noise_matrix[]) != 0
            noise_C1 = sum(file_C1.*(noise_matrix[])')/sum(noise_matrix[])
            noise_C2 = sum(file_C2.*(noise_matrix[])')/sum(noise_matrix[])
            noise_C3 = sum(file_C3.*(noise_matrix[])')/sum(noise_matrix[])
        else
            noise_C1 = 0
            noise_C2 = 0
            noise_C3 = 0
        end

        for n in 1:length(ROI_intensity_C1)
            push!(ROI_intensity_C1[n][], sum(file_C1.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C1)
            push!(ROI_real_time_C1[n][], time)
            push!(ROI_intensity_C2[n][], sum(file_C2.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C2)
            push!(ROI_real_time_C2[n][], time)
            push!(ROI_intensity_C3[n][], sum(file_C3.*(ROI_matrix[][n+1,:,:])')/sum(ROI_matrix[][n+1,:,:])[1,1,1] - noise_C3)
            push!(ROI_real_time_C3[n][], time)

            notify(ROI_intensity_C1[n])
            notify(ROI_real_time_C1[n])
            notify(ROI_intensity_C2[n])
            notify(ROI_real_time_C2[n])
            notify(ROI_intensity_C3[n])
            notify(ROI_real_time_C3[n])
        end

        push!(intensity_C1[], sum(file_C1)/sum(ROI_matrix[][1,:,:]) - noise_C1)
        push!(real_time_C1[], time)
        push!(intensity_C2[], sum(file_C2)/sum(ROI_matrix[][1,:,:]) - noise_C2)
        push!(real_time_C2[], time)
        push!(intensity_C3[], sum(file_C3)/sum(ROI_matrix[][1,:,:]) - noise_C3)
        push!(real_time_C3[], time)

        notify(array_C1)
        notify(intensity_C1)
        notify(real_time_C1)
        notify(array_C2)
        notify(intensity_C2)
        notify(real_time_C2)
        notify(array_C3)
        notify(intensity_C3)
        notify(real_time_C3)

        push!(real_time_ratio[], time)

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

        notify(intensity_ratio)
        notify(real_time_ratio)

        # limits!(axis_C1_plot, max(0, present_time-10), max(10, present_time), minimum(intensity_C1[])*0.9, maximum(intensity_C1[])*1.1)
        autolimits!(axis_C1_plot)
        xlims!(axis_C1_plot, max(0, time - display_x_limits), max(display_x_limits, time))
        autolimits!(axis_C2_plot)
        xlims!(axis_C2_plot, max(0, time - display_x_limits), max(display_x_limits, time))
        autolimits!(axis_C3_plot)
        xlims!(axis_C3_plot, max(0, time - display_x_limits), max(display_x_limits, time))
        autolimits!(axis_ratio_plot)
        xlims!(axis_ratio_plot, max(0, time - display_x_limits), max(display_x_limits, time))

        i += 1

        current_time = now()

        # println(round(1/((current_time - previous_time)/Millisecond(1)/1000), sigdigits=2), " Hz")

        label_frequency.text = "Frequency = $(string(round(1/((current_time - previous_time)/Millisecond(1)/1000), sigdigits=2))) Hz"

        previous_time = current_time

        sleep(0.001)
    end
end

image_colors = [RGB(1,0.25,0.25), RGB(0.25,1,0.25), RGB(0.25,0.25,1), RGB(0.25,1,1), RGB(1,0.25,1), RGB(1,1,0.25), RGB(1, 1, 1)]
ROI_colors = [RGB(1,0,0), RGB(1,0.5,0), RGB(1,1,0), RGB(0,1,0), RGB(0,1,1), RGB(0,0,1), RGB(0.5,0,1), RGB(1,0,1), RGB(1,1,1)]
ROI_colors_dimmed = [RGBA(1,0,0, 0.1), RGBA(1,0.5,0, 0.1), RGBA(1,1,0, 0.1), RGBA(0,1,0, 0.1), RGBA(0,1,1, 0.1), RGBA(0,0,1, 0.1), RGBA(0.5,0,1, 0.1), RGBA(1,0,1, 0.1), RGBA(1,1,1, 0.1)]

set_theme!(backgroundcolor = :gray95)

display_x_limits = 10
display_averaging_limits = 50

# folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/VMS_test"
# folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-18_TestSoftware/test save online-001/Bliq VMS"
folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-18_TestSoftware/test save online-001/Test VMS"

C1_path = joinpath(folderpath, "C1")
C2_path = joinpath(folderpath, "C2")
C3_path = joinpath(folderpath, "C3")

array_C1 = Observable(zeros(1024,512))
array_C1[][1,1] = 1
array_C2 = Observable(zeros(1024,512))
array_C2[][1,1] = 1
array_C3 = Observable(zeros(1024,512))
array_C3[][1,1] = 1

fig = Figure(size = (1440, 847))

fig[1, 1] = buttongrid = GLMakie.GridLayout(tellheight=false)
fig[1, 2] = figgrid = GLMakie.GridLayout()
fig[1, 3] = plotgrid = GLMakie.GridLayout()

colsize!(fig.layout, 1, GLMakie.Fixed(300))

label_frequency = buttongrid[0, 1] = Label(fig, "Frequency : 0 Hz", halign=:left, fontsize=12)

button_start = buttongrid[1, 1] = Button(fig, label="START", font=:bold, height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)
button_reset = buttongrid[1, 2] = Button(fig, label="RESET", font=:bold, height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)

button_ROI = buttongrid[3, 1] = Button(fig, label="Set ROI", height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)
button_noise = buttongrid[3, 2] = Button(fig, label="Set noise", height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)

label_color_C1 = buttongrid[5, 1] = Label(fig, "Color C1 :", halign=:right, fontsize=20)
label_color_C2 = buttongrid[6, 1] = Label(fig, "Color C2 :", halign=:right, fontsize=20)
label_color_C3 = buttongrid[7, 1] = Label(fig, "Color C3 :", halign=:right, fontsize=20)

menu_color_C1 = buttongrid[5, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "White"], default = "Red", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)
menu_color_C2 = buttongrid[6, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "White"], default = "Green", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)
menu_color_C3 = buttongrid[7, 2] = Menu(fig, options=["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "White"], default = "White", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)

label_ratio = buttongrid[9, 1] = Label(fig, "Ratio : ", halign=:right, fontsize=20)
menu_ratio = buttongrid[9, 2] = Menu(fig, options=["C1/C2", "C1/C3", "C2/C1", "C2/C3", "C3/C1", "C3/C2"], default = "C1/C2", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)

label_averaging = buttongrid[11, 1:2] = Label(fig, "Display averaging : 1", fontsize=20)
slider_averaging = buttongrid[12, 1:2] = Slider(fig, range = 1:1:display_averaging_limits, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_limit_C1 = buttongrid[13, 1:2] = Label(fig, "Display limit C1 : 50", fontsize=20)
slider_limit_C1 = buttongrid[14, 1:2] = Slider(fig, range = 3:1:99, startvalue = 50, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_limit_C2 = buttongrid[15, 1:2] = Label(fig, "Display limit C2 : 50", fontsize=20)
slider_limit_C2 = buttongrid[16, 1:2] = Slider(fig, range = 3:1:99, startvalue = 50, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

label_limit_C3 = buttongrid[17, 1:2] = Label(fig, "Display limit C3 : 50", fontsize=20)
slider_limit_C3 = buttongrid[18, 1:2] = Slider(fig, range = 6:1:99, startvalue = 50, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8))

axis_C1_fig = figgrid[1, 1][1, 1] = GLMakie.Axis(fig, yreversed = true, xrectzoom=false, yrectzoom=false, title="Channel 1", titlesize=16)
axis_C2_fig = figgrid[2, 1][1, 1] = GLMakie.Axis(fig, yreversed = true, xrectzoom=false, yrectzoom=false, title="Channel 2", titlesize=16)
axis_C3_fig = figgrid[3, 1][1, 1] = GLMakie.Axis(fig, yreversed = true, xrectzoom=false, yrectzoom=false, title="Channel 3", titlesize=16)

global axis_C1_plot = plotgrid[1, 1] = GLMakie.Axis(fig, title="Channel 1", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))
global axis_C2_plot = plotgrid[2, 1] = GLMakie.Axis(fig, title="Channel 2", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))
global axis_C3_plot = plotgrid[3, 1] = GLMakie.Axis(fig, title="Channel 3", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))
global axis_ratio_plot = plotgrid[4, 1] = GLMakie.Axis(fig, title="Ratio of C1/C2", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16, xlabelsize=14, ylabelsize=14, xticks = LinearTicks(5))

limits!(axis_C1_fig, 0, 1024, 0, 512)
limits!(axis_C2_fig, 0, 1024, 0, 512)
limits!(axis_C3_fig, 0, 1024, 0, 512)

hidedecorations!(axis_C1_fig)
hidedecorations!(axis_C2_fig)
hidedecorations!(axis_C3_fig)

colorbar_C1 = figgrid[1, 1][1, 2] = Colorbar(fig, halign=:left, limits=(0, 50), colormap=[:black, image_colors[menu_color_C1.i_selected[]]])
colorbar_C2 = figgrid[2, 1][1, 2] = Colorbar(fig, halign=:left, limits=(0, 50), colormap=[:black, image_colors[menu_color_C2.i_selected[]]])
colorbar_C3 = figgrid[3, 1][1, 2] = Colorbar(fig, halign=:left, limits=(0, 50), colormap=[:black, image_colors[menu_color_C3.i_selected[]]])

img_C1 = image!(axis_C1_fig, array_C1, colorrange=[0,50], colormap=[:black, image_colors[menu_color_C1.i_selected[]]])
img_C2 = image!(axis_C2_fig, array_C2, colorrange=[0,50], colormap=[:black, image_colors[menu_color_C2.i_selected[]]])
img_C3 = image!(axis_C3_fig, array_C3, colorrange=[0,50], colormap=[:black, image_colors[menu_color_C3.i_selected[]]])

dragging = false
idx = 1

running = Observable(false)
current_ROI = Observable(1)
current_noise = Observable(1)
region_points = Observable(Point2f[])

new_ROI = Observable(0)

ROI_matrix = Observable(ones(1, 512,1024))
noise_matrix = Observable(zeros(512,1024))

region_plot_C1 = poly!(axis_C1_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
region_plot_C2 = poly!(axis_C2_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
region_plot_C3 = poly!(axis_C3_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)

region_scatter = scatter!(axis_C1_fig, region_points, markersize = 12, overdraw=true, color=strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)])

hideydecorations!(axis_C1_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
hideydecorations!(axis_C2_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
hideydecorations!(axis_C3_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
hideydecorations!(axis_ratio_plot, label = false, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)

on(slider_averaging.value) do value
    label_averaging.text = "Display averaging : $value"
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
    label_limit_C3.text = "Display limit C3 : $value"
    axis_C3_fig.scene.plots[1].colorrange = [0, value]
    colorbar_C3.limits = (0, round(Int, value))
end

on(button_start.clicks) do _
    if !running[]
        running[] = true
        # @async start_real_time(folderpath, display_x_limits)
        @async start_playback(folderpath, display_x_limits)
    end
end

on(button_reset.clicks) do _
    running[] = false
    current_ROI[] = 1
    current_noise[] = 1
    region_points[] = Point2f[]
    noise_matrix[] = zeros(512,1024)

    array_C1[] = zeros(1024,512)
    array_C1[][1,1] = 1
    array_C2[] = zeros(1024,512)
    array_C2[][1,1] = 1
    array_C3[] = zeros(1024,512)
    array_C3[][1,1] = 1
    
    a_slice = ROI_matrix[][1, :, :]
    ROI_matrix[] = reshape(a_slice, 1, size(a_slice)...)

    empty!(axis_C1_fig)
    empty!(axis_C2_fig)
    empty!(axis_C3_fig)
    empty!(axis_C1_plot)
    empty!(axis_C2_plot)
    empty!(axis_C3_plot)
    empty!(axis_ratio_plot)

    image!(axis_C1_fig, array_C1, colormap=[:black, image_colors[menu_color_C1.i_selected[]]], colorrange = [0, slider_limit_C1.value[]])
    image!(axis_C2_fig, array_C2, colormap=[:black, image_colors[menu_color_C2.i_selected[]]], colorrange = [0, slider_limit_C2.value[]])
    image!(axis_C3_fig, array_C3, colormap=[:black, image_colors[menu_color_C3.i_selected[]]], colorrange = [0, slider_limit_C3.value[]])

    poly!(axis_C1_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
    poly!(axis_C2_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)
    poly!(axis_C3_fig, region_points, color=ROI_colors_dimmed[mod1(current_ROI[], length(ROI_colors_dimmed)-1)], strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)], strokewidth=2)

    scatter!(axis_C1_fig, region_points, markersize = 12, overdraw=true, color=strokecolor=ROI_colors[mod1(current_ROI[], length(ROI_colors)-1)])
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
