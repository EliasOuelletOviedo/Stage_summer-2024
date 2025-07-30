using GLMakie
using Statistics
using NativeFileDialog
using ZipFile
using Observables
using Dates
using DelimitedFiles
using LibSerialPort
using CSV
using DataFrames
using SavitzkyGolay
using Lowess
using Optim
using LsqFit

include("vec_to_lifetime.jl")

function get_folder_path(str)
    if isfile(str*"_folder_path.txt")
        folderpath = open(f->read(f, String), str*"_folder_path.txt")

        if isempty(folderpath)
            folderpath = pick_folder()

            open(str*"_folder_path.txt", "w") do io
                write(io, folderpath)
            end
        end
    else
        folderpath = pick_folder()

        open(str*"_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    end

    if isempty(folderpath)
        println("Folder does not exist")

        folderpath = get_folder_path(str)
    end

    return folderpath
end

function get_new_folder_path(str)
    folderpath = pick_folder()

    if !isempty(folderpath)
        open(str*"_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    else
        folderpath = open(f->read(f, String), str*"_folder_path.txt")
    end

    return folderpath
end

function reshape_to_vec(file, histogram_resolution)
    return sum(reshape(file, (histogram_resolution, :)), dims = 2)[:,1]
end

function open_SDT_file(filepath)
    open(filepath, "r") do io
        seek(io, 14)
        header = read!(io, Vector{UInt8}(undef, 12))

        seek(io, header[12]*0x100 + header[11] + 82)
        infos = read!(io, Vector{UInt8}(undef, 2))
        histogram_resolution = infos[2]*0x100 + infos[1]
    
        seek(io, header[12]*0x100 + header[11] + 215)
        infos_2 = read!(io, Vector{UInt8}(undef, 4))
        time = reinterpret(Float32, infos_2)[1]

        if (header[8]*0x100 + header[7]) in (512, 128)
            seek(io, header[2]*0x100 + header[1] + 22)
            vector = read(io)  # lit tout le reste du fichier
            # Conversion rapide en UInt16 sans boucle
            n = div(length(vector), 2)
            new_vector = reinterpret(UInt16, vector[1:2*n])
            return new_vector, histogram_resolution, time
        else
            seek(io, header[2]*0x100 + header[1] + 2)
            file_info = read!(io, Vector{UInt8}(undef, 20))
            shift_1 = file_info[8]*0x1000000 + file_info[7]*0x10000 + file_info[6]*0x100 + file_info[5] -
                      (file_info[2]*0x100 + file_info[1])
            shift_2 = file_info[20]*0x1000000 + file_info[19]*0x10000 + file_info[18]*0x100 + file_info[17]

            buffer = IOBuffer(read(io, shift_1))
            zip_file = first(ZipFile.Reader(buffer).files)

            file_data = Vector{UInt8}(undef, shift_2)
            read!(zip_file, file_data)

            vector = reshape_to_vec(reinterpret(UInt16, file_data), histogram_resolution)
            return vector, histogram_resolution, time
        end
    end
end

function add_to_monitor(msg::String)
    currentText = monitor.text[]
    lines = split(currentText, "\n")
    push!(lines, split(msg, "\r\n")[1])
    lines = filter(!isempty, lines)
    MAX_LINES = 7
    if length(lines) > MAX_LINES
        monitor.text[] = join(lines[end - (MAX_LINES - 1):end], "\n")
    else
        monitor.text[] = join(lines, "\n")
    end
end

function connect_to_port(portName)
    try
        # ser = LibSerialPort.open("/dev/tty." * portName, 115200)
        ser = LibSerialPort.open(portName, 115200)
        LibSerialPort.set_read_timeout(ser, 3)
        sleep(0.1)
        add_to_monitor("Connection successful")
        println("Connection successful")
        connected[] = true
        button_connect.label[] = "Disconnect"
        return ser
    catch e
        add_to_monitor("Connection error")
        println("Connection error")
    end
end

function send_command(serialConn, commandStr)
    # println("sending signal")
    # println(commandStr)
    write(serialConn, commandStr)
    # sleep(0.01)
    # response = strip(read(serialConn, String))
    # println(response)
    # add_to_monitor("$response")
    # return response
end

function clear(serial)
    send_command(serial, "C Z\nC R\n")
end

function ewma!(smoothed, raw, a)
    smoothed[1] = raw[1]

    for i in eachindex(raw)[2:end]
        smoothed[i] = a*raw[i] + (1-a)*smoothed[i-1]
    end
end

function exponential_model(x, p)
    y = p[1] .* exp.(-p[2] .* (x .- p[3])) .+ p[4]
    y .= ifelse.(x .< 0, p[1], y)
    return y
end

function loess_smooth(x,
                      y;
                      f::Real = 0.3)
    @assert length(x) == length(y) "x and y must have the same length"
    @assert issorted(x)        "x must be sorted in ascending order"

    model = lowess_model(x, y, f)

    return model(x)
end

function update_smooth!(smooth, times, lifetime, smoothing)
    N = length(times)
    
    if N > 100
        w = (N-100):N

        tview = @view times[w]
        lview = @view lifetime[w]
        sview = @view smooth[w]

        fit = lowess(tview, lview, smoothing)

        blend = range(0.0, 1.0, length=length(fit))

        @. sview = blend * fit + (1 - blend) * sview
    else
        result = lowess(times, lifetime, smoothing)
        @. smooth = result
    end

    return smooth
end

function analyse_time(times)
    try
        println("Analysing")

        # deltas = zeros(length(times)-1)

        # for (i, time) in enumerate(times[2:end])
        #     deltas[i] = (time - times[i])*1000
        # end

        deltas = 1000*times

        stdev = std(deltas)
        meann = mean(deltas)
        mediann = median(deltas)

        for (i, delta) in enumerate(deltas)
            if delta > meann + 3*stdev
                deltas[i] = meann
            end
        end

        stdev = std(deltas)
        meann = mean(deltas)
        mediann = median(deltas)

        hist!(axis_counts, deltas, color=Makie.wong_colors()[1], bins=50)
        vlines!(axis_counts, meann, color=Makie.wong_colors()[2])
        vlines!(axis_counts, mediann, color=Makie.wong_colors()[3])
        autolimits!(axis_counts)

        println("Analyse done")
        println("Mean : $(round(meann, digits = 2))")
        println("Median : $(round(mediann, digits = 2))")
        println("Stdev : $(round(stdev, digits = 2))")
        println("Max : $(round(maximum(deltas), digits = 2))")
        println("Min : $(round(minimum(deltas), digits = 2))")
    catch e
        println(e)
    end
end

function start_real_time()
    GC.gc()
    println(folderpath)

    hist = Observable(zeros(256))
    fit = Observable(zeros(256))
    hist_time = Observable(zeros(256))
    
    lifetime = Observable(Float64[])
    smooth = Observable(Float64[])
    counts = Observable(Float64[])
    global signal = Observable(Float64[])
    times = Observable(Float64[])
    binning_list = []
    smooth_list = []
    prop_list = []
    integ_list = []
    deriv_list = []
    setpoint_list = []
    
    vectors = zeros(100, 256)
    n_vectors = size(vectors, 1)
    histogram_resolution = 256
    
    barplot!(axis_hist, hist_time, hist, gap=0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_life, times, smooth, color=Makie.wong_colors()[3])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    lifetimes = parse(Int, lifetimes_textbox.stored_string[])

    if lifetimes == 1
        params = [5.0, -0.25, 1.0e-5]
    elseif lifetimes == 2
        params = [6.0, 0.75, 1.0, 0.0, 1.0e-5]
    else
        params = [6.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    end

    error_I = 0
    error_P = 0
    old_error = 0
    display_period = 1/10

    i = 1
    n = 0
    t_0 = now()
    previous_files = readdir(folderpath)
    dt = 0

    println("Starting real time")
    add_to_monitor("Starting real time")

    while running[]
        try
            current_files = readdir(folderpath)
            nb_files = length(current_files)
            new_files = setdiff(current_files, previous_files)
            sleep(0.0001)
            
            if !isempty(new_files)
                if i == 1
                    println("Files detected")
                    add_to_monitor("Files detected")
                    t_0 = now()
                end

                t_now = (now() - t_0)/Millisecond(1)/1000

                time_range = parse(Int, range_textbox.stored_string[])
                averaging = parse(Int, averaging_textbox.stored_string[])
                # smoothing = exp(-parse(Int, smoothing_textbox.stored_string[])/1.5)
                smoothing = (parse(Int, smoothing_textbox.stored_string[])/10)^2
                prop_val = parse(Int, proportionnal_textbox.stored_string[])*10
                integ_val = parse(Int, integral_textbox.stored_string[])/10
                deriv_val = parse(Int, derivative_textbox.stored_string[])*1000
                setpoint_val = parse(Float64, setpoint_textbox.stored_string[])
                inverted = checkbox_inverted.active[]

                push!(binning_list, averaging)
                push!(smooth_list, smoothing)
                push!(prop_list, prop_val)
                push!(integ_list, integ_val)
                push!(deriv_list, deriv_val)
                push!(setpoint_list, setpoint_val)

                vector, histogram_resolution, duration = open_SDT_file(joinpath(folderpath, new_files[end]))

                # if i == 1
                #     for n in 1:100
                #         vectors[n, 1:histogram_resolution] = vector
                #     end
                # else
                #     for n in 1:99
                #         vectors[n, 1:histogram_resolution] = vectors[n+1, 1:histogram_resolution]
                #     end
                #     vectors[end, 1:histogram_resolution] = vector
                # end
                
                if i == 1
                    vectors[:, 1:histogram_resolution] .= reshape(vector, (1, :))
                else
                    vectors[1:end-1, 1:histogram_resolution] .= vectors[2:end, 1:histogram_resolution]
                    vectors[end, 1:histogram_resolution] .= vector
                end

                idx_start = max(1, n_vectors - averaging + 1)
                final_vector = @views sum(vectors[idx_start:end, 1:histogram_resolution]; dims=1)[1, :] / averaging

                params_raw, data = vec_to_lifetime(final_vector, guess=params, histogram_resolution=histogram_resolution)

                nb_photons = sum(final_vector)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)*nb_photons

                push!(times[], t_now)
                push!(counts[], nb_photons)
                push!(lifetime[], params[1])
                push!(smooth[], params[1])

                if i > 3 && active[]
                    if i > 100
                        smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                        smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                    else
                        smooth[] = lowess(times[], lifetime[], smoothing)
                    end
                end

                # error = setpoint_val - params[1]

                dt = duration

                if active[] 
                    error_P = smooth[][end] - setpoint_val

                    if i > 100
                        error_I += error_P * duration
                        error_D = (smooth[][end] - smooth[][end-1])/(1*duration)
                    else
                        error_I = 0
                        error_D = 0
                    end

                    old_error = error_P

                    if inverted
                        signal_value = clamp((prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    else
                        signal_value = clamp(-(prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    end
                else
                    signal_value = 0
                    error_I = 0
                end

                # try
                #     dt = t_now - times[][end]
                # catch e
                #     dt = t_now
                #     # println(e)
                # end

                # if reset_I[]
                #     error_I = 0
                #     reset_I[] = false
                # end

                # error_I += error * dt

                # if i == 1
                #     error_I = 0
                # end

                # if inverted
                #     signal_value = clamp(-(prop_val * error + integ_val * error_I), 0, 100)
                # else
                #     signal_value = clamp(prop_val * error + integ_val * error_I, 0, 100)
                # end

                push!(signal[], signal_value)

                instant_freq =  i / (t_now + 1e-6)
                nb_i = ceil(instant_freq*display_period)

                # if mod1(i, nb_i) == nb_i
                # if mod1(n, 3) == 3
                #     label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                #     label_number.text = "Files : $i"
                # end

                # x_min = max(0, times[][end] - time_range + 0.1)
                # x_max = max(time_range, times[][end] + 0.1)

                # notify(hist_time)
                # notify(times)

                # autolimits!(axis_hist)
                # autolimits!(axis_counts)
                # autolimits!(axis_life)

                # xlims!(axis_counts, x_min, x_max)
                # xlims!(axis_life, x_min, x_max)
                # xlims!(axis_signal, x_min, x_max)
                # # end

                if mod1(i, nb_i) == nb_i
                    if mod1(n, 3) == 3
                        label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                        label_number.text = "Files : $i/$nb_files"
                    end

                    if i > 3 && !active[]
                        if i > 100
                            smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                            smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                        else
                            smooth[] = lowess(times[], lifetime[], smoothing)
                        end
                    end

                    x_min = max(0, times[][end] - time_range + 0.1)
                    x_max = max(time_range, times[][end] + 0.1)

                    notify(hist_time)
                    notify(times)

                    autolimits!(axis_hist)
                    autolimits!(axis_counts)
                    autolimits!(axis_life)

                    xlims!(axis_counts, x_min, x_max)
                    xlims!(axis_life, x_min, x_max)
                    xlims!(axis_signal, x_min, x_max)

                    n += 1
                end

                previous_files = current_files
                i += 1
            end

            while pause[]
                sleep(0.1)
            end
        catch e
            println(e)
        end
                
        if !running[]
            add_to_monitor("Aborted")
            break
        end
    end

    filepath = save_file()

    if !isempty(filepath)
        if filepath[end-4:end] != ".csv"
            filepath *= ".csv"
        end

        println("Saving as $(split(filepath, "/")[end])")
        add_to_monitor("Saving as $(split(filepath, "/")[end])")

        df = DataFrame(
            Time_s = times[],
            Number_of_photons = counts[],
            Lifetime_ns = lifetime[],
            Smooth = smooth[],          
            Signal = signal[],
            Binning = binning_list,
            Smoothing = smooth_list,
            Proportionnal = prop_list,
            Integral = integ_list,
            Derivative = deriv_list,
            Setpoint = setpoint_list
        )
        
        open(filepath, "w") do io
            println(io, "Timestamp: ", now())
            CSV.write(io, df)
        end
    end

    filename = string(now())*"2.csv"

    open(filename, "w") do io
        writedlm(io, [string(now())])
        writedlm(io, [""])
        writedlm(io, ["Time [s]" "Number of photons" "Lifetime [ns]" "Smooth" "Signal [%]" "Binning" "Smoothing" "Proportionnal" "Integral" "Derivative" "Setpoint"], ",")
        writedlm(io, [times[] counts[] lifetime[] smooth[] signal[] binning_list smooth_list prop_list integ_list deriv_list setpoint_list], ",")
    end
end

function start_protocol_real_time()
    GC.gc()
    println(folderpath)

    hist = Observable(zeros(256))
    fit = Observable(zeros(256))
    hist_time = Observable(zeros(256))
    
    lifetime = Observable(Float64[])
    smooth = Observable(Float64[])
    counts = Observable(Float64[])
    global signal = Observable(Float64[])
    times = Observable(Float64[])
    binning_list = []
    smooth_list = []
    prop_list = []
    integ_list = []
    deriv_list = []
    setpoint_list = []
    
    vectors = zeros(100, 256)
    n_vectors = size(vectors, 1)
    histogram_resolution = 256
    
    barplot!(axis_hist, hist_time, hist, gap=0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_life, times, smooth, color=Makie.wong_colors()[3])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    lifetimes = parse(Int, lifetimes_textbox.stored_string[])

    if lifetimes == 1
        params = [5.0, -0.25, 1.0e-5]
    elseif lifetimes == 2
        params = [6.0, 0.75, 1.0, 0.0, 1.0e-5]
    else
        params = [6.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    end

    error_I = 0
    error_P = 0
    old_error = 0
    display_period = 1/10
    old_setpoint = 5
    protocol_copy = deepcopy(protocol)
    local_active = false

    i = 1
    n = 0
    t_0 = now()
    previous_files = readdir(folderpath)
    dt = 0

    println("Starting real time")
    add_to_monitor("Starting real time")

    while running[]
        try
            current_files = readdir(folderpath)
            nb_files = length(current_files)
            new_files = setdiff(current_files, previous_files)
            sleep(0.0001)
            
            if !isempty(new_files)
                if i == 1
                    println("Files detected")
                    add_to_monitor("Files detected")
                    t_0 = now()
                end

                t_now = (now() - t_0)/Millisecond(1)/1000

                time_range = parse(Int, range_textbox.stored_string[])
                averaging = parse(Int, averaging_textbox.stored_string[])
                # smoothing = exp(-parse(Int, smoothing_textbox.stored_string[])/1.5)
                smoothing = (parse(Int, smoothing_textbox.stored_string[])/10)^2
                prop_val = parse(Int, proportionnal_textbox.stored_string[])*10
                integ_val = parse(Int, integral_textbox.stored_string[])/10
                deriv_val = parse(Int, derivative_textbox.stored_string[])*1000
                # setpoint_val = parse(Float64, setpoint_textbox.stored_string[])
                inverted = checkbox_inverted.active[]

                setpoint_val = protocol_copy[round(Int, t_now + 1)]
                if isnan(setpoint_val)
                    local_active = false
                else
                    local_active = true
                end

                if (isnan(setpoint_val) && !isnan(old_setpoint)) || (!isnan(setpoint_val) && isnan(old_setpoint))
                    if setpoint_val != old_setpoint
                        println("Setpoint changed")
                        setpoint_textbox.stored_string[] = string(setpoint_val)
                        setpoint_textbox.displayed_string[] = string(setpoint_val)
                        old_setpoint = setpoint_val
                    end
                end

                push!(binning_list, averaging)
                push!(smooth_list, smoothing)
                push!(prop_list, prop_val)
                push!(integ_list, integ_val)
                push!(deriv_list, deriv_val)
                push!(setpoint_list, setpoint_val)

                vector, histogram_resolution, duration = open_SDT_file(joinpath(folderpath, new_files[end]))
                
                if i == 1
                    vectors[:, 1:histogram_resolution] .= reshape(vector, (1, :))
                else
                    vectors[1:end-1, 1:histogram_resolution] .= vectors[2:end, 1:histogram_resolution]
                    vectors[end, 1:histogram_resolution] .= vector
                end

                idx_start = max(1, n_vectors - averaging + 1)
                final_vector = @views sum(vectors[idx_start:end, 1:histogram_resolution]; dims=1)[1, :] / averaging

                params_raw, data = vec_to_lifetime(final_vector, guess=params, histogram_resolution=histogram_resolution)

                nb_photons = sum(final_vector)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)*nb_photons

                push!(times[], t_now)
                push!(counts[], nb_photons)
                push!(lifetime[], params[1])
                push!(smooth[], params[1])

                if i > 3 && active[]
                    if i > 100
                        smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                        smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                    else
                        smooth[] = lowess(times[], lifetime[], smoothing)
                    end
                end

                dt = duration

                if active[] && local_active
                    error_P = smooth[][end] - setpoint_val

                    if i > 100
                        error_I += error_P * duration
                        error_D = (smooth[][end] - smooth[][end-1])/(1*duration)
                    else
                        error_I = 0
                        error_D = 0
                    end

                    old_error = error_P

                    if inverted
                        signal_value = clamp((prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    else
                        signal_value = clamp(-(prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    end
                else
                    signal_value = 0
                    error_I = 0
                end

                push!(signal[], signal_value)

                instant_freq =  i / (t_now + 1e-6)
                nb_i = ceil(instant_freq*display_period)

                if mod1(i, nb_i) == nb_i
                    if mod1(n, 3) == 3
                        label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                        label_number.text = "Files : $i/$nb_files"
                    end

                    if i > 3 && !active[]
                        if i > 100
                            smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                            smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                        else
                            smooth[] = lowess(times[], lifetime[], smoothing)
                        end
                    end

                    x_min = max(0, times[][end] - time_range + 0.1)
                    x_max = max(time_range, times[][end] + 0.1)

                    notify(hist_time)
                    notify(times)

                    autolimits!(axis_hist)
                    autolimits!(axis_counts)
                    autolimits!(axis_life)

                    xlims!(axis_counts, x_min, x_max)
                    xlims!(axis_life, x_min, x_max)
                    xlims!(axis_signal, x_min, x_max)

                    n += 1
                end

                previous_files = current_files
                i += 1
            end

            while pause[]
                sleep(0.1)
            end
        catch e
            println(typeof(e))
            add_to_monitor("Protocol finished")
            running[] = false
        end
                
        if !running[]
            add_to_monitor("Aborted")
            break
        end
    end

    filepath = save_file()

    if !isempty(filepath)
        if filepath[end-4:end] != ".csv"
            filepath *= ".csv"
        end

        println("Saving as $(split(filepath, "/")[end])")
        add_to_monitor("Saving as $(split(filepath, "/")[end])")

        df = DataFrame(
            Time_s = times[],
            Number_of_photons = counts[],
            Lifetime_ns = lifetime[],
            Smooth = smooth[],          
            Signal = signal[],
            Binning = binning_list,
            Smoothing = smooth_list,
            Proportionnal = prop_list,
            Integral = integ_list,
            Derivative = deriv_list,
            Setpoint = setpoint_list
        )
        
        open(filepath, "w") do io
            println(io, "Timestamp: ", now())
            CSV.write(io, df)
        end
    end

    filename = string(now())*"2.csv"

    open(filename, "w") do io
        writedlm(io, [string(now())])
        writedlm(io, [""])
        writedlm(io, ["Time [s]" "Number of photons" "Lifetime [ns]" "Smooth" "Signal [%]" "Binning" "Smoothing" "Proportionnal" "Integral" "Derivative" "Setpoint"], ",")
        writedlm(io, [times[] counts[] lifetime[] smooth[] signal[] binning_list smooth_list prop_list integ_list deriv_list setpoint_list], ",")
    end
end

function start_playback()
    GC.gc()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)
    nb_files = length(files)

    hist = Observable(zeros(256))
    fit = Observable(zeros(256))
    hist_time = Observable(zeros(256))

    lifetime = Observable(Float64[])
    smooth = Observable(Float64[])
    counts = Observable(Float64[])
    global signal = Observable(Float64[])
    times = Observable(Float64[])
    test_times = []

    vectors = zeros(100, 256)
    n_vectors = size(vectors, 1)
    histogram_resolution = 256

    barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_life, times, smooth, color=Makie.wong_colors()[3])
    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    lifetimes = parse(Int, lifetimes_textbox.stored_string[])

    if lifetimes == 1
        params = [5.0, 0, 5.0e-5]
    elseif lifetimes == 2
        params = [6.0, 0.75, 1.0, 0.0, 1.0e-5]
    else
        params = [6.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    end

    error_I = 0
    error_P = 0
    old_error = 0
    display_period = 1/10

    i = 1
    n = 1
    t_0 = now()
    cum_time = 0

    println("Starting playback")
    add_to_monitor("Starting playback")

    sleep(0.1)

    while running[]
        for file in files
            try 
                t = time()

                time_range = parse(Int, range_textbox.stored_string[])
                averaging = parse(Int, averaging_textbox.stored_string[])
                # smoothing = exp(-parse(Int, smoothing_textbox.stored_string[])/1.5)
                smoothing = (parse(Int, smoothing_textbox.stored_string[])/10)^2
                prop_val = parse(Int, proportionnal_textbox.stored_string[])*10
                integ_val = parse(Int, integral_textbox.stored_string[])/10
                deriv_val = parse(Int, derivative_textbox.stored_string[])*1000
                setpoint_val = parse(Float64, setpoint_textbox.stored_string[])
                inverted = checkbox_inverted.active[]

                # if smoothing == 1 smoothing = 0.95 end
                
                if !running[]
                    add_to_monitor("Aborted")
                    break
                end
                
                vector, histogram_resolution, duration = open_SDT_file(joinpath(folderpath, file))
                cum_time += duration
                
                if i == 1
                    vectors[:, 1:histogram_resolution] .= reshape(vector, (1, :))
                else
                    vectors[1:end-1, 1:histogram_resolution] .= vectors[2:end, 1:histogram_resolution]
                    vectors[end, 1:histogram_resolution] .= vector
                end

                idx_start = max(1, n_vectors - averaging + 1)
                final_vector = @views sum(vectors[idx_start:end, 1:histogram_resolution]; dims=1)[1, :] / averaging

                params_raw, data = vec_to_lifetime(final_vector, guess=params, histogram_resolution=histogram_resolution)

                nb_photons = sum(final_vector)

                # if nb_photons > 4e4
                #     println(i)
                # end

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)*nb_photons

                push!(times[], cum_time)
                push!(counts[], nb_photons)
                push!(lifetime[], params[1])
                push!(smooth[], params[1])
                # ewma!(smooth[], lifetime[], smoothing)
                # ewma!(smooth[], smooth[], smoothing)

                if i > 3 && active[]
                    # if i > 100
                    #     smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                    #     smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                    # else
                    #     smooth[] = lowess(times[], lifetime[], smoothing)
                    # end

                    update_smooth!(smooth[], times[], lifetime[], smoothing)
                end

                if i == 1
                    t_0 = now()
                end

                delta_t = (now() - t_0) / Microsecond(1) / 1e6

                if active[] 
                    error_P = smooth[][end] - setpoint_val

                    if i > 100
                        error_I += error_P * duration
                        error_D = (smooth[][end] - smooth[][end-50])/(50*duration)
                    else
                        error_I = 0
                        error_D = 0
                    end

                    old_error = error_P

                    if inverted
                        signal_value = clamp((prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    else
                        signal_value = clamp(-(prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    end
                else
                    signal_value = 0
                    error_I = 0
                end

                push!(signal[], signal_value)

                instant_freq =  i / (delta_t + 1e-6)
                nb_i = ceil(instant_freq*display_period)

                if mod1(i, nb_i) == nb_i
                    if i > 3 && !active[]
                        # if i > 100
                        #     smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                        #     smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                        # else
                        #     smooth[] = lowess(times[], lifetime[], smoothing)
                        # end

                        update_smooth!(smooth[], times[], lifetime[], smoothing)
                    end

                    x_min = max(0, times[][end] - time_range + 0.1)
                    x_max = max(time_range, times[][end] + 0.1)

                    notify(hist_time)
                    notify(times)

                    autolimits!(axis_hist)
                    autolimits!(axis_counts)
                    autolimits!(axis_life)

                    xlims!(axis_counts, x_min, x_max)
                    xlims!(axis_life, x_min, x_max)
                    xlims!(axis_signal, x_min, x_max)

                    if mod1(n, 3) == 3
                        label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                        label_number.text = "Files : $i/$nb_files"

                        push!(test_times, time() - t)
                    end

                    n += 1
                end

                i += 1
                sleep(0.0001)

                while pause[]
                    sleep(0.05)
                end
            catch e
                println(e)
            end
        end

        sleep(0.0001)
    end
    # analyse_time(test_times)
    # lines!(axis_signal, test_times[2:end], color=Makie.wong_colors()[1])
    # autolimits!(axis_signal)
end

function start_protocol_playback()
    GC.gc()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)
    nb_files = length(files)

    hist = Observable(zeros(256))
    fit = Observable(zeros(256))
    hist_time = Observable(zeros(256))

    lifetime = Observable(Float64[])
    smooth = Observable(Float64[])
    counts = Observable(Float64[])
    global signal = Observable(Float64[])
    times = Observable(Float64[])
    test_times = []

    vectors = zeros(100, 256)
    n_vectors = size(vectors, 1)
    histogram_resolution = 256

    barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_life, times, smooth, color=Makie.wong_colors()[3])
    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    lifetimes = parse(Int, lifetimes_textbox.stored_string[])

    if lifetimes == 1
        params = [5.0, -0.25, 1.0e-5]
    elseif lifetimes == 2
        params = [6.0, 0.75, 1.0, 0.0, 1.0e-5]
    else
        params = [6.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    end

    error_I = 0
    error_P = 0
    old_error = 0
    display_period = 1/10
    old_setpoint = NaN
    local_active = false

    i = 1
    n = 1
    t_0 = now()
    cum_time = 0

    setpoint_textbox.stored_string[] = string(NaN)
    setpoint_textbox.displayed_string[] = string(NaN)

    println("Starting playback")
    add_to_monitor("Starting playback")
    # println(protocol)

    protocol_copy = deepcopy(protocol)

    sleep(0.1)

    while running[]
        for file in files
            try 
                t = time()

                time_range = parse(Int, range_textbox.stored_string[])
                averaging = parse(Int, averaging_textbox.stored_string[])
                # smoothing = exp(-parse(Int, smoothing_textbox.stored_string[])/1.5)
                smoothing = (parse(Int, smoothing_textbox.stored_string[])/10)^2
                prop_val = parse(Int, proportionnal_textbox.stored_string[])*10
                integ_val = parse(Int, integral_textbox.stored_string[])/10
                deriv_val = parse(Int, derivative_textbox.stored_string[])*1000
                # setpoint_val = parse(Float64, setpoint_textbox.stored_string[])
                inverted = checkbox_inverted.active[]

                # if smoothing == 1 smoothing = 0.95 end
                
                if !running[]
                    add_to_monitor("Aborted")
                    break
                end
                
                vector, histogram_resolution, duration = open_SDT_file(joinpath(folderpath, file))
                cum_time += duration

                setpoint_val = protocol_copy[round(Int, cum_time + 1)]
                if isnan(setpoint_val)
                    local_active = false
                else
                    local_active = true
                end

                if (isnan(setpoint_val) && !isnan(old_setpoint)) || (!isnan(setpoint_val) && isnan(old_setpoint))
                    if setpoint_val != old_setpoint
                        println("Setpoint changed")
                        setpoint_textbox.stored_string[] = string(setpoint_val)
                        setpoint_textbox.displayed_string[] = string(setpoint_val)
                        old_setpoint = setpoint_val
                    end
                end

                if i == 1
                    vectors[:, 1:histogram_resolution] .= reshape(vector, (1, :))
                else
                    vectors[1:end-1, 1:histogram_resolution] .= vectors[2:end, 1:histogram_resolution]
                    vectors[end, 1:histogram_resolution] .= vector
                end

                idx_start = max(1, n_vectors - averaging + 1)
                final_vector = @views sum(vectors[idx_start:end, 1:histogram_resolution]; dims=1)[1, :] / averaging

                params_raw, data = vec_to_lifetime(final_vector, guess=params, histogram_resolution=histogram_resolution)

                nb_photons = sum(final_vector)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)*nb_photons

                push!(times[], cum_time)
                push!(counts[], nb_photons)
                push!(lifetime[], params[1])
                push!(smooth[], params[1])
                # ewma!(smooth[], lifetime[], smoothing)
                # ewma!(smooth[], smooth[], smoothing)

                if i > 3 && active[]
                    # if i > 100
                    #     smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                    #     smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                    # else
                    #     smooth[] = lowess(times[], lifetime[], smoothing)
                    # end

                    update_smooth!(smooth[], times[], lifetime[], smoothing)
                end

                if i == 1
                    t_0 = now()
                end

                delta_t = (now() - t_0) / Microsecond(1) / 1e6

                if active[] && local_active
                    error_P = smooth[][end] - setpoint_val

                    if i > 100
                        error_I += error_P * duration
                        error_D = (smooth[][end] - smooth[][end-50])/(50*duration)
                    else
                        error_I = 0
                        error_D = 0
                    end

                    old_error = error_P

                    if inverted
                        signal_value = clamp((prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    else
                        signal_value = clamp(-(prop_val * error_P + integ_val * error_I + deriv_val * error_D) + 0, 0, 100)
                    end
                else
                    signal_value = 0
                    error_I = 0
                end

                push!(signal[], signal_value)

                instant_freq =  i / (delta_t + 1e-6)
                nb_i = ceil(instant_freq*display_period)

                if mod1(i, nb_i) == nb_i
                    if i > 3 && !active[]
                        # if i > 100
                        #     smooth_fit = lowess(times[][end-100:end], lifetime[][end-100:end], smoothing)
                        #     smooth[][end-100:end] = collect(0:0.01:1) .* smooth_fit .+ collect(1:-0.01:0) .* smooth[][end-100:end]
                        # else
                        #     smooth[] = lowess(times[], lifetime[], smoothing)
                        # end

                        update_smooth!(smooth[], times[], lifetime[], smoothing)
                    end

                    x_min = max(0, times[][end] - time_range + 0.1)
                    x_max = max(time_range, times[][end] + 0.1)

                    notify(hist_time)
                    notify(times)

                    autolimits!(axis_hist)
                    autolimits!(axis_counts)
                    autolimits!(axis_life)

                    xlims!(axis_counts, x_min, x_max)
                    xlims!(axis_life, x_min, x_max)
                    xlims!(axis_signal, x_min, x_max)

                    if mod1(n, 3) == 3
                        label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                        label_number.text = "Files : $i/$nb_files"

                        push!(test_times, time() - t)
                    end

                    n += 1
                end

                i += 1
                sleep(0.0001)

                while pause[]
                    sleep(0.05)
                end
            catch e
                println(typeof(e))
                add_to_monitor("Protocol finished")
                running[] = false
            end
        end

        sleep(0.0001)
    end
    # analyse_time(test_times)
    # lines!(axis_signal, test_times[2:end], color=Makie.wong_colors()[1])
    # autolimits!(axis_signal)
end

function start_save()
    GC.gc()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)
    nb_files = length(files)
    checkpoints = round.(Int, nb_files .* (0:1:100) ./ 100)

    lifetime = Float64[]
    lifetime_i = Float64[]
    lifetime_m = Float64[]
    lifetime_2 = Float64[]
    counts = Float64[]
    times = Float64[]

    lifetimes = parse(Int, lifetimes_textbox.stored_string[])

    if lifetimes == 1
        params = [5.0, 0, 5.0e-5]
    elseif lifetimes == 2
        params = [6.0, 0.75, 1.0, 0.0, 1.0e-5]
    # else
    #     params = [6.0, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    end

    i = 1
    cum_time = 0

    println("Starting analysis")
    add_to_monitor("Starting analysis")

    histogram_resolution = 256

    box_progress.width = 10

    beginning = time()

    sleep(0.1)

    for file in files
        println(file)
        try
            if !running[]
                add_to_monitor("Aborted")
                break
            end

            vector, histogram_resolution, delta_time = open_SDT_file(joinpath(folderpath, file))
            cum_time += delta_time

            params_raw, data = vec_to_lifetime(convert(Vector{Float64}, vector), guess=params, histogram_resolution=histogram_resolution)

            if !isnan(params_raw[1])
                params = params_raw
            end

            push!(lifetime, params[1])
            push!(counts, sum(vector))
            push!(times, cum_time)

            if lifetimes == 2
                push!(lifetime_i, params[1] * params[2] + params[3] * (1 - params[2]))
                push!(lifetime_m, (params[1]^2 * params[2] + params[3]^2 * (1 - params[2]))/(params[1] * params[2] + params[3] * (1 - params[2])))
                push!(lifetime_2, params[3])
            end

            i += 1

            sleep(0.0001)

            if i in checkpoints
                box_progress.width = 10 + 1290*i/nb_files
                label_progress.text = "$(round(Int, i/nb_files*100)) %"
            end

            while pause[]
                sleep(0.05)
            end
        catch e
            println(e)
        end
    end

    total_time = time() - beginning
    println("Time : $(round(total_time, sigdigits=3)) s")

    label_frequency.text = "Frequency : $(round(i/total_time, sigdigits=3)) Hz"
    label_number.text = "Files : $i"

    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])

    if lifetimes == 2
        lines!(axis_life, times, lifetime_i, color=Makie.wong_colors()[2])
        lines!(axis_life, times, lifetime_m, color=Makie.wong_colors()[3])
        lines!(axis_life, times, lifetime_2, color=Makie.wong_colors()[4])
    end

    autolimits!(axis_counts)
    autolimits!(axis_life)

    filepath = save_file()

    if !isempty(filepath)
        if filepath[end-4:end] != ".csv"
            filepath *= ".csv"
        end

        println("Saving as $(split(filepath, "/")[end])")
        add_to_monitor("Saving as $(split(filepath, "/")[end])")

        df = DataFrame(
            Time_s = times,
            Number_of_photons = counts,
            Lifetime_ns = lifetime,
            Lifetime_i_ns = lifetime_i,
            Lifetime_m_ns = lifetime_m,
            Lifetime_2_ns = lifetime_2
        )
        
        open(filepath, "w") do io
            println(io, "Timestamp: ", now())
            CSV.write(io, df)
        end
    end

    # filename = string(now())*".csv"

    # df = DataFrame(
    #     Time_s = times,
    #     Number_of_photons = counts,
    #     Lifetime_ns = lifetime
    # )
    
    # open(filename, "w") do io
    #     println(io, "Timestamp: ", now())
    #     CSV.write(io, df)
    # end

    add_to_monitor("Done")
    
    running[] = false
    pause[] = false

    button_start.label = "START"
    button_stop.label = "CLEAR"
end

function start_analysis()
    println("Starting analysis")

    filepath = pick_file()

    if isempty(filepath)
        println("No file selected")
        return
    end

    println("Opening file : ", filepath)

    df = CSV.read(filepath, DataFrame)

    times    = df.Time_s
    counts   = df.Number_of_photons
    lifetime = df.Lifetime_ns
    signal   = df.Signal
    setpoint = df.Setpoint

    # setpoint[signal .== 0] .= NaN

    starts = Int[]
    ends   = Int[]
    in_region = false

    for i in eachindex(setpoint)
        if !isnan(setpoint[i]) && !in_region
            push!(starts, i)
            in_region = true
        elseif (isnan(setpoint[i]) || i == lastindex(setpoint)) && in_region
            push!(ends, isnan(setpoint[i]) ? i-1 : i)
            in_region = false
        end
    end

    if in_region
        push!(ends, lastindex(a))
    end

    times_i = times[starts]
    times_f = times[ends]

    # println(starts, ends)
    # println(times_i, times_f)

    vspan!(axis_life, times_i, times_f, color=Makie.wong_colors()[2], alpha=0.2)

    # exponential_model(x, p) = p[1] .* exp.(-p[2] .* (x .- p[3])) .+ p[4]
    fits = []

    for i in eachindex(starts)
        i0 = starts[i]
        i1 = ends[i]

        xseg = collect(times[i0:i1])
        yseg = lifetime[i0:i1]

        p0 = [
            maximum(yseg) - minimum(yseg),
            0.1,
            xseg[1],
            minimum(yseg),
        ]

        fit = curve_fit(exponential_model, xseg, yseg, p0)
        push!(fits, fit)
        # println(fit.param)
        println("Descente $i")
        # println("Tau = $(round(fit.param[2], digits=3)) s")
        println("Rise time = $(round(Int, 1/fit.param[2])) s")
        # println("x offset = $(round(fit.param[3] - times_i[i], digits=1)) s")
        println("y offset = $(round(fit.param[4], digits=1)) ns\n")

        i0 = ends[i]
        try
            i1 = starts[i+1]
        catch
            # println("end")
            i1 = lastindex(times)
        end

        xseg = collect(times[i0:i1])
        yseg = lifetime[i0:i1]

        p0 = [
            minimum(yseg) - maximum(yseg),
            0.1,
            xseg[1],
            maximum(yseg),
        ]

        fit = curve_fit(exponential_model, xseg, yseg, p0)
        push!(fits, fit)
        # println(fit.param)
        println("Montée $i")
        # println("Tau = $(round(fit.param[2], digits=3)) s")
        println("Rise time = $(round(Int, 1/fit.param[2])) s")
        # println("x offset = $(round(fit.param[3] - times_f[i], digits=1)) s")
        println("y offset = $(round(fit.param[4], digits=1)) ns\n")
    end

    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_life, times, setpoint, color=Makie.wong_colors()[2])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    for i in eachindex(starts)
        i0 = starts[i]
        i1 = ends[i]

        xseg = collect(times[i0:i1])

        xfit = range(xseg[1], xseg[end], length=200)
        yfit = exponential_model(xfit, coef(fits[2*i-1]))

        lines!(axis_life, xfit, yfit, color=:red)

        i0 = ends[i]
        try
            i1 = starts[i+1]
        catch
            i1 = lastindex(times)
        end

        xseg = collect(times[i0:i1])

        xfit = range(xseg[1], xseg[end], length=200)
        yfit = exponential_model(xfit, coef(fits[2*i]))

        lines!(axis_life, xfit, yfit, color=:red)
    end

    autolimits!(axis_counts)
    autolimits!(axis_life)

    add_to_monitor("Done")
    
    running[] = false
    pause[] = false

    button_start.label = "START"
    button_stop.label = "CLEAR"
end

function start_signal(frequency=100, refresh=1, min=3)
    println("sending signals")

    t_0 = now()
    delta_t = refresh*1000/frequency
    min_val = min*frequency/10
    max_val = 100 - min_val

    i = 0

    while active[]
        if running[] && Dates.value(now()) - Dates.value(t_0) > 1000*i/parse(Float64, textbox_refresh.stored_string[])
            try
                signal_val = signal[][end]
                # signal_val = 50
                # println(signal_val)

                # command = "C 1\n"
                # command *= "S 1 S 1\n"
                # command *= "A 1 DP 2 $frequency $(signal_val + 0.4*delay*frequency)\n"
                # command *= "A 1 WT $(delay*1000)\n"
                # command *= "A 1 DP 1 $frequency $signal_val\n"
                # command *= "A 0 TG 1\n"
                # println(command)

                # send_command(serial, command)

                if menu_output.selection[] == "Digital"
                    if signal_val < min_val/2
                        send_command(serial, "A 0 AO 2 0\n")

                    elseif signal_val > 100 - min_val/2
                        send_command(serial, "A 0 AO 2 5000\n")

                    else
                        signal_val = round(Int, clamp(signal_val, min_val, max_val))

                        send_command(serial, "A 0 AP 2 $(Int(frequency)) $signal_val 0 5000\n")
                        println("A 0 AP 2 $frequency $signal_val 0 5000\n")
                    end
                else
                    send_command(serial, "A 0 AO 1 $(0.5*signal_val*parse(Float64, textbox_min.stored_string[]))\n")
                    send_command(serial, "A 0 AO 2 $(0.5*signal_val*parse(Float64, textbox_fps.stored_string[]))\n")
                end

                i += 1
            catch e 
                println(e)
            end
        end

        sleep(0.000001)
    end
end
