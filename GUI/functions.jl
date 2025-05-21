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

include("vec_to_lifetime.jl")

function get_folder_path(str)
    if isfile("GUI/"*str*"_folder_path.txt")
        folderpath = open(f->read(f, String), "GUI/"*str*"_folder_path.txt")

        if isempty(folderpath)
            folderpath = pick_folder()

            open("GUI/"*str*"_folder_path.txt", "w") do io
                write(io, folderpath)
            end
        end
    else
        folderpath = pick_folder()

        open("GUI/"*str*"_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    end

    if isempty(folderpath)
        println("Folder does not exist")

        folderpath = get_folder_path("GUI/"*str)
    end

    return folderpath
end

function get_new_folder_path(str)
    folderpath = pick_folder()

    if !isempty(folderpath)
        open("GUI/"*str*"_folder_path.txt", "w") do io
            write(io, folderpath)
        end
    else
        folderpath = open(f->read(f, String), "GUI/"*str*"_folder_path.txt")
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
        ser = LibSerialPort.open("/dev/tty." * portName, 115200)
        # ser = LibSerialPort.open(portName, 115200)
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
    for i in 2:length(raw)
        smoothed[i] = a*raw[i] + (1-a)*smoothed[i-1]
    end
end

function analyse_time(times)
    try
        println("Analysing")

        # deltas = zeros(length(times)-1)

        # for (i, time) in enumerate(times[2:end])
        #     deltas[i] = (time - times[i])*1000
        # end

        deltas = times

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
    println(folderpath)
    println("Initializing")

    hist = Observable(zeros(256))
    fit = Observable(zeros(256))
    hist_time = Observable(zeros(256))
    
    lifetime = Observable(Float64[])
    smooth = Observable(Float64[])
    counts = Observable(Float64[])
    global signal = Observable(Float64[])
    times = Observable(Float64[])
    averaging_list = []
    prop_list = []
    integ_list = []
    setpoint_list = []
    
    vectors = zeros(100, 256)
    n_vectors = size(vectors, 1)
    
    barplot!(axis_hist, hist_time, hist, gap=0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_life, times, smooth, color=Makie.wong_colors()[3])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    params = [4.0, 2.0, 1.0e-4]
    error_I = 0

    i = 1
    n = 0
    t_0 = now()
    previous_files = readdir(folderpath)
    dt = 0
    display_period = 1/10

    println("Starting real time")
    add_to_monitor("Starting real time")

    while running[]
        try
            current_files = readdir(folderpath)
            new_files = setdiff(current_files, previous_files)
            sleep(0.001)
            
            if !isempty(new_files)
                if i == 1
                    t_0 = now()
                end

                time_range = parse(Int, range_textbox.stored_string[])
                averaging = parse(Int, averaging_textbox.stored_string[])
                smoothing = parse(Int, smoothing_textbox.stored_string[])/10
                prop_val = parse(Int, proportionnal_textbox.stored_string[])
                integ_val = parse(Int, integral_textbox.stored_string[])/100
                setpoint_val = parse(Float64, setpoint_textbox.stored_string[])
                inverted = checkbox_inverted.active[]

                push!(averaging_list, averaging)
                push!(prop_list, prop_val)
                push!(integ_list, integ_val)
                push!(setpoint_list, setpoint_val)

                vector, histogram_resolution, time = open_SDT_file(joinpath(folderpath, new_files[end]))

                if i == 1
                    for n in 1:100
                        vectors[n, 1:histogram_resolution] = vector
                    end
                else
                    for n in 1:99
                        vectors[n, 1:histogram_resolution] = vectors[n+1, 1:histogram_resolution]
                    end
                    vectors[end, 1:histogram_resolution] = vector
                end

                idx_start = max(1, n_vectors - averaging + 1)
                final_vector = @views sum(vectors[idx_start:end, 1:histogram_resolution]; dims=1)[1, :] / averaging

                params_raw, data = vec_to_lifetime(final_vector, guess=params, histogram_resolution=histogram_resolution)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)

                error = setpoint_val - params[1]

                t_now = (now() - t_0)/Millisecond(1)/1000

                try
                    dt = t_now - times[][end]
                catch e
                    dt = t_now
                    # println(e)
                end

                if reset_I[]
                    error_I = 0
                    reset_I[] = false
                end

                error_I += error * dt

                if i == 1
                    error_I = 0
                end

                if inverted
                    signal_value = clamp(-(prop_val * error + integ_val * error_I), 0, 100)
                else
                    signal_value = clamp(prop_val * error + integ_val * error_I, 0, 100)
                end

                push!(lifetime[], params[1])
                push!(counts[], sum(@view vectors[idx_start:end, :]) / averaging)
                push!(signal[], signal_value)
                push!(times[], t_now)

                instant_freq =  i / (t_now + 1e-6)
                # nb_i = ceil(instant_freq*display_period)

                # if mod1(i, nb_i) == nb_i
                if mod1(n, 3) == 3
                    label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                    label_number.text = "Files : $i"
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
                # end

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
            Signal_sent = signal[],
            Averaging = averaging_list,
            Proportionnal = prop_list,
            Integral = integ_list,
            Setpoint = setpoint_list
        )
        
        open(filepath, "w") do io
            println(io, "Timestamp: ", now())
            CSV.write(io, df)
        end
    end

    filename = string(now())*".csv"

    open(filename, "w") do io
        writedlm(io, [string(now())])
        writedlm(io, [""])
        writedlm(io, ["Time [s]" "Number of photons" "Lifetime [ns]" "Signal sent [%]" "Averaging" "Proportionnal" "Integral" "Setpoint"], ",")
        writedlm(io, [times[] counts[] lifetime[] signal[] averaging_list prop_list integ_list setpoint_list], ",")
    end
end

function start_playback()
    GC.gc()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)

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
                time_range = parse(Int, range_textbox.stored_string[])
                averaging = parse(Int, averaging_textbox.stored_string[])
                smoothing = exp(-parse(Int, smoothing_textbox.stored_string[])/1.5)
                prop_val = parse(Int, proportionnal_textbox.stored_string[])*10
                integ_val = parse(Int, integral_textbox.stored_string[])/100
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

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)*nb_photons

                push!(lifetime[], params[1])
                push!(smooth[], params[1])
                ewma!(smooth[], lifetime[], smoothing)
                ewma!(smooth[], smooth[], smoothing)

                if i == 1
                    t_0 = now()
                end

                delta_t = (now() - t_0) / Microsecond(1) / 1e6

                if active[] 
                    error_P = smooth[][end] - setpoint_val

                    if i < 4
                        error_I = 0
                        error_D = 0
                    else
                        error_I += error_P * duration
                        error_D = ((smooth[][end] - smooth[][end-1])/duration + (smooth[][end-1] - smooth[][end-2])/duration + (smooth[][end-2] - smooth[][end-3])/duration)/3
                    end

                    old_error = error_P

                    if inverted
                        signal_value = clamp((prop_val * error_P + integ_val * error_I + deriv_val * error_D), 0, 100)
                    else
                        signal_value = clamp(-(prop_val * error_P + integ_val * error_I + deriv_val * error_D), 0, 100)
                    end
                else
                    signal_value = 0
                    error_I = 0
                end

                push!(counts[], nb_photons)
                push!(signal[], signal_value)
                push!(times[], cum_time)

                instant_freq =  i / (delta_t + 1e-6)
                nb_i = ceil(instant_freq*display_period)

                if mod1(i, nb_i) == nb_i
                    if mod1(n, 3) == 3
                        label_frequency.text = "Frequency : $(round(instant_freq, sigdigits=2)) Hz"
                        label_number.text = "Files : $i"
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
end

function start_save()
    GC.gc()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)
    nb_files = length(files)
    checkpoints = round.(Int, nb_files .* (0:1:100) ./ 100)

    lifetime = Float64[]
    counts = Float64[]
    times = Float64[]

    params = [4.0, 2.0, 1.0e-4]

    i = 1
    cum_time = 0

    println("Starting analysis")
    add_to_monitor("Starting analysis")

    histogram_resolution = 256

    box_progress.width = 10

    beginning = time()

    sleep(0.1)

    for file in files
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
            Lifetime_ns = lifetime
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

    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])

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
        if running[] && Dates.value(now()) - Dates.value(t_0) > i*delta_t
            try
                signal_val = signal[][end]
                # println(signal_val)

                if menu_output.selection[] == "Digital"
                    if signal_val < min_val/2
                        send_command(serial, "A 0 DO 1 0\n")

                    elseif signal_val > 100 - min_val/2
                        send_command(serial, "A 0 DO 1 1\n")

                    else
                        signal_val = round(Int, clamp(signal_val, min_val, max_val))

                        send_command(serial, "A 0 DP 1 $frequency $signal_val\n")
                    end
                else
                    send_command(serial, "A 0 AO 1 $(signal_val*10)\n")
                end

                i += 1
            catch e 
                println(e)
            end
        end

        sleep(0.000001)
    end
end
