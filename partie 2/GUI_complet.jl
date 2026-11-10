using GLMakie
using PyCall
using FFTW
using Statistics
using NativeFileDialog
using ZipFile
using Observables
using Base.Threads
using Dates
using DelimitedFiles
using LibSerialPort
using DataStructures
using BenchmarkTools
using Optim
using LineSearches
using CSV
using DataFrames

dark_gray  = RGBf(0.5, 0.5, 0.5)
mid_gray   = RGBf(0.7, 0.7, 0.7)
light_gray = RGBf(0.85, 0.85, 0.85)
off_white  = RGBf(0.9, 0.9, 0.9)

try
    LibSerialPort.close(serial)
catch
end

isnotnan(x) = !isnan(x)
smaller_or_eq_zero(x) = x <= 0

function get_irf(; channel=1)
    """
    Loads an .sdt IRF file. If irf_filepath.txt exists, load the IRF from this filepath, otherwise, open a file dialog create irf_filepath.txt
    """
    if isfile("irf_filepath.txt")
        filepath = open(f->read(f, String), "irf_filepath.txt")
        if !ispath(filepath)
            println("IRF filepath does not exist. Please select valid .sdt file.")
            filepath = open_dialog("IRF: Choose .sdt file to open.")
            open("irf_filepath.txt", "w") do io
                write(io, filepath)
            end
        end
    else
        filepath = open_dialog("IRF: Choose .sdt file to open.")
        open("irf_filepath.txt", "w") do io
            write(io, filepath)
        end
    end

    sdt = pyimport("sdtfile")
    file = sdt.SdtFile(filepath)
    counts = convert.(Int64, file.data[channel])
    median_irf = round(median(counts))
    counts = counts.-round(Int, median_irf)
    counts[counts .<= 0] .= 0
    times = convert.(Float64, vec(file.times[1]).*1e9)
    data = zeros(Float64, length(times), 2)::Matrix{Float64}
    data[:, 1] = times
    data[:, 2] = vec(counts)
    return data
end

function get_new_irf(; channel=1)
    println("Getting IRF file...")

    filepath = pick_file()

    open("irf_filepath.txt", "w") do io
        write(io, filepath)
    end

    sdt = pyimport("sdtfile")
    file = sdt.SdtFile(filepath)
    data = convert.(Int64, file.data[channel])
    median_irf = round(median(data))
    data = data.-round(Int, 1.25*median_irf)
    data[data .<= 0] .= 0
    times = convert.(Float64, vec(file.times[1]).*1e9)
    data = [times, vec(data)]

    println("Done")
    println()

    return data
end

function get_irf_bin_size()
    """
    Gets the bin size in nanoseconds from the IRF.
    """
    irf = get_irf()
    h = Inf::Float64
    for i in 1:1:length(irf[:, 1])
        if i == 1
            continue
        end
        if irf[i, 1]-irf[i-1, 1] < h
            h = irf[i, 1]-irf[i-1, 1]
        end
    end
    return h
end

function convolve(irf::Vector{Float64}, decay::Vector{Float64}; histogram_resolution=256::Int64)
    """
    Convolves the IRF with fluorescence decay
    """
    y = real.(ifft_plan*((fft_plan*irf) .* (fft_plan*decay)))

    return (y/sum(y[13:histogram_resolution-12]))::Vector{Float64}
end

function irf_shift(data_irf::Matrix{Float64}, shift)
    """
    Translation of the IRF to account for systematic delays in instrument
    REF: Enderlein, 1997
    """
    if isnan(shift)
        shift = 0
    end

    n = length(data_irf[:, 1])
    channel = collect(1:1:length(data_irf[:, 1]))
    irf_counts = data_irf[:, 2]
    
    index_1 = vec( ((channel .- (floor(Int, shift) - 1)).% n .% n).% n) .+ 1
    index_1[index_1 .< 1] .+= length(index_1)
    index_2 = vec(((channel .- (ceil(Int, shift) - 1)).% n .+ n).% n) .+ 1
    index_2[index_2 .< 1] .+= length(index_2)
    
    return vec((1 - shift + floor(shift)) .* irf_counts[index_1] + (shift - floor(shift)) .* irf_counts[index_2])::Vector{Float64}
end

function conv_irf_data(x_data::Vector{Float64}, params::Tuple{Float64, Float64, Float64}, irf::Matrix{Float64}; histogram_resolution=256::Int64)
    return convolve(irf_shift(irf, params[2]), exp.(-x_data./params[1]), histogram_resolution=histogram_resolution)
end

function MLE_model_func(params::Vector{Float64}, x_data::Vector{Float64}, y_data::Vector{Float64}, irf::Matrix{Float64}, gating_function::Vector{UInt8}, histogram_resolution::Int64)
    """
    Model function for MLE optimizer.
    # REF: Bajzer, et al. 1991
    # REF: Maus, et al. 2001
    """
    poisson_deviance = 0.0::Float64

    exp_val_c = gating_function.*((1-params[3]).*conv_irf_data(x_data, (params[1], params[2], params[3]), irf, histogram_resolution=histogram_resolution).+params[3]).*float(sum(y_data))

    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, exp_val_c)
    ln_data = log.(y_data./exp_val_c)
    poisson_deviance = 2*sum(y_data.*ln_data.-y_data.+exp_val_c)/(histogram_resolution-length(params))

    if isnan(poisson_deviance)
        println(params, " ", poisson_deviance)
        println(exp_val_c)
        println(y_data)
    end

    return poisson_deviance::Float64
end

function find_mean_arrival_time(counts)
    num = 0.0
    for (i, c) in enumerate(counts)
        num += i*c
    end
    return num/sum(counts)
end

function lifetime_estimate(counts)
    mean_irf_arrival_time = find_mean_arrival_time(irf[:, 2])
    mean_data_arrival_time = find_mean_arrival_time(counts)
    return (mean_data_arrival_time-mean_irf_arrival_time)*irf_bin_size
end

function MLE_iterative_reconvolution_jl(irf::Matrix{Float64}, data_xy::Vector{Vector{Float64}}; params, gating_function=ones(UInt8, 320), histogram_resolution=256, laser_pulse_period=12.5, fixed_parameters=Float64[NaN, NaN, NaN])
    """
    Curve fitting of fluorescence decay using Maximum Likelihood Estimation
    """
    x_data = vec(data_xy[1])
    y_data = vec(data_xy[2]).*gating_function
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, y_data)
    irf = deepcopy(irf)

    params_copy = deepcopy(params)

    if all(isnotnan.(fixed_parameters))
        return fixed_parameters
    end

    params_copy[1] = lifetime_estimate(y_data)
    lower_bounds = Float64[2*irf_bin_size, -16.0, 0.0]
    upper_bounds = Float64[laser_pulse_period*2, 32.0, 1.0]

    if any(isnotnan.(fixed_parameters))
        for i in 1:1:length(fixed_parameters)
            if !isnan(fixed_parameters[i])
                if fixed_parameters[i] < lower_bounds[i]
                    fixed_parameters[i] = lower_bounds[i]
                elseif fixed_parameters[i] > upper_bounds[i]
                    fixed_parameters[i] = upper[i]
                end
                params_copy[i] = fixed_parameters[i]
                if fixed_parameters[i]-1e-6 >= lower_bounds[i]
                    lower_bounds[i] = fixed_parameters[i]-1e-6
                end
                if fixed_parameters[i]+1e-6 <= upper_bounds[i]
                    upper_bounds[i] = fixed_parameters[i]+1e-6
                end
            end
        end
    end

    fit = optimize(x->MLE_model_func(x, x_data, y_data, irf, gating_function, histogram_resolution), lower_bounds, upper_bounds, params_copy, Fminbox(LBFGS(linesearch = LineSearches.BackTracking())), Optim.Options(outer_iterations=5, x_abstol=5e-16, outer_x_abstol=5e-16, allow_f_increases=false))
    res = Optim.minimizer(fit)
    
    if res[1] == lower_bounds[1] || res[1] == upper_bounds[1] || !Optim.converged(fit)
        return Float64[NaN, NaN, NaN]
    end

    if Optim.minimum(fit) > 5 && sum(y_data) < 50000
        println("Ignoring bad fit with optimization result: ", Optim.minimum(fit), " found values: ", res)
        return Float64[NaN, NaN, NaN]
    end

    return res::Vector{Float64}
end

function vec_to_lifetime(x; guess=[3.0, 1.0, 1e-6], laser_pulse_period=12.5, histogram_resolution=256, tac_low_cut=5.0980392, tac_high_cut=94.901962, fixed_parameters=Float64[NaN, NaN, NaN])
    """
    Calculate fluorescence lifetime from vector of counts
    """
    total_channels = round(Int, laser_pulse_period*histogram_resolution/tcspc_window_size)
    x_data = vec(collect(irf_bin_size:irf_bin_size:total_channels*irf_bin_size))

    # if total_channels != histogram_resolution
    #     println("Mismatch")
    #     global fft_plan = plan_fft(zeros(Float64, total_channels))
    #     global ifft_plan = plan_ifft(zeros(Float64, total_channels))
    #     new_irf = zeros(Float64, total_channels, 2)
    #     new_irf[:, 1] = collect(irf_bin_size:irf_bin_size:irf_bin_size*total_channels)
    #     new_irf[1:length(irf[:, 2]), 2] = irf[:, 2]
    #     global irf = new_irf
    #     append!(vec(x), zeros(total_channels-histogram_resolution))
    # end
    
    if sum(x) < 100
        return Float64[NaN], [x_data, vec(x)]::Vector{Vector{Float64}}
    else
        gating_function = ones(UInt8, total_channels)
        gating_function[1:round(Int, tac_low_cut/100*histogram_resolution)] .= 0

        # tac_low_cut_index = round(Int, tac_low_cut/100*histogram_resolution)
        gating_function[round(Int, tac_high_cut/100*histogram_resolution):total_channels] .= 0
        # tcspc_high_cut_index = total_channels-round(Int, tac_high_cut/100*histogram_resolution)

        data_xy = [x_data, vec(x)]
        tau_fit = MLE_iterative_reconvolution_jl(irf, data_xy, params=guess, gating_function=gating_function, histogram_resolution=total_channels, fixed_parameters=fixed_parameters, laser_pulse_period=laser_pulse_period)

        return tau_fit::Vector{Float64}, data_xy::Vector{Vector{Float64}}
    end
end

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
    MAX_LINES = 9
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
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    params = [4.0, 2.0, 1.0e-4]
    erreur_I = 0

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

                time_range = params_sliders.sliders[1].value[]
                averaging = params_sliders.sliders[2].value[]
                prop_val = params_sliders.sliders[3].value[]
                integ_val = params_sliders.sliders[4].value[]
                setpoint_val = params_sliders.sliders[5].value[]
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

                erreur = setpoint_val - params[1]

                t_now = (now() - t_0)/Millisecond(1)/1000

                try
                    dt = t_now - times[][end]
                catch e
                    dt = t_now
                    # println(e)
                end

                if reset_I[]
                    erreur_I = 0
                    reset_I[] = false
                end

                erreur_I += erreur * dt

                if i == 1
                    erreur_I = 0
                end

                if inverted
                    signal_value = clamp(-(prop_val * erreur + integ_val * erreur_I), 0, 100)
                else
                    signal_value = clamp(prop_val * erreur + integ_val * erreur_I, 0, 100)
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
    counts = Observable(Float64[])
    global signal = Observable(Float64[])
    times = Observable(Float64[])
    test_times = []

    vectors = zeros(100, 256)
    n_vectors = size(vectors, 1)
    histogram_resolution = 256

    barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])
    lines!(axis_signal, times, signal, color=Makie.wong_colors()[1])

    params = [5.0, 2.0, 1.0e-4]
    erreur_I = 0
    display_period = 1/10

    i = 1
    n = 1
    t_0 = now()
    cum_time = 0

    println("Starting playback")
    add_to_monitor("Starting playback")

    while running[]
        for file in files
            try 
                averaging = params_sliders.sliders[2].value[]
                setpoint_val = params_sliders.sliders[5].value[]
                prop_val = params_sliders.sliders[3].value[]
                integ_val = params_sliders.sliders[4].value[]
                time_range = params_sliders.sliders[1].value[]
                inverted = checkbox_inverted.active[]
                
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

                t_test = time()
                params_raw, data = vec_to_lifetime(final_vector, guess=params, histogram_resolution=histogram_resolution)
                push!(test_times, (time() - t_test)*1000)

                nb_photons = sum(final_vector)

                # println(params_raw)
                # println(data)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=histogram_resolution)*nb_photons

                if i == 1
                    t_0 = now()
                end

                erreur = setpoint_val - params[1]

                delta_t = (now() - t_0) / Microsecond(1) / 1e6

                if reset_I[]
                    erreur_I = 0
                    reset_I[] = false
                end

                erreur_I += erreur * duration

                if i == 1
                    erreur_I = 0
                end

                if inverted
                    signal_value = clamp(-(prop_val * erreur + integ_val * erreur_I), 0, 100)
                else
                    signal_value = clamp(prop_val * erreur + integ_val * erreur_I, 0, 100)
                end

                push!(lifetime[], params[1])
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
    analyse_time(test_times)
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

    box_progress.width = 8

    beginning = time()

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
                box_progress.width = 8 + 1345*i/nb_files
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

    filename = string(now())*".csv"

    df = DataFrame(
        Time_s = times,
        Number_of_photons = counts,
        Lifetime_ns = lifetime
    )
    
    open(filename, "w") do io
        println(io, "Timestamp: ", now())
        CSV.write(io, df)
    end

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

set_theme!(backgroundcolor = :gray95)

fig = Figure(size = (1440, 847), msaa = 1)

fig[1, 1]          = plotgrid     = GridLayout()
fig[2, 1]          = buttongrid   = GridLayout(tellwidth=true)
fig[3, 1]          = progressgrid = GridLayout(tellwidth=true)
buttongrid[1:3, 8] = paramsgrid   = GridLayout(tellwidth=true)

# Box(buttongrid[1:3, 1])
# Box(buttongrid[1, 1])
# Box(buttongrid[2, 1])
# Box(buttongrid[3, 1])
# Box(buttongrid[3, 1][1, 1])
# Box(buttongrid[3, 1][1, 1])

global axis_hist   = Axis(plotgrid[1, 1], alignmode = Outside(), title="Histogram",       xlabel = "Time [ns]", ylabel = "Photon counts", titlesize=16, width = 600)
global axis_counts = Axis(plotgrid[2, 1], alignmode = Outside(), title="Photon counts",   xlabel = "Time [s]",  ylabel = "Photon counts", titlesize=16, width = 600)
global axis_life   = Axis(plotgrid[1, 3], alignmode = Outside(), title="Lifetime",        xlabel = "Time [s]",  ylabel = "Lifetime [ns]", titlesize=16, width = 600)
global axis_signal = Axis(plotgrid[2, 3], alignmode = Outside(), title="Light intensity", xlabel = "Time [s]",  ylabel = "[%]",           titlesize=16, width = 600)

ylims!(axis_signal, -5, 105)

label_frequency    = Label(buttongrid[1, 1],       justification=:left, text="Frequency : 0 Hz", width = Auto())
label_number       = Label(buttongrid[2, 1],       justification=:left, text="Files : 0",        width = nothing)
label_inverted     = Label(buttongrid[3, 1][1, 1], justification=:left, text="Inverted",         width = nothing)

checkbox_inverted  = Toggle(buttongrid[3, 1][1, 2], active = true, rimfraction=0.24, buttoncolor=dark_gray, framecolor_active=mid_gray, framecolor_inactive=light_gray)

box_attributes = Dict(
    :color          => :black,
    :cornerradius   => 4,
    :height         => 40,
    :strokewidth    => 0.5,
    :width          => 120,
)

button_attributes = Dict(
    :buttoncolor         => light_gray,
    :buttoncolor_active  => dark_gray,
    :buttoncolor_hover   => mid_gray,
    :cornerradius        => 4,
    :font                => :bold,
    :fontsize            => 16,
    :height              => 40,
    :labelcolor_active   => :black,
    :width               => 120,
)

slider_attributes = Dict(
    :color_active        => dark_gray,
    :color_active_dimmed => mid_gray,
    :color_inactive      => light_gray,
    :linewidth           => 16,
    :tellheight          => true,
)

textbox_attributes = Dict(
    :bordercolor               => :black,
    :bordercolor_focused       => :black,
    :bordercolor_hover         => :black,
    :borderwidth               => 0.5,
    :boxcolor                  => :white,
    :boxcolor_focused          => :white,
    :boxcolor_focused_invalid  => RGBf(1, 0.9, 0.9),
    :boxcolor_hover            => off_white,
    :cornerradius              => 3,
    :halign                    => :right,
    :height                    => 26,
    :placeholder               => " ",
    :textpadding               => (10, 8, 8, 4.5),
    :width                     => 60,
)

spinnerbox_attributes = Dict(
    :bordercolor                 => :black,
    :bordercolor_focused         => :black,
    :bordercolor_focused_invalid => :black,
    :bordercolor_hover           => :black,
    :borderwidth                 => 0,
    :boxcolor                    => :white,
    :boxcolor_focused            => :white,
    :boxcolor_focused_invalid    => RGBf(1, 0.9, 0.9),
    :boxcolor_hover              => off_white,
    :cornerradius                => 0,
    :halign                      => :left,
    :height                      => 24,
    :placeholder                 => " ",
    :reset_on_defocus            => true,
    :textcolor_placeholder       => :black,
    :textpadding                 => (8, 0, 0, 4),
    :width                       => 64
)

spinnerup_attributes = Dict(
    :buttoncolor        => :gray95,
    :buttoncolor_active => dark_gray,
    :buttoncolor_hover  => light_gray,
    :cornerradius       => 0,
    :font               => :bold,
    :fontsize           => 6,
    :halign             => :right,
    :height             => 12,
    :label              => "▲",
    :labelcolor_active  => :black,
    :strokecolor        => :black,
    :strokewidth        => 1,
    :valign             => :top,
    :width              => 32
)

spinnerdown_attributes = Dict(
    :buttoncolor        => :gray95,
    :buttoncolor_active => dark_gray,
    :buttoncolor_hover  => light_gray,
    :cornerradius       => 0,
    :font               => :bold,
    :fontsize           => 6,
    :halign             => :right,
    :height             => 12,
    :label              => "▼",
    :labelcolor_active  => :black,
    :strokecolor        => :black,
    :strokewidth        => 1,
    :valign             => :bottom,
    :width              => 32
)

Box(buttongrid[1, 2]; box_attributes...)
Box(buttongrid[2, 2]; box_attributes...)
Box(buttongrid[3, 2]; box_attributes...)
Box(buttongrid[1, 3]; box_attributes...)
Box(buttongrid[2, 3]; box_attributes...)
Box(buttongrid[3, 3]; merge(box_attributes, Dict(:cornerradius => -0.1))...)
Box(buttongrid[1, 6]; box_attributes...)
Box(buttongrid[2, 6]; box_attributes...)
Box(buttongrid[3, 6]; merge(box_attributes, Dict(:cornerradius => -0.1))...)

button_start    = Button(buttongrid[1, 2]; merge(button_attributes, Dict(:label => "START"))...)
button_stop     = Button(buttongrid[1, 3]; merge(button_attributes, Dict(:label => "CLEAR"))...)
button_load     = Button(buttongrid[2, 2]; merge(button_attributes, Dict(:label => "LOAD"))...)
button_load_IRF = Button(buttongrid[2, 3]; merge(button_attributes, Dict(:label => "LOAD IRF"))...)
button_reset    = Button(buttongrid[3, 2]; merge(button_attributes, Dict(:label => "RESET"))...)

menu_mode = Menu(buttongrid[3, 3],
    cell_color_active=RGBf(0.8, 0.8, 0.8),
    cell_color_inactive_even=RGBf(0.9, 0.9, 0.9),
    cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9),
    cell_color_hover=RGBf(0.6, 0.6, 0.6),
    default = "Real-time",
    fontsize=20,
    height = 40,
    options=["Analysis", "Playback", "Real-time", "Save"],
    selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8),
    textcolor=:black,
    width = 120
)

portname_textbox = Textbox(buttongrid[1:3, 4:5][1,1:2]; merge(textbox_attributes, Dict(:cornerradius => 3, :height => 26, :tellwidth => false, :width => nothing))...)
portname_textbox.displayed_string[] = "usbmodem2302"
portname_textbox.stored_string[]    = "usbmodem2302"

label_fps   = Label(buttongrid[1:3, 4:5][2, 1], text=L"\text{fps [Hz] :}", halign=:right, width => nothing)
textbox_fps = Textbox(buttongrid[1:3, 4:5][2, 2]; textbox_attributes...)
textbox_fps.displayed_string[] = "50"
textbox_fps.stored_string[]    = "50"

label_min   = Label(buttongrid[1:3, 4:5][3, 1], text=L"\text{Minimum [ms] :}", halign=:right, width => nothing)
textbox_min = Textbox(buttongrid[1:3, 4:5][3, 2]; textbox_attributes...)
textbox_min.displayed_string[] = "3"
textbox_min.stored_string[]    = "3"

label_max   = Label(buttongrid[1:3, 4:5][4, 1], text=L"\text{Maximum [ms] :}", halign=:right, width => nothing)
textbox_max = Textbox(buttongrid[1:3, 4:5][4, 2]; textbox_attributes...)
textbox_max.displayed_string[] = "14"
textbox_max.stored_string[]    = "14"

button_active  = Button(buttongrid[1, 6]; merge(button_attributes, Dict(:label => "Activate", :buttoncolor => mid_gray, :buttoncolor_hover => dark_gray))...)
button_connect = Button(buttongrid[2, 6]; merge(button_attributes, Dict(:label => "Connect",  :buttoncolor => mid_gray, :buttoncolor_hover => dark_gray))...)

menu_output = Menu(buttongrid[3, 6],
    cell_color_active = RGBf(0.8, 0.8, 0.8),
    cell_color_inactive_even = RGBf(0.9, 0.9, 0.9),
    cell_color_inactive_odd = RGBf(0.9, 0.9, 0.9),
    cell_color_hover = RGBf(0.6, 0.6, 0.6),
    default = "Digital",
    fontsize = 20,
    height = 40,
    options = ["Analog", "Digital"],
    selection_cell_color_inactive = RGBf(0.8, 0.8, 0.8),
    textcolor = :black,
    width = 120
)

monitor_box = Box(buttongrid[1:3, 7],
    color = :white,
    cornerradius = 0,
    height = 160,
    width = 320
)

global monitor = Label(buttongrid[1:3, 7],
    "\n",
    height = 160,
    justification = :left,
    padding = (8, 8, 8, 8),
    tellheight = false,
    tellwidth = false,
    valign = :bottom,
    width = 320
)

# println(monitor_box.width[])

# params_sliders = SliderGrid(
#     buttongrid[1:3, 8],
#     NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Time range",   :range => 10:10:3600, :startvalue => 30))),
#     NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Averaging",    :range => 1:1:100,    :startvalue => 1))),
#     NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Proportional", :range => 0:1:100,    :startvalue => 0))),
#     NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Integral",     :range => 0:1:100,    :startvalue => 0))),
#     NamedTuple(merge(slider_attributes, Dict(:format => "{:.2f}", :label => "Setpoint",     :range => 3:0.01:7,   :startvalue => 5))),
#     width = 360
# )

# fps_label = Label(paramsgrid[1, 1],
#     fontsize = 16,
#     justification = :right,
#     text = "Max fps:")

Box(paramsgrid[3, 1], strokewidth = 1, width = 96)
Box(paramsgrid[3, 2], strokewidth = 1, width = 96)
Box(paramsgrid[3, 3], strokewidth = 1, width = 96)
Box(paramsgrid[6, 1], strokewidth = 1, width = 96)
Box(paramsgrid[6, 2], strokewidth = 1, width = 96)
Box(paramsgrid[6, 3], strokewidth = 1, width = 96)


function range_validator(a::Integer, b::Integer)
    @assert a ≤ b "Lower bound must not exceed upper bound"
    return function(s::String)
        s_trim = strip(s)
        x = tryparse(Int, s_trim)
        return x !== nothing && a ≤ x ≤ b
    end
end

Label(paramsgrid[1, 1:3], text="General settings", halign=:center, height = nothing, tellwidth = true, width = Auto())

Label(paramsgrid[2, 1], text="Time range", halign=:center, height = nothing, tellwidth = true, width = Auto())
range_textbox = Textbox(paramsgrid[3, 1]; merge(spinnerbox_attributes, Dict(:displayed_string => "60", :stored_string => "60", :validator => UInt32))...)
range_up      = Button(paramsgrid[3, 1];  spinnerup_attributes...)
range_down    = Button(paramsgrid[3, 1];  spinnerdown_attributes...)

Label(paramsgrid[2, 2], text="Averaging", halign=:center, height = nothing, tellwidth = true, width = Auto())
averaging_textbox = Textbox(paramsgrid[3, 2]; merge(spinnerbox_attributes, Dict(:displayed_string => "1", :stored_string => "1", :validator => range_validator(1, 100)))...)
averaging_up      = Button(paramsgrid[3, 2];  spinnerup_attributes...)
averaging_down    = Button(paramsgrid[3, 2];  spinnerdown_attributes...)

Label(paramsgrid[2, 3], text="Lifetimes", halign=:center, height = nothing, tellwidth = true, width = Auto())
lifetimes_textbox = Textbox(paramsgrid[3, 3]; merge(spinnerbox_attributes, Dict(:displayed_string => "1", :stored_string => "1", :validator => range_validator(1, 3)))...)
lifetimes_up      = Button(paramsgrid[3, 3];  spinnerup_attributes...)
lifetimes_down    = Button(paramsgrid[3, 3];  spinnerdown_attributes...)

Label(paramsgrid[4, 1:3], text="Controller settings", halign=:center, height = nothing, tellwidth = true, width = Auto())

Label(paramsgrid[5, 1], text="Proportionnal", halign=:center, height = nothing, tellwidth = true, width = Auto())
proportionnal_textbox = Textbox(paramsgrid[6, 1]; merge(spinnerbox_attributes, Dict(:displayed_string => "0", :stored_string => "0", :validator => UInt32))...)
proportionnal_up      = Button(paramsgrid[6, 1];  spinnerup_attributes...)
proportionnal_down    = Button(paramsgrid[6, 1];  spinnerdown_attributes...)

Label(paramsgrid[5, 2], text="Integral", halign=:center, height = nothing, tellwidth = true, width = Auto())
integral_textbox = Textbox(paramsgrid[6, 2]; merge(spinnerbox_attributes, Dict(:displayed_string => "0", :stored_string => "0", :validator => UInt32))...)
integral_up      = Button(paramsgrid[6, 2];  spinnerup_attributes...)
integral_down    = Button(paramsgrid[6, 2];  spinnerdown_attributes...)

Label(paramsgrid[5, 3], text="Setpoint", halign=:center, height = nothing, tellwidth = true, width = Auto())
setpoint_textbox = Textbox(paramsgrid[6, 3]; merge(spinnerbox_attributes, Dict(:displayed_string => "5.0", :stored_string => "5.0", :validator => Float64))...)
setpoint_up      = Button(paramsgrid[6, 3];  spinnerup_attributes...)
setpoint_down    = Button(paramsgrid[6, 3];  spinnerdown_attributes...)

rowgap!(paramsgrid, 4)
rowgap!(paramsgrid, 3, 24)

# Box(paramsgrid[1:6, 1:2])
# Box(paramsgrid[1, 2])

# rowgap!(paramsgrid, 1, 4)
# rowgap!(paramsgrid, 3, 4)
# rowgap!(paramsgrid, 5, 4)

# Label(paramsgrid[1, 2], text="Time range", halign=:center, height = nothing, tellwidth = true, width = Auto())
# Box(paramsgrid[2, 2]; spinnerbox_attributes...)
# Box(paramsgrid[2, 2]; spinnerup_attributes...)
# Box(paramsgrid[2, 2]; spinnerdown_attributes...)

# Label(paramsgrid[3, 2], text="Time range", halign=:center, height = nothing, tellwidth = true, width = Auto())
# Box(paramsgrid[4, 2]; spinnerbox_attributes...)
# Box(paramsgrid[4, 2]; spinnerup_attributes...)
# Box(paramsgrid[4, 2]; spinnerdown_attributes...)

# Label(paramsgrid[5, 2], text="Time range", halign=:center, height = nothing, tellwidth = true, width = Auto())
# Box(paramsgrid[6, 2]; spinnerbox_attributes...)
# Box(paramsgrid[6, 2]; spinnerup_attributes...)
# Box(paramsgrid[6, 2]; spinnerdown_attributes...)


Box(progressgrid[1, 1],
    color = light_gray,
    cornerradius = 4,
    halign = :left,
    height = 8,
    strokevisible=false,
    width = nothing,
)

box_progress = Box(progressgrid[1, 1],
    color = mid_gray,
    cornerradius = 4,
    halign = :left,
    height = 8,
    strokevisible=false,
    tellwidth = false,
    width = 0
)

label_progress = Label(progressgrid[1, 2], text="100 %", halign=:right, tellwidth = true, width = Auto())

colsize!(plotgrid, 2, 20)
colsize!(buttongrid, 4, 120)
colsize!(buttongrid, 5, 60)
# colsize!(buttongrid, 8, 0)
colgap!(buttongrid, 4, 0)
rowgap!(buttongrid, 1, 20)
rowgap!(buttongrid, 2, 20)

global irf = get_irf()
global irf_bin_size = get_irf_bin_size()
global tcspc_window_size = round(irf[end, 1]+irf[2, 1], sigdigits=4)
global fft_plan = plan_fft(zeros(Float64, 256))
global ifft_plan = plan_ifft(zeros(Float64, 256))
global folderpath = get_folder_path("FLIM")
global serial = nothing
global running = Observable(false)
global pause = Observable(false)
global active = Observable(false)
global connected = Observable(false)
global reset_I = Observable(false)

setpoint = hlines!(axis_life, parse(Float64, setpoint_textbox.stored_string[]), color=Makie.wong_colors()[2])
# setpoint = hlines!(axis_life, params_sliders.sliders[5].value[], color=Makie.wong_colors()[2])

on(button_start.clicks) do _
    if !running[]
        running[] = true

        button_start.label = "PAUSE"
        button_stop.label = "STOP"

        box_progress.width = 0
        label_progress.text = "0 %"
        label_frequency.text = "Frequency : 0 Hz"
        label_number.text = "Files : 0"

        empty!(axis_hist)
        empty!(axis_counts)
        empty!(axis_life)
        empty!(axis_signal)

        global setpoint = hlines!(axis_life, params_sliders.sliders[5].value[], color=Makie.wong_colors()[2])

        if menu_mode.selection[] == "Playback"
            @async start_playback()
        elseif menu_mode.selection[] == "Analysis"
            start_analysis()
        elseif menu_mode.selection[] == "Save"
            @async start_save()
        elseif menu_mode.selection[] == "Real-time"
            @async start_real_time()
        end
    elseif pause[]
        pause[] = false

        button_start.label = "PAUSE"
        button_stop.label = "STOP"
    else
        pause[] = true

        button_start.label = "START"
        button_stop.label = "CLEAR"
    end
end

on(button_stop.clicks) do _
    running[] = false
    pause[] = false

    if connected[]
        clear(serial)
    end

    button_start.label = "START"
    button_stop.label = "CLEAR"

    box_progress.width = 0
    label_progress.text = "0 %"
    label_frequency.text = "Frequency : 0 Hz"
    label_number.text = "Files : 0"

    empty!(axis_hist)
    empty!(axis_counts)
    empty!(axis_life)
    empty!(axis_signal)

    global setpoint = hlines!(axis_life, params_sliders.sliders[5].value[], color=Makie.wong_colors()[2])
end

on(button_load.clicks) do _
    global folderpath = get_new_folder_path("FLIM")
end

on(button_load_IRF.clicks) do _
    global irf = get_new_irf()
    global irf_bin_size = get_irf_bin_size()
end

on(button_reset.clicks) do _
    reset_I[] = true
end

on(button_connect.clicks) do _
    if connected[]
        clear(serial)
        LibSerialPort.close(serial)
        add_to_monitor("Disconnection successful")
        add_to_monitor("Deactivated")
        println("Disconnection successful")

        button_connect.label[] = "Connect"
        button_active.label = "Activate"
        connected[] = false
        active[] = false

        button_active.buttoncolor[] = mid_gray
        button_active.buttoncolor_hover[] = dark_gray

        button_connect.buttoncolor[] = mid_gray
        button_connect.buttoncolor_hover[] = dark_gray
    else
        global serial = connect_to_port(portname_textbox.stored_string[])

        sleep(0.1)
        
        try
            send_command(serial, "C Z\nC R\n")
            send_command(serial, "S 1 S 1\n")
        catch e
            println(e)
        end

        button_connect.label[] = "Disconnect"

        button_connect.buttoncolor[] = light_gray
        button_connect.buttoncolor_hover[] = mid_gray
    end
end

on(button_active.clicks) do _
    if connected[]
        if active[]
            active[] = false
            button_active.label = "Activate"

            button_active.buttoncolor[] = mid_gray
            button_active.buttoncolor_hover[] = dark_gray

            add_to_monitor("Deactivated")
            clear(serial)
        else
            active[] = true
            button_active.label = "Deactivate"

            button_active.buttoncolor[] = light_gray
            button_active.buttoncolor_hover[] = mid_gray

            add_to_monitor("Activated")
            @async start_signal(parse(Float64, textbox_fps.stored_string[]), parse(Float64, textbox_fps.stored_string[])/10, parse(Float64, textbox_min.stored_string[]))
        end
    end
end

on(range_up.clicks) do _
    if !isnothing(tryparse(UInt32, range_textbox.displayed_string[]))
        value = parse(UInt32, range_textbox.displayed_string[])
        multiplier = floor(Int, round(log(10, value), digits = 4))
        res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

        range_textbox.stored_string[] = string(res)
        range_textbox.displayed_string[] = string(res)
    else
        range_textbox.stored_string[] = "60"
        range_textbox.displayed_string[] = "60"
    end
end

on(range_down.clicks) do _
    if !isnothing(tryparse(UInt32, range_textbox.displayed_string[]))
        value = parse(UInt32, range_textbox.displayed_string[])

        if value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            range_textbox.stored_string[] = string(res)
            range_textbox.displayed_string[] = string(res)
        end
    else
        range_textbox.stored_string[] = "60"
        range_textbox.displayed_string[] = "60"
    end
end

on(averaging_up.clicks) do _
    if !isnothing(tryparse(UInt32, averaging_textbox.displayed_string[]))
        value = parse(UInt32, averaging_textbox.displayed_string[])

        if value < 1 | value > 100
            averaging_textbox.stored_string[] = "1"
            averaging_textbox.displayed_string[] = "1"
        end
        
        if value != 100
            averaging_textbox.stored_string[] = string(value + 1)
            averaging_textbox.displayed_string[] = string(value + 1)
        end
    else
        averaging_textbox.stored_string[] = "1"
        averaging_textbox.displayed_string[] = "1"
    end
end

on(averaging_down.clicks) do _
    if !isnothing(tryparse(UInt32, averaging_textbox.displayed_string[]))
        value = parse(UInt32, averaging_textbox.displayed_string[])

        if value < 1 | value > 100
            averaging_textbox.stored_string[] = "1"
            averaging_textbox.displayed_string[] = "1"
        end

        if value != 1
            averaging_textbox.stored_string[] = string(value - 1)
            averaging_textbox.displayed_string[] = string(value - 1)
        end
    else
        averaging_textbox.stored_string[] = "1"
        averaging_textbox.displayed_string[] = "1"
    end
end

on(lifetimes_up.clicks) do _
    if !isnothing(tryparse(UInt32, lifetimes_textbox.displayed_string[]))
        value = parse(UInt32, lifetimes_textbox.displayed_string[])

        if value == 1
            lifetimes_textbox.stored_string[] = "2"
            lifetimes_textbox.displayed_string[] = "2"
        elseif value == 2
            lifetimes_textbox.stored_string[] = "3"
            lifetimes_textbox.displayed_string[] = "3"
        elseif value == 3
        else 
            lifetimes_textbox.stored_string[] = "1"
            lifetimes_textbox.displayed_string[] = "1"
        end
    else
        lifetimes_textbox.stored_string[] = "1"
        lifetimes_textbox.displayed_string[] = "1"
    end
end

on(lifetimes_down.clicks) do _
    if !isnothing(tryparse(UInt32, lifetimes_textbox.displayed_string[]))
        value = parse(UInt32, lifetimes_textbox.displayed_string[])

        if value == 1
        elseif value == 2
            lifetimes_textbox.stored_string[] = "1"
            lifetimes_textbox.displayed_string[] = "1"
        elseif value == 3
            lifetimes_textbox.stored_string[] = "2"
            lifetimes_textbox.displayed_string[] = "2"
        else 
            lifetimes_textbox.stored_string[] = "1"
            lifetimes_textbox.displayed_string[] = "1"
        end
    else
        lifetimes_textbox.stored_string[] = "1"
        lifetimes_textbox.displayed_string[] = "1"
    end
end

# on(proportionnal_up.clicks) do _
#     if !isnothing(tryparse(Float64, proportionnal_textbox.displayed_string[]))
#         value = parse(Float64, proportionnal_textbox.displayed_string[])
#         if value < 0
#             proportionnal_textbox.stored_string[] = "0"
#             proportionnal_textbox.displayed_string[] = "0"
#         elseif value == 0
#             proportionnal_textbox.stored_string[] = "0.1"
#             proportionnal_textbox.displayed_string[] = "0.1"
#         else
#             multiplier = max(floor(Int, round(log(10, value) - 1, digits = 4)), -1)

#             proportionnal_textbox.stored_string[] = string(round(floor(Int, value*0.1^multiplier)*10.0^multiplier + 1*10.0^multiplier, digits=1))
#             proportionnal_textbox.displayed_string[] = string(round(floor(Int, value*0.1^multiplier)*10.0^multiplier + 1*10.0^multiplier, digits=1))
#         end
#     else
#         proportionnal_textbox.stored_string[] = "0"
#         proportionnal_textbox.displayed_string[] = "0"
#     end
# end

# on(proportionnal_down.clicks) do _
#     if !isnothing(tryparse(Float64, proportionnal_textbox.displayed_string[]))
#         value = parse(Float64, proportionnal_textbox.displayed_string[])
#         if value <= 0
#             proportionnal_textbox.stored_string[] = "0"
#             proportionnal_textbox.displayed_string[] = "0"
#         # elseif value != 1
#         else
#             multiplier = max(floor(Int, round(log(10, value - .1) - 1, digits = 8)), -1)

#             proportionnal_textbox.stored_string[] = string(round(ceil(Int, round(value*0.1^multiplier, digits = 8))*10.0^multiplier - 1*10.0^multiplier, digits=1))
#             proportionnal_textbox.displayed_string[] = string(round(ceil(Int, round(value*0.1^multiplier, digits = 8))*10.0^multiplier - 1*10.0^multiplier, digits=1))
#         end
#     else
#         proportionnal_textbox.stored_string[] = "0"
#         proportionnal_textbox.displayed_string[] = "0"
#     end
# end

on(proportionnal_up.clicks) do _
    if !isnothing(tryparse(UInt32, proportionnal_textbox.displayed_string[]))
        value = parse(UInt32, proportionnal_textbox.displayed_string[])
        if value != 0
            multiplier = floor(Int, round(log(10, value), digits = 4))
            res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

            proportionnal_textbox.stored_string[] = string(res)
            proportionnal_textbox.displayed_string[] = string(res)
        else
            proportionnal_textbox.stored_string[] = "1"
            proportionnal_textbox.displayed_string[] = "1"
        end
    else
        proportionnal_textbox.stored_string[] = "0"
        proportionnal_textbox.displayed_string[] = "0"
    end
end

on(proportionnal_down.clicks) do _
    if !isnothing(tryparse(UInt32, proportionnal_textbox.displayed_string[]))
        value = parse(UInt32, proportionnal_textbox.displayed_string[])

        if value == 0
        elseif value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            proportionnal_textbox.stored_string[] = string(res)
            proportionnal_textbox.displayed_string[] = string(res)
        else
            proportionnal_textbox.stored_string[] = "0"
            proportionnal_textbox.displayed_string[] = "0"
        end
    else
        proportionnal_textbox.stored_string[] = "0"
        proportionnal_textbox.displayed_string[] = "0"
    end
end

on(integral_up.clicks) do _
    if !isnothing(tryparse(UInt32, integral_textbox.displayed_string[]))
        value = parse(UInt32, integral_textbox.displayed_string[])
        if value != 0
            multiplier = floor(Int, round(log(10, value), digits = 4))
            res = floor(Int, value*0.1^multiplier)*10^multiplier + 1*10^multiplier

            integral_textbox.stored_string[] = string(res)
            integral_textbox.displayed_string[] = string(res)
        else
            integral_textbox.stored_string[] = "1"
            integral_textbox.displayed_string[] = "1"
        end
    else
        integral_textbox.stored_string[] = "0"
        integral_textbox.displayed_string[] = "0"
    end
end

on(integral_down.clicks) do _
    if !isnothing(tryparse(UInt32, integral_textbox.displayed_string[]))
        value = parse(UInt32, integral_textbox.displayed_string[])

        if value == 0
        elseif value != 1
            multiplier = floor(Int, round(log(10, value - 1), digits = 8))
            res = ceil(Int, round(value*0.1^multiplier, digits = 8))*10^multiplier - 1*10^multiplier

            integral_textbox.stored_string[] = string(res)
            integral_textbox.displayed_string[] = string(res)
        else
            integral_textbox.stored_string[] = "0"
            integral_textbox.displayed_string[] = "0"
        end
    else
        integral_textbox.stored_string[] = "0"
        integral_textbox.displayed_string[] = "0"
    end
end

on(setpoint_up.clicks) do _
    if !isnothing(tryparse(Float64, setpoint_textbox.displayed_string[]))
        value = parse(Float64, setpoint_textbox.displayed_string[])
        if value < 0
            setpoint_textbox.stored_string[] = "0"
            setpoint_textbox.displayed_string[] = "0"
        elseif value == 0
            setpoint_textbox.stored_string[] = "0.1"
            setpoint_textbox.displayed_string[] = "0.1"
        else
            multiplier = max(floor(Int, round(log(10, value) - 1, digits = 4)), -1)

            setpoint_textbox.stored_string[] = string(round(floor(Int, value*0.1^multiplier)*10.0^multiplier + 1*10.0^multiplier, digits=1))
            setpoint_textbox.displayed_string[] = string(round(floor(Int, value*0.1^multiplier)*10.0^multiplier + 1*10.0^multiplier, digits=1))
        end
    else
        setpoint_textbox.stored_string[] = "5"
        setpoint_textbox.displayed_string[] = "5"
    end
end

on(setpoint_down.clicks) do _
    if !isnothing(tryparse(Float64, setpoint_textbox.displayed_string[]))
        value = parse(Float64, setpoint_textbox.displayed_string[])
        if value <= 0
            setpoint_textbox.stored_string[] = "0"
            setpoint_textbox.displayed_string[] = "0"
        # elseif value != 1
        else
            multiplier = max(floor(Int, round(log(10, value - .1) - 1, digits = 8)), -1)

            setpoint_textbox.stored_string[] = string(round(ceil(Int, round(value*0.1^multiplier, digits = 8))*10.0^multiplier - 1*10.0^multiplier, digits=1))
            setpoint_textbox.displayed_string[] = string(round(ceil(Int, round(value*0.1^multiplier, digits = 8))*10.0^multiplier - 1*10.0^multiplier, digits=1))
        end
    else
        setpoint_textbox.stored_string[] = "5"
        setpoint_textbox.displayed_string[] = "5"
    end
end

on(setpoint_textbox.stored_string) do value
    setpoint.args[1][] = parse(Float64, value)
end

# on(params_sliders.sliders[5].value) do value
#     setpoint.args[1][] = value
# end

fig
