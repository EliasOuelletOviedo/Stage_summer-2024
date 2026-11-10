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

dark_gray = RGBf(0.4, 0.4, 0.4)
mid_gray = RGBf(0.6, 0.6, 0.6)
light_gray = RGBf(0.8, 0.8, 0.8)
off_white = RGBf(0.9, 0.9, 0.9)

try
    LibSerialPort.close(serial)
catch
end

sp_o = pyimport("scipy.optimize")
sdt = pyimport("sdtfile")

nanmean(x) = mean(filter(!isnan, x))
nanmean(x, dims) = mapslices(nanmean, x, dims=dims)

nanstd(x) = std(filter(!isnan, x))
nanstd(x, dims) = mapslices(nanstd, x, dims=dims)

smaller_or_eq_zero(x) = x <= 0

function get_irf(; channel=1)
    println("Getting IRF file...")

    if isfile("irf_filepath.txt")
        filepath = open(f->read(f, String), "irf_filepath.txt")
    else
        filepath = pick_file()

        open("irf_filepath.txt", "w") do io
            write(io, filepath)
        end
    end

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

function get_irf_bin_size(irf)
    h = Inf
    for (i, t) in enumerate(irf[1])
        if i == 1
            continue
        end
        if irf[1][i]-irf[1][i-1] < h
            h = irf[1][i]-irf[1][i-1]
        end
    end
    return float(h)
end

function convolve(irf, decay, y_data; histogram_resolution=256)
    irf_data = deepcopy(irf)
    size_data = size(decay)[1]

    if size_data > length(irf_data)
        append!(irf_data, zeros(size_data-length(irf_data), 1))
    else
        irf_data = irf_data[1:size_data]
    end

    exp_fft = fft(decay)
    irf_data_fft = fft(irf_data)
    y = real.(ifft(irf_data_fft .* exp_fft))
    y = y*(sum(y_data)/sum(y[1:histogram_resolution]))

    return y
end

function irf_shift(data_irf, shift)
    if isnan(shift)
        shift = 0
    end

    n = length(data_irf[1])
    channel = collect(1:1:length(data_irf[1]))
    irf = convert.(Int, data_irf[2])

    index_1 = reshape(convert(Vector{Int64}, ((channel .- (floor(shift) - 1)).% n .% n).% n), 1, :) .+ 1
    index_1[index_1 .< 1] .+= length(index_1)
    index_2 = reshape(convert(Vector{Int64}, ((channel .- (ceil(shift) - 1)).% n .+ n).% n), 1, :) .+ 1
    index_2[index_2 .< 1] .+= length(index_2)
    terme_1 = (1 - shift + floor(shift)) .* irf[index_1]
    terme_2 = (shift - floor(shift)) .* irf[index_2]

    shifted_irf = terme_1 + terme_2

    return vec(shifted_irf)
end

function MLE_iterative_reconvolution(irf, data_xy; params, gating_function=ones(UInt8, 256), histogram_resolution=256)
    x_data = vec(data_xy[1])
    y_data = vec(data_xy[2]).*gating_function
    irf = deepcopy(irf)
    params_copy = deepcopy(params)
    min_offset = floor(mean(y_data[1:argmax(y_data)-10]))+ceil(std(y_data[1:argmax(y_data)-10]))
    
    if min_offset < 5
        min_offset = 5
    end

    push!(params_copy, min_offset)
    
    bounds = [(0.1, 12.5), (-16, 64), (0, min_offset)]
    res = sp_o.minimize(MLE_model_func, params_copy, args=(x_data, y_data, irf, gating_function, histogram_resolution), bounds=bounds, method="SLSQP")
    # res = sp_o.minimize(MLE_model_func, params_copy, args=(x_data, y_data, irf, gating_function, histogram_resolution), bounds=bounds, method="L-BFGS-B")
    # res = let
    #     lb = [b[1] for b in bounds]
    #     ub = [b[2] for b in bounds]
    #     δ = 1e-5
    #     params_copy = clamp.(params_copy, lb .+ δ, ub .- δ)      # <— Ajustement ici
    
    #     o = optimize(
    #         p -> MLE_model_func(p, x_data, y_data, irf, gating_function, histogram_resolution),
    #         lb, ub, params_copy,
    #         Fminbox(LBFGS()),                                      # boîtes + LBFGS
    #         Optim.Options()                                        # options externes
    #     )
    
    #     Dict(
    #         "x"       => Optim.minimizer(o),
    #         "fun"     => Optim.minimum(o),
    #         "success" => Optim.converged(o),
    #     )
    # end
    

    if res["x"] == bounds[1][1] || res["x"] == bounds[1][2] || !res["success"] || (res["fun"] > 5 && sum(y_data) < 10000)
        return [NaN, NaN, NaN]
    end

    return res["x"]
end

function MLE_model_func(params, x_data, y_data, irf, gating_function, histogram_resolution)
    exp_val_c = conv_irf_data(x_data, params, irf, y_data, histogram_resolution=histogram_resolution).*gating_function
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, y_data)
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, exp_val_c)
    ln_data = log.(y_data./exp_val_c)
    poisson_deviance = 2*sum(y_data.*ln_data.-y_data.+exp_val_c)/(histogram_resolution-length(params))

    return poisson_deviance
end

function conv_irf_data(x_data, params, irf, y_data; histogram_resolution=256)
    t_1 = params[1]
    d_0 = params[2]
    y_offset = params[3]
    data = deepcopy(x_data)
    irf_data = deepcopy(irf)
    irf_y_data = irf_shift(irf_data, d_0)
    exp_1 = exp.(-data./t_1)

    for previous_pulse in 1:1:5
        exp_1 += exp.(-(data.+12.5.*previous_pulse)./t_1)
    end
    
    return convolve(irf_y_data, exp_1, y_data, histogram_resolution=histogram_resolution).+y_offset
end

function vec_to_lifetime(x; bin_size=3.9088729e-2, guess=[2.0, 1.0], histogram_resolution=256)
    x_data = vec(collect(bin_size:bin_size:histogram_resolution*bin_size))

    if sum(x) < 100
        return sum(x)
    else
        global gating_function = ones(UInt8, histogram_resolution)
        gating_function[1:round(Int, 0.050980392*histogram_resolution)] .= 0
        gating_function[round(Int, 0.94901962*histogram_resolution):histogram_resolution] .= 0

        data_xy = [x_data, vec(x)]
        tau_fit = MLE_iterative_reconvolution(irf, data_xy, params=guess, gating_function=gating_function, histogram_resolution=histogram_resolution)

        return tau_fit, data_xy
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

    params = ones(3)
    erreur_I = 0

    i = 1
    n = 0
    t_0 = now()
    previous_files = readdir(folderpath)
    dt = 0
    display_period = 1/3

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
                
                params_raw, data = vec_to_lifetime(final_vector, bin_size=irf_bin_size, guess=params[1:2], histogram_resolution=histogram_resolution)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

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

    params = ones(3)
    erreur_I = 0
    display_period = 1/3

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
                params_raw, data = vec_to_lifetime(final_vector, bin_size=irf_bin_size, guess=params[1:2], histogram_resolution=histogram_resolution)
                push!(test_times, (time() - t_test)*1000)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                @inbounds hist[][1:length(data[2])] = data[2]
                @inbounds hist_time[][1:length(data[2])] = data[1]
                @inbounds fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

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
                push!(counts[], sum(@view vectors[idx_start:end, :]) / averaging)
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

function start_analysis()
    GC.gc()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)
    nb_files = length(files)
    checkpoints = round.(Int, nb_files .* (0:0.5:100) ./ 100)

    lifetime = Float64[]
    counts = Float64[]
    times = Float64[]

    params = ones(3)

    i = 1
    cum_time = 0

    println("Starting analysis")
    add_to_monitor("Starting analysis")

    histogram_resolution = 256

    box_progress.width = 16

    beginning = time()

    for file in files
        try
            if !running[]
                add_to_monitor("Aborted")
                break
            end

            vector, histogram_resolution, time = open_SDT_file(joinpath(folderpath, file))
            vector = convert(Vector{Float64}, vector)
            cum_time += time

            params_raw, data = vec_to_lifetime(vector, bin_size=irf_bin_size, guess=params[1:2], histogram_resolution=histogram_resolution)

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

    lines!(axis_counts, times, counts, color=Makie.wong_colors()[1])
    lines!(axis_life, times, lifetime, color=Makie.wong_colors()[1])

    autolimits!(axis_counts)
    autolimits!(axis_life)

    filename = string(now())*".csv"

    open(filename, "w") do io
        writedlm(io, [string(now())])
        writedlm(io, [""])
        writedlm(io, ["Time [s]" "Number of photons" "Lifetime [ns]"], ",")
        writedlm(io, [times counts lifetime], ",")
    end

    println(time() - beginning)
    add_to_monitor("Done")
    
    running[] = false
    pause[] = false

    button_start.label = "START"
    button_stop.label = "CLEAR"
end

# function start_save()
#     running[] = true

#     # rowsize!(fig.layout, 3, Auto())

#     axis_progress = progressgrid[1, 2] = Axis(fig, height = 32, xticksvisible = false, yticksvisible = false, backgroundcolor=:transparent)
#     xlims!(axis_progress, -0.01, 1.01)
#     hidespines!(axis_progress)
#     hidedecorations!(axis_progress)

#     progress = Observable(0.0)
#     line = @lift([Point2f(0, 0), Point2f($progress, 0)])

#     lines!(axis_progress, [Point2f(0, 0), Point2f(1, 0)], linecap=:round, linewidth=14, color=RGBf(0.8, 0.8, 0.8))
#     lines!(axis_progress, line, linecap=:round, linewidth=14, color=RGBf(0.6, 0.6, 0.6))

#     files_label = progressgrid[1, 1] = Label(fig, text="Files : 0/0", width = 120, justification=:left, halign=:left)
#     percentage_label = progressgrid[1, 3] = Label(fig, @lift("$(Int(round(100 * $progress)))%"), width = 120, justification=:right, halign=:right)

#     function saving(lifetimes)
#         println(folderpath)
#         all_files = readdir(folderpath)
#         files = filter(f -> occursin(".sdt", f), all_files)
    
#         lifetime_1 = Any[]
#         lifetime_2 = Any[]
#         lifetime_3 = Any[]
#         counts = Any[]
#         time = Any[]
        
#         params = ones(2*lifetimes+1)
    
#         i = 1
    
#         println("Starting")
#         add_to_monitor("Starting")
        
#         for file in files
#             if running[]
#                 vector, histogram_resolution = open_SDT_file(joinpath(folderpath, file))

#                 params, data = vec_to_lifetime(vector, histogram_resolution=histogram_resolution, guess=params[1:2*lifetimes], tcspc_window_size=12.5, bin_size=irf_bin_size)

#                 push!(lifetime_1, params[1])
#                 push!(counts, sum(vector))
#                 push!(time, i)

#                 if lifetimes == 2
#                     push!(lifetime_2, params[3])
#                     # push!(lifetime_i, (params[2]*params[1]^2 + (1 - params[2])*params[3]^2)/(params[2]*params[1] + (1 - params[2])*params[3]))
#                     # push!(lifetime_m, (params[2]*params[1] + (1 - params[2])*params[3]))
#                 elseif lifetimes == 3
#                     push!(lifetime_2, params[3])
#                     push!(lifetime_3, params[5])
#                     # push!(lifetime_i, (params[2]*params[1]^2 + params[4]*params[3]^2 + (1 - params[2] - params[4])*params[5]^2)/(params[2]*params[1] + params[4]*params[3] + (1 - params[2] - params[4])*params[5]))
#                     # push!(lifetime_m, (params[2]*params[1] + params[4]*params[3] + (1 - params[2] - params[4])*params[5]))
#                 end

#                 progress[] = i/length(files)

#                 files_label.text[] = "Files : $i/$(length(files))"

#                 i += 1
#                 sleep(0.001)
#             end
#         end

#         rowsize!(fig.layout, 3, Fixed(0))

#         delete!(axis_progress)
#         delete!(files_label)
#         delete!(percentage_label)

#         filepath = save_file()

#         if !isempty(filepath)
#             open(filepath, "w") do io
#                 writedlm(io, [string(now())])
#                 writedlm(io, [""])
#                 writedlm(io, [time counts lifetime_1 concentration], ",")
#             end
#         end
#     end

#     # if running[]
#     #     Threads.@spawn saving(menu_lifetimes.i_selected[])
#     # end
# end

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

fig = Figure(size = (1440, 847))

fig[1, 1] = plotgrid = GridLayout()
fig[2, 1] = buttongrid = GridLayout(tellwidth=true)
fig[3, 1] = progressgrid = GridLayout(tellwidth=false)

global axis_hist    = plotgrid[1, 1] = Axis(fig, alignmode = Outside(), title="Histogram",       xlabel = "Time [ns]", ylabel = "Photon counts", titlesize=16)
global axis_counts  = plotgrid[2, 1] = Axis(fig, alignmode = Outside(), title="Photon counts",   xlabel = "Time [s]",  ylabel = "Photon counts", titlesize=16)
global axis_life    = plotgrid[1, 3] = Axis(fig, alignmode = Outside(), title="Lifetime",        xlabel = "Time [s]",  ylabel = "Lifetime [ns]", titlesize=16)
global axis_signal  = plotgrid[2, 3] = Axis(fig, alignmode = Outside(), title="Light intensity", xlabel = "Time [s]",  ylabel = "[%]",           titlesize=16)

ylims!(axis_signal, -5, 105)

label_frequency     = Label(buttongrid[1, 1],       justification=:left, text="Frequency : 0 Hz", width = 120)
label_number        = Label(buttongrid[2, 1],       justification=:left, text="Files : 0",        width = 120)
label_inverted      = Label(buttongrid[3, 1][1, 1], justification=:left, text="Inverted",         width = 70)

checkbox_inverted   = Toggle(buttongrid[3, 1][1, 2], active = true, rimfraction=0.24, buttoncolor=dark_gray, framecolor_active=mid_gray, framecolor_inactive=light_gray)

box_attributes = Dict(
    :color          => :black,
    :cornerradius   => 6,
    :height         => 40,
    :strokewidth    => 0.5,
    :width          => 120,
)

button_attributes = Dict(
    :buttoncolor         => light_gray,
    :buttoncolor_active  => dark_gray,
    :buttoncolor_hover   => mid_gray,
    :cornerradius        => 6,
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
    :cornerradius              => 4,
    :halign              => :right,
    :height                    => 26,
    :placeholder               => " ",
    :textpadding               => (10, 8, 8, 4.5),
    :width                     => 60,
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
    default = "Playback",
    fontsize=20,
    height = 40,
    options=["Real-time", "Playback", "Analysis"],
    selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8),
    textcolor=:black,
    width = 120
)

portname_textbox = Textbox(buttongrid[1:3, 4:5][1,1:2]; merge(textbox_attributes, Dict(:cornerradius => 6, :height => 26, :width => 180))...)
portname_textbox.displayed_string[] = "usbmodem2302"
portname_textbox.stored_string[]    = "usbmodem2302"

label_fps   = Label(buttongrid[1:3, 4:5][2, 1], text=L"\text{fps [Hz] :}", halign=:right)
textbox_fps = Textbox(buttongrid[1:3, 4:5][2, 2]; textbox_attributes...)
textbox_fps.displayed_string[] = "50"
textbox_fps.stored_string[]    = "50"

label_min   = Label(buttongrid[1:3, 4:5][3, 1], text=L"\text{Minimum [ms] :}", halign=:right)
textbox_min = Textbox(buttongrid[1:3, 4:5][3, 2]; textbox_attributes...)
textbox_min.displayed_string[] = "3"
textbox_min.stored_string[]    = "3"

label_max   = Label(buttongrid[1:3, 4:5][4, 1], text=L"\text{Maximum [ms] :}", halign=:right)
textbox_max = Textbox(buttongrid[1:3, 4:5][4, 2]; textbox_attributes...)
textbox_max.displayed_string[] = "14"
textbox_max.stored_string[]    = "14"

button_active  = Button(buttongrid[1, 6]; merge(button_attributes, Dict(:label => "Activate"), Dict(:buttoncolor => mid_gray), Dict(:buttoncolor_hover => dark_gray))...)
button_connect = Button(buttongrid[2, 6]; merge(button_attributes, Dict(:label => "Connect"),  Dict(:buttoncolor => mid_gray), Dict(:buttoncolor_hover => dark_gray))...)

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
    width = 260
)

global monitor = Label(buttongrid[1:3, 7],
    "\n",
    height = 160,
    justification = :left,
    padding = (8, 8, 8, 8),
    tellheight = false,
    tellwidth = false,
    valign = :bottom,
    width = 260
)

params_sliders = SliderGrid(
    buttongrid[1:3, 8],
    NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Time range",   :range => 10:10:3600, :startvalue => 30))),
    NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Averaging",    :range => 1:1:100,    :startvalue => 1))),
    NamedTuple(merge(slider_attributes, Dict(:format => "{:d}",   :label => "Proportional", :range => 0:1:100,    :startvalue => 0))),
    NamedTuple(merge(slider_attributes, Dict(:format => "{:.1f}", :label => "Integral",     :range => 0:1:100,    :startvalue => 0))),
    NamedTuple(merge(slider_attributes, Dict(:format => "{:.2f}", :label => "Setpoint",     :range => 0:0.01:10,  :startvalue => 2.50))),
    width = 360
)

Box(progressgrid[1, 1],
    color = light_gray,
    cornerradius = 4,
    halign = :left,
    height = 8,
    strokevisible=false,
    width = 1353
)

box_progress = Box(progressgrid[1, 1],
    color = mid_gray,
    cornerradius = 4,
    halign = :left,
    height = 8,
    strokevisible=false,
    width = 0
)

label_progress = Label(progressgrid[1, 2], text="0 %", halign=:right, tellwidth = false)

colgap!(buttongrid, 3, 50)
colgap!(buttongrid, 5, 50)
colsize!(plotgrid, 2, 10)
rowgap!(buttongrid, 1, 20)
rowgap!(buttongrid, 2, 20)

global irf = get_irf()
global irf_bin_size = get_irf_bin_size(irf)
global folderpath = get_folder_path("FLIM")
global serial = nothing
global running = Observable(false)
global pause = Observable(false)
global active = Observable(false)
global connected = Observable(false)
global reset_I = Observable(false)

setpoint = hlines!(axis_life, params_sliders.sliders[5].value[], color=Makie.wong_colors()[2])

on(button_start.clicks) do _
    if !running[]
        running[] = true

        button_start.label = "PAUSE"
        button_stop.label = "STOP"

        if menu_mode.selection[] == "Playback"
            @async start_playback()
        elseif menu_mode.selection[] == "Analysis"
            @async start_analysis()
        else
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
    global irf_bin_size = get_irf_bin_size(irf)
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

on(params_sliders.sliders[5].value) do value
    setpoint.args[1][] = value
end

fig
