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

dark_gray = RGBf(0.4, 0.4, 0.4)
mid_gray = RGBf(0.6, 0.6, 0.6)
light_gray = RGBf(0.8, 0.8, 0.8)
off_white = RGBf(0.9, 0.9, 0.9)

try
    LibSerialPort.close(serial)
catch
    # Nothing to close
    println("No serial connection to close.")
end

sp_o = pyimport("scipy.optimize")

# functions ignoring NaNs
nanmean(x) = mean(filter(!isnan, x))
nanmean(x, dims) = mapslices(nanmean, x, dims=dims)

nanstd(x) = std(filter(!isnan, x))
nanstd(x, dims) = mapslices(nanstd, x, dims=dims)

nanmax(x) = maximum(filter(!isnan, x))
nanmin(x) = minimum(filter(!isnan, x))

is_zero(x) = x == 0
smaller_or_eq_zero(x) = x <= 0

function get_irf(; channel=1)
    """
    Loads an .sdt IRF file. If irf_filepath.txt exists, load the IRF from this filepath, otherwise, open a file dialog create irf_filepath.txt
    """

    println("Getting IRF file...")

    if isfile("irf_filepath.txt")
        filepath = open(f->read(f, String), "irf_filepath.txt")
    else
        filepath = pick_file()

        open("irf_filepath.txt", "w") do io
            write(io, filepath)
        end
    end
    sdt = pyimport("sdtfile") # use python package sdtfile
    file = sdt.SdtFile(filepath)
    data = convert.(Int64, file.data[channel])
    median_irf = round(median(data)) # get median
    data = data.-round(Int, 1.25*median_irf) # remove noise from the IRF
    data[data .<= 0] .= 0 # set negative values to 0
    times = convert.(Float64, vec(file.times[1]).*1e9) # convert from seconds to nanoseconds
    data = [times, vec(data)]

    println("Done")
    println()

    return data
end

function get_new_irf(; channel=1)
    """
    Loads an .sdt IRF file. If irf_filepath.txt exists, load the IRF from this filepath, otherwise, open a file dialog create irf_filepath.txt
    """

    println("Getting IRF file...")

    filepath = pick_file()

    open("irf_filepath.txt", "w") do io
        write(io, filepath)
    end

    sdt = pyimport("sdtfile") # use python package sdtfile
    file = sdt.SdtFile(filepath)
    data = convert.(Int64, file.data[channel])
    median_irf = round(median(data)) # get median
    data = data.-round(Int, 1.25*median_irf) # remove noise from the IRF
    data[data .<= 0] .= 0 # set negative values to 0
    times = convert.(Float64, vec(file.times[1]).*1e9) # convert from seconds to nanoseconds
    data = [times, vec(data)]

    println("Done")
    println()

    return data
end

function get_irf_bin_size(irf)
    """
    Gets the bin size in nanoseconds from the IRF.
    """
    # irf = get_irf()
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
    """
    Convolves the IRF with fluorescence decay
    """
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
    """
    Shifts the IRF by a number of bins
    REF: Enderlein, 1997
    """
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

function MLE_iterative_reconvolution(irf, data_xy; params, gating_function=ones(UInt8, 256), histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
    """
    Curve fitting of fluorescence decay using Maximum Likelihood Estimation
    """
    #global x_data = vec(deepcopy(data_xy[1]))
    #global y_data = vec(deepcopy(data_xy[2])).*gating_function
    #global irf = deepcopy(irf)
    #global histogram_resolution_fit = histogram_resolution
    x_data = vec(data_xy[1])
    y_data = vec(data_xy[2]).*gating_function
    irf = deepcopy(irf)
    #argmax_y_data = argmax(y_data)

    number_of_lifetimes = convert(Int, round((length(params))/2))
    params_copy = deepcopy(params)
    min_offset = floor(mean(y_data[1:argmax(y_data)-10]))+ceil(std(y_data[1:argmax(y_data)-10])) # get height before the peak of the decay
    if min_offset < 5
        min_offset = 5
    end
    push!(params_copy, min_offset)
    if number_of_lifetimes == 1
        bounds = [(0.1, 12.5), (-16, 64), (0, min_offset)]
        res = sp_o.minimize(MLE_model_func, params_copy, args=(x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period), bounds=bounds, method="L-BFGS-B")
        if res["x"] == bounds[1][1] || res["x"] == bounds[1][2]
            return [NaN, NaN, NaN]
        end
        #res = sp_o.minimize(MLE_model_func, params_copy, tol=1e-3)
    elseif number_of_lifetimes == 2
        bounds = [(0.1, 12.5), (0.0, 1.0), (0.1, 12.5), (-30, histogram_resolution), (0, min_offset)]
        res = sp_o.minimize(MLE_model_func, params_copy, args=(x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period), bounds=bounds, method="L-BFGS-B", options=Dict("maxiter" => 100000))
    elseif number_of_lifetimes == 3
        bounds = [(0.1, 12.5), (0.0, 1.0), (0.1, 12.5), (0.0, 1.0), (0.1, 12.5), (-30, histogram_resolution), (0, min_offset)]
        res = sp_o.minimize(MLE_model_func, params_copy, args=(x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period), bounds=bounds, method="L-BFGS-B", options=Dict("maxiter" => 100000))
    end
    if !res["success"]
        if number_of_lifetimes == 1
            return [NaN, NaN, NaN]
        elseif number_of_lifetimes == 2
            return [NaN, NaN, NaN, NaN, NaN]
        elseif number_of_lifetimes == 3
            return [NaN, NaN, NaN, NaN, NaN, NaN, NaN]
        end
    end
    #fit = conv_irf_data(x_data, res["x"])
    #chi2 = sum((y_data.-fit).^2 ./ (histogram_resolution.*y_data.+1))
    #println(chi2)
    #if chi2 > 3
    if res["fun"] > 5 && sum(y_data) < 10000
        println("Ignoring bad fit with optimization result: ", res["fun"], " found values: ", res["x"])
        if number_of_lifetimes == 1
            return [NaN, NaN, NaN]
        elseif number_of_lifetimes == 2
            return [NaN, NaN, NaN, NaN, NaN]
        elseif number_of_lifetimes == 3
            return [NaN, NaN, NaN, NaN, NaN, NaN, NaN]
        end
    end
    if number_of_lifetimes == 2 && res["x"][3] > res["x"][1]
        t_1 = res["x"][1]
        t_2 = res["x"][3]
        res["x"][1] = t_2
        res["x"][3] = t_1
        res["x"][2] = 1-res["x"][2]
    end
    return res["x"]
end

function MLE_model_func(params, x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period)
    """
    Model function for MLE optimizer.
    # REF: Bajzer, et al. 1991
    # REF: Maus, et al. 2001
    """
    number_of_lifetimes = convert(Int, round((length(params)-1)/2))
    if number_of_lifetimes == 1 || number_of_lifetimes == 0
        exp_val_c = conv_irf_data(x_data, (params[1], params[2], params[3]), irf, y_data, histogram_resolution=histogram_resolution, number_of_previous_pulses=number_of_previous_pulses, laser_pulse_period=laser_pulse_period).*gating_function
    elseif number_of_lifetimes == 2
        exp_val_c = conv_irf_data(x_data, (params[1], params[2], params[3], params[4], params[5]), irf, y_data, histogram_resolution=histogram_resolution, number_of_previous_pulses=number_of_previous_pulses, laser_pulse_period=laser_pulse_period).*gating_function
    elseif number_of_lifetimes == 3
        exp_val_c = conv_irf_data(x_data, (params[1], params[2], params[3], params[4], params[5], params[6], params[7]), irf, y_data, histogram_resolution=histogram_resolution, number_of_previous_pulses=number_of_previous_pulses, laser_pulse_period=laser_pulse_period).*gating_function
    end
 
    output_resolution = histogram_resolution
    #downsample_ratio = Int(histogram_resolution/output_resolution)
    #new_y_data = zeros(output_resolution)
    #new_decay = zeros(output_resolution)
    #for (downsampled_i, i) in enumerate(1:downsample_ratio:histogram_resolution)
    #    new_y_data[downsampled_i] = sum(y_data[i:i+downsample_ratio-1])
    #    new_decay[downsampled_i] = sum(exp_val_c[i:i+downsample_ratio-1])
    #    #new_data[i÷downsample_ratio] = Int(floor(sum(data[i-downsample_ratio+1:i])/downsample_ratio))
    #end

    #Nbins = 64
    #new_y = zeros(Nbins)
    #new_decay = zeros(Nbins)
    #for i in 1:1:Nbins
    #    new_y[i] = y_data[floor(Int, 2*(i-1)+j)]+y_data[floor(Int, 2*i)]
    #    new_decay[i] = exp_val_c[floor(Int, 2*i-1)]+exp_val_c[floor(Int, 2*i)]
    #end
    #y_data = new_y_data
    #exp_val_c = new_decay
    #output_resolution=256

    #y_data[y_data .<= 0.0] .= 1e-256
    #exp_val_c[exp_val_c .<= 0.0] .= 1e-256
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, y_data) # replace values that are smaller or equal to zero by 1e-256 (very small but not zero)
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, exp_val_c) # same as above. This is because log(0) is undefined
    ln_data = log.(y_data./exp_val_c)
    poisson_deviance = 2*sum(y_data.*ln_data.-y_data.+exp_val_c)/(output_resolution-length(params))
    #poisson_deviance = sum((y_data.*ln_data))/(256-length(params))
    #println("Params: ", params, " poisson_deviance: ", poisson_deviance)
    #fig = CairoMakie.Figure()
    #ax = CairoMakie.Axis(fig[1, 1])
    #lines!(ax, 1:1:output_resolution, y_data[1:output_resolution], color=:grey)
    #lines!(ax, 1:1:output_resolution, exp_val_c[1:output_resolution], color=:red)
    #display(fig)
    return poisson_deviance
end

function conv_irf_data(x_data, params, irf, y_data; histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
    """
    Creates decay function from parameters
    """
    number_of_params = length(params)
    if number_of_params == 3
        t_1 = params[1]
        d_0 = params[2]
        y_offset = params[3]
        data = deepcopy(x_data)
        irf_data = deepcopy(irf)
        irf_y_data = irf_shift(irf_data, d_0)
        exp_1 = exp.(-data./t_1)
        for previous_pulse in 1:1:number_of_previous_pulses
            exp_1 += exp.(-(data.+laser_pulse_period.*previous_pulse)./t_1)
        end
        decay = convolve(irf_y_data, exp_1, y_data, histogram_resolution=histogram_resolution).+y_offset
        return decay
    elseif number_of_params == 5
        t_1 = params[1]
        A_1 = params[2]
        t_2 = params[3]
        d_0 = params[4]
        y_offset = params[5]
        data = deepcopy(x_data)
        irf_data = deepcopy(irf)
        irf_y_data = irf_shift(irf_data, d_0)
        exp_1 = A_1.*exp.(-data./t_1)
        exp_2 = (1.0-A_1).*exp.(-data./t_2)
        for previous_pulse in 1:1:number_of_previous_pulses
            exp_1 += A_1.*exp.(-(data.+laser_pulse_period.*previous_pulse)./t_1)
            exp_2 += (1.0-A_1).*exp.(-(data.+laser_pulse_period.*previous_pulse)./t_2)
        end
        decay = convolve(irf_y_data, (exp_1.+exp_2), y_data, histogram_resolution=histogram_resolution).+y_offset
        return decay
    elseif number_of_params == 7
        t_1 = params[1]
        A_1 = params[2]
        t_2 = params[3]
        A_2 = params[4]
        t_3 = params[5]
        d_0 = params[6]
        y_offset = params[7]
        data = deepcopy(x_data)
        irf_data = deepcopy(irf)
        irf_y_data = irf_shift(irf_data, d_0)
        exp_1 = A_1.*exp.(-data./t_1)
        exp_2 = A_2.*exp.(-data./t_2)
        exp_3 = (1.0-A_1-A_2).*exp.(-data./t_3)
        for previous_pulse in 1:1:number_of_previous_pulses
            exp_1 += A_1.*exp.(-(data.+laser_pulse_period.*previous_pulse)./t_1)
            exp_2 += A_2.*exp.(-(data.+laser_pulse_period.*previous_pulse)./t_1)
            exp_3 += (1.0-A_1-A_2).*exp.(-(data.+laser_pulse_period.*previous_pulse)./t_3)
        end
        decay = convolve(irf_y_data, (exp_1.+exp_2.+exp_3), y_data, histogram_resolution=histogram_resolution).+y_offset
        return decay
    end
end

function LS_model_func_1_lt(data, t_1, d_0, y_offset)
    """
    Model function for LS optimizer for 1 lifetime
    """
    x_data = data[1]
    irf = data[2]
    y_data = data[3]
    A = conv_irf_data(x_data, [t_1, d_0, y_offset], irf, y_data, histogram_resolution=histogram_resolution_fit).*gating_function
    return A
end

function LS_model_func_2_lt(data, t_1, A_1, t_2, d_0, y_offset)
    """
    Model function for LS optimizer for 2 lifetimes
    """
    x_data = data[1]
    irf = data[2]
    y_data = data[3]
    A = conv_irf_data(x_data, [t_1, A_1, t_2, d_0, y_offset], irf, y_data, histogram_resolution=histogram_resolution_fit).*gating_function
    return A
end

function LS_model_func_3_lt(data, t_1, A_1, t_2, A_2, t_3, d_0, y_offset)
    """
    Model function for LS optimizer for 3 lifetimes
    """
    x_data = x_data[1]
    irf = x_data[2]
    y_data = data[3]
    A = conv_irf_data(x_data, [t_1, A_1, t_2, A_2, t_3, d_0, y_offset], irf, y_data, histogram_resolution=histogram_resolution_fit).*gating_function
    return A
end

function LS_iterative_reconvolution(irf, data_xy; guess=[2.0, 1.0], bounds=([0.0, -30], [12.5, 256]), irf_shift=0, gating_function=ones(UInt8, 256), histogram_resolution=256)
    """
    Curve fitting of fluorescence decay using Least Squares
    """
    x_data = vec(deepcopy(data_xy[1]))
    y_data = vec(deepcopy(data_xy[2])).*gating_function
    #global y_data = vec(deepcopy(data_xy[2])).*gating_function
    #global irf = deepcopy(irf)
    #global histogram_resolution_fit = histogram_resolution
    number_of_lifetimes = convert(Int, round((length(guess))/2))
    data = [x_data, irf, y_data]

    min_offset = floor(mean(y_data[1:argmax(y_data)-10]))+ceil(std(y_data[1:argmax(y_data)-10]))
    guess_copy = deepcopy(guess)
    if min_offset < 5
        min_offset = 5
    end
    push!(guess_copy, min_offset)
    if number_of_lifetimes == 1
        try
            res, cov = sp_o.curve_fit(LS_model_func_1_lt, data, y_data, guess_copy, method="lm")
        catch e
            res = [NaN, 0.0, 0.0]
        end
        #res, cov = sp_o.curve_fit(LS_model_func_1_lt, x_data, y_data, guess_copy, method="lm")
    elseif number_of_lifetimes == 2
        bounds = ([0.0, 0.0, 0.0, -30, 0.0], [12.5, 1.0, 12.5, histogram_resolution, min_offset])
        try
            res, cov = sp_o.curve_fit(LS_model_func_2_lt, data, y_data, guess_copy, bounds=bounds) 
        catch
            res = [NaN, NaN, NaN, 0.0, 0.0]
        end
        if res[3] > res[1]
            tau_1 = res[1]
            tau_2 = res[3]
            res[1] = tau_2
            res[3] = tau_1
            res[2] = 1-res[2]
        end
    elseif number_of_lifetimes == 3
        bounds = ([0.0, 0.0, 0.0, 0.0, 0.0, -30, 0.0], [12.5, 1.0, 12.5, 1.0, 12.5, histogram_resolution, min_offset])
        res, cov = sp_o.curve_fit(LS_model_func_3_lt, data, y_data, guess_copy, bounds=bounds)
    end
    return res
end

function vec_to_lifetime(x; bin_size=3.9088729e-2, hist_size_threshold=100, method="MLE", broadcast_call=false, guess=[2.0, 1.0], laser_pulse_period=12.5, histogram_resolution=256, number_of_previous_pulses=5, tcspc_window_size=10, tac_low_cut=5.0980392, tac_high_cut=94.901962)
    """
    Calculate fluorescence lifetime from vector of counts
    """
    if broadcast_call
        x = x[2]
    end
    total_channels = Int(round(laser_pulse_period*histogram_resolution/tcspc_window_size))
    append!(vec(x), zeros(total_channels-histogram_resolution))
    x_data = vec(collect(bin_size:bin_size:total_channels*bin_size))
    if sum(x) < hist_size_threshold
        if method == "Intensity"
            # Instead of fluorescence lifetime, calculate number of photons in histogram
            return sum(x)
        #elseif method == "Estimate"
        #    return ([0], [x_data, vec(x)])
        end
        if broadcast_call
            return NaN
        else
            return NaN, [x_data, vec(x)]
        end
    else
        if method == "Intensity"
            return sum(x)
        #elseif method == "Estimate"
        #    return ([lifetime_estimate(x[1:histogram_resolution], bin_size=bin_size)], [x_data[1:histogram_resolution], vec(x)[1:histogram_resolution]])
        #elseif method == "CNN"
        #    return CNN(reshape(x[1:histogram_resolution], 256, 1, 1))
        end
        global gating_function = ones(UInt8, total_channels)
        gating_function[1:round(Int, tac_low_cut/100*histogram_resolution)] .= 0
        gating_function[round(Int, tac_high_cut/100*histogram_resolution):total_channels] .= 0
        if method == "MLE"
            data_xy = [x_data, vec(x)]
            tau_fit = MLE_iterative_reconvolution(irf, data_xy, params=guess, gating_function=gating_function, histogram_resolution=histogram_resolution, number_of_previous_pulses=number_of_previous_pulses, laser_pulse_period=laser_pulse_period)
            if broadcast_call
                return tau_fit[1]
            end
        elseif method == "LS"
            data_xy = [x_data, vec(x)]
            tau_fit = LS_iterative_reconvolution(irf, data_xy, guess=guess, gating_function=gating_function, histogram_resolution=histogram_resolution)
            if broadcast_call
                return tau_fit[1]
            end
        end
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

function optimized_sum(file, histogram_resolution)
    len = length(file)
    result = zeros(Int, histogram_resolution)

    @threads for i in 1:histogram_resolution
        idx_start = (i-1) * (len ÷ histogram_resolution) + 1
        idx_end = i * (len ÷ histogram_resolution)
        result[i] = sum(file[idx_start:idx_end])
    end

    return result
end

function open_SDT_file(filepath)
    open(filepath,"r") do io
        seek(io, 14)
        header = read(io, 12)

        seek(io, header[12]*16^2 + header[11]+82)
        infos = read(io, 2)
        histogram_resolution = infos[2]*16^2 + infos[1]

        if header[8]*16^2+header[7] == 512 || header[8]*16^2+header[7] == 128
            seek(io, header[2]*16^2+header[1]+22)

            vector = read(io)
            new_vector = UInt16[]

            for i in 1:2:length(vector)
                push!(new_vector, vector[i+1]*16^2 + vector[i])
            end

            return new_vector, histogram_resolution
        else
            seek(io, header[2]*16^2+header[1]+2)

            file_info = read(io, 20)

            shift_1 = file_info[8]*16^6 + file_info[7]*16^4 + file_info[6]*16^2 + file_info[5] - file_info[2]*16^2 - file_info[1]
            shift_2 = file_info[20]*16^6 + file_info[19]*16^4 + file_info[18]*16^2 + file_info[17]

            buffer = IOBuffer(read(io, shift_1))
            zip_file = ZipFile.Reader(buffer).files[1]

            file_data = Vector{UInt8}(undef, shift_2)
            read!(zip_file, file_data)

            vector = reshape_to_vec(reinterpret(UInt16, file_data), histogram_resolution)
    
            return vector, histogram_resolution
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
        ser = LibSerialPort.open("/dev/tty." * portName, 115200)
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

function clear()
    send_command(serial, "C Z\nC R\n")
end

function start_real_time()
    println(folderpath)
    if !initialized[]
        println("Initializing")

        hist = Observable(zeros(256))
        fit = Observable(zeros(256))
        hist_time = Observable(zeros(256))
        
        lifetime = Observable(Float64[])
        counts = Observable(Float64[])
        concentration = Observable(Float64[])
        global signal = Observable(Float64[])
        time = Observable(Float64[])
        
        vectors = zeros(100, 256)
        
        barplot!(axis_hist, hist_time, hist, gap=0, color=Makie.wong_colors()[1])
        lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
        lines!(axis_counts, time, counts, color=Makie.wong_colors()[1])
        lines!(axis_conc, time, concentration, color=Makie.wong_colors()[1])
        lines!(axis_signal, time, signal, color=Makie.wong_colors()[1])
        
        previous_time = now()
        current_time = now()
        initialized[] = true
        params = ones(3)
        frequency = 1
        erreur_I = 0
        
        i = 1
        t_0 = now()
        previous_files = readdir(folderpath)
    end

    println("Starting real time")
    println(folderpath)

    while running[]
        try
            current_files = readdir(folderpath)
            new_files = setdiff(current_files, previous_files)
            sleep(0.001)
            
            if !isempty(new_files)
                if i == 1
                    println("Files detected")
                end

                println(new_files[end])

                vector, histogram_resolution = open_SDT_file(joinpath(folderpath, new_files[end]))
                
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

                avg_vector = sum(vectors[end-params_sliders.sliders[2].value[]+1:end, 1:histogram_resolution], dims=1)[1, :] / params_sliders.sliders[2].value[]
                
                params_raw, data = vec_to_lifetime(avg_vector, histogram_resolution=histogram_resolution, guess=params[1:2], tcspc_window_size=12.5, bin_size=irf_bin_size)
                if !isnan(params_raw[1])
                    params = params_raw
                end

                hist[][1:length(data[2])] = data[2]
                hist_time[][1:length(data[2])] = data[1]
                fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

                erreur = params_sliders.sliders[5].value[] - (parse(Float64, textbox_tau_0.stored_string[])/params[1] - 1) / parse(Float64, textbox_K.stored_string[])*1000
                dt = 0

                try
                    dt = (now() - t_0)/Millisecond(1)/1000 - time[][end]
                catch e
                    println(e)
                end

                erreur_I += erreur * dt

                if checkbox_inverted.active[]
                    signal_value = clamp(-(params_sliders.sliders[3].value[]*erreur + params_sliders.sliders[4].value[]*erreur_I), 0, 100)
                else
                    signal_value = clamp(params_sliders.sliders[3].value[]*erreur + params_sliders.sliders[4].value[]*erreur_I, 0, 100)
                end

                push!(lifetime[], params[1])
                push!(counts[], sum(vectors[end-params_sliders.sliders[2].value[]+1:end, :]) / params_sliders.sliders[2].value[])
                push!(concentration[], (parse(Float64, textbox_tau_0.stored_string[])/params[1] - 1) / parse(Float64, textbox_K.stored_string[])*1000)
                push!(signal[], signal_value)
                push!(time[], (now() - t_0)/Millisecond(1)/1000)

                # send_command(serialConn, "A 0 AO 1 $(signal[][end]*1000)\n")
                
                if i == 1
                    t_0 = now()
                end

                if mod1(i, 8) == 8
                    current_time = now()
                    frequency = round(1/((current_time - previous_time)/Millisecond(1))*1000*8, sigdigits=2)
                    label_frequency.text = "Frequency : $(string(frequency)) Hz"
                    previous_time = current_time
                end

                if mod1(i, ceil(frequency/16)) == ceil(frequency/16)
                    notify(hist_time)
                    notify(time)
                    
                    autolimits!(axis_hist)
                    autolimits!(axis_counts)
                    autolimits!(axis_conc)
                    # autolimits!(axis_signal)
                    
                    xlims!(axis_counts, max(0, time[][end] - params_sliders.sliders[1].value[] + 0.1), max(params_sliders.sliders[1].value[], time[][end] + 0.1))
                    xlims!(axis_conc, max(0, time[][end] - params_sliders.sliders[1].value[] + 0.1), max(params_sliders.sliders[1].value[], time[][end] + 0.1))
                    xlims!(axis_signal, max(0, time[][end] - params_sliders.sliders[1].value[] + 0.1), max(params_sliders.sliders[1].value[], time[][end] + 0.1))
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
end

# function start_real_time(lifetimes)
#     if !initialized[]
#         println("Initializing")

#         hist = Observable(zeros(256))
#         fit = Observable(zeros(256))
#         hist_time = Observable(zeros(256))

#         lifetime_1 = Observable(Float64[])
#         lifetime_2 = Observable(Float64[])
#         lifetime_3 = Observable(Float64[])
#         lifetime_i = Observable(Float64[])
#         lifetime_m = Observable(Float64[])
#         counts = Observable(Float64[])
#         concentration = Observable(Float64[])

#         time = Observable(Float64[])

#         vectors = zeros(100, 256)

#         barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
#         lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
#         # lines!(axis_life, time, lifetime_1, color=Makie.wong_colors()[1], label = L"\tau_1")
#         lines!(axis_counts, time, counts, color=Makie.wong_colors()[1])
#         lines!(axis_conc, time, concentration, color=Makie.wong_colors()[menu_main_lifetime.i_selected[]])

#         if lifetimes == 2
#             # lines!(axis_life, time, lifetime_2, color=Makie.wong_colors()[2], label = L"\tau_2")
#             # lines!(axis_life, time, lifetime_i, color=Makie.wong_colors()[3], label = L"\tau_i")
#             # lines!(axis_life, time, lifetime_m, color=Makie.wong_colors()[4], label = L"\tau_m")
#         elseif lifetimes == 3
#             # lines!(axis_life, time, lifetime_2, color=Makie.wong_colors()[2], label = L"\tau_2")
#             # lines!(axis_life, time, lifetime_3, color=Makie.wong_colors()[3], label = L"\tau_3")
#             # lines!(axis_life, time, lifetime_i, color=Makie.wong_colors()[4], label = L"\tau_i")
#             # lines!(axis_life, time, lifetime_m, color=Makie.wong_colors()[5], label = L"\tau_m")
#         end

#         # global legend = axislegend(axis_life)

#         previous_time = now()
#         current_time = now()
#         initialized[] = true
#         params = ones(2*lifetimes+1)
#         frequency = 1

#         previous_files = readdir(folderpath)

#         i = 1
#         t_0 = now()
#     end

#     println(folderpath)
#     println("Waiting")

#     while running[]
#         current_files = readdir(folderpath)

#         new_files = setdiff(current_files, previous_files)

#         sleep(0.001)

#         if !isempty(new_files)
#             if i == 1
#                 println("Files detected")
#             end

#             vector, histogram_resolution = open_SDT_file(joinpath(folderpath, new_files[end]))

#             if i == 1
#                 for n in 1:100
#                     vectors[n,begin:histogram_resolution] = vector
#                 end
#             else
#                 for n in 1:99
#                     vectors[n,begin:histogram_resolution] = vectors[n+1,begin:histogram_resolution]
#                 end

#                 vectors[end,begin:histogram_resolution] = vector
#             end

#             # vector = sum(vectors, dims=1)[1,:]/100

#             params_raw, data = vec_to_lifetime(sum(vectors[end-slider_averaging.value[]+1:end,begin:histogram_resolution], dims=1)[1,:]/slider_averaging.value[], histogram_resolution=histogram_resolution, guess=params[1:2*lifetimes], tcspc_window_size=12.5, bin_size=irf_bin_size)

#             if !isnan(params_raw[1])
#                 params = params_raw
#             end

#             hist[][1:length(data[2])] = data[2]
#             hist_time[][1:length(data[2])] = data[1]
#             fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

#             push!(lifetime_1[], params[1])
#             push!(counts[], sum(vectors[end-slider_averaging.value[]+1:end,:])/slider_averaging.value[])
#             push!(concentration[], 1/(1/parse(Float64, textbox_tau_0.stored_string[]) + parse(Float64, textbox_K.stored_string[])*params[2*menu_main_lifetime.i_selected[] - 1]))

#             if lifetimes == 2
#                 push!(lifetime_2[], params[3])
#             elseif lifetimes == 3
#                 push!(lifetime_2[], params[3])
#                 push!(lifetime_3[], params[5])
#             end

#             push!(time[], (now() - t_0)/Millisecond(1)/1000)

#             if i == 1
#                 t_0 = now()
#             end

#             if mod1(i, 8) == 8
#                 current_time = now()
#                 label_frequency.text = "Frequency : $(string(round(1/((current_time - previous_time)/Millisecond(1))*1000*8, sigdigits=2))) Hz"
#                 previous_time = current_time

#                 # label_number.text = "Files : $(round(UInt64, i))"
#             end

#             if mod1(i, ceil(frequency/16)) == ceil(frequency/16)
#                 notify(hist_time)
#                 notify(time)

#                 autolimits!(axis_hist)
#                 # autolimits!(axis_life)
#                 autolimits!(axis_counts)
#                 autolimits!(axis_conc)
                
#                 # xlims!(axis_life, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
#                 xlims!(axis_counts, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
#                 xlims!(axis_conc, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
#             end

#             previous_files = current_files

#             i += 1
#         end

#         while pause[]
#             sleep(0.1)
#         end
#     end
# end

function start_playback()
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)

    # if !initialized[]
    println("Initializing")

    hist = Observable(zeros(256))
    fit = Observable(zeros(256))
    hist_time = Observable(zeros(256))

    lifetime = Observable(Float64[])
    counts = Observable(Float64[])
    concentration = Observable(Float64[])
    global signal = Observable(Float64[])

    time = Observable(Float64[])

    vectors = zeros(100, 256)

    barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
    lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
    lines!(axis_counts, time, counts, color=Makie.wong_colors()[1])
    lines!(axis_conc, time, concentration, color=Makie.wong_colors()[1])
    lines!(axis_signal, time, signal, color=Makie.wong_colors()[1])

    previous_time = now()
    current_time = now()
    initialized[] = true
    params = ones(3)
    frequency = 1
    erreur_I = 0

    i = 1
    t_0 = now()
    # end

    println("Starting")
    add_to_monitor("Starting")

    while running[]
        sleep(0.001)
        for file in files
            if running[]
                vector, histogram_resolution = open_SDT_file(joinpath(folderpath, file))

                if i == 1
                    for n in 1:100
                        vectors[n,begin:histogram_resolution] = vector
                    end
                else
                    for n in 1:99
                        vectors[n,begin:histogram_resolution] = vectors[n+1,begin:histogram_resolution]
                    end

                    vectors[end,begin:histogram_resolution] = vector
                end
                
                final_vector = sum(vectors[end-params_sliders.sliders[2].value[]+1:end,begin:histogram_resolution], dims=1)[1,:]/params_sliders.sliders[2].value[]

                params_raw, data = vec_to_lifetime(final_vector, histogram_resolution=histogram_resolution, guess=params[1:2], tcspc_window_size=12.5, bin_size=irf_bin_size)

                if !isnan(params_raw[1])
                    params = params_raw
                end

                hist[][1:length(data[2])] = data[2]
                hist_time[][1:length(data[2])] = data[1]
                fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

                if i == 1
                    t_0 = now()
                end

                erreur = params_sliders.sliders[5].value[] - (parse(Float64, textbox_tau_0.stored_string[])/params[1] - 1)/parse(Float64, textbox_K.stored_string[])*1000
                dt = 0

                try
                    dt = (now() - t_0)/Millisecond(1)/1000 - time[][end]
                catch e
                    println(e)
                end

                erreur_I += erreur*dt

                if i == 1
                    erreur_I = 0
                end

                if checkbox_inverted.active[]
                    signal_value = clamp(-(params_sliders.sliders[3].value[]*erreur + params_sliders.sliders[4].value[]*erreur_I), 0, 100)
                else
                    signal_value = clamp(params_sliders.sliders[3].value[]*erreur + params_sliders.sliders[4].value[]*erreur_I, 0, 100)
                end

                push!(lifetime[], params[1])
                push!(counts[], sum(vectors[end-params_sliders.sliders[2].value[]+1:end,:])/params_sliders.sliders[2].value[])
                push!(concentration[], (parse(Float64, textbox_tau_0.stored_string[])/params[1] - 1)/parse(Float64, textbox_K.stored_string[])*1000)
                push!(signal[], signal_value)
                push!(time[], (now() - t_0)/Millisecond(1)/1000)

                # send_command(serial, "A 0 AO 1 $(signal_value*50)\n")

                # Box_function(signal_value)

                if mod1(i, 16) == 16
                    current_time = now()
                    frequency = round(1/((current_time - previous_time)/Millisecond(1))*1000*16, sigdigits=2)
                    label_frequency.text = "Frequency : $frequency Hz"
                    previous_time = current_time
        
                    label_number.text = "Files : $(round(UInt64, i))"
                end

                if mod1(i, ceil(frequency/16)) == ceil(frequency/16)
                    notify(hist_time)
                    notify(time)

                    autolimits!(axis_hist)
                    autolimits!(axis_counts)
                    autolimits!(axis_conc)
                    # autolimits!(axis_signal)
                    
                    xlims!(axis_counts, max(0, time[][end] - params_sliders.sliders[1].value[] + 0.1), max(params_sliders.sliders[1].value[], time[][end] + 0.1))
                    xlims!(axis_conc, max(0, time[][end] - params_sliders.sliders[1].value[] + 0.1), max(params_sliders.sliders[1].value[], time[][end] + 0.1))
                    xlims!(axis_signal, max(0, time[][end] - params_sliders.sliders[1].value[] + 0.1), max(params_sliders.sliders[1].value[], time[][end] + 0.1))
                end

                i += 1
                sleep(0.001)

                while pause[]
                    sleep(0.1)
                end
            end
        end
    end
end

function start_save()
    running[] = true

    # rowsize!(fig.layout, 3, Auto())

    axis_progress = progressgrid[1, 2] = Axis(fig, height = 32, xticksvisible = false, yticksvisible = false, backgroundcolor=:transparent)
    xlims!(axis_progress, -0.01, 1.01)
    hidespines!(axis_progress)
    hidedecorations!(axis_progress)

    progress = Observable(0.0)
    line = @lift([Point2f(0, 0), Point2f($progress, 0)])

    lines!(axis_progress, [Point2f(0, 0), Point2f(1, 0)], linecap=:round, linewidth=14, color=RGBf(0.8, 0.8, 0.8))
    lines!(axis_progress, line, linecap=:round, linewidth=14, color=RGBf(0.6, 0.6, 0.6))

    files_label = progressgrid[1, 1] = Label(fig, text="Files : 0/0", width = 120, justification=:left, halign=:left)
    percentage_label = progressgrid[1, 3] = Label(fig, @lift("$(Int(round(100 * $progress)))%"), width = 120, justification=:right, halign=:right)

    function saving(lifetimes)
        println(folderpath)
        all_files = readdir(folderpath)
        files = filter(f -> occursin(".sdt", f), all_files)
    
        lifetime_1 = Any[]
        lifetime_2 = Any[]
        lifetime_3 = Any[]
        counts = Any[]
        time = Any[]
        
        params = ones(2*lifetimes+1)
    
        i = 1
    
        println("Starting")
        add_to_monitor("Starting")
        
        for file in files
            if running[]
                vector, histogram_resolution = open_SDT_file(joinpath(folderpath, file))

                params, data = vec_to_lifetime(vector, histogram_resolution=histogram_resolution, guess=params[1:2*lifetimes], tcspc_window_size=12.5, bin_size=irf_bin_size)

                push!(lifetime_1, params[1])
                push!(counts, sum(vector))
                push!(time, i)

                if lifetimes == 2
                    push!(lifetime_2, params[3])
                    # push!(lifetime_i, (params[2]*params[1]^2 + (1 - params[2])*params[3]^2)/(params[2]*params[1] + (1 - params[2])*params[3]))
                    # push!(lifetime_m, (params[2]*params[1] + (1 - params[2])*params[3]))
                elseif lifetimes == 3
                    push!(lifetime_2, params[3])
                    push!(lifetime_3, params[5])
                    # push!(lifetime_i, (params[2]*params[1]^2 + params[4]*params[3]^2 + (1 - params[2] - params[4])*params[5]^2)/(params[2]*params[1] + params[4]*params[3] + (1 - params[2] - params[4])*params[5]))
                    # push!(lifetime_m, (params[2]*params[1] + params[4]*params[3] + (1 - params[2] - params[4])*params[5]))
                end

                progress[] = i/length(files)

                files_label.text[] = "Files : $i/$(length(files))"

                i += 1
                sleep(0.001)
            end
        end

        rowsize!(fig.layout, 3, Fixed(0))

        delete!(axis_progress)
        delete!(files_label)
        delete!(percentage_label)

        filepath = save_file()

        if !isempty(filepath)
            open(filepath, "w") do io
                writedlm(io, [string(now())])
                writedlm(io, [""])
                writedlm(io, [time counts lifetime_1 concentration], ",")
            end
        end
    end

    # if running[]
    #     Threads.@spawn saving(menu_lifetimes.i_selected[])
    # end
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

                if signal_val < min_val/2
                    send_command(serial, "A 0 DO 1 0\n")

                elseif signal_val > 100 - min_val/2
                    send_command(serial, "A 0 DO 1 1\n")

                else
                    signal_val = round(Int, clamp(signal_val, min_val, max_val))

                    send_command(serial, "A 0 DP 1 $frequency $signal_val\n")
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

global axis_hist    = plotgrid[1, 1] = Axis(fig, alignmode = Outside(), title="Histogram", xlabel = "Time [ns]", ylabel = "Photon counts", titlesize=16)
global axis_counts  = plotgrid[2, 1] = Axis(fig, alignmode = Outside(), title="Photon counts", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16)
global axis_conc    = plotgrid[1, 3] = Axis(fig, alignmode = Outside(), title="Ion concentration", xlabel = "Time [s]", ylabel = "Ion concentration [mM]", titlesize=16)
global axis_signal  = plotgrid[2, 3] = Axis(fig, alignmode = Outside(), title="Light intensity", xlabel = "Time [s]", ylabel = "[%]", titlesize=16)

ylims!(axis_signal, -5, 105)

label_frequency     = Label(buttongrid[1, 1], justification=:left, text="Frequency : 0 Hz", width = 120)
label_number        = Label(buttongrid[2, 1], justification=:left, text="Files : 0", width = 120)
label_inverted      = Label(buttongrid[3, 1][1, 1], justification=:left, text="Inverted", width = 70)
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
Box(buttongrid[3, 6]; box_attributes...)

button_start    = Button(buttongrid[1, 2]; merge(button_attributes, Dict(:label => "START"))...)
button_stop     = Button(buttongrid[1, 3]; merge(button_attributes, Dict(:label => "CLEAR"))...)
button_load     = Button(buttongrid[2, 2]; merge(button_attributes, Dict(:label => "LOAD"))...)
button_load_IRF = Button(buttongrid[2, 3]; merge(button_attributes, Dict(:label => "LOAD IRF"))...)
button_save     = Button(buttongrid[3, 2]; merge(button_attributes, Dict(:label => "SAVE"))...)

menu_mode = Menu(buttongrid[3, 3],
    cell_color_active=RGBf(0.8, 0.8, 0.8),
    cell_color_inactive_even=RGBf(0.9, 0.9, 0.9),
    cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9),
    cell_color_hover=RGBf(0.6, 0.6, 0.6),
    default = "Playback",
    fontsize=20,
    height = 40,
    options=["Real-time", "Playback"],
    selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8),
    textcolor=:black,
    width = 120
)

label_K = Label(buttongrid[1:3, 4:5][1, 1], text=L"$K_q$ [$\text{M}^{-1}$] :", halign=:right)
textbox_K = Textbox(buttongrid[1:3, 4:5][1, 2]; textbox_attributes...)
textbox_K.displayed_string[] = "204"
textbox_K.stored_string[] = "204"

label_tau_0 = Label(buttongrid[1:3, 4:5][2, 1], text=L"$\tau_0$ [ns] :", halign=:right)
textbox_tau_0 = Textbox(buttongrid[1:3, 4:5][2, 2]; textbox_attributes...)
textbox_tau_0.displayed_string[] = "18"
textbox_tau_0.stored_string[] = "18"

label_fps = Label(buttongrid[1:3, 4:5][3, 1], text=L"\text{fps [Hz] :}", halign=:right)
textbox_fps = Textbox(buttongrid[1:3, 4:5][3, 2]; textbox_attributes...)
textbox_fps.displayed_string[] = "50"
textbox_fps.stored_string[] = "50"

label_min = Label(buttongrid[1:3, 4:5][4, 1], text=L"\text{puff min. [ms] :}", halign=:right)
textbox_min = Textbox(buttongrid[1:3, 4:5][4, 2]; textbox_attributes...)
textbox_min.displayed_string[] = "3"
textbox_min.stored_string[] = "3"

button_active = Button(buttongrid[1, 6],; merge(button_attributes, Dict(:label => "INACTIVE"))...)
button_connect = Button(buttongrid[2, 6]; merge(button_attributes, Dict(:label => "Connect"))...)

portname_textbox = Textbox(buttongrid[3, 6]; merge(textbox_attributes, Dict(:borderwidth => 0, :cornerradius => 6, :height => 40, :textpadding => (12, 8, 8, 12), :width => 120))...)
portname_textbox.displayed_string[] = "usbmodem2302"
portname_textbox.stored_string[] = "usbmodem2302"

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
    NamedTuple(merge(slider_attributes, Dict( :format => "{:d}", :label => "Time range", :range => 1:1:600, :startvalue => 30))),
    NamedTuple(merge(slider_attributes, Dict( :format => "{:d}", :label => "Averaging", :range => 1:1:100, :startvalue => 1))),
    NamedTuple(merge(slider_attributes, Dict( :format => "{:.2f}", :label => "Proportional", :range => 0:0.01:1, :startvalue => 0))),
    NamedTuple(merge(slider_attributes, Dict( :format => "{:.2f}", :label => "Integral", :range => 0:0.01:1, :startvalue => 0))),
    NamedTuple(merge(slider_attributes, Dict( :format => "{:.1f}", :label => "Setpoint", :range => 0:0.1:100, :startvalue => 25))),
    width = 360
)

colgap!(buttongrid, 3, 50)
colgap!(buttongrid, 5, 50)
colsize!(plotgrid, 2, 10)
rowgap!(buttongrid, 1, 20)
rowgap!(buttongrid, 2, 20)

global irf = get_irf()
global irf_bin_size = get_irf_bin_size(irf)
global folderpath = get_folder_path("FLIM")
global serial = nothing
global initialized = Observable(false)
global running = Observable(false)
global pause = Observable(false)
global active = Observable(false)
global connected = Observable(false)
# global inverted = Observable(false)
global signal_intensity = 0

setpoint = hlines!(axis_conc, params_sliders.sliders[5].value[], color=Makie.wong_colors()[2])

on(button_start.clicks) do _
    if !running[]
        running[] = true

        button_start.label = "PAUSE"
        button_stop.label = "STOP"

        if menu_mode.selection[] == "Playback"
            @async start_playback()
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
    initialized[] = false

    button_start.label = "START"
    button_stop.label = "CLEAR"

    empty!(axis_hist)
    empty!(axis_counts)
    empty!(axis_conc)
    empty!(axis_signal)

    global setpoint = hlines!(axis_conc, params_sliders.sliders[5].value[], color=Makie.wong_colors()[2])
end

on(button_load.clicks) do _
    global folderpath = get_new_folder_path("FLIM")
end

on(button_load_IRF.clicks) do _
    global irf = get_new_irf()
    global irf_bin_size = get_irf_bin_size(irf)
end

on(button_save.clicks) do _
    start_save()
end

on(button_connect.clicks) do _
    if !connected[]
        global serial = connect_to_port(portname_textbox.stored_string[])

        sleep(0.1)
        
        try
            send_command(serial, "C Z\nC R\n")
            send_command(serial, "S 1 S 1\n")
        catch e
            println(e)
        end

        button_connect.label[] = "Disconnect"
    else
        LibSerialPort.close(serial)
        add_to_monitor("Disconnection successful")
        println("Disconnection successful")
        button_connect.label[] = "Connect"
        connected[] = false
    end
end

on(button_active.clicks) do _
    if active[]
        active[] = false
        button_active.label = "INACTIVE"
        add_to_monitor("Deactivated")
    else
        active[] = true
        button_active.label = "ACTIVE"

        add_to_monitor("Activated")
        @async start_signal(parse(Float64, textbox_fps.stored_string[]), parse(Float64, textbox_fps.stored_string[])/10, parse(Float64, textbox_min.stored_string[]))
    end
end

on(params_sliders.sliders[5].value) do value
    setpoint.args[1][] = value
end

fig
