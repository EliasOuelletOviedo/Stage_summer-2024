using GLMakie
using PyCall
using FFTW
using Statistics
using NativeFileDialog
using ZipFile
using Observables
using Base.Threads
using Dates

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
        filepath = open_dialog("IRF: Choose .sdt file to open.")
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

function MLE_iterative_reconvolution(irf, data_xy; params, gating_function=ones(UInt8, 320), histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
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

function LS_iterative_reconvolution(irf, data_xy; guess=[2.0, 1.0], bounds=([0.0, -30], [12.5, 256]), irf_shift=0, gating_function=ones(UInt8, 320), histogram_resolution=256)
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

function start_real_time(lifetimes)
    if !initialized[]
        println("Initializing")

        hist = Observable(zeros(320))
        fit = Observable(zeros(320))
        hist_time = Observable(zeros(320))

        lifetime_1 = Observable(Float64[])
        lifetime_2 = Observable(Float64[])
        lifetime_3 = Observable(Float64[])
        lifetime_i = Observable(Float64[])
        lifetime_m = Observable(Float64[])
        counts = Observable(Float64[])
        concentration = Observable(Float64[])

        time = Observable(Float64[])

        vectors = zeros(100, 256)

        barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
        lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
        lines!(axis_life, time, lifetime_1, color=Makie.wong_colors()[1], label = L"\tau_1")
        lines!(axis_counts, time, counts, color=Makie.wong_colors()[1])
        lines!(axis_conc, time, concentration, color=Makie.wong_colors()[menu_main_lifetime.i_selected[]])

        if lifetimes == 2
            lines!(axis_life, time, lifetime_2, color=Makie.wong_colors()[2], label = L"\tau_2")
            lines!(axis_life, time, lifetime_i, color=Makie.wong_colors()[3], label = L"\tau_i")
            lines!(axis_life, time, lifetime_m, color=Makie.wong_colors()[4], label = L"\tau_m")
        elseif lifetimes == 3
            lines!(axis_life, time, lifetime_2, color=Makie.wong_colors()[2], label = L"\tau_2")
            lines!(axis_life, time, lifetime_3, color=Makie.wong_colors()[3], label = L"\tau_3")
            lines!(axis_life, time, lifetime_i, color=Makie.wong_colors()[4], label = L"\tau_i")
            lines!(axis_life, time, lifetime_m, color=Makie.wong_colors()[5], label = L"\tau_m")
        end

        global legend = axislegend(axis_life)

        previous_time = now()
        current_time = now()
        initialized[] = true
        params = ones(2*lifetimes+1)
        frequency = 1

        i = 1
        t_0 = now()
    end

    println("Waiting")

    while running[]
        current_files = readdir(folderpath)

        new_files = setdiff(current_files, previous_files)

        sleep(0.001)

        if !isempty(new_files)
            if i == 1
                println("Files detected")
            end

            vector, histogram_resolution = open_SDT_file(joinpath(folderpath, new_files[end]))

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

            # vector = sum(vectors, dims=1)[1,:]/100

            params, data = vec_to_lifetime(sum(vectors[end-slider_averaging.value[]+1:end,begin:histogram_resolution], dims=1)[1,:]/slider_averaging.value[], histogram_resolution=histogram_resolution, guess=params[1:2*lifetimes])
            # params = vec_to_lifetime(vector, histogram_resolution=histogram_resolution)

            hist[][1:length(data[2])] = data[2]
            hist_time[][1:length(data[2])] = data[1]
            fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

            if i == 1
                t_0 = now()
            end

            push!(lifetime_1[], params[1])
            push!(counts[], sum(vectors[end-slider_averaging.value[]+1:end,:])/slider_averaging.value[])
            push!(concentration[], 1/(1/parse(Float64, textbox_tau_0.stored_string[]) + parse(Float64, textbox_K.stored_string[])*params[2*menu_main_lifetime.i_selected[] - 1]))
            # push!(concentration[], 1/(1/parse(Int8, textbox_tau_0.stored_string[]) + parse(Int8, textbox_tau_0.stored_string[])*params[2*menu_main_lifetime.i_selected[] - 1]))
            push!(time[], (now() - t_0)/Millisecond(1)/1000)
            # push!(time[], i)

            if lifetimes == 2
                push!(lifetime_2[], params[3])
            elseif lifetimes == 3
                push!(lifetime_2[], params[3])
                push!(lifetime_3[], params[5])
            end

            if mod1(i, 8) == 8
                current_time = now()
                label_frequency.text = "Frequency : $(string(round(1/((current_time - previous_time)/Millisecond(1))*1000*8, sigdigits=2))) Hz"
                previous_time = current_time

                # label_number.text = "Files : $(round(UInt64, i))"
            end

            if mod1(i, ceil(frequency/16)) == ceil(frequency/16)
                notify(hist_time)
                notify(time)

                autolimits!(axis_hist)
                autolimits!(axis_life)
                autolimits!(axis_counts)
                autolimits!(axis_conc)
                
                xlims!(axis_life, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
                xlims!(axis_counts, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
                xlims!(axis_conc, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
            end

            previous_files = current_files

            i += 1
        end
    end
end

function start_playback(lifetimes)
    println(folderpath)
    all_files = readdir(folderpath)
    files = filter(f -> occursin(".sdt", f), all_files)

    if !initialized[]
        println("Initializing")

        hist = Observable(zeros(320))
        fit = Observable(zeros(320))
        hist_time = Observable(zeros(320))

        lifetime_1 = Observable(Float64[])
        lifetime_2 = Observable(Float64[])
        lifetime_3 = Observable(Float64[])
        lifetime_i = Observable(Float64[])
        lifetime_m = Observable(Float64[])
        counts = Observable(Float64[])
        concentration = Observable(Float64[])

        time = Observable(Float64[])

        vectors = zeros(100, 256)

        barplot!(axis_hist, hist_time, hist, gap = 0, color=Makie.wong_colors()[1])
        lines!(axis_hist, hist_time, fit, color=Makie.wong_colors()[6])
        lines!(axis_life, time, lifetime_1, color=Makie.wong_colors()[1], label = L"\tau_1")
        lines!(axis_counts, time, counts, color=Makie.wong_colors()[1])
        lines!(axis_conc, time, concentration, color=Makie.wong_colors()[menu_main_lifetime.i_selected[]])

        if lifetimes == 2
            lines!(axis_life, time, lifetime_2, color=Makie.wong_colors()[2], label = L"\tau_2")
            lines!(axis_life, time, lifetime_i, color=Makie.wong_colors()[3], label = L"\tau_i")
            lines!(axis_life, time, lifetime_m, color=Makie.wong_colors()[4], label = L"\tau_m")
        elseif lifetimes == 3
            lines!(axis_life, time, lifetime_2, color=Makie.wong_colors()[2], label = L"\tau_2")
            lines!(axis_life, time, lifetime_3, color=Makie.wong_colors()[3], label = L"\tau_3")
            lines!(axis_life, time, lifetime_i, color=Makie.wong_colors()[4], label = L"\tau_i")
            lines!(axis_life, time, lifetime_m, color=Makie.wong_colors()[5], label = L"\tau_m")
        end

        global legend = axislegend(axis_life)

        previous_time = now()
        current_time = now()
        initialized[] = true
        params = ones(2*lifetimes+1)
        frequency = 1

        i = 1
        t_0 = now()
    end

    println("Starting")

    while running[]
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

                # vector = sum(vectors, dims=1)[1,:]/100

                params, data = vec_to_lifetime(sum(vectors[end-slider_averaging.value[]+1:end,begin:histogram_resolution], dims=1)[1,:]/slider_averaging.value[], histogram_resolution=histogram_resolution, guess=params[1:2*lifetimes])
                # params = vec_to_lifetime(vector, histogram_resolution=histogram_resolution)

                hist[][1:length(data[2])] = data[2]
                hist_time[][1:length(data[2])] = data[1]
                fit[][1:length(data[2])] = conv_irf_data(data[1], params, irf, data[2], histogram_resolution=histogram_resolution)

                if i == 1
                    t_0 = now()
                end

                push!(lifetime_1[], params[1])
                push!(counts[], sum(vectors[end-slider_averaging.value[]+1:end,:])/slider_averaging.value[])
                push!(concentration[], 1/(1/parse(Float64, textbox_tau_0.stored_string[]) + parse(Float64, textbox_K.stored_string[])*params[2*menu_main_lifetime.i_selected[] - 1]))
                # push!(concentration[], 1/(1/parse(Int8, textbox_tau_0.stored_string[]) + parse(Int8, textbox_tau_0.stored_string[])*params[2*menu_main_lifetime.i_selected[] - 1]))
                push!(time[], (now() - t_0)/Millisecond(1)/1000)
                # push!(time[], i)

                if lifetimes == 2
                    push!(lifetime_2[], params[3])
                    push!(lifetime_i[], (params[2]*params[1]^2 + (1 - params[2])*params[3]^2)/(params[2]*params[1] + (1 - params[2])*params[3]))
                    push!(lifetime_m[], (params[2]*params[1] + (1 - params[2])*params[3]))
                elseif lifetimes == 3
                    push!(lifetime_2[], params[3])
                    push!(lifetime_3[], params[5])
                    push!(lifetime_i[], (params[2]*params[1]^2 + params[4]*params[3]^2 + (1 - params[2] - params[4])*params[5]^2)/(params[2]*params[1] + params[4]*params[3] + (1 - params[2] - params[4])*params[5]))
                    push!(lifetime_m[], (params[2]*params[1] + params[4]*params[3] + (1 - params[2] - params[4])*params[5]))
                end

                if mod1(i, 8) == 8
                    current_time = now()
                    frequency = round(1/((current_time - previous_time)/Millisecond(1))*1000*8, sigdigits=2)
                    label_frequency.text = "Frequency : $frequency Hz"
                    previous_time = current_time
        
                    label_number.text = "Files : $(round(UInt64, i))"
                end

                if mod1(i, ceil(frequency/16)) == ceil(frequency/16)
                    notify(hist_time)
                    notify(time)

                    autolimits!(axis_hist)
                    autolimits!(axis_life)
                    autolimits!(axis_counts)
                    autolimits!(axis_conc)
                    
                    xlims!(axis_life, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
                    xlims!(axis_counts, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
                    xlims!(axis_conc, max(0, time[][end] - slider_time.value[]+0.1), max(slider_time.value[], time[][end]+0.1))
                end

                i += 1
                sleep(0.001)
            end
        end
    end
end

fig = Figure(size = (1440, 847))

set_theme!(backgroundcolor = :gray95)

fig[2, 1] = buttongrid = GridLayout(tellwidth=false)
fig[1, 1] = plotgrid = GridLayout()

global axis_hist = plotgrid[1, 1] = Axis(fig, title="Histogram", xlabel = "Time [ns]", ylabel = "Photon counts", titlesize=16)
global axis_life = plotgrid[1, 3] = Axis(fig, title="Lifetime", xlabel = "Time [s]", ylabel = "Lifetime [ns]", titlesize=16)
global axis_counts = plotgrid[2, 1] = Axis(fig, title="Photon counts", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=16)
global axis_conc = plotgrid[2, 3] = Axis(fig, title="Ion concentration", xlabel = "Time [s]", ylabel = "Ion concentration [mM]", titlesize=16)

colsize!(plotgrid, 2, 10)

label_frequency = buttongrid[1, 1] = Label(fig, text="Frquency : 0 Hz")
label_number = buttongrid[2, 1] = Label(fig, text="Files : 0")

button_start = buttongrid[1, 2] = Button(fig, label="START", font=:bold, width = 120, height = 40, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)
button_stop = buttongrid[1, 3] = Button(fig, label="CLEAR", font=:bold, width = 120, height = 40, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)

button_load = buttongrid[2, 2] = Button(fig, label="LOAD", font=:bold, width = 120, height = 40, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)
button_load_IRF = buttongrid[2, 3] = Button(fig, label="LOAD IRF", font=:bold, width = 120, height = 40, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)

button_save = buttongrid[1, 4] = Button(fig, label="SAVE", font=:bold, width = 120, height = 40, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=16, cornerradius=8, labelcolor_active=:black)
# menu_mode = buttongrid[2, 4] = Menu(fig, options=["Real-time", "Playback"], default = "Real-time", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)
menu_mode = buttongrid[2, 4] = Menu(fig, options=["Real-time", "Playback"], default = "Playback", height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)

label_lifetimes = buttongrid[1:2, 5][1, 1] = Label(fig, text="Number of lifetimes :", halign=:right)
menu_lifetimes = buttongrid[1:2, 5][1, 2] = Menu(fig, options=["1", "2", "3"], default = "1", height = 32, width = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black, dropdown_arrow_size=0)

label_main_lifetime = buttongrid[1:2, 5][2, 1] = Label(fig, text="Principal lifetime :", halign=:right)
menu_main_lifetime = buttongrid[1:2, 5][2, 2] = Menu(fig, options=["1"], default = "1", height = 32, width = 32, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black, dropdown_arrow_size=0)

label_K = buttongrid[1:2, 6][1, 1] = Label(fig, text=L"$K_q$ [-] :", halign=:right)
textbox_K = buttongrid[1:2, 6][1, 2] = Textbox(fig, placeholder = " ", height = 32, width = 64, cornerradius=8, boxcolor=:white, boxcolor_focused=:white, boxcolor_hover=:white, bordercolor=:transparent, bordercolor_hover=RGBf(0.8, 0.8, 0.8), bordercolor_focused=RGBf(0.6, 0.6, 0.6))
textbox_K.displayed_string[] = "1"
textbox_K.stored_string[] = "1"

label_tau_0 = buttongrid[1:2, 6][2, 1] = Label(fig, text=L"$\tau_0$ [ns] :", halign=:right)
textbox_tau_0 = buttongrid[1:2, 6][2, 2] = Textbox(fig, placeholder = " ", height = 32, width = 64, cornerradius=8, boxcolor=:white, boxcolor_focused=:white, boxcolor_hover=:white, bordercolor=:transparent, bordercolor_hover=RGBf(0.8, 0.8, 0.8), bordercolor_focused=RGBf(0.6, 0.6, 0.6))
textbox_tau_0.displayed_string[] = "1"
textbox_tau_0.stored_string[] = "1"

slider_time_label = buttongrid[1, 7][1, 1] = Label(fig, text="Time range : 10", tellheight=false, valign=:top)
slider_time = buttongrid[1, 7][2, 1] = Slider(fig, range = 1:1:60, startvalue = 10, width = 240, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false)
textbox_time = buttongrid[1, 7][2, 2] = Textbox(fig, placeholder = " ", height = 30, width = 40, textpadding=(8,6,6,7), halign=:center, valign=:center, bordercolor=:transparent, bordercolor_focused=:transparent, bordercolor_focused_invalid=:transparent, bordercolor_hover=:transparent, textcolor_placeholder=RGBf(0.4, 0.4, 0.4), textcolor=RGBf(0.4, 0.4, 0.4), validator=UInt16)
textbox_time.displayed_string[] = "60"
textbox_time.stored_string[] = "60"

slider_averaging_label = buttongrid[2, 7][1, 1] = Label(fig, text="Averaging : 1", tellheight=false, valign=:top)
slider_averaging = buttongrid[2, 7][2, 1] = Slider(fig, range = 1:1:100, startvalue = 1, width = 240, linewidth = 14, color_active=RGBf(0.4, 0.4, 0.4), color_active_dimmed=RGBf(0.6, 0.6, 0.6), color_inactive=RGBf(0.8, 0.8, 0.8), tellheight=false)
textbox_averaging = buttongrid[2, 7][2, 2] = Textbox(fig, placeholder = " ", height = 30, width = 40, textpadding=(8,6,6,7), halign=:center, valign=:center, bordercolor=:transparent, bordercolor_focused=:transparent, bordercolor_focused_invalid=:transparent, bordercolor_hover=:transparent, textcolor_placeholder=RGBf(0.4, 0.4, 0.4), textcolor=RGBf(0.4, 0.4, 0.4), validator=UInt16)
textbox_averaging.displayed_string[] = "100"
textbox_averaging.stored_string[] = "100"

colsize!(buttongrid, 4, 200)
colsize!(buttongrid, 6, 200)

global irf = get_irf()
global irf_bin_size = get_irf_bin_size(irf)
global folderpath = get_folder_path("FLIM")

global initialized = Observable(false)
running = Observable(false)
ran = Observable(false)
saved = false

on(button_start.clicks) do _
    if !running[]
        running[] = true
        ran[] = true
        saved = false

        button_start.label = "PAUSE"
        button_stop.label = "STOP"

        if menu_mode.selection[] == "Playback"
            @async start_playback(menu_lifetimes.i_selected[])
        else
            @async start_real_time(menu_lifetimes.i_selected[])
        end
    else
        running[] = false

        button_start.label = "START"
        button_stop.label = "CLEAR"
    end
end

on(button_stop.clicks) do _
    running[] = false
    ran[] = false
    saved = false
    initialized[] = false

    button_start.label = "START"
    button_stop.label = "CLEAR"

    empty!(axis_hist)
    empty!(axis_life)
    empty!(axis_counts)
    empty!(axis_conc)

    delete!(legend)
end

on(button_load.clicks) do _
    global folderpath = get_new_folder_path("FLIM")
end

on(menu_lifetimes.i_selected) do index
    menu_main_lifetime.options[] = menu_lifetimes.options[][1:index]

    if menu_main_lifetime.i_selected[] == 0
        menu_main_lifetime.i_selected[] = index
    end
end

on(slider_time.value) do value
    slider_time_label.text = "Time range : $value"
end

on(textbox_time.displayed_string) do string
    if length(string) > 4
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

    println(string)

    slider_time.range = 1:1:parse(UInt16, string)
end

on(slider_averaging.value) do value
    slider_averaging_label.text = "Averaging : $value"
end

on(textbox_averaging.displayed_string) do string
    if length(string) > 4
        textbox_averaging.displayed_string[] = string[begin:end-1]
        string = string[begin:end-1]
    end
end

on(textbox_averaging.stored_string) do string
    if parse(UInt16, string) < 2
        textbox_averaging.stored_string[] = "2"
        textbox_averaging.displayed_string[] = "2"
        string = "2"
    elseif parse(UInt16, string) > 100
        textbox_averaging.stored_string[] = "100"
        textbox_averaging.displayed_string[] = "100"
        string = "100"
    end

    println(string)

    slider_averaging.range = 1:1:parse(UInt16, string)
end

fig
