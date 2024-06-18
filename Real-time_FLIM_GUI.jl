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
using GLMakie
using Observables
using Dates
using Colors

@everywhere sp_o = pyimport("scipy.optimize")

# functions ignoring NaNs
@everywhere nanmean(x) = mean(filter(!isnan, x))
@everywhere nanmean(x, dims) = mapslices(nanmean, x, dims=dims)

@everywhere nanstd(x) = std(filter(!isnan, x))
@everywhere nanstd(x, dims) = mapslices(nanstd, x, dims=dims)

@everywhere nanmax(x) = maximum(filter(!isnan, x))
@everywhere nanmin(x) = minimum(filter(!isnan, x))

@everywhere is_zero(x) = x == 0
@everywhere smaller_or_eq_zero(x) = x <= 0

@everywhere function get_irf(; channel=1)
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

@everywhere function get_irf_bin_size(irf)
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

@everywhere function convolve(irf, decay, y_data; histogram_resolution=256)
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

@everywhere function irf_shift(data_irf, shift)
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

@everywhere function MLE_iterative_reconvolution(irf, data_xy; params, gating_function=ones(UInt8, 320), histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
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

@everywhere function MLE_model_func(params, x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period)
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

@everywhere function conv_irf_data(x_data, params, irf, y_data; histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
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

@everywhere function LS_model_func_1_lt(data, t_1, d_0, y_offset)
    """
    Model function for LS optimizer for 1 lifetime
    """
    x_data = data[1]
    irf = data[2]
    y_data = data[3]
    A = conv_irf_data(x_data, [t_1, d_0, y_offset], irf, y_data, histogram_resolution=histogram_resolution_fit).*gating_function
    return A
end

@everywhere function LS_model_func_2_lt(data, t_1, A_1, t_2, d_0, y_offset)
    """
    Model function for LS optimizer for 2 lifetimes
    """
    x_data = data[1]
    irf = data[2]
    y_data = data[3]
    A = conv_irf_data(x_data, [t_1, A_1, t_2, d_0, y_offset], irf, y_data, histogram_resolution=histogram_resolution_fit).*gating_function
    return A
end

@everywhere function LS_model_func_3_lt(data, t_1, A_1, t_2, A_2, t_3, d_0, y_offset)
    """
    Model function for LS optimizer for 3 lifetimes
    """
    x_data = x_data[1]
    irf = x_data[2]
    y_data = data[3]
    A = conv_irf_data(x_data, [t_1, A_1, t_2, A_2, t_3, d_0, y_offset], irf, y_data, histogram_resolution=histogram_resolution_fit).*gating_function
    return A
end

@everywhere function LS_iterative_reconvolution(irf, data_xy; guess=[2.0, 1.0], bounds=([0.0, -30], [12.5, 256]), irf_shift=0, gating_function=ones(UInt8, 320), histogram_resolution=256)
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

@everywhere function vec_to_lifetime(x; bin_size=3.9088729e-2, hist_size_threshold=100, method="MLE", broadcast_call=false, guess=[4.0, 3.0, 2.0, 1.0], laser_pulse_period=12.5, histogram_resolution=256, number_of_previous_pulses=5, tcspc_window_size=10, tac_low_cut=5.0980392, tac_high_cut=94.901962)
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

@everywhere function get_folder_path()
    if isfile("folder_path.txt")
        folderpath = open(f->read(f, String), "folder_path.txt")

        if isempty(folderpath)
            py"""
            from tkinter.filedialog import askdirectory

            path = askdirectory()
            """

            folderpath = py"path"

            open("folder_path.txt", "w") do io
                write(io, folderpath)
            end
            end
    else
        py"""
        from tkinter.filedialog import askdirectory

        path = askdirectory()
        """

        folderpath = py"path"

        open("folder_path.txt", "w") do io
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
    py"""
    from tkinter.filedialog import askdirectory

    path = askdirectory()
    """

    folderpath = py"path"

    open("folder_path.txt", "w") do io
        write(io, folderpath)
    end

    if isempty(folderpath)
        println("Folder does not exist")

        folderpath = get_folder_path()
    end

    return folderpath
end

@everywhere function SDTFile(filepath)
    open(filepath,"r") do io
        seek(io, 14)
        header = read(io, 2)
        seek(io, header[2]*16^2+header[1])
    
        file_info = read(io, 22)
    
        zip = IOBuffer(read(io, file_info[10]*16^6 + file_info[9]*16^4 + file_info[8]*16^2 + file_info[7] - (file_info[6]*10^6 + file_info[5]*16^4 + file_info[4]*16^2 + file_info[3])))

        return reinterpret(UInt16, read(ZipFile.Reader(zip).files[1], file_info[22]*16^6 + file_info[21]*16^4 + file_info[20]*16^2 + file_info[19]))
    end
end

@everywhere function SDTFile_single(filepath)
    open(filepath,"r") do io
        seek(io, 14)
        header = read(io, 2)
        seek(io, header[2]*16^2+header[1]+22)

        vector = read(io)
        new_vector = UInt16[]

        for i in 1:2:length(vector)
            push!(new_vector, vector[i+1]*16^2 + vector[i])
        end

        return new_vector
    end
end

@everywhere function get_file_lifetime(filepath)
    file = SDTFile(filepath)
    
    array = reshape(file, (256, :))

    vector = sum(array, dims = 2)[:,1]

    return vec_to_lifetime(vector)
end

@everywhere function display_lifetime_folder(folderpath)

    for i in readdir(folderpath)
        if endswith(i, ".sdt") || endswith(i[begin:end-6], ".sdt") || endswith(i[begin:end-8], ".sdt")
            filepath = joinpath(folderpath, i)

            @time region_lifetime = get_file_lifetime(filepath)[1][1]

            push!(time, length(lifetime))
            push!(lifetime, region_lifetime)
        end
    end
end

@everywhere function start_real_time(folderpath, type, hist_res, axis_histogram, axis_lifetime, axis_photons, axis_concentration)
    channels = Int(round(laser_pulse_period[]*hist_res/tcspc_window_size[]))

    hist_time = LinRange(0, laser_pulse_period[], channels)
    hist = Observable(zeros(channels))
    curve_fit = Observable(zeros(channels))
    real_time = Observable(Float64[])
    lifetime = Observable(Float64[])
    photons = Observable(Int64[])
    concentration = Observable(Float64[])

    line_histogram = barplot!(axis_histogram, hist_time, hist, gap = 0, color = Colors.JULIA_LOGO_COLORS.blue)
    line_curve_fit = lines!(axis_histogram, hist_time, curve_fit, color = Colors.JULIA_LOGO_COLORS.red, linewidth=2)
    line_lifetime = lines!(axis_lifetime, real_time, lifetime, color = Colors.JULIA_LOGO_COLORS.blue)
    line_photons = lines!(axis_photons, real_time, photons, color = Colors.JULIA_LOGO_COLORS.blue)
    line_concentration = lines!(axis_concentration, real_time, concentration, color = Colors.JULIA_LOGO_COLORS.blue)

    print("Folder path : ")
    println(folderpath)
    
    if type == "Single"
        println(size(hist[]))
        println("Single")
        @async new_file_lifetime_single(folderpath, hist_res, axis_histogram, axis_lifetime, axis_photons, axis_concentration, hist, curve_fit, real_time, lifetime, photons, concentration)
        # @async new_file_lifetime_single(folderpath, hist_res, axis_histogram, axis_lifetime, axis_photons, axis_concentration)
    else
        println(size(hist[]))
        println("Image")
        @async new_file_lifetime(folderpath, hist_res, axis_histogram, axis_lifetime, axis_photons, axis_concentration, hist, curve_fit, real_time, lifetime, photons, concentration)
    end
end

@everywhere function new_file_lifetime(path, hist_res, axis_histogram, axis_lifetime, axis_photons, axis_concentration, hist, curve_fit, real_time, lifetime, photons, concentration)
    previous_files = readdir(path)

    t_0 = now()

    i = 0

    println("Search started")

    while running[]
        current_files = readdir(path)

        new_files = setdiff(current_files, previous_files)

        if !isempty(new_files)
            
            println("File detected!")

            for new_file in new_files
                # if endswith(new_file, ".sdt") || endswith(new_file[begin:end-6], ".sdt") || endswith(new_file[begin:end-8], ".sdt") || endswith(new_file[begin:end-9], ".sdt")
                println(joinpath(path, new_file))
                file = SDTFile(joinpath(path, new_file))
                array = reshape(file, (hist_res, :))
                file_params = vec_to_lifetime(sum(array, dims = 2)[:,1], histogram_resolution=hist_res, laser_pulse_period=laser_pulse_period[], tcspc_window_size=tcspc_window_size[])

                file_hist = file_params[2][2]
                file_curve_fit = conv_irf_data(file_params[2][1], file_params[1], irf, file_params[2][2], histogram_resolution=hist_res)
                file_lifetime = file_params[1][1]

                if isnan(file_lifetime)
                    if i == 0
                        file_lifetime = 1
                    else
                        file_lifetime = lifetime[][end]
                    end
                end

                file_photons = sum(file_hist)[1]
                file_concentration = 1/(1/tau_0[] + K_q[]*file_lifetime)
                # file_concentration = file_params[1][3]

                if i == 0
                    t_0 = now()
                    i += 1
                end

                present_time = (now() - t_0)/Millisecond(1)/1000

                hist[] = file_hist
                curve_fit[] = file_curve_fit
                push!(real_time[], present_time)
                push!(lifetime[], file_lifetime)
                push!(photons[], file_photons)
                push!(concentration[], file_concentration)

                notify(hist)
                notify(curve_fit)
                notify(real_time)
                notify(lifetime)
                notify(photons)
                notify(concentration)

                limits!(axis_histogram, 0, laser_pulse_period[], 0, maximum(hist[])*1.1)
                limits!(axis_lifetime, max(0, present_time-10), max(10, present_time), minimum(lifetime[])*0.99, maximum(lifetime[])*1.01)
                limits!(axis_photons, max(0, present_time-10), max(10, present_time), minimum(photons[])*0.99, maximum(photons[])*1.01)
                limits!(axis_concentration, max(0, present_time-10), max(10, present_time), minimum(concentration[])*0.99, maximum(concentration[])*1.01)
                
                # end

                # println(i)
                sleep(0.01)
            end
        end

        previous_files = current_files

        if save[] == true
            filename = tb_file_name.displayed_string[]*".csv"

            if filename == "Enter a file name.csv"
                filename = split(path, "/")[end]*".csv"
            end

            open(filename, "w") do io
                writedlm(io, [string(now())])
                writedlm(io, [""])
                writedlm(io, [real_time[] lifetime[] photons[] concentration[]], ",")
            end

            save[] = false
        end

        sleep(0.01)
    end
end

@everywhere function new_file_lifetime_single(path, hist_res, axis_histogram, axis_lifetime, axis_photons, axis_concentration, hist, curve_fit, real_time, lifetime, photons, concentration)
    previous_files = readdir(path)

    t_0 = now()

    i = 0

    println("Search started")

    while running[]
        current_files = readdir(path)

        new_files = setdiff(current_files, previous_files)

        if !isempty(new_files)

            println("File detected!")

            for new_file in new_files
                # if endswith(new_file, ".sdt") || endswith(new_file[begin:end-6], ".sdt") || endswith(new_file[begin:end-8], ".sdt") || endswith(new_file[begin:end-9], ".sdt")
                println(joinpath(path, new_file))
                
                vector = SDTFile_single(joinpath(path, new_file))

                file_params = vec_to_lifetime(vector, histogram_resolution=hist_res, laser_pulse_period=laser_pulse_period[], tcspc_window_size=tcspc_window_size[])
                file_hist = file_params[2][2]
                file_curve_fit = conv_irf_data(file_params[2][1], file_params[1], irf, file_params[2][2], histogram_resolution=hist_res)
                file_lifetime = file_params[1][1]

                if isnan(file_lifetime)
                    if i == 0
                        file_lifetime = 1
                    else
                        file_lifetime = lifetime[][end]
                    end
                end

                file_photons = sum(file_hist)[1]
                file_concentration = 1/(1/tau_0[] + K_q[]*file_lifetime)

                if i == 0
                    t_0 = now()
                    i += 1
                end

                present_time = (now() - t_0)/Millisecond(1)/1000

                hist[] = file_hist
                curve_fit[] = file_curve_fit
                push!(real_time[], present_time)
                push!(lifetime[], file_lifetime)
                push!(photons[], file_photons)
                push!(concentration[], file_concentration)

                notify(hist)
                notify(curve_fit)
                notify(real_time)
                notify(lifetime)
                notify(photons)
                notify(concentration)

                limits!(axis_histogram, 0, laser_pulse_period[], 0, maximum(hist[])*1.1)
                limits!(axis_lifetime, max(0, present_time-10), max(10, present_time), minimum(lifetime[])*0.99, maximum(lifetime[])*1.01)
                limits!(axis_photons, max(0, present_time-10), max(10, present_time), minimum(photons[])*0.99, maximum(photons[])*1.01)
                limits!(axis_concentration, max(0, present_time-10), max(10, present_time), minimum(concentration[])*0.99, maximum(concentration[])*1.01)
                
                # end

                # println(i)
                sleep(0.01)
            end

            # new_file = new_files[end]
            # # if endswith(new_file, ".sdt") || endswith(new_file[begin:end-6], ".sdt") || endswith(new_file[begin:end-8], ".sdt") || endswith(new_file[begin:end-9], ".sdt")
            # vector = SDTFile_single(joinpath(path, new_file))

            # file_params = vec_to_lifetime(vector, histogram_resolution=hist_res, laser_pulse_period=laser_pulse_period[], tcspc_window_size=tcspc_window_size[])
            # file_hist = file_params[2][2]
            # file_curve_fit = conv_irf_data(file_params[2][1], file_params[1], irf, file_params[2][2], histogram_resolution=hist_res)
            # file_lifetime = file_params[1][1]

            # if isnan(file_lifetime)
            #     if i == 0
            #         file_lifetime = 1
            #     else
            #         file_lifetime = lifetime[][end]
            #     end
            # end

            # file_photons = sum(file_hist)[1]
            # file_concentration = 1/(1/tau_0[] + K_q[]*file_lifetime)

            # if i == 0
            #     t_0 = now()
            #     i += 1
            # end

            # present_time = (now() - t_0)/Millisecond(1)/1000

            # hist[] = file_hist
            # curve_fit[] = file_curve_fit
            # push!(real_time[], present_time)
            # push!(lifetime[], file_lifetime)
            # push!(photons[], file_photons)
            # push!(concentration[], file_concentration)

            # notify(hist)
            # notify(curve_fit)
            # notify(real_time)
            # notify(lifetime)
            # notify(photons)
            # notify(concentration)

            # limits!(axis_histogram, 0, laser_pulse_period[], 0, maximum(hist[])*1.1)
            # limits!(axis_lifetime, max(0, present_time-10), max(10, present_time), minimum(lifetime[])*0.99, maximum(lifetime[])*1.01)
            # limits!(axis_photons, max(0, present_time-10), max(10, present_time), minimum(photons[])*0.99, maximum(photons[])*1.01)
            # limits!(axis_concentration, max(0, present_time-10), max(10, present_time), minimum(concentration[])*0.99, maximum(concentration[])*1.01)
            
            # end

            # println(i)
            sleep(0.01)
        end

        previous_files = current_files

        if save[] == true
            filename = tb_file_name.displayed_string[]*".csv"

            if filename == "Enter a file name.csv"
                filename = split(path, "/")[end]*".csv"
            end

            open(filename, "w") do io
                writedlm(io, [string(now())])
                writedlm(io, [""])
                writedlm(io, [real_time[] lifetime[] photons[] concentration[]], ",")
            end

            save[] = false
        end

        sleep(0.01)
    end
end

@everywhere global irf = get_irf()
@everywhere global irf_bin_size = get_irf_bin_size(irf)
@everywhere global folderpath = get_folder_path()

set_theme!(backgroundcolor = :gray95)

tau_0 = Observable(1)
K_q = Observable(1)

laser_pulse_period = Observable(10)
tcspc_window_size = Observable(10)

running = Observable(false)
ran = Observable(false)
save = Observable(false)

fig = Figure(size = (1440, 847))

Label(fig[0, 1:4], text = "Title", font=:bold, fontsize = 32)

fig[1, 1:4] = axisgrid = GridLayout()

axis_histogram = axisgrid[1, 1] = Axis(fig, title="Histogram of the photon counts", xlabel = "Time [ns]", ylabel = "Photon counts", titlesize=20, xlabelsize=18, ylabelsize=18)
axis_lifetime = axisgrid[1, 3] = Axis(fig, title="Real-time plot of the lifetime", xlabel = "Time [s]", ylabel = "Lifetime [ns]", titlesize=20, xlabelsize=18, ylabelsize=18)
axis_photons = axisgrid[2, 1] = Axis(fig, title="Real-time plot of the total photon counts", xlabel = "Time [s]", ylabel = "Photon counts", titlesize=20, xlabelsize=18, ylabelsize=18)
axis_concentration = axisgrid[2, 3] = Axis(fig, title="Real-time plot of the chlorine concentration", xlabel = "Time [s]", ylabel = "Chlorine concentration [mM]", titlesize=20, xlabelsize=18, ylabelsize=18)

limits!(axis_histogram, 0, 12.5, 0, 10)

colsize!(axisgrid, 1, Auto(16))
colsize!(axisgrid, 3, Auto(16))

fig[2, 1:4] = buttongrid = GridLayout()

button_start = buttongrid[1, 1] = Button(fig, label="START", font=:bold, height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)
button_reset = buttongrid[1, 2] = Button(fig, label="RESET", font=:bold, height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)
menu_type = buttongrid[1, 3] = Menu(fig, options=["Image", "Single"], height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)
menu_time_res = buttongrid[1, 4] = Menu(fig, options=["256", "64"], height = 40, width = 120, selection_cell_color_inactive=RGBf(0.8, 0.8, 0.8), cell_color_hover=RGBf(0.6, 0.6, 0.6), cell_color_active=RGBf(0.8, 0.8, 0.8), cell_color_inactive_even=RGBf(0.9, 0.9, 0.9), cell_color_inactive_odd=RGBf(0.9, 0.9, 0.9), fontsize=20, textcolor=:black)

tb_file_name = buttongrid[2, 1:3] = Textbox(fig, placeholder="Enter a file name", textcolor=RGBf(0.4, 0.4, 0.4), textcolor_placeholder=RGBf(0.6, 0.6, 0.6),height = 40, width = 397, borderwidth=2, bordercolor=RGBf(0.8, 0.8, 0.8), bordercolor_focused=RGBf(0.4, 0.4, 0.4), bordercolor_hover=RGBf(0.6, 0.6, 0.6), boxcolor=:white, boxcolor_focused=:white, boxcolor_hover=:white, fontsize=20, cornerradius=8)
button_save = buttongrid[2, 4] = Button(fig, label="Save", height = 40, width = 120, buttoncolor=RGBf(0.8, 0.8, 0.8), buttoncolor_hover=RGBf(0.6, 0.6, 0.6), buttoncolor_active=RGBf(0.4, 0.4, 0.4), fontsize=20, cornerradius=8, labelcolor_active=:black)

on(button_start.clicks) do _
    ran[] = true

    if running[] == false
        running[] = true
        println("Button pressed")

        start_real_time(folderpath, menu_type.selection[], parse(Int64, menu_time_res.selection[]), axis_histogram, axis_lifetime, axis_photons, axis_concentration)
    end
end

on(button_reset.clicks) do _
    running[] = false

    empty!(axis_histogram)
    empty!(axis_lifetime)
    empty!(axis_photons)
    empty!(axis_concentration)
end

on(button_save.clicks) do _
    # if running[] == false && ran[] == true
    #     filename = tb_file_name.displayed_string[]*".csv"

    #     if filename == "Enter a file name.csv"
    #         filename = split(folderpath, "/")[end]*".csv"
    #     end

    #     open(filename, "w") do io
    #         writedlm(io, [real_time[] lifetime[] photons[] concentration[]], ",")
    #     end
    # end

    save[] = true
end

fig
