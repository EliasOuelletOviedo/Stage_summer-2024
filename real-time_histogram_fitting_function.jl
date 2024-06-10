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
using Plots
# using Gtk4
# using FilePathsBase
# using Gadfly

# gr()

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

@everywhere function vec_to_lifetime(x; bin_size=3.9088729e-2, hist_size_threshold=100, method="MLE", broadcast_call=false, guess=[2.0, 1.0], laser_pulse_period=12.5, histogram_resolution=256, number_of_previous_pulses=5, tcspc_window_size=10, tac_low_cut=5.0980392, tac_high_cut=94.901962)
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

@everywhere function ROI_coordinates(x_region, y_region)
    ROI = [(0, 0)]
    for x in x_region
        for y in y_region
            push!(ROI, (x, y))
        end
    end
    return ROI[2:end]
end

@everywhere function get_ROI_vector(array, ROI)
    ROI_arr = zeros(256)

    for (x, y) in ROI
        ROI_arr += array[:, x + (y-1)*2048]
    end

    return ROI_arr
end

@everywhere function get_ROI_lifetime(filepath, x_region, y_region)
    ROI = ROI_coordinates(x_region, y_region)

    file = SDTFile(filepath)

    array = reshape(file, (256, :))

    ROI_vector = get_ROI_vector(array, ROI)

    return vec_to_lifetime(ROI_vector)
end
    

@everywhere function get_file_lifetime(filepath)
    file = SDTFile(filepath)
    
    array = reshape(file, (256, :))

    vector = sum(array, dims = 2)[:,1]

    return vec_to_lifetime(vector)
end

@everywhere function display_realtime_lifetime(folderpath)
    t = zeros(0)
    lifetime = zeros(0)

    while true
        println("Iteration #"*string(length(lifetime) + 1))

        filepath = folderpath*readdir(folderpath)[end]
        region_lifetime = get_file_lifetime(filepath)[1][1]
        # region_lifetime = get_ROI_lifetime(filepath, 462:562, 206:306)[1][1]

        println("   "*string(region_lifetime)*" ns")
        println()

        push!(t, time())
        push!(t, time()-t[1])
        push!(lifetime, region_lifetime)

        p = plot(t[2:2:end],
            lifetime,
            xlabel = "Time [s]",
            ylabel = "Lifetime [ns]",
            linewidth = 2,
            linecolor = :black,
            legend = false)
        display(p)
    end
end

@everywhere function display_lifetime_folder(folderpath)
    t = zeros(0)
    lifetime = zeros(0)

    for i in tqdm(readdir(folderpath))
        println("Iteration #"*string(length(lifetime) + 1))
        
        if endswith(i, ".sdt") || endswith(i[begin:end-6], ".sdt") || endswith(i[begin:end-8], ".sdt")
        # if endswith(i, ".sdt")
            filepath = joinpath(folderpath, i)

            @time region_lifetime = get_file_lifetime(filepath)[1][1]
            # region_lifetime = get_ROI_lifetime(filepath, 462:562, 206:306)[1][1]

            println("  "*string(region_lifetime)*" ns")
            println()

            push!(t, length(lifetime))
            push!(lifetime, region_lifetime)

            p = plot(t,
                lifetime,
                xlabel = "Time [s]",
                ylabel = "Lifetime [ns]",
                linewidth = 2,
                linecolor = :black,
                legend = false)
            display(p)
        else
            println("Not a .sdt file")
            println()
        end
    end
end

@everywhere function new_file_lifetime(path)
    previous_files = readdir(path)

    t = zeros(0)

    lifetime = zeros(0)

    t_0 = time()

    println("Waiting for file")

    while true
        current_files = readdir(path)

        new_files = setdiff(current_files, previous_files)

        if !isempty(new_files)
            # println("New file detected")
            for new_file in new_files
                if endswith(new_file, ".sdt") || endswith(new_file[begin:end-6], ".sdt") || endswith(new_file[begin:end-8], ".sdt") || endswith(new_file[begin:end-9], ".sdt")

                    file_lifetime = get_file_lifetime(joinpath(path, new_file))[1][1]

                    # println("File lifetime : "*string(file_lifetime)*" ns")
                    # println()

                    push!(t, time()-t_0)
                    push!(lifetime, file_lifetime)

                    p = plot(t,
                        lifetime,
                        xlabel = "Time [s]",
                        ylabel = "Lifetime [ns]",
                        linewidth = 2,
                        linecolor = :black,
                        legend = false)
            
                    display(p)
                    display(p)
                end
            end
        end

        previous_files = current_files
    end
end


@everywhere global irf = get_irf()
@everywhere global irf_bin_size = get_irf_bin_size(irf)
@everywhere global folderpath = get_folder_path()
# @everywhere global folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichier_test/"
# @everywhere global folderpath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/20240607_test_taille_image/"
# @everywhere global filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichiers_sdt/4.sdt"
# @everywhere global filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/Fichier_test/test_60X60X256_2_c001.sdt" 

# file = get_file_lifetime(filepath)

# @everywhere global BenchmarkTools.DEFAULT_PARAMETERS.seconds = 12

# @benchmark region_lifetime = get_ROI_lifetime(filepath, 462:562, 206:306) samples = 10
# @benchmark region_lifetime = get_file_lifetime(filepath) samples = 10

# display_realtime_lifetime(folderpath)
# display_lifetime_folder(folderpath)
new_file_lifetime(folderpath)

# t = zeros(0)
# lifetime = zeros(0)

# while true
    # filepath = folderpath*readdir(folderpath)[end]
    # println(filepath)
    # @time region_lifetime = get_file_lifetime(filepath)
    # # @time region_lifetime = get_ROI_lifetime(filepath, 462:562, 206:306)
    # println(region_lifetime[1][1])

    # push!(t, time())
    # push!(t, time()-t[1])
    # push!(lifetime, region_lifetime[1][1])

    # p = plot(t[2:2:end],
    #     lifetime,
    #     xlabel = "Time [s]",
    #     ylabel = "Lifetime [nm]",
    #     linewidth = 2,
    #     linecolor = :black,
    #     legend = false)
    # display(p)
# end


# for i in tqdm(readdir(folderpath))
#     if i != ".DS_Store"
#         filepath = folderpath*i
        
#         @time region_lifetime = get_ROI_lifetime(filepath, 462:562, 206:306)
#         # @time complete_lifetime = get_file_lifetime(filepath)

#         println(region_lifetime[1][1])
#         # println(complete_lifetime[1][1])
#         println()
#     else
#         println("Not a .sdt file")
#         println()
#     end
# end

# GC.gc()