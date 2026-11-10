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

function MLE_iterative_reconvolution(irf, data_xy; params, gating_function=ones(UInt8, 256), histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
    x_data = vec(data_xy[1])
    y_data = vec(data_xy[2]).*gating_function
    irf = deepcopy(irf)
    params_copy = deepcopy(params)
    min_offset = floor(mean(y_data[1:argmax(y_data)-10]))+ceil(std(y_data[1:argmax(y_data)-10])) # get height before the peak of the decay
    
    if min_offset < 5
        min_offset = 5
    end

    push!(params_copy, min_offset)
    
    bounds = [(0.1, 12.5), (-16, 64), (0, min_offset)]
    res = sp_o.minimize(MLE_model_func, params_copy, args=(x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period), bounds=bounds, method="L-BFGS-B")
    
    if res["x"] == bounds[1][1] || res["x"] == bounds[1][2] || !res["success"] || (res["fun"] > 5 && sum(y_data) < 10000)
        return [NaN, NaN, NaN]
    end

    return res["x"]
end

function MLE_model_func(params, x_data, y_data, irf, gating_function, histogram_resolution, number_of_previous_pulses, laser_pulse_period)
    exp_val_c = conv_irf_data(x_data, (params[1], params[2], params[3]), irf, y_data, histogram_resolution=histogram_resolution, number_of_previous_pulses=number_of_previous_pulses, laser_pulse_period=laser_pulse_period).*gating_function
    output_resolution = histogram_resolution
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, y_data)
    replace!(x->smaller_or_eq_zero(x) ? 1e-256 : x, exp_val_c)
    ln_data = log.(y_data./exp_val_c)
    poisson_deviance = 2*sum(y_data.*ln_data.-y_data.+exp_val_c)/(output_resolution-length(params))

    return poisson_deviance
end

function conv_irf_data(x_data, params, irf, y_data; histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
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
    
    return convolve(irf_y_data, exp_1, y_data, histogram_resolution=histogram_resolution).+y_offset
end

function vec_to_lifetime(x; bin_size=3.9088729e-2, hist_size_threshold=100, method="MLE", guess=[2.0, 1.0], laser_pulse_period=12.5, histogram_resolution=256, number_of_previous_pulses=5, tcspc_window_size=10, tac_low_cut=5.0980392, tac_high_cut=94.901962)
    total_channels = Int(round(laser_pulse_period*histogram_resolution/tcspc_window_size))
    append!(vec(x), zeros(total_channels-histogram_resolution))
    x_data = vec(collect(bin_size:bin_size:total_channels*bin_size))

    if sum(x) < hist_size_threshold
        return sum(x)
    else
        global gating_function = ones(UInt8, total_channels)
        gating_function[1:round(Int, tac_low_cut/100*histogram_resolution)] .= 0
        gating_function[round(Int, tac_high_cut/100*histogram_resolution):total_channels] .= 0

        data_xy = [x_data, vec(x)]
        tau_fit = MLE_iterative_reconvolution(irf, data_xy, params=guess, gating_function=gating_function, histogram_resolution=histogram_resolution, number_of_previous_pulses=number_of_previous_pulses, laser_pulse_period=laser_pulse_period)

        return tau_fit, data_xy
    end
end