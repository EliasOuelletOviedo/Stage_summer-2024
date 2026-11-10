using PyCall
using FFTW
using Optim
using Statistics
using LineSearches
using Plots
using BenchmarkTools

sp_o = pyimport("scipy.optimize")

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

function vec_to_lifetime(x; guess=[3.0, 1.0, 1e-6], laser_pulse_period=10.0, histogram_resolution=256, tac_low_cut=5.0980392, tac_high_cut=94.901962, fixed_parameters=Float64[NaN, NaN, NaN])
    """
    Calculate fluorescence lifetime from vector of counts
    """
    total_channels = round(Int, laser_pulse_period*histogram_resolution/tcspc_window_size)
    x_data = vec(collect(irf_bin_size:irf_bin_size:total_channels*irf_bin_size))

    if total_channels != histogram_resolution
        println("Mismatch")
        global fft_plan = plan_fft(zeros(Float64, total_channels))
        global ifft_plan = plan_ifft(zeros(Float64, total_channels))
        new_irf = zeros(Float64, total_channels, 2)
        new_irf[:, 1] = collect(irf_bin_size:irf_bin_size:irf_bin_size*total_channels)
        new_irf[1:length(irf[:, 2]), 2] = irf[:, 2]
        global irf = new_irf
        append!(vec(x), zeros(total_channels-histogram_resolution))
    end
    
    if sum(x) < 10
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

global irf = get_irf()
global irf_bin_size = get_irf_bin_size()
global tcspc_window_size = round(irf[end, 1]+irf[2, 1], sigdigits=4)
global fft_plan = plan_fft(zeros(Float64, 256))
global ifft_plan = plan_ifft(zeros(Float64, 256))

println(tcspc_window_size)

vector = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 142.0, 129.0, 156.0, 121.0, 118.0, 122.0, 122.0, 150.0, 374.0, 740.0, 920.0, 1035.0, 1126.0, 1103.0, 1178.0, 1124.0, 1139.0, 1115.0, 1126.0, 1003.0, 1051.0, 1027.0, 1042.0, 1000.0, 1033.0, 940.0, 987.0, 963.0, 937.0, 891.0, 862.0, 894.0, 895.0, 886.0, 873.0, 856.0, 843.0, 831.0, 822.0, 848.0, 802.0, 782.0, 774.0, 797.0, 757.0, 837.0, 777.0, 739.0, 710.0, 782.0, 766.0, 700.0, 713.0, 716.0, 698.0, 724.0, 671.0, 662.0, 752.0, 645.0, 634.0, 688.0, 701.0, 632.0, 628.0, 625.0, 616.0, 638.0, 622.0, 560.0, 588.0, 585.0, 652.0, 598.0, 570.0, 564.0, 577.0, 579.0, 547.0, 544.0, 536.0, 551.0, 568.0, 523.0, 500.0, 505.0, 522.0, 523.0, 519.0, 521.0, 517.0, 497.0, 500.0, 484.0, 494.0, 452.0, 514.0, 484.0, 476.0, 484.0, 491.0, 423.0, 461.0, 450.0, 433.0, 460.0, 477.0, 452.0, 456.0, 469.0, 409.0, 429.0, 393.0, 375.0, 397.0, 351.0, 411.0, 401.0, 374.0, 371.0, 377.0, 382.0, 380.0, 379.0, 363.0, 411.0, 361.0, 398.0, 366.0, 324.0, 363.0, 360.0, 342.0, 305.0, 334.0, 355.0, 308.0, 329.0, 275.0, 323.0, 322.0, 344.0, 308.0, 320.0, 294.0, 322.0, 338.0, 294.0, 298.0, 336.0, 290.0, 270.0, 313.0, 291.0, 288.0, 273.0, 298.0, 268.0, 262.0, 250.0, 269.0, 283.0, 237.0, 255.0, 251.0, 264.0, 268.0, 263.0, 254.0, 254.0, 251.0, 253.0, 266.0, 268.0, 237.0, 230.0, 251.0, 229.0, 206.0, 218.0, 229.0, 243.0, 243.0, 216.0, 231.0, 227.0, 234.0, 202.0, 199.0, 213.0, 214.0, 234.0, 212.0, 212.0, 216.0, 211.0, 215.0, 201.0, 212.0, 178.0, 201.0, 220.0, 195.0, 200.0, 205.0, 206.0, 199.0, 183.0, 186.0, 171.0, 175.0, 197.0, 186.0, 187.0, 169.0, 174.0, 169.0, 157.0, 159.0, 162.0, 179.0, 153.0, 168.0, 174.0, 161.0, 146.0, 179.0, 153.0, 156.0, 153.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
println(length(vector))
data = vec_to_lifetime(vector)
data = vec_to_lifetime(vector)

# println(data)

params, xy_data = data
println(params)

test_params = Tuple(params)

println(test_params)

fit  = conv_irf_data(xy_data[1], test_params, irf, histogram_resolution=256)
println(fit)

x_data = xy_data[1]
y_data = xy_data[2]

plot(x_data, fit)

# println(typeof(irf))

# @benchmark vec_to_lifetime(vector)