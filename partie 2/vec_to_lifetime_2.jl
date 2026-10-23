nanmean(x) = mean(filter(!isnan, x))
nanmean(x, dims) = mapslices(nanmean, x, dims=dims)

function vector_moving_average(vector; operator_size=1)
    n = length(vector)
    result = zeros(n)
    for i in 1:n
        lo = max(1, i - operator_size)
        hi = min(n, i + operator_size)
        result[i] = nanmean(vector[lo:hi])
    end
    return result
end

function lifetime_estimate(counts; bin_size=0.039, tol=0.5)
    i0 = argmax(counts)
    mean_tail = mean(counts[i0:end])
    sm = vector_moving_average(counts, operator_size=2)
    for (j, v) in enumerate(sm[i0:end])
        if abs(v - mean_tail) <= tol
            return (i0 + j - 1) * bin_size
        end
    end
    return lifetime_estimate(counts; bin_size=bin_size, tol=tol*2)
end

function irf_shift(irf_data, shift)
    if isnan(shift)
        shift = 0
    end
    n = length(irf_data)
    f = floor(Int, shift)
    c = ceil(Int, shift)
    w1, w2 = 1 - (shift - f), (shift - f)
    # circular shift indices
    idx1 = mod1.(collect(1:n) .- f, n)
    idx2 = mod1.(collect(1:n) .- c, n)
    return w1 .* irf_data[idx1] .+ w2 .* irf_data[idx2]
end

function convolve(irf::Vector{<:Real}, decay::Vector{<:Real}, y_data; histogram_resolution=256)
    N = max(length(irf), length(decay))
    pad_irf = length(irf) < N ? vcat(irf, zeros(N-length(irf))) : irf[1:N]
    pad_dec = length(decay) < N ? vcat(decay, zeros(N-length(decay))) : decay[1:N]
    y = real.(ifft(fft(pad_irf) .* fft(pad_dec)))
    # normalize to data
    y .*= (sum(y_data) / sum(y[1:histogram_resolution]))
    return y
end

function conv_irf_data(x, params::NTuple{3,Real}, irf, y_data; histogram_resolution=256, number_of_previous_pulses=5, laser_pulse_period=12.5)
    t1, d0, yoff = params
    irf_y, irf_x = irf_shift(irf, d0)

    expv = exp.(-x ./ t1)
    for k in 1:number_of_previous_pulses
        expv .+= exp.(-(x .+ k*laser_pulse_period) ./ t1)
    end
    # println(irf_y)
    return convolve(irf_y, expv, y_data; histogram_resolution=histogram_resolution) .+ yoff
end

function MLE_model_func(params, x, y, irf, gating, hist_res, nprev, period)
    decay = conv_irf_data(x, (params[1], params[2], params[3]), irf, y; histogram_resolution=hist_res,
                           number_of_previous_pulses=nprev, laser_pulse_period=period) .* gating
    replace!((v)->v<=0 ? 1e-256 : v, y)
    replace!((v)->v<=0 ? 1e-256 : v, decay)
    ln = log.(y ./ decay)
    return 2*sum(y .* ln .- y .+ decay) / (hist_res - length(params))
end

function MLE_iterative_reconvolution(irf, data_xy; params=[2.0,1.0], gating=ones(UInt8,256), histogram_resolution=256,
                                     number_of_previous_pulses=5, laser_pulse_period=12.5)
    x, y = data_xy
    y = y .* gating
    m = argmax(y)
    offs = max(floor(mean(y[1:m-10]))+ceil(std(y[1:m-10])), 5)
    p = vcat(params, offs)
    bounds = [(2e-3,12.5),(-16,64),(0,offs)]
    res = sp_o.minimize(MLE_model_func, p, args=(x,y,irf,gating,histogram_resolution,number_of_previous_pulses,laser_pulse_period), bounds=bounds, method="L-BFGS-B")
    return res["x"]
end


function vec_to_lifetime(x, irf; bin_size=0.5, hist_size_threshold=8, method="MLE",
                         guess=[2.0,1.0], laser_pulse_period=12.5, histogram_resolution=256,
                         number_of_previous_pulses=5, tcspc_window_size=12.5,
                         tac_low_cut=5.1, tac_high_cut=94.9)
    v = copy(x)

    total_ch = Int(round(laser_pulse_period*histogram_resolution/tcspc_window_size))
    resize!(v, total_ch)

    x_data = collect(bin_size:bin_size:total_ch*bin_size)

    if sum(v) < hist_size_threshold
        return method=="Intensity" ? sum(v) : (method=="Estimate" ? ([0], [x_data, v]) : NaN)
    end

    if method=="Intensity"
        return sum(v)
    end
    
    gate = ones(UInt8, total_ch)
    gate[1:round(Int,tac_low_cut/100*histogram_resolution)] .= 0
    gate[round(Int,tac_high_cut/100*histogram_resolution):end] .= 0
    data_xy = ([x_data, v])

    return MLE_iterative_reconvolution(irf, data_xy; params=guess, gating=gate,
                                          histogram_resolution=histogram_resolution,
                                          number_of_previous_pulses=number_of_previous_pulses,
                                          laser_pulse_period=laser_pulse_period)[1], data_xy
end

