using ZipFile
using Plots
# using ProgressBars
using NaNStatistics
using ProgressMeter
# using ProgressMeter
# using Term.Progress

include("vec_to_lifetime.jl")

global irf               = get_irf()
global irf_bin_size      = get_irf_bin_size()
global tcspc_window_size = round(irf[end, 1] + irf[2, 1], sigdigits=4)
global fft_plan          = plan_fft(zeros(Float64, 256))
global ifft_plan         = plan_ifft(zeros(Float64, 256))

sdt = pyimport("sdtfile")

function sdt_to_vec(path)
    py"""
    import numpy as np
    import sdtfile

    file = sdtfile.SdtFile($path)
    vector = np.array(file.data)
    """

    vector = Array(py"vector")
    # println(array)
    return vector
end

function sdt_to_vec_2(path)
    py"""
    import numpy as np
    import sdtfile

    file = sdtfile.SdtFile($path)

    # print(file.header)

    array = np.squeeze(np.array(file.data)).reshape((-1, 1024, 256))[::2,:,:]
    image = np.sum(array, axis = 2)
    vector = np.sum(array, axis = (0,1))
    """

    vector, array, image = Array(py"vector"), Array(py"array"), Array(py"image")
    # println(array)
    return (vector, array, image)
end

function pool2d(A::AbstractArray, ph::Int, pw::Int; op=sum)
    H,W = size(A)
    @assert H % ph == 0 && W % pw == 0
    A2 = reshape(A, ph, div(H,ph), pw, div(W,pw))
    B = op(A2, dims=(1,3))
    return dropdims(B, dims=(1,3))
end

function pool3d(A::AbstractArray, ph::Int, pw::Int; op=sum)
    H,W,T = size(A)
    @assert H % ph == 0 && W % pw == 0
    A2 = reshape(A, ph, div(H,ph), pw, div(W,pw), T)
    B = op(A2, dims=(1,3))
    return dropdims(B, dims=(1,3))
end

function compute_lifetime!(lifetime_image, pooled_array, pooled_image)
    for I in ProgressBar(CartesianIndices(pooled_image))
        params, data = vec_to_lifetime(pooled_array[I,:])
        lifetime_image[I] = params[1]
    end

    print("\e[A\e[2K")

    return nanmean(lifetime_image)
end

# function compute_lifetime!(lifetime_image, pooled_array, pooled_image; update_inner=nothing)
#     for I in CartesianIndices(pooled_image)
#         params, data = vec_to_lifetime(pooled_array[I,:])
#         lifetime_image[I] = params[1]
#         update_inner !== nothing && update_inner()
#     end

#     return nanmean(lifetime_image)
# end


path = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie2/Données/2026-03-12/"

ROI_paths = []
cells = 2
repeats = 5
acquisitions = ["Single_1"]
lines = [512, 256, 128, 60]
# lines = ["512_1", "256_1", "128_1", "512_2", "256_2", "128_2", "512_3", "256_3", "128_3"]
files = 5

lifetimes = []
mean_lifetimes = []

# total_iters = length(concentrations) * length(acquisitions) * length(lines) * length(files)
total_iters = cells * repeats * length(lines) * files
p = Progress(total_iters, desc="Processing")

for cell in 1:cells
    # for acquisition in acquisitions
    for repeat in 1:repeats
        for line in lines
            full_vec = zeros(UInt32, 256)

            for file in 1:files
                # filepath = joinpath(path, acquisition, string(line), "single_c$file.sdt")
                # filepath = joinpath(path, string(concentration), acquisition, string(line), "single_c$file.sdt")
                filepath = joinpath(path, "cell$cell", "$(line)_$repeat", "single_c$file.sdt")
                if !isfile(filepath)
                    println("error")
                    println(filepath)
                end
                # println("Processing file: $filepath")
                vector = sdt_to_vec(filepath)
                params, data = vec_to_lifetime(vector)

                full_vec .+= vec(vector)
                push!(lifetimes, params[1])

                next!(p)
            end

            full_params, full_data = vec_to_lifetime(full_vec)

            push!(mean_lifetimes, full_params[1])
        end
    end
    # end
end



plot(mean_lifetimes)

# x1 = [1, 2, 3, 4]
# x2 = [5, 6, 7, 8]
# x3 = [9, 10, 11, 12]

# # x1 = [512, 256, 128, 60]
# # x2 = [512, 256, 128, 60]
# # x3 = [512, 256, 128, 60]

# # println(length(mean_lifetimes))

# plot(0,0, xlabel="Temps écoulé", ylabel="Temps de vie [ns]")
# # plot(0,0, xlabel="Nombre de lignes", ylabel="Temps de vie [ns]")
# scatter!(repeat(x1, 3), lifetimes[1:12],  markersize=3, label="Série 1", color=:blue)
# scatter!(repeat(x2, 3), lifetimes[13:24], markersize=3, label="Série 2", color=:red)
# scatter!(repeat(x3, 3), lifetimes[25:36], markersize=3, label="Série 3", color=:green)
# plot!(x1, mean_lifetimes[1:4],  marker=:circle, label="Série 1 (valeur moyenne)", color=:blue)
# plot!(x2, mean_lifetimes[5:8],  marker=:circle, label="Série 2 (valeur moyenne)", color=:red)
# plot!(x3, mean_lifetimes[9:12], marker=:circle, label="Série 3 (valeur moyenne)", color=:green)

# plot(0,0)
# cell1_sweep1 = scatter!(x, mean_lifetimes[1:3])
# cell1_sweep2 = scatter!(x, mean_lifetimes[4:6])
# cell1_sweep3 = scatter!(x, mean_lifetimes[7:9])

# cell2_sweep1 = scatter!(x, mean_lifetimes[10:12])
# cell2_sweep2 = scatter!(x, mean_lifetimes[13:15])
# cell2_sweep3 = scatter!(x, mean_lifetimes[16:18])

# scatter(mean_lifetimes[1:12])

# Full_0  = scatter!(lines, mean_lifetimes[1                            : 1*length(mean_lifetimes)÷4], label="Concentration : 0,  Acquisition : Full")
# Iris_0  = scatter!(lines, mean_lifetimes[1*length(mean_lifetimes)÷4+1 : 2*length(mean_lifetimes)÷4], label="Concentration : 0,  Acquisition : Iris")
# Full_80 = scatter!(lines, mean_lifetimes[2*length(mean_lifetimes)÷4+1 : 3*length(mean_lifetimes)÷4], label="Concentration : 80, Acquisition : Full")
# Iris_80 = scatter!(lines, mean_lifetimes[3*length(mean_lifetimes)÷4+1 : end],                        label="Concentration : 80, Acquisition : Iris")

# Full_0_2  = scatter!(repeat(lines, inner=length(files)), lifetimes[1                       : 1*length(lifetimes)÷4], label="Concentration : 0,  Acquisition : Full")
# Iris_0_2  = scatter!(repeat(lines, inner=length(files)), lifetimes[1*length(lifetimes)÷4+1 : 2*length(lifetimes)÷4], label="Concentration : 0,  Acquisition : Iris")
# Full_80_2 = scatter!(repeat(lines, inner=length(files)), lifetimes[2*length(lifetimes)÷4+1 : 3*length(lifetimes)÷4], label="Concentration : 80, Acquisition : Full")
# Iris_80_2 = scatter!(repeat(lines, inner=length(files)), lifetimes[3*length(lifetimes)÷4+1 : end],                   label="Concentration : 80, Acquisition : Iris")

# plot(0, 0)

# scatter(repeat(lines[1:4], inner=length(files)), lifetimes[1 : 4*length(lifetimes)÷10], label="Concentration : 0,  Acquisition : Full")
# scatter!(repeat([lines[5]], inner=length(files)), lifetimes[4*length(lifetimes)÷10+1 : 1*length(lifetimes)÷2], label="Concentration : 0,  Acquisition : Full")

# scatter!(repeat(lines[1:4], inner=length(files)), lifetimes[1*length(lifetimes)÷2+1 : 1*length(lifetimes)÷2+2*length(lifetimes)÷5], label="Concentration : 0,  Acquisition : Iris")
# scatter!(repeat([lines[5]], inner=length(files)), lifetimes[1*length(lifetimes)÷2+2*length(lifetimes)÷5+1 : end], label="Concentration : 0,  Acquisition : Iris")

# scatter(repeat(lines, inner=length(files)), lifetimes[1 : 1*length(lifetimes)÷2], label="Concentration : 0,  Acquisition : Full")
# scatter!(repeat(lines, inner=length(files)), lifetimes[1*length(lifetimes)÷2 : end], label="Concentration : 0,  Acquisition : Iris")

# mean_0 = mean(lifetimes[1 : 1*length(lifetimes)÷2])
# mean_80 = mean(lifetimes[1*length(lifetimes)÷2 : end])

# b = 1/mean_0
# a = (1/mean_80 - b)/80

# concentrations = (lifetimes.^-1 .- b) ./ a

