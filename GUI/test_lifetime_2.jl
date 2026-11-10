using ZipFile
using Plots
using ProgressBars
using NaNStatistics
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
cells = 1
acquisitions = ["Iris", "Full"]
# lines = [60, 128, 256, 512]
# lines = [512, 256, 128, 60]
lines = [512, 256, 128, 60]
powers = [1, 5, 10, 50]
repeats = 5
hdf_paths = []
sdt_paths = []
single_paths = []

for cell in 1:cells
    for repeat in 1:repeats
        for line in lines
            # for power in powers
            for i in 1:5
                push!(sdt_paths, joinpath(path, "cell$cell", "$(line)_$repeat", "single_c$i.sdt"))
            end
            # end

            # push!(ROI_paths, joinpath(path, "ROI_$line.roi"))
        end
    end
end

for sdt_path in sdt_paths
    if !isfile(sdt_path)
        println("error")
        println(sdt_path)
    end
end

# path = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie2/Données/2025-11-25/10/Full/512/50/image_c1.sdt"

# vector, array, image = sdt_to_vec_2(path)

# heatmap(image, aspect_ratio=:equal, color = :greys)

# n = 4

# pooled_image = pool2d(image, n, n)
# pooled_array = pool3d(array, n, n)

# lifetime_image = zeros(size(pooled_image))

# println(compute_lifetime!(lifetime_image, pooled_array, pooled_image))
# println(vec_to_lifetime(vector)[1][1])

# heatmap(lifetime_image, aspect_ratio=:equal, color = :greys, background_color_inside=:red)
# println(sum(image)/(80000000*10)/(1024*512))
# heatmap(image, aspect_ratio=:equal, color = :greys, background_color_inside=:red)
# /(80000000*10)*1000000
# plot(dropdims(nanmean(lifetime_image, dims=1), dims=1))
# plot!(dropdims(nanmean(lifetime_image, dims=2), dims=2))

# image_lifetimes = Float64[]
# image_mean_lifetimes = Float64[]
# single_lifetimes = Float64[]

# wrapper compute_lifetime! pour accepter callback (voir plus bas)
# compute_lifetime!(lifetime_image, pooled_array, pooled_image, guess; update_inner=nothing)

####################################################

# for i in sdt_paths
#     println(i)
# end

println("Starting processing...")

values = []
single = []
plots_Cl0 = []

for path in ProgressBar(sdt_paths)
    vector, array, image = sdt_to_vec_2(path)

    # n = 8
    # pooled_image = pool2d(image, n, n)
    # pooled_array = pool3d(array, n, n)
    # lifetime_image = zeros(size(pooled_image))

    # push!(values, compute_lifetime!(lifetime_image, pooled_array, pooled_image))
    params, data = vec_to_lifetime(vector)

    fit = conv_irf_data(data[1], Tuple(params), irf, histogram_resolution=256)*sum(data[2])

    push!(single, params[1])
    push!(plots_Cl0, plot(data[1], [data[2], fit]))
end

plot(0,0)
# plot!(values[begin:length(sdt_paths)÷2], label="Pooled Full")
plot!(single[begin:length(sdt_paths)÷2], label="Single Full")
# plot!(values[length(sdt_paths)÷2+1:end], label="Pooled Iris")
plot!(single[length(sdt_paths)÷2+1:end], label="Single Iris")

Cl0 = plot([single[1:length(sdt_paths)÷2], single[length(sdt_paths)÷2+1:end]], label=["Full" "Iris"])

Cl0_power = [plot([single[1:3:length(sdt_paths)÷2], single[length(sdt_paths)÷2+1:3:end]], title="Lifetime vs. number of lines\nPower = 1%", label=["Full" "Iris"], xticks = ([1, 2, 3, 4], ["60", "128", "256", "512"]), xlabel="Number of lines", ylabel="Lifetime [ns]"),
             plot([single[2:3:length(sdt_paths)÷2], single[length(sdt_paths)÷2+2:3:end]], title="Lifetime vs. number of lines\nPower = 5%", label=["Full" "Iris"], xticks = ([1, 2, 3, 4], ["60", "128", "256", "512"]), xlabel="Number of lines", ylabel="Lifetime [ns]"),
             plot([single[3:3:length(sdt_paths)÷2], single[length(sdt_paths)÷2+3:3:end]], title="Lifetime vs. number of lines\nPower = 10%", label=["Full" "Iris"], xticks = ([1, 2, 3, 4], ["60", "128", "256", "512"]), xlabel="Number of lines", ylabel="Lifetime [ns]")
]

Cl0_lines = [plot([single[1:3], single[length(sdt_paths)÷2+1:length(sdt_paths)÷2+3]], title="Lifetime vs. power\nNumber of lines = 60", label=["Full" "Iris"], xticks = ([1, 2, 3], ["1", "5", "10"]), xlabel="Power [%]", ylabel="Lifetime [ns]"),
             plot([single[4:6], single[length(sdt_paths)÷2+4:length(sdt_paths)÷2+6]], title="Lifetime vs. power\nNumber of lines = 128", label=["Full" "Iris"], xticks = ([1, 2, 3], ["1", "5", "10"]), xlabel="Power [%]", ylabel="Lifetime [ns]"),
             plot([single[7:9], single[length(sdt_paths)÷2+7:length(sdt_paths)÷2+9]], title="Lifetime vs. power\nNumber of lines = 256", label=["Full" "Iris"], xticks = ([1, 2, 3], ["1", "5", "10"]), xlabel="Power [%]", ylabel="Lifetime [ns]"),
             plot([single[10:12], single[length(sdt_paths)÷2+10:length(sdt_paths)÷2+12]], title="Lifetime vs. power\nNumber of lines = 512", label=["Full" "Iris"], xticks = ([1, 2, 3], ["1", "5", "10"]), xlabel="Power [%]", ylabel="Lifetime [ns]")
]

####################################################

# println(vec_to_lifetime(sdt_to_vec(path))[1][1])

# println("Starting processing...")

# values_full = []
# values_iris = []

# for idx in ProgressBar(1:length(single_paths)÷2)
#     # vector_full, array_full, image_full = sdt_to_vec_2(sdt_paths[idx])
#     # vector_iris, array_iris, image_iris = sdt_to_vec_2(sdt_paths[idx + length(sdt_paths)÷2])

#     # n = 32
#     # pooled_image_full = pool2d(image_full, n, n)
#     # pooled_array_full = pool3d(array_full, n, n)
#     # lifetime_image_full = zeros(size(pooled_image_full))
#     # pooled_image_iris = pool2d(image_iris, n, n)
#     # pooled_array_iris = pool3d(array_iris, n, n)
#     # lifetime_image_iris = zeros(size(pooled_image_iris))

#     push!(values_full, vec_to_lifetime(sdt_to_vec(single_paths[idx]))[1][1])
#     push!(values_iris, vec_to_lifetime(sdt_to_vec(single_paths[idx + length(single_paths)÷2]))[1][1])
# end

# # plot(values_full)
# plot!(values_iris)

####################################################

# for path in ProgressBar(sdt_paths)
#     vector, array, image = sdt_to_vec_2(path)

#     # heatmap(image, aspect_ratio=:equal, color = :greys)

#     n = 16

#     pooled_image = pool2d(image, n, n)
#     pooled_array = pool3d(array, n, n)

#     lifetime_image = zeros(size(pooled_image))

#     push!(image_lifetimes, compute_lifetime!(lifetime_image, pooled_array, pooled_image))
#     push!(image_mean_lifetimes, vec_to_lifetime(vector)[1][1])
# end

# plot(image_lifetimes)
# plot!(image_mean_lifetimes)

