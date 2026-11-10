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

files = readdir(path)

# for file in files
#     println(file)
#     filepath = joinpath(path, file)
#     vector, array, image = sdt_to_vec_2(filepath)

# end


# filepath = joinpath(path, files[25])
# vector, array, image = sdt_to_vec_2(filepath)

circ_mask_1 = [ (i - 140)^2 + (j - 205)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_2 = [ (i - 110)^2 + (j - 535)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_3 = [ (i - 100)^2 + (j - 875)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_4 = [ (i - 265)^2 + (j - 165)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_5 = [ (i - 265)^2 + (j - 545)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_6 = [ (i - 250)^2 + (j - 895)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_7 = [ (i - 430)^2 + (j - 205)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_8 = [ (i - 410)^2 + (j - 540)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]
circ_mask_9 = [ (i - 395)^2 + (j - 895)^2 <= 50^2
             for i in 1:512, j in 1:1024 ]


full_image = zeros(512, 1024)

for i in 1:9
    filepath = joinpath(path, files[i*3])
    vector, array, image = sdt_to_vec_2(filepath)
    circ_mask = eval(Symbol("circ_mask_$i"))
    full_image .+= image .* circ_mask

    # calculate lifetime for masked region
    masked_array = array .* circ_mask
    region_vector = sum(masked_array, dims=(1,2))[:]
    params, data = vec_to_lifetime(region_vector)
    println("Lifetime for region $i: ", params[1])
end

heatmap(full_image)