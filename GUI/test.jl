using Lowess, Plots
using SavitzkyGolay

function ewma!(smoothed, raw, a)
    smoothed[1] = raw[1]

    for i in eachindex(raw)[2:end]
        smoothed[i] = a*raw[i] + (1-a)*smoothed[i-1]
    end
end

alpha = 0.1

xs = 10 .* rand(100)
xs = sort(xs)
ys = sin.(xs) .+ 0.5 * rand(100)

model = lowess_model(xs, ys, 0.1)

us = range(extrema(xs)...; step = 0.1)
vs_1 = model(us)

vs_2 = []

for i in eachindex(ys)
    push!(vs_2, ys[i])

    ewma!(vs_2, ys[1:i], 0.05)
end

vs_3 = savitzky_golay(ys, 21, 1).y

using Lowess

"""
    loess_smooth(x, y; f=0.3, iter=2)

Compute a LOWESS (LOESS) smoothing of the data `(x, y)`.

# Arguments
- `x::AbstractVector{<:Real}` — sorted independent–variable samples  
- `y::AbstractVector{<:Real}` — dependent–variable samples (same length as `x`)  
- `f::Real=0.3` — span (fraction of data used in each local regression)  
- `iter::Integer=2` — number of robustness (re-weighting) iterations  

# Returns
- `y_smooth::Vector{Float64}` — smoothed `y` at the original `x` points  
"""
function loess_smooth(x::AbstractVector{<:Real},
                      y::AbstractVector{<:Real};
                      f::Real = 0.3)
    @assert length(x) == length(y) "x and y must have the same length"
    @assert issorted(x)        "x must be sorted in ascending order"

    # 1. Fit the LOWESS model
    model = lowess_model(x, y, f)

    # 2. Predict smoothed values at the original x
    return model(x)
end

vs_4 = loess_smooth(xs, ys, f=0.1)


scatter(xs, ys)
plot!(us, vs_1, legend=true)
plot!(xs, vs_2, legend=true)
plot!(xs, vs_3, legend=true)
plot!(xs, vs_4, legend=true)
