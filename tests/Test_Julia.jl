using Random
using Plots
using LsqFit

function weighted_random_exponential(tau::Float64, min_val::Int, max_val::Int)
    u = rand() 
    x = -tau * log(1 - u)
    x = round(Int, x)
    if min_val <= x <= max_val
        return x
    else
        return weighted_random_exponential(tau, min_val, max_val)
    end
end

nb_bins = 5

model(t, p) = p[1] * exp.(-p[2] * t)
p_0 = [1.0, 2.0]

tau = nb_bins/2
min_val, max_val = 1, nb_bins

x = range(1, nb_bins, length = 100)
y = zeros(nb_bins)

for i in 1:20*nb_bins
    value = weighted_random_exponential(tau, min_val, max_val)

    if min_val <= value <= max_val
        y[value] += 1
    end
end

fit_params = curve_fit(model, 1:nb_bins, y, p_0)
params = fit_params.param

fit = params[1] * exp.(-params[2] .* x)

bar(y, grid = false, legend=false, label=false, ticks = false, thickness_scaling=3, ylim=(0.0, maximum(fit)*1.1), size = (400,300))
plot!(x, fit, linewidth=10)