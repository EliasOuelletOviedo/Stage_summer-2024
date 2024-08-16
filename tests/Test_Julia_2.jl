using GLMakie
using Observables
using Random
using LsqFit

fig = Figure(size = (600, 400))
set_theme!(backgroundcolor = :white)

# time = Observable(0.0)
# x = range(0, 10, length = 1000)
# y = @lift(exp.(-(0.25+$time/10).*x))
# y_0 = exp.(-0.25.*x)
# tau = @lift([0, 1/(0.25+$time/10), 1/(0.25+$time/10)])
# tau_0 = [4, 4]

# xticks = @lift(([0, 1/(0.25+$time/10), 4, 5, 10], ["0", L"\tau", L"\tau_0", "5", "10"]))

# ax = Axis(fig[1, 1], title="Effects of an augmentation of chloride concentration on MQAE lifetime", xticks=xticks, yticks=([0, exp(-1), 0.5, 1], ["0", L"1/e", "0.5", "1"]), xgridvisible=false, ygridvisible=false, xlabel="Time [ns]", ylabel="Normalized fluorescence intensity")

# lines!(ax, x, y_0, color=:black, linestyle=:dash, linewidth=2)
# lines!(ax, tau_0, [exp(-1), 0], color=:black, linestyle=:dot, linewidth=1)
# lines!(ax, tau, [exp(-1), exp(-1), 0], color=:black, linewidth=2)
# hlines!(ax, [exp(-1)], color=:black, linestyle=:dot, linewidth=1)
# lines!(ax, x, y, color=:red, linewidth=3)
# limits!(ax, 0, 10, 0, 1)

# framerate = 30
# timestamps = range(0, 2.5, step=1/framerate)

# record(fig, "time_animation.mp4", timestamps;
#         framerate = framerate) do t
#     time[] = t
# end

function weighted_random_exponential(tau, min_val, max_val)
    u = rand()
    x = -tau * log(1 - u)
    x = round(Int, x)

    if min_val <= x <= max_val
        return x
    else
        return weighted_random_exponential(tau, min_val, max_val)
    end
end

model(t, p) = p[1] * exp.(-p[2] * t)
p_0 = [0.0, 20.0]

tau_0 = 20.0
min_val, max_val = 1, 64
x_i = range(0, 12.5, length = 65)
x = x_i[2:end]
y = Observable(zeros(64))
fit = Observable(p_0[1] * exp.(-p_0[2] .* x))
tau_x = Observable([0.0, 0.0, 0.0])
tau_y = Observable([0.0, 0.0, 0.0])
tau = Observable(4.0)

xticks = @lift(([0, 5, 10, $tau], ["0", "5", "10", L"\tau"]))

ax = Axis(fig[1, 1], title="Histogram of delay between excitation and emission", xgridvisible=false, xticks=xticks, ygridvisible=false, xlabel="Time [ns]", ylabel="Number of photons")

barplot!(ax, x, y, gap = 0.15, color=Makie.wong_colors()[1])
lines!(ax, tau_x, tau_y, color=:black, linewidth=2)
lines!(ax, x, fit, color=Makie.wong_colors()[6], linewidth=3)
limits!(ax, 0, 12.75, 0, 1)

time = Observable(0.0)
framerate = 20
timestamps = range(0, 5, step=1/framerate)

record(fig, "time_animation_2.mp4", timestamps;
        framerate = framerate) do t
    if t != 0 
        time[] = t

        for i in 1:(round(Int, t^4/2))+1
            value = weighted_random_exponential(tau_0, min_val, max_val)
            y[][value] += 1

            fit_params = curve_fit(model, x, y[], p_0)
            params = fit_params.param

            fit[] = params[1] * exp.(-params[2] .* x)

            tau[] = 1/params[2]
            tau_x[] = [0, 1/params[2], 1/params[2]]
            tau_y[] = [exp(-1)*params[1], exp(-1)*params[1], 0]            
        end

        notify(y)
        notify(tau_y)
        limits!(ax, 0, 12.75, 0, maximum(y[])*1.05)
    end
end
