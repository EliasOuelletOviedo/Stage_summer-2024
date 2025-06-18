using Plots
using CSV
using DataFrames
using LsqFit
using Statistics
using StatsBase
# using CurveFit

function exponential_model(x, p)
    y = p[1] .* exp.(-p[2] .* (x .- p[3])) .+ p[4]
    y .= ifelse.(x .< 0, p[1], y)
    return y
end

function normal_model(x, p)
    y = 1/(p[2]*sqrt(2*pi)) .* exp.(-0.5 .* ((x .- p[1]) ./ p[2]).^2)
    return y
end

println("Starting analysis")

filepath = "/Users/eliasouellet-oviedo/Desktop/test3.csv.csv"
# filepath = "/Users/eliasouellet-oviedo/Desktop/analyse1.csv.csv"
# filepath = "/Users/eliasouellet-oviedo/Desktop/test2/analyse2.csv"

println("Opening file : ", filepath)

df = CSV.read(filepath, DataFrame)

# times    = df.Time_s[begin:3492]
# counts   = df.Number_of_photons[begin:3492]
# lifetime = df.Lifetime_ns[begin:3492]
# signal   = df.Signal[begin:3492]
# setpoint = df.Setpoint[begin:3492]

times    = df.Time_s
counts   = df.Number_of_photons
lifetime = df.Lifetime_ns
signal   = df.Signal
setpoint = df.Setpoint

setpoint[signal .== 0] .= NaN

starts = Int[]
ends   = Int[]
global in_region = false

for i in eachindex(setpoint)
    if !isnan(setpoint[i]) && !in_region
        push!(starts, i)
        global in_region = true
    elseif (isnan(setpoint[i]) || i == lastindex(setpoint)) && in_region
        push!(ends, isnan(setpoint[i]) ? i-1 : i)
        global in_region = false
    end
end

if in_region
    push!(ends, lastindex(a))
end

times_i = times[starts]
times_f = times[ends]

println(starts)

fits = []
delay = []
offset = []

for i in eachindex(starts)
    i0 = starts[i]
    i1 = ends[i]

    xseg = collect(times[i0:i1])
    yseg = lifetime[i0:i1]

    p0 = [
        maximum(yseg) - minimum(yseg),
        0.1,
        xseg[1],
        minimum(yseg),
    ]

    fit = curve_fit(exponential_model, xseg, yseg, p0)
    push!(fits, fit)
    push!(delay, 1/fit.param[2])
    push!(offset, fit.param[4])
    # println(fit.param)
    println("Descente $i")
    # println("Tau = $(round(fit.param[2], digits=3)) s")
    println("Rise time = $(round(Int, 1/fit.param[2])) s")
    # println("x offset = $(round(fit.param[3] - times_i[i], digits=1)) s")
    println("y offset = $(round(fit.param[4], digits=1)) ns\n")

    i0 = ends[i]
    try
        i1 = starts[i+1]
    catch
        # println("end")
        i1 = lastindex(times)
    end

    xseg = collect(times[i0:i1])
    yseg = lifetime[i0:i1]

    p0 = [
        minimum(yseg) - maximum(yseg),
        0.1,
        xseg[1],
        maximum(yseg),
    ]

    fit = curve_fit(exponential_model, xseg, yseg, p0)
    push!(fits, fit)
    # println(fit.param)
    println("Montée $i")
    # println("Tau = $(round(fit.param[2], digits=3)) s")
    println("Rise time = $(round(Int, 1/fit.param[2])) s")
    # println("x offset = $(round(fit.param[3] - times_f[i], digits=1)) s")
    println("y offset = $(round(fit.param[4], digits=1)) ns\n")
end

for i in eachindex(starts)
    i0 = starts[i]
    i1 = ends[i]

    xseg = collect(times[i0:i1])

    xfit = range(xseg[1], xseg[end], length=200)
    yfit = exponential_model(xfit, coef(fits[2*i-1]))

    i0 = ends[i]
    try
        i1 = starts[i+1]
    catch
        i1 = lastindex(times)
    end

    xseg = collect(times[i0:i1])

    xfit = range(xseg[1], xseg[end], length=200)
    yfit = exponential_model(xfit, coef(fits[2*i]))
end

plot(0, 0, background_color=:gray8, background_color_outside=:gray12, xlabel = "Time [s]", ylabel = "Lifetime [ns]", title = "Lifetime", legend = false, linewidth = 1, color = :royalblue)

span_list = collect(Iterators.flatten(zip(times_i, times_f)))
vspan!(span_list, color = :gold, alpha = 0.2)

plot!(times, lifetime, background_color=:gray8, background_color_outside=:gray12, xlabel = "Time [s]", ylabel = "Lifetime [ns]", title = "Lifetime", legend = false, linewidth = 1, color = :royalblue)
# plot!(times, lifetime, background_color=:gray8, background_color_outside=:gray12, xlabel = "Time [s]", ylabel = "Lifetime [ns]", title = "Lifetime", legend = false, linewidth = 1, color = :transparent)

region = []
# normal_fit = []
# xsegs = []
# ysegs = []

for i in eachindex(starts)
    i0 = starts[i]
    i2 = ends[i]

    # if i == 1
    #     t1 = times[i0] + 3.5*delay[i]
    # elseif i == 4
    #     t1 = times[i0] + 5*delay[i]
    # else
    #     t1 = times[i0] + 4*delay[i]
    # end
    t1 = times[i0] + 3*delay[i]
    i1 = findfirst(x -> x >= t1, times)

    if i1 === nothing
        i1 = lastindex(times)
    end

    if i1 > lastindex(times)
        i1 = lastindex(times)
    end

    xseg = collect(times[i1:i2])
    yseg = lifetime[i1:i2]

    # push!(xsegs, xseg)

    # fit = curve_fit(normal_model, xseg, yseg, [maximum(yseg), mean(xseg), std(xseg)])
    # push!(normal_fit, fit)
    # println(fit.param)

    # push!(ysegs, normal_model(xseg, fit.param))

    push!(region, yseg)

    # println(xseg)

    plot!(xseg, yseg, background_color=:gray8, background_color_outside=:gray12, xlabel = "Time [s]", ylabel = "Lifetime [ns]", title = "Lifetime", legend = false, linewidth = 1, color = :royalblue)
end

# x = range(0, 10, length=100)
# y = sin.(x)
plot!(0, 0, legend = false, title = "Lifetime")

# println(starts)

setpoints = [4.8, 4.9, 5.1, 4.7]
# setpoints = [4.9, 4.9, 4.9, 4.9]

i = 2

mean_value = mean(region[i])
std_value = std(region[i])
RMS = rmsd(region[i], setpoints[i].*ones(1, length(region[i])))
SEM = sem(region[i])
# SEM = std_value/sqrt(length(region[i]))

println("RMS error : $(round(RMS, digits = 3))")

println("SEM : $(round(SEM, digits = 4))")

println("Setpoint: $(setpoints[i]), Mean: $(round(mean_value, digits = 2)), Std: $(round(std_value, digits = 3))")

x = range(mean_value - 5*std_value, mean_value + 5*std_value, length=length(region[i]))
bins = range(mean_value - 5*std_value, mean_value + 5*std_value, length=16)
# println(bins)

histogram(region[i], bins = bins, normalize=:pdf, background_color=:gray8, background_color_outside=:gray12, xlabel = "Lifetime", ylabel = "Probability Density", label = "Region $i", title = "Lifetime Distribution", legend = :topright, linewidth = 1, color = :royalblue)
vline!([mean(region[i])], linewidth = 4, label = "Mean", color = :darkgreen)
vline!([setpoints[i]], linewidth = 4, label = "Setpoint", color = :brown)
plot!(x, normal_model(x, [mean_value, std_value]), label = "Normal fit", linewidth = 4, color = :chocolate)
