using LsqFit
using Plots
using Statistics
using CSV
using DataFrames
using FileIO
using Images
using Colors

folderpath = "/Users/eliasouellet-oviedo/Desktop/20250718_SCLM_prot/Analyse/Fichiers"
files = readdir(folderpath)

println(files)
# files = files[3:end]

hill(x, p) = (p[1] * x.^p[3]) ./ (exp(p[2])^p[3] .+ x.^p[3]) .+ p[4]
# p = [Vmax, Kd, n, offset]

plot(0, 0, xscale=:log10)

# file = files[3]

# for file in files
#     filepath = joinpath(folderpath, file)
#     df = CSV.read(filepath, DataFrame)
#     lifetime = df.Lifetime_m_ns
#     # println(lifetime)


#     groups = [lifetime[i:i+2] for i in 1:3:length(lifetime)]  # divise en groupes de 3
#     means = [mean(g) for g in groups]
#     stds  = [std(g) for g in groups]
#     x = [0.01, 20, 40, 146, 100000]
#     x_vals = [0.01,0.01,0.01, 20,20,20, 40,40,40, 146,146,146, 100000, 100000, 100000]

#     # println(means)

#     initial_params = [5.0, 1.0, 1.0, 0]  # Vmax, Kd, n (choix de départ)
#     fit = curve_fit(hill, x_vals, lifetime, initial_params; lower = [0.0, 0.0, 0.0, 0.0], upper = [Inf, Inf, Inf, Inf])
#     params = fit.param  # les paramètres ajustés

#     println(exp(params[2]))

#     xfit = 10 .^ range(-2, stop=5, length=200)
#     yfit = hill(xfit, params)

#     # plot!(x, means .- means[1], yerror=stds, seriestype=:scatter, title = "Lifetime of $(replace(file,".csv" => "" ))", label=nothing, xlabel="[Cl⁻] [mM]", ylabel="τ [ns]")
#     plot!(x, means, yerror=stds, seriestype=:scatter, title = "Lifetime variation of different SCLM variant", label="$(replace(file,".csv" => "" ))", xlabel="[Cl⁻] [mM]", ylabel="Δτ [ns]",
#     xticks = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "10⁵"]), legend_position = :topleft)

#     plot!(xfit, yfit, label = "Fit Hill", lw=2)
# end

filepath = joinpath(folderpath, files[1])
df = CSV.read(filepath, DataFrame)
lifetime = df.Lifetime_m_ns
# println(lifetime)


groups = [lifetime[i:i+2] for i in 1:3:length(lifetime)]  # divise en groupes de 3
means = [mean(g) for g in groups]
stds  = [std(g) for g in groups]
x = [0.01, 20, 40, 146, 100000]
x_vals = [0.01,0.01,0.01, 20,20,20, 40,40,40, 146,146,146, 100000, 100000, 100000]

# println(means)

initial_params = [5.0, 1.0, 1.0, 0]  # Vmax, Kd, n (choix de départ)
fit = curve_fit(hill, x_vals, lifetime, initial_params; lower = [0.0, 0.0, 0.0, 0.0], upper = [Inf, Inf, Inf, Inf])
params = fit.param  # les paramètres ajustés

# println(exp(params[2]))
# println(minimum(means), " ", maximum(means))

println("SCLM")
println("Δτ     = ", round(params[1], digits = 2))
println("Δτ [%] = ", round(params[1]/(params[1] + params[4])*100, digits = 2), " %")
println("τ_0    = ", round(params[4], digits = 2))
println("EC50      = ", round(exp(params[2]), digits = 2))

xfit = 10 .^ range(-2, stop=5, length=200)
yfit = hill(xfit, params)

# plot!(x, means .- means[1], yerror=stds, seriestype=:scatter, title = "Lifetime of $(replace(files[1],".csv" => "" ))", label=nothing, xlabel="[Cl⁻] [mM]", ylabel="τ [ns]")
plot!(x, means, yerror=stds, seriestype=:scatter, title = "Lifetime variation of different SCLM variant", label="$(replace(files[1],".csv" => "" ))", xlabel="[Cl⁻] [mM]", ylabel="Δτ [ns]",
xticks = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "10⁵"]), legend_position = :topleft, color = 1)

plot!(xfit, yfit, label = "Fit Hill", lw=2, color = 1)

filepath = joinpath(folderpath, files[3])
df = CSV.read(filepath, DataFrame)
lifetime = df.Lifetime_m_ns
# println(lifetime)


groups = [lifetime[i:i+2] for i in 1:3:length(lifetime)]  # divise en groupes de 3
means = [mean(g) for g in groups]
stds  = [std(g) for g in groups]
x = [0.01, 20, 40, 146, 100000]
x_vals = [0.01,0.01,0.01, 20,20,20, 40,40,40, 146,146,146, 100000, 100000, 100000]

# println(means)

initial_params = [5.0, 1.0, 1.0, 0]  # Vmax, Kd, n (choix de départ)
fit = curve_fit(hill, x_vals, lifetime, initial_params; lower = [0.0, 0.0, 0.0, 0.0], upper = [Inf, Inf, Inf, Inf])
params = fit.param  # les paramètres ajustés

# println(exp(params[2]))
# println(minimum(means), " ", maximum(means))

println("SCLM_mCer3")
println("Δτ     = ", round(params[1], digits = 2))
println("Δτ [%] = ", round(params[1]/(params[1] + params[4])*100, digits = 2), " %")
println("τ_0    = ", round(params[4], digits = 2))
println("EC50      = ", round(exp(params[2]), digits = 2))

xfit = 10 .^ range(-2, stop=5, length=200)
yfit = hill(xfit, params)

# plot!(x, means .- means[1], yerror=stds, seriestype=:scatter, title = "Lifetime of $(replace(files[3],".csv" => "" ))", label=nothing, xlabel="[Cl⁻] [mM]", ylabel="τ [ns]")
plot!(x, means, yerror=stds, seriestype=:scatter, title = "Lifetime of different SCLM variant", label="$(replace(files[3],".csv" => "" ))", xlabel="[Cl⁻] [mM]", ylabel="τ [ns]",
xticks = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "10⁵"]), legend_position = :topleft, color = 2)

plot!(xfit, yfit, label = "Fit Hill", lw=2, color = 2)

plot(0,0, xscale=:log10)

folderpath = "/Users/eliasouellet-oviedo/Desktop/20250718_SCLM_prot"
subfolders = ["SCLM", "SCLM_mCer3", "SCLM_A206K", "mCer3"]
# subfolders = ["SCLM", "SCLM_mCer3"]
concentrations = ["0Cl-00", "20Cl-00", "40Cl-00", "146Cl-00", "146F-00"]

mask = [ (i - 270)^2 + (j - 567)^2 <= 160^2
         for i in 1:512, j in 1:1024 ]

for (i, subfolder) in enumerate(subfolders)
    subfolderpath = joinpath(folderpath, subfolder)
    files = readdir(subfolderpath)

    data = Float64[]

    for concentration in concentrations
        concentration_files = filter(f -> occursin(concentration, f), files)

        # ratios = Float64[]

        for i in 1:3
            # println(concentration_files[i])
            filepath_1 = joinpath(folderpath, subfolder, concentration_files[i], "Bliq VMS/C1")
            filepath_2 = joinpath(folderpath, subfolder, concentration_files[i], "Bliq VMS/C2")
            file_1 = readdir(filepath_1)
            file_2 = readdir(filepath_2)
            filepath_1 = joinpath(filepath_1, file_1[1])
            filepath_2 = joinpath(filepath_2, file_2[1])

            # println(filepath_1)
            data_1 = load(filepath_1)
            data_2 = load(filepath_2)

            # data_1[.!mask] .= zero(eltype(data_1))
            # data_2[.!mask] .= zero(eltype(data_2))

            # sum_1 = sum(data_1)
            # sum_2 = sum(data_2)

            # ratio = sum_1 / sum_2

            in_1 = data_1[mask]
            out_1 = data_1[.!mask]
            in_2 = data_2[mask]
            out_2 = data_2[.!mask]
            ratio = (mean(in_1) - mean(out_1)) / (mean(in_2) - mean(out_2))

            push!(data, ratio)
        end

        # push!(data, ratios)

        # println("Pour une concentration de $concentration, moyenne de :", mean(ratios))
    end

    groups = [data[i:i+2] for i in 1:3:length(data)]  # divise en groupes de 3
    means = [mean(g) for g in groups]
    stds  = [std(g) for g in groups]
    x = [0.01, 20, 40, 146, 100000]
    x_vals = [0.01,0.01,0.01, 20,20,20, 40,40,40, 146,146,146, 100000, 100000, 100000]

    initial_params = [5.0, 1.0, 1.0, 0]  # Vmax, Kd, n (choix de départ)
    # fit = curve_fit(hill, x_vals, data/0.3147782107623443, initial_params)
    fit = curve_fit(hill, x_vals, data, initial_params)
    params = fit.param  # les paramètres ajustés

    println("FLIM $i ", params)
    # println(exp(params[2]))

    println(subfolder)
    println("ΔFRET     = ", round(params[1], digits = 2))
    println("ΔFRET [%] = ", round(params[1]/(params[1] + params[4])*100, digits = 2), " %")
    println("FRET_0    = ", round(params[4], digits = 2))
    println("EC50      = ", round(exp(params[2]), digits = 2))

    xfit = 10 .^ range(-2, stop=5, length=200)
    yfit = hill(xfit, params)

    plot!(x, means, yerror=stds, seriestype=:scatter, title = "Normalized FRET ratio of different variations of SCLM", label="$subfolder", xlabel="[Cl⁻] [mM]", ylabel="Normalized FRET ratio",
    xticks = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "10⁵"]), legend_position = :topright, color = i)
    plot!(xfit, yfit, label = "Fit Hill", lw=2, color = i)
    # println(data)
end

# ylims!(0, 7)
plot!(0,0)

# filepath = "/Users/eliasouellet-oviedo/Desktop/20250718_SCLM_prot/SCLM_mCer3/20Cl-003/Bliq VMS/C1/20Cl-003-C1-T001.tif"

# data = load(filepath)*100

# # data[data .< 1] .= zero(eltype(data))

# mask = [ (i - 270)^2 + (j - 567)^2 <= 160^2
#          for i in 1:512, j in 1:1024 ]

# data[.!mask] .= zero(eltype(data))

# println(typeof(data))

# plot(data)
