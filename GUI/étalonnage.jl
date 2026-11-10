using Plots
using CSV
using DataFrames
using LsqFit
using Statistics
using StatsBase
using ProgressBars
using LsqFit
using LaTeXStrings


include("functions.jl")

global irf               = get_irf()
global irf_bin_size      = get_irf_bin_size()
global tcspc_window_size = round(irf[end, 1] + irf[2, 1], sigdigits=4)
global fft_plan          = plan_fft(zeros(Float64, 256))
global ifft_plan         = plan_ifft(zeros(Float64, 256))
global folderpath        = get_folder_path("FLIM")
global serial            = nothing
global running           = Observable(false)
global pause             = Observable(false)
global active            = Observable(false)
global connected         = Observable(false)
global reset_I           = Observable(false)
global protocol_active   = Observable(false)

concentrations = Dict(
    "mCer3"=>Dict("0Cl"=>[], "2p5Cl"=>[], "5Cl"=>[], "10Cl"=>[], "20Cl"=>[], "40Cl"=>[], "80Cl"=>[], "146Cl"=>[], "146F"=>[]),
    "SCLM" =>Dict("0Cl"=>[], "2p5Cl"=>[], "5Cl"=>[], "10Cl"=>[], "20Cl"=>[], "40Cl"=>[], "80Cl"=>[], "146Cl"=>[], "146F"=>[])
)

main_path = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie2/Données/2025-10-17/Calibration"

folders = readdir(main_path)
filter!(f -> !startswith(f, "."), folders)
# pop!(folders, ".DS_Store")

for folder in folders
    folderpath = joinpath(main_path, folder)

    files = readdir(folderpath)
    filter!(f -> !startswith(f, "."), files)
    # pop!(files, ".DS_Store")

    cell = Float64[]

    full_vector = zeros(Float64, 256)

    rx = r"^(.*?)[_.]"
    protein = match(rx, folder)[1]

    for file in files
        if endswith(file, ".sdt")
            filepath = joinpath(folderpath, file)

            # println(filepath)

            vector, histogram_resolution, duration = open_SDT_file(filepath)

            rx = r"^(.*?)[_.]"
            concentration = match(rx, file)[1]

            push!(concentrations[protein][concentration], vec_to_lifetime(Float64.(vector))[1][1])

            # index = findfirst(isequal("_"), file)
            # println(file[1:index[1]-1])

            # full_vector .+= vector

            # params_raw, data = vec_to_lifetime(Float64.(vector))

            # println(params_raw[1])

            # println(index[1])

            # println(file[1:index[1]-1])

            # push!(cell, params_raw[1])
        end
    end

    # params_raw, data = vec_to_lifetime(Float64.(full_vector))

    # index = findfirst("m", folder)

    # push!(concentrations[string(folder[1:index[1]-1])], mean(filter(!isnan, cell)))
    # push!(concentrations[string(folder[1:index[1]-1])], params_raw[1])
end

# println(concentrations)

# plot(0,0; layout=(2, 1))

hill(x, p) = (p[1] * x.^p[3]) ./ (exp(p[2])^p[3] .+ x.^p[3]) .+ p[4]

plot(0,0)

# n = 1
scatters = []

for (i, (protein, conc_dict)) in enumerate(concentrations)
    println("Protein: $protein")

    x = []
    y = []

    for (i, (concentration, values)) in enumerate(conc_dict)
        # println("  Concentration: $concentration, Values: $values")
        # println(concentration)
        if concentration != "146F"
            append!(x, fill(.01 + parse(Float64, replace(replace(replace(concentration, "Cl"=>""), "p"=>"."), "146F"=>"300")), length(filter(!isnan, values))))
            append!(y, filter(!isnan, values))
            # println(x)
        else
            # println(concentration)

            append!(x, fill(10000, length(filter(!isnan, values))))
            append!(y, filter(!isnan, values) .+ 0.13)

            # println(mean(filter(!isnan, values)))
        end
    end

    xticks = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "", "10⁰", "10¹", "10²", "", "∞"])

    push!(scatters, scatter(x, y; label="$protein", markersize=5, color = i, legend_position = :topleft, xlabel="[Cl⁻] [mM]", ylabel="τ [ns]", xticks=xticks))

    initial_params = [0.5, 2.5, 1, 3]  # Vmax, Kd, n (choix de départ)
    
    fit = curve_fit(hill, x, y, initial_params; lower = [0.0, 0.0, 0.0, 0.0], upper = [Inf, Inf, Inf, Inf])
    params = fit.param  # les paramètres ajustés
    println(params)

    xfit = 10 .^ range(-2, stop=4, length=200)
    yfit = hill(xfit, params)

    # println(yfit[end])

    push!(scatters, plot!(xfit, yfit, label = "Fit Hill", lw=2, color = i))
end

# println

plot(scatters[1], scatters[3], layout=(2,1), xscale=:log10)
plot!(scatters[2], scatters[4], layout=(2,1), xscale=:log10)

main_path = "/Users/eliasouellet-oviedo/Documents/Stage2/Partie2/Données/2025-10-17/Limits"

folders = readdir(main_path)
filter!(f -> !startswith(f, "."), folders)
# pop!(folders, ".DS_Store")

concentrations = Dict(
    "mCer3"=>Dict(
        "0Cl_750_25"=>[], "0Cl_750_50"=>[], "0Cl_750_100"=>[],
        "0Cl_800_15"=>[], "0Cl_800_31"=>[], "0Cl_800_62"=>[],
        "0Cl_850_20"=>[], "0Cl_850_41"=>[], "0Cl_850_82"=>[],
        "146Cl_750_25"=>[], "146Cl_750_50"=>[], "146Cl_750_100"=>[],
        "146Cl_800_15"=>[], "146Cl_800_31"=>[], "146Cl_800_62"=>[],
        "146Cl_850_20"=>[], "146Cl_850_41"=>[], "146Cl_850_82"=>[]),
    "SCLM" =>Dict(
        "0Cl_750_25"=>[], "0Cl_750_50"=>[], "0Cl_750_100"=>[],
        "0Cl_800_15"=>[], "0Cl_800_31"=>[], "0Cl_800_62"=>[],
        "0Cl_850_20"=>[], "0Cl_850_41"=>[], "0Cl_850_82"=>[],
        "146Cl_750_25"=>[], "146Cl_750_50"=>[], "146Cl_750_100"=>[],
        "146Cl_800_15"=>[], "146Cl_800_31"=>[], "146Cl_800_62"=>[],
        "146Cl_850_20"=>[], "146Cl_850_41"=>[], "146Cl_850_82"=>[])
)

function extract_label(filename::AbstractString)
    # nom de fichier sans chemin ni extension
    name = splitext(basename(filename))[1]

    # regex strict (1: préfixe avant premier "_", 2: digits avant "nmexc", 3: digits après "em_" ou "-em_")
    pat = Regex("^([^_]+).*?(\\d+)nmexc.*?em[_-]?(\\d+)", "i")
    m = match(pat, name)
    if m === nothing
        error("Nom de fichier invalide — ne correspond pas au format attendu : $filename")
    end

    prefix, exc, em = m.captures
    return string(prefix, "_", exc, "_", em)
end

for folder in folders
    folderpath = joinpath(main_path, folder)

    files = readdir(folderpath)
    filter!(f -> !startswith(f, "."), files)
    # pop!(files, ".DS_Store")

    cell = Float64[]

    full_vector = zeros(Float64, 256)

    rx = r"^(.*?)[_.]"
    protein = match(rx, folder)[1]

    for file in files
        if endswith(file, ".sdt")
            filepath = joinpath(folderpath, file)

            # println(filepath)

            vector, histogram_resolution, duration = open_SDT_file(filepath)

            label = extract_label(file)

            push!(concentrations[protein][label], vec_to_lifetime(Float64.(vector))[1][1])

            # index = findfirst(isequal("_"), file)
            # println(file[1:index[1]-1])

            # full_vector .+= vector

            # params_raw, data = vec_to_lifetime(Float64.(vector))

            # println(params_raw[1])

            # println(index[1])

            # println(file[1:index[1]-1])

            # push!(cell, params_raw[1])
        end
    end

    # params_raw, data = vec_to_lifetime(Float64.(full_vector))

    # index = findfirst("m", folder)

    # push!(concentrations[string(folder[1:index[1]-1])], mean(filter(!isnan, cell)))
    # push!(concentrations[string(folder[1:index[1]-1])], params_raw[1])
end

# println(concentrations)

scatters = []

plot(0,0)

for (i, (protein, conc_dict)) in enumerate(concentrations)
    println("Protein: $protein")

    x = []
    y = []

    wavelength_dict = Dict(
        "750" => [[],[], [],[]],
        "800" => [[],[], [],[]],
        "850" => [[],[], [],[]]
    )

    power_dict = Dict(
        "25"  => [[],[], [],[]],
        "50"  => [[],[], [],[]],
        "100" => [[],[], [],[]]
    )

    for (i, (key, values)) in enumerate(conc_dict)
        concentration = split(key, "_")[1]
        wavelength = split(key, "_")[2]
        power = split(key, "_")[3]

        if power == "15" || power == "20"
            power = "25"
        elseif power == "31" || power == "41"
            power = "50"
        elseif power == "62" || power == "82"
            power = "100"
        end

        # println(split(concentration, "_"))
        append!(x, fill(.01 + parse(Float64, replace(concentration, "Cl"=>"")), length(filter(!isnan, values))))
        append!(y, filter(!isnan, values))

        if parse(Float64, replace(concentration, "Cl"=>"")) == 0
            append!(wavelength_dict[wavelength][1], fill(.01 + parse(Float64, power), length(filter(!isnan, values))))
            append!(power_dict[power][1], fill(.01 + parse(Float64, wavelength), length(filter(!isnan, values))))
            append!(wavelength_dict[wavelength][2], filter(!isnan, values))
            append!(power_dict[power][2], filter(!isnan, values))
        else
            append!(wavelength_dict[wavelength][3], fill(.01 + parse(Float64, power), length(filter(!isnan, values))))
            append!(power_dict[power][3], fill(.01 + parse(Float64, wavelength), length(filter(!isnan, values))))
            append!(wavelength_dict[wavelength][4], filter(!isnan, values))
            append!(power_dict[power][4], filter(!isnan, values))
        end
    end

    ylims = [(3.3, 4.1), (2.4, 3.2)]

    push!(scatters, scatter(wavelength_dict["750"][1], wavelength_dict["750"][2]; title="$protein 750 nm", label = "0 mM", markersize=5, color = 1, xlabel="Power [%]", ylabel="τ [ns]", xlims=(0, 110), ylims=ylims[i], legend = false))
    push!(scatters, scatter!(wavelength_dict["750"][3], wavelength_dict["750"][4]; title="$protein 750 nm", label = "146 mM", markersize=5, color = 2, xlabel="Power [%]", ylabel="τ [ns]", xlims=(0, 110), legend = false))
    push!(scatters, scatter(wavelength_dict["800"][1], wavelength_dict["800"][2]; title="$protein 800 nm", label = "0 mM", markersize=5, color = 1, xlabel="Power [%]", ylabel="τ [ns]", xlims=(0, 110), ylims=ylims[i], legend = false))
    push!(scatters, scatter!(wavelength_dict["800"][3], wavelength_dict["800"][4]; title="$protein 800 nm", label = "146 mM", markersize=5, color = 2, xlabel="Power [%]", ylabel="τ [ns]", xlims=(0, 110), legend = false))
    push!(scatters, scatter(wavelength_dict["850"][1], wavelength_dict["850"][2]; title="$protein 850 nm", label = "0 mM", markersize=5, color = 1, xlabel="Power [%]", ylabel="τ [ns]", xlims=(0, 110), ylims=ylims[i], legend = false))
    push!(scatters, scatter!(wavelength_dict["850"][3], wavelength_dict["850"][4]; title="$protein 850 nm", label = "146 mM", markersize=5, color = 2, xlabel="Power [%]", ylabel="τ [ns]", xlims=(0, 110), legend = false))

    push!(scatters, scatter(power_dict["25"][1], power_dict["25"][2]; title="$protein 12 mM", label = "0 mM", markersize=5, color = 1, xlabel="Wavelength [nm]", ylabel="τ [ns]", xlims=(700, 900), ylims=ylims[i], legend = false))
    push!(scatters, scatter!(power_dict["25"][3], power_dict["25"][4]; title="$protein 12 mM", label = "146 mM", markersize=5, color = 2, xlabel="Wavelength [nm]", ylabel="τ [ns]", xlims=(700, 900), legend = false))
    push!(scatters, scatter(power_dict["50"][1], power_dict["50"][2]; title="$protein 18 mM", label = "0 mM", markersize=5, color = 1, xlabel="Wavelength [nm]", ylabel="τ [ns]", xlims=(700, 900), ylims=ylims[i], legend = false))
    push!(scatters, scatter!(power_dict["50"][3], power_dict["50"][4]; title="$protein 18 mM", label = "146 mM", markersize=5, color = 2, xlabel="Wavelength [nm]", ylabel="τ [ns]", xlims=(700, 900), legend = false))
    push!(scatters, scatter(power_dict["100"][1], power_dict["100"][2]; title="$protein 27 mM", label = "0 mM", markersize=5, color = 1, xlabel="Wavelength [nm]", ylabel="τ [ns]", xlims=(700, 900), ylims=ylims[i], legend = false))
    push!(scatters, scatter!(power_dict["100"][3], power_dict["100"][4]; title="$protein 27 mM", label = "146 mM", markersize=5, color = 2, xlabel="Wavelength [nm]", ylabel="τ [ns]", xlims=(700, 900), legend = false))
end

plot(scatters[2:2:12]..., layout=(2,3))

plot!(size=(1200,800))