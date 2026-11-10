using Plots
using LsqFit
using Statistics

include("SdtFile.jl")
using .SdtFile
include("lifetime_analysis.jl")

# lifetime_analysis.jl is `include`d directly (not loaded as a package), so
# its __init__() is never invoked automatically by Julia -- call it here to
# populate RUNTIME with default FFT plans before we override the IRF fields.
# warmup_lifetime_fitting! (below) also references this constant.
const DEFAULT_HISTOGRAM_RESOLUTION = 256
__init__()

# -----------------------------------------------------------------------------
# Dossier de données et IRF
# -----------------------------------------------------------------------------

base_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-07-22_GECO"

# Aucun fichier IRF n'est fourni avec ce nouveau jeu de données -- l'IRF
# caractérise l'électronique TCSPC du détecteur, pas l'échantillon, et cette
# expérience a été acquise deux jours seulement après celle de SCLM sur le
# même instrument ; on réutilise donc son IRF.
irf_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-07-20_SCLM/irf_20260710_16x_sclm.sdt"

irf_data = load_irf_from_sdt(irf_path)
RUNTIME[].irf           = irf_data
RUNTIME[].irf_bin_size   = compute_irf_bin_size(irf_data)
RUNTIME[].tcspc_window_size = round(irf_data[end, 1] + irf_data[2, 1], sigdigits=4)

# Compile la chaîne vec_to_lifetime -> Optim -> convolve une première fois :
# sans cet appel, le premier vrai fit ci-dessous a le temps du JIT (plusieurs
# secondes) qui dépasse le time_limit de l'optimiseur et renvoie NaN.
warmup_lifetime_fitting!()

# -----------------------------------------------------------------------------
# Les dossiers sont nommés d'après la concentration de CaEGTA (0-10 mM), mais
# c'est le [Ca²⁺] libre résultant (μM, table de calibration CaEGTA fournie)
# qui a un sens biologique pour FR-GECO (un indicateur calcique) et qui
# alimente le fit de Hill -- gamme non-linéaire (0.017 à 39 μM, presque 4
# ordres de grandeur), d'où un axe log comme l'expérience SCLM plutôt qu'un
# axe linéaire. "0 mM" (0 μM libre, un vrai zéro) doit être placé *très* près
# de 0 sur l'axe log, pas juste "un peu moins que 0.017" : à 0.005 (à peine 3x
# sous le premier point réel), le fit de Hill confondait carrément la montée
# rapide entre "0 mM" et "1 mM" avec la position de K_D lui-même (K_D et τ_0
# sortaient tous les deux faux, ex. τ_0 ≈ 0.02 ns alors que les données à
# 0 mM donnent ≈ 0.43 ns). Vérifié empiriquement : le fit se stabilise (même
# résultat à la précision numérique) pour tout placeholder ≤ 1e-4, donc 1e-5
# ci-dessous a une bonne marge de sécurité.
# Une seule protéine ici (FR-GECO) -- la variable qui distingue les séries
# dans les graphiques est la longueur d'onde, pas la protéine.
# -----------------------------------------------------------------------------

flim_wavelengths = ["750", "800"]        # "1064" existe comme dossier mais est vide (aucun .sdt)

concentration_folders  = ["0 mM", "1 mM", "2 mM", "3 mM", "4 mM", "5 mM", "6 mM", "7 mM", "8 mM", "9 mM", "10 mM"]
free_ca_by_caegta_conc = Dict("0 mM" => 0.0, "1 mM" => 0.017, "2 mM" => 0.038, "3 mM" => 0.065, "4 mM" => 0.100,
    "5 mM" => 0.150, "6 mM" => 0.225, "7 mM" => 0.351, "8 mM" => 0.602, "9 mM" => 1.35, "10 mM" => 39.0)
concentration_values = [c == 0.0 ? 1e-5 : c for c in getindex.(Ref(free_ca_by_caegta_conc), concentration_folders)]

hill(x, p) = (p[1] * x .^ p[3]) ./ (exp(p[2])^p[3] .+ x .^ p[3]) .+ p[4]
# p = [Δτ, ln(K_D), n, τ_0]

# vec_to_lifetime's default (first_fit=false) uses the tight ~20ms time_limit
# tuned for live acquisition display; under real system load a single fit can
# occasionally miss that budget and return NaN even after warmup_lifetime_fitting!.
# first_fit=true gives it the more generous ~60ms budget, and retrying on a NaN
# result (cheap here -- a few dozen fits total, one-off batch analysis) makes
# this robust to the rare case where even that isn't enough.
function fit_decay(counts::Vector{Float64}, histogram_resolution::Int, guess::Vector{Float64}; max_attempts::Int=10)
    for _ in 1:max_attempts
        params, data_xy = vec_to_lifetime(counts; guess=copy(guess), histogram_resolution=histogram_resolution, first_fit=true)
        isnan(params[1]) || return params, data_xy
    end
    error("Fit did not converge after $max_attempts attempts (guess length $(length(guess)))")
end

# -----------------------------------------------------------------------------
# Sauvegarde des histogrammes sommés : pour chaque longueur d'onde x
# concentration (22 dossiers), les 3 réplicats sont sommés en un seul "gros"
# histogramme, ajusté avec les deux modèles (1 et 2 exponentielles) et
# sauvegardé -- 22 x 2 = 44 fichiers au total (même principe que l'expérience
# SCLM, où c'était 27 x 2 = 54).
# -----------------------------------------------------------------------------

histogram_plots_folder = joinpath(@__DIR__, "histogram_fits_GECO")
mkpath(histogram_plots_folder)

sanitize_for_filename(s::String) = replace(s, " " => "")

function summed_counts_of_folder(folderpath::String)
    files = sort(filter(f -> endswith(f, ".sdt"), readdir(folderpath)))
    frames = [read_sdt_frame(joinpath(folderpath, f)) for f in files]
    histogram_resolution = frames[1][3]
    summed_counts = sum(Float64.(frame[1]) for frame in frames)
    return summed_counts, histogram_resolution
end

# data_xy[2] (les comptes bruts) contient de vrais zéros avant l'arrivée de
# l'impulsion ; on plafonne à 0.5 uniquement pour l'affichage en échelle log
# (reconstruct_model_counts plafonne déjà le modèle, mais à 1e-256 -- bien
# trop petit pour l'axe -- d'où le second plafonnage ici).
function save_decay_fit_plot(data_xy::Vector{Vector{Float64}}, params::Vector{Float64}, chi2::Float64, title_str::String, filepath::String)
    x_data, y_data_raw = data_xy
    model_counts, _, _ = reconstruct_model_counts(params, data_xy)

    p = plot(x_data, max.(y_data_raw, 0.5), seriestype=:scatter, markersize=2.5, markeralpha=0.4,
        markerstrokewidth=0, color=:black, label="data", yscale=:log10, framestyle=:box, grid=false,
        legend=:topright, xlabel="Time [ns]", ylabel="Counts", title=title_str, titlefontsize=10,
        size=(600, 450), dpi=150)
    plot!(p, x_data, max.(model_counts, 0.5), lw=1.5, color=:red, label="fit (χ² = $(round(chi2, digits=2)))")
    ylims!(p, (minimum(filter(>(100), y_data_raw)), ylims(p)[2]))

    savefig(p, filepath)
end

function save_summed_histogram_fits(wavelength::String, conc_folder::String, folderpath::String)
    summed_counts, histogram_resolution = summed_counts_of_folder(folderpath)
    conc_label = sanitize_for_filename(conc_folder)

    params_1exp, data_xy_1exp = fit_decay(summed_counts, histogram_resolution, [3.0, 0.0, 1e-6])
    params_2exp, data_xy_2exp = fit_decay(summed_counts, histogram_resolution, [3.0, 0.5, 0.5, 0.0, 1e-6])

    chi2_1 = fit_reduced_chi_squared(params_1exp, data_xy_1exp)
    chi2_2 = fit_reduced_chi_squared(params_2exp, data_xy_2exp)

    save_decay_fit_plot(data_xy_1exp, params_1exp, chi2_1, "$(wavelength) nm $conc_folder (sum of 3) -- 1-exponential",
        joinpath(histogram_plots_folder, "$(wavelength)nm_$(conc_label)_sum_1exp.png"))
    save_decay_fit_plot(data_xy_2exp, params_2exp, chi2_2, "$(wavelength) nm $conc_folder (sum of 3) -- 2-exponential",
        joinpath(histogram_plots_folder, "$(wavelength)nm_$(conc_label)_sum_2exp.png"))
end

# Fits both a 1- and a 2-exponential decay model for one histogram (no PNG
# saved here -- see save_summed_histogram_fits above for that). Returns both
# tau (the 1-exponential lifetime) and tau_mean (the 2-exponential fit's
# intensity-weighted mean ⟨τ⟩ = (a₁τ₁² + (1-a₁)τ₂²)/(a₁τ₁ + (1-a₁)τ₂),
# params_2exp = [τ₁, a₁, τ₂, shift, offset]) -- the two quantities plotted
# side by side below. chi2_1/chi2_2 (see `fit_reduced_chi_squared` in
# lifetime_analysis.jl) let the per-wavelength model-selection step further
# down report which model actually describes the decay at that wavelength.
function lifetime_of_file(filepath::String)
    counts, _, histogram_resolution, _ = read_sdt_frame(filepath)
    x = Float64.(counts)

    params_1exp, data_xy_1exp = fit_decay(x, histogram_resolution, [3.0, 0.0, 1e-6])
    params_2exp, data_xy_2exp = fit_decay(x, histogram_resolution, [3.0, 0.5, 0.5, 0.0, 1e-6])

    chi2_1 = fit_reduced_chi_squared(params_1exp, data_xy_1exp)
    chi2_2 = fit_reduced_chi_squared(params_2exp, data_xy_2exp)

    tau1, a1, tau2 = params_2exp[1], params_2exp[2], params_2exp[3]
    tau_mean = (a1 * tau1^2 + (1 - a1) * tau2^2) / (a1 * tau1 + (1 - a1) * tau2)

    return (tau=params_1exp[1], tau_mean=tau_mean, chi2_1=chi2_1, chi2_2=chi2_2)
end

function lifetimes_of_folder(folderpath::String)
    files = sort(filter(f -> endswith(f, ".sdt"), readdir(folderpath)))
    return [lifetime_of_file(joinpath(folderpath, f)) for f in files]
end

# -----------------------------------------------------------------------------
# Calcul des durées de vie pour chaque longueur d'onde / concentration
# -----------------------------------------------------------------------------

results          = Dict{String, Dict{String, Vector{Float64}}}()
results_tau_mean = Dict{String, Dict{String, Vector{Float64}}}()
chi2_1_by_wavelength = Dict{String, Vector{Float64}}()
chi2_2_by_wavelength = Dict{String, Vector{Float64}}()

for wavelength in flim_wavelengths
    results[wavelength]          = Dict{String, Vector{Float64}}()
    results_tau_mean[wavelength] = Dict{String, Vector{Float64}}()
    chi2_1_by_wavelength[wavelength] = Float64[]
    chi2_2_by_wavelength[wavelength] = Float64[]

    for conc_folder in concentration_folders
        folderpath = joinpath(base_path, conc_folder, wavelength)
        file_results = lifetimes_of_folder(folderpath)
        save_summed_histogram_fits(wavelength, conc_folder, folderpath)

        results[wavelength][conc_folder]          = [r.tau for r in file_results]
        results_tau_mean[wavelength][conc_folder] = [r.tau_mean for r in file_results]
        append!(chi2_1_by_wavelength[wavelength], (r.chi2_1 for r in file_results))
        append!(chi2_2_by_wavelength[wavelength], (r.chi2_2 for r in file_results))

        println("$(wavelength)nm / $conc_folder : τ=", results[wavelength][conc_folder], " τᵢ=", results_tau_mean[wavelength][conc_folder])
    end
end

# -----------------------------------------------------------------------------
# Choix du modèle de décroissance (1 vs 2 exponentielles) par longueur d'onde,
# à partir du χ² moyen de chaque modèle sur tous ses histogrammes (toutes
# concentrations et réplicats confondus) -- le modèle au χ² moyen le plus bas
# est retenu. Purement diagnostique (imprimé en console) ; n'affecte pas la
# durée de vie (1-exponentielle) qui alimente le graphique ci-dessous.
# -----------------------------------------------------------------------------

for wavelength in flim_wavelengths
    mean_chi2_1 = mean(chi2_1_by_wavelength[wavelength])
    mean_chi2_2 = mean(chi2_2_by_wavelength[wavelength])
    n_exp = mean_chi2_2 < mean_chi2_1 ? 2 : 1

    println("$(wavelength)nm : ⟨χ²⟩ 1-exp = $(round(mean_chi2_1, digits=2)), ⟨χ²⟩ 2-exp = $(round(mean_chi2_2, digits=2)) -> $n_exp-exponential")
end

# -----------------------------------------------------------------------------
# Fit de Hill (partagé entre τ et τᵢ). `series_labels` remplace la liste fixe
# de protéines de l'expérience SCLM -- ici les deux longueurs d'onde FLIM.
# -----------------------------------------------------------------------------

function gather_xy(data_dict::Dict{String, Dict{String, Vector{Float64}}}, series_label::String)
    x_vals = Float64[]
    y_vals = Float64[]
    for (conc_value, conc_folder) in zip(concentration_values, concentration_folders)
        vals = data_dict[series_label][conc_folder]
        append!(x_vals, fill(conc_value, length(vals)))
        append!(y_vals, vals)
    end
    return x_vals, y_vals
end

function fit_hill(x_vals::Vector{Float64}, y_vals::Vector{Float64})
    # ln(1.0) μM comme point de départ pour K_D (plutôt que ln(20.0) mM dans
    # l'expérience SCLM) : la plage de [Ca²⁺] libre testée ici est 0.017-39 μM.
    initial_params = [1.0, log(1.0), 1.0, minimum(y_vals)]
    fit = curve_fit(hill, x_vals, y_vals, initial_params; lower=[0.0, -16.0, 0.0, 0.0], upper=[Inf, 16.0, Inf, Inf])
    return fit.param
end

# `sub` distingue τᵢ de τ sans inventer un nouveau "symbol" par appelant : ""
# reproduit l'étiquette plate d'origine ("τ", "τ₀", "τ_∞"), alors que sub="i"
# donne l'étiquette composée "τ_i", "τ_i₀", "τ_i_∞".
function hill_summary_lines(params::Vector{Float64}; base::String, sub::String="", unit::String)
    amplitude, ln_kd, hill_n, baseline = params
    label = isempty(sub) ? base : "$(base)_$(sub)"
    return [
        "Δ$(label) = $(round(amplitude, digits=3)) $(unit)",
        "Δ$(label) [%] = $(round(amplitude / (amplitude + baseline) * 100, digits=1)) %",
        "$(label)₀ = $(round(baseline, digits=3)) $(unit)",
        "$(label)_∞ = $(round(amplitude + baseline, digits=3)) $(unit)",
        "n = $(round(hill_n, digits=2))",
        "K_D = $(round(exp(ln_kd), digits=2)) μM",
    ]
end

function compute_metric_stats(data_dict::Dict{String, Dict{String, Vector{Float64}}}, series_labels::Vector{String}; base::String, sub::String="", unit::String)
    stats = Dict{String, NamedTuple}()
    for series_label in series_labels
        x_vals, y_vals = gather_xy(data_dict, series_label)
        params = fit_hill(x_vals, y_vals)

        println(series_label)
        foreach(println, hill_summary_lines(params; base=base, sub=sub, unit=unit))

        stats[series_label] = (x=x_vals, y=y_vals, yfit=hill(xfit, params), params=params)
    end
    return stats
end

# -----------------------------------------------------------------------------
# Graphique : un panneau par mesure (durée de vie, ratio d'intensité), les
# longueurs d'onde superposées dans chacun. Axe log (comme l'expérience SCLM) :
# le [Ca²⁺] libre testé couvre 0.017-39 μM, presque 4 ordres de grandeur.
# -----------------------------------------------------------------------------

marker_shapes = [:circle, :rect, :diamond]
xticks_log = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ["0", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹", "10²"])
xfit = 10 .^ range(-5.3, stop=2, length=200)

base_panel(xlabel, ylabel) = plot(xscale=:log10, framestyle=:box, grid=false,
    tickfontsize=11, guidefontsize=13, legendfontsize=10, titlefontsize=13,
    xlabel=xlabel, ylabel=ylabel, xticks=xticks_log, xlims=(6e-6, 150), palette=:Set1_3,
    left_margin=4Plots.mm, right_margin=3Plots.mm, top_margin=2Plots.mm, bottom_margin=2Plots.mm)

function make_combined_panel(stats::Dict{String, NamedTuple}, series_labels::Vector{String}; title::String, ylabel::String)
    panel = base_panel("[Ca²⁺] [μM]", ylabel)
    plot!(panel, title=title, legend=:topleft)

    for (i, series_label) in enumerate(series_labels)
        d = stats[series_label]
        scatter!(panel, d.x, d.y, label="$(series_label) nm", color=i, markershape=marker_shapes[i],
            markersize=5, markeralpha=0.5, markerstrokewidth=0.75, markerstrokecolor=:black)
        plot!(panel, xfit, d.yfit, label=nothing, lw=2, color=i)
    end

    # Un peu plus de marge que le padding automatique de Plots au-dessus/en
    # dessous des données, pour que les courbes ne touchent pas le cadre.
    all_y = reduce(vcat, (stats[s].y for s in series_labels))
    ymin, ymax = extrema(all_y)
    yspan = ymax - ymin
    ylims!(panel, (ymin - 0.12 * yspan, ymax + 0.12 * yspan))

    return panel
end

tau_stats      = compute_metric_stats(results, flim_wavelengths; base="τ", unit="ns")
tau_mean_stats = compute_metric_stats(results_tau_mean, flim_wavelengths; base="τ", sub="i", unit="ns")

tau_combined      = make_combined_panel(tau_stats, flim_wavelengths; title="τ (1-exponential) vs [Ca²⁺]", ylabel="τ [ns]")
tau_mean_combined = make_combined_panel(tau_mean_stats, flim_wavelengths; title="τᵢ (2-exp mean) vs [Ca²⁺]", ylabel="τᵢ [ns]")

# Mêmes limites d'axe y pour les deux panneaux (calculées sur l'union des
# données de τ et τᵢ, avec le même padding de 12% que make_combined_panel
# applique individuellement) -- pour rester directement comparables visuellement.
shared_y = vcat((tau_stats[w].y for w in flim_wavelengths)..., (tau_mean_stats[w].y for w in flim_wavelengths)...)
shared_ymin, shared_ymax = extrema(shared_y)
shared_yspan = shared_ymax - shared_ymin
shared_ylims = (shared_ymin - 0.12 * shared_yspan, shared_ymax + 0.12 * shared_yspan)
ylims!(tau_combined, shared_ylims)
ylims!(tau_mean_combined, shared_ylims)

# -----------------------------------------------------------------------------
# Tableaux récapitulatifs (mêmes valeurs que hill_summary_lines, une ligne par
# longueur d'onde), un pour τ et un pour τᵢ, affichés sous leur graphique
# respectif ("les descriptions sous les graphiques").
# -----------------------------------------------------------------------------

function table_row(series_label::String, params::Vector{Float64})
    amplitude, ln_kd, hill_n, baseline = params
    return [
        "$(series_label) nm",
        string(round(amplitude, digits=2)),
        string(round(baseline, digits=2)),
        string(round(amplitude + baseline, digits=2)),
        string(round(hill_n, digits=2)),
        string(round(exp(ln_kd), digits=3)),
    ]
end

# Pas de "τ₀"/"τ∞"/"τ_i₀"/"τ_i∞" ici (contrairement aux axes des graphiques,
# qui utilisent la police par défaut de GR) : la police "Courier" ci-dessous
# n'a pas ces glyphes indice Unicode ("GKS: glyph missing from current
# font") -- d'où "Baseline"/"Asymptote" en toutes lettres.
tau_table_headers      = ["Wavelength", "Δτ", "τ_0", "τ_∞", "n", "K_D"]
tau_mean_table_headers = ["Wavelength", "Δτ_i", "τ_i0", "τ_i∞", "n", "K_D"]

tau_table_rows      = [table_row(wavelength, tau_stats[wavelength].params) for wavelength in flim_wavelengths]
tau_mean_table_rows = [table_row(wavelength, tau_mean_stats[wavelength].params) for wavelength in flim_wavelengths]

# Un panneau sans axes contenant le texte formaté en colonnes alignées --
# Plots.jl n'a pas de type de série "tableau" ; "Courier" est la seule police
# à la fois monospace ET reconnue par GR sans avertissement (testé :
# "monospace" et "DejaVu Sans Mono" déclenchent tous deux un "GKS: could not
# find font"). annotate! n'est pas rogné aux limites de son panneau, donc un
# texte trop large déborde carrément dans le panneau voisin plutôt que d'être
# coupé -- 12pt convient ici (2-3 lignes de données seulement).
function render_table_panel(headers::Vector{String}, rows::Vector{Vector{String}})
    all_rows = vcat([headers], rows)
    ncols = length(headers)
    widths = [maximum(length(r[j]) for r in all_rows) for j in 1:ncols]
    fmt_row(r) = join([rpad(r[j], widths[j]) for j in 1:ncols], "  ")

    lines = [fmt_row(headers), fmt_row(["-"^widths[j] for j in 1:ncols])]
    append!(lines, [fmt_row(r) for r in rows])
    table_text = join(lines, "\n")

    panel = plot(axis=false, grid=false, framestyle=:none, xlims=(0, 1), ylims=(0, 1))
    annotate!(panel, 0.0, 1.0, text(table_text, 12, :left, :top, "Courier"))
    return panel
end

tau_table_panel      = render_table_panel(tau_table_headers, tau_table_rows)
tau_mean_table_panel = render_table_panel(tau_mean_table_headers, tau_mean_table_rows)

# -----------------------------------------------------------------------------
# Figure finale : τ et τᵢ côte à côte (1x2), une description sous chacun
# -----------------------------------------------------------------------------

layout = @layout [a b; c{0.2h} d{0.2h}]
plot(tau_combined, tau_mean_combined, tau_table_panel, tau_mean_table_panel, layout=layout, size=(1400, 900), dpi=200)
