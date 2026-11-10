using Plots
using LsqFit
using FileIO
using Images
using LaTeXStrings
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

base_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-07-20_SCLM"
irf_path  = joinpath(base_path, "irf_20260710_16x_sclm.sdt")

irf_data = load_irf_from_sdt(irf_path)
RUNTIME[].irf           = irf_data
RUNTIME[].irf_bin_size   = compute_irf_bin_size(irf_data)
RUNTIME[].tcspc_window_size = round(irf_data[end, 1] + irf_data[2, 1], sigdigits=4)

# Compile la chaîne vec_to_lifetime -> Optim -> convolve une première fois :
# sans cet appel, le premier vrai fit ci-dessous a le temps du JIT (plusieurs
# secondes) qui dépasse le time_limit de l'optimiseur et renvoie NaN.
warmup_lifetime_fitting!()

# -----------------------------------------------------------------------------
# Protéines et concentrations en Cl- à analyser (3 réplicats par dossier)
# -----------------------------------------------------------------------------

proteins = ["SCLM", "M2-SCLM", "mCer3_Topaz"]

# "F" représente une concentration infinie en Cl- ; "0 mM" est placé à 0.01
# pour rester affichable sur un axe log, comme dans analyse_SCLM_1.jl.
concentration_folders = ["0 mM", "2.5 mM", "5 mM", "10 mM", "20 mM", "40 mM", "80 mM", "146 mM", "F"]
concentration_values  = [0.01, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 146.0, 100000.0]

hill(x, p) = (p[1] * x .^ p[3]) ./ (exp(p[2])^p[3] .+ x .^ p[3]) .+ p[4]
# p = [Δτ, ln(K_D), n, τ_0]

# vec_to_lifetime's default (first_fit=false) uses the tight ~20ms time_limit
# tuned for live acquisition display; under real system load a single fit can
# occasionally miss that budget and return NaN even after warmup_lifetime_fitting!.
# first_fit=true gives it the more generous ~60ms budget, and retrying on a NaN
# result (cheap here -- a few hundred fits total, one-off batch analysis) makes
# this robust to the rare case where even that isn't enough.
function fit_decay(counts::Vector{Float64}, histogram_resolution::Int, guess::Vector{Float64}; max_attempts::Int=5)
    for _ in 1:max_attempts
        params, data_xy = vec_to_lifetime(counts; guess=copy(guess), histogram_resolution=histogram_resolution, first_fit=true)
        isnan(params[1]) || return params, data_xy
    end
    error("Fit did not converge after $max_attempts attempts (guess length $(length(guess)))")
end

# Fits both a 1- and a 2-exponential decay model for one histogram (no PNG
# saved here -- see "Sauvegarde des histogrammes sommés" below for that).
# Returns all four "lifetime-like" values the final figure plots as separate
# columns: tau_1exp (the 1-exponential τ), tau_mean (the 2-exponential fit's
# intensity-weighted mean ⟨τ⟩ = (a₁τ₁² + (1-a₁)τ₂²)/(a₁τ₁ + (1-a₁)τ₂)), and
# the individual 2-exponential components tau_1/tau_2 (params_2exp = [τ₁, a₁,
# τ₂, shift, offset]; τ₁ is always the longer of the two, by
# mle_reconvolution_fit's own sorting convention). chi2_1/chi2_2 (see
# `fit_reduced_chi_squared` in lifetime_analysis.jl) let the per-protein
# model-selection step further down decide which model actually describes
# each protein's decay.
function lifetime_of_file(filepath::String)
    counts, _, histogram_resolution, _ = read_sdt_frame(filepath)
    x = Float64.(counts)

    params_1exp, data_xy_1exp = fit_decay(x, histogram_resolution, [3.0, 0.0, 1e-6])
    params_2exp, data_xy_2exp = fit_decay(x, histogram_resolution, [3.0, 0.5, 0.5, 0.0, 1e-6])

    chi2_1 = fit_reduced_chi_squared(params_1exp, data_xy_1exp)
    chi2_2 = fit_reduced_chi_squared(params_2exp, data_xy_2exp)

    tau1, a1, tau2 = params_2exp[1], params_2exp[2], params_2exp[3]
    tau_mean = (a1 * tau1^2 + (1 - a1) * tau2^2)/(a1 * tau1 + (1 - a1) * tau2)

    return (tau_1exp=params_1exp[1], tau_mean=tau_mean, tau_1=tau1, tau_2=tau2, a1=a1, chi2_1=chi2_1, chi2_2=chi2_2)
end

function lifetimes_of_folder(folderpath::String)
    files = sort(filter(f -> endswith(f, ".sdt"), readdir(folderpath)))
    return [lifetime_of_file(joinpath(folderpath, f)) for f in files]
end

# -----------------------------------------------------------------------------
# Sauvegarde des histogrammes sommés : pour chaque protéine x concentration
# (27 dossiers), les 3 réplicats sont sommés en un seul "gros" histogramme,
# ajusté avec les deux modèles (1 et 2 exponentielles) et sauvegardé -- 27 x 2
# = 54 fichiers au total.
# -----------------------------------------------------------------------------

histogram_plots_folder = joinpath(@__DIR__, "histogram_fits")
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
    ylims!(p, (minimum(filter(>(100), y_data_raw))*0.85, ylims(p)[2]))

    savefig(p, filepath)
end

function save_summed_histogram_fits(protein::String, conc_folder::String, folderpath::String)
    summed_counts, histogram_resolution = summed_counts_of_folder(folderpath)
    conc_label = sanitize_for_filename(conc_folder)

    params_1exp, data_xy_1exp = fit_decay(summed_counts, histogram_resolution, [3.0, 0.0, 1e-6])
    params_2exp, data_xy_2exp = fit_decay(summed_counts, histogram_resolution, [3.0, 0.5, 0.5, 0.0, 1e-6])

    chi2_1 = fit_reduced_chi_squared(params_1exp, data_xy_1exp)
    chi2_2 = fit_reduced_chi_squared(params_2exp, data_xy_2exp)

    save_decay_fit_plot(data_xy_1exp, params_1exp, chi2_1, "$protein $conc_folder -- 1-exponential",
        joinpath(histogram_plots_folder, "$(protein)_$(conc_label)_sum_1exp.png"))
    save_decay_fit_plot(data_xy_2exp, params_2exp, chi2_2, "$protein $conc_folder -- 2-exponential",
        joinpath(histogram_plots_folder, "$(protein)_$(conc_label)_sum_2exp.png"))
end

# -----------------------------------------------------------------------------
# Calcul des durées de vie pour chaque protéine / concentration
# -----------------------------------------------------------------------------

# Cinq représentations de "la" durée de vie / du fit à 2 exponentielles, une
# par colonne de la figure finale : τ (1-exponentielle), τₘ (moyenne pondérée
# par intensité du fit à 2 exponentielles), les deux composantes individuelles
# τ₁/τ₂ de ce même fit, et l'amplitude a₁ de la composante τ₁.
results_tau      = Dict{String, Dict{String, Vector{Float64}}}()
results_tau_mean = Dict{String, Dict{String, Vector{Float64}}}()
results_tau1     = Dict{String, Dict{String, Vector{Float64}}}()
results_tau2     = Dict{String, Dict{String, Vector{Float64}}}()
results_a1       = Dict{String, Dict{String, Vector{Float64}}}()
results_chi2_1   = Dict{String, Dict{String, Vector{Float64}}}()
results_chi2_2   = Dict{String, Dict{String, Vector{Float64}}}()
chi2_1_by_protein = Dict{String, Vector{Float64}}()
chi2_2_by_protein = Dict{String, Vector{Float64}}()

for protein in proteins
    results_tau[protein]      = Dict{String, Vector{Float64}}()
    results_tau_mean[protein] = Dict{String, Vector{Float64}}()
    results_tau1[protein]     = Dict{String, Vector{Float64}}()
    results_tau2[protein]     = Dict{String, Vector{Float64}}()
    results_a1[protein]       = Dict{String, Vector{Float64}}()
    results_chi2_1[protein]   = Dict{String, Vector{Float64}}()
    results_chi2_2[protein]   = Dict{String, Vector{Float64}}()
    chi2_1_by_protein[protein] = Float64[]
    chi2_2_by_protein[protein] = Float64[]

    for conc_folder in concentration_folders
        folderpath = joinpath(base_path, protein, conc_folder)
        file_results = lifetimes_of_folder(folderpath)
        save_summed_histogram_fits(protein, conc_folder, folderpath)

        results_tau[protein][conc_folder]      = [r.tau_1exp for r in file_results]
        results_tau_mean[protein][conc_folder] = [r.tau_mean for r in file_results]
        results_tau1[protein][conc_folder]     = [r.tau_1 for r in file_results]
        results_tau2[protein][conc_folder]     = [r.tau_2 for r in file_results]
        results_a1[protein][conc_folder]       = [r.a1 for r in file_results]
        results_chi2_1[protein][conc_folder]   = [r.chi2_1 for r in file_results]
        results_chi2_2[protein][conc_folder]   = [r.chi2_2 for r in file_results]
        append!(chi2_1_by_protein[protein], (r.chi2_1 for r in file_results))
        append!(chi2_2_by_protein[protein], (r.chi2_2 for r in file_results))

        println("$protein / $conc_folder : τ=", results_tau[protein][conc_folder], " τₘ=", results_tau_mean[protein][conc_folder])
    end
end

# -----------------------------------------------------------------------------
# Choix du modèle de décroissance (1 vs 2 exponentielles) par protéine, à
# partir du χ² moyen de chaque modèle sur tous ses histogrammes (toutes
# concentrations et réplicats confondus) -- le modèle au χ² moyen le plus
# proche de 1 (donc le plus bas, les deux moyennes étant > 1 en pratique) est
# retenu.
# -----------------------------------------------------------------------------

decay_model_line = Dict{String, String}()
decay_model_rows = Vector{String}[]

for protein in proteins
    mean_chi2_1 = mean(chi2_1_by_protein[protein])
    mean_chi2_2 = mean(chi2_2_by_protein[protein])
    n_exp = mean_chi2_2 < mean_chi2_1 ? 2 : 1
    best_chi2 = n_exp == 2 ? mean_chi2_2 : mean_chi2_1

    println("$protein : ⟨χ²⟩ 1-exp = $(round(mean_chi2_1, digits=2)), ⟨χ²⟩ 2-exp = $(round(mean_chi2_2, digits=2)) -> $n_exp-exponential")

    # decay_model_line[protein] = string(L"\mathrm{Decay:}\ %$(n_exp)\textrm{-exponential}\ (\langle\chi^2\rangle=%$(round(best_chi2, digits=2)))")
    decay_model_line[protein] = ""
    push!(decay_model_rows, [protein, string(mean_chi2_1), string(mean_chi2_2), string(n_exp)])
end

# -----------------------------------------------------------------------------
# Ratio FRET à partir des images TIFF (2 canaux, Bliq VMS/C1 et C2)
# -----------------------------------------------------------------------------

fret_base_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-07-20"

# Le nom du sous-dossier d'acquisition (ex. "850nm_Av200_65mW_...KIRIN1-016-001")
# est le même gabarit pour toutes les protéines/concentrations, mais on le
# retrouve dynamiquement plutôt que de le coder en dur.
function experiment_subfolder(protein::String, conc_folder::String)
    conc_path = joinpath(fret_base_path, protein, conc_folder)
    name = first(filter(f -> !startswith(f, "."), readdir(conc_path)))
    return joinpath(conc_path, name)
end

function channel_image_files(protein::String, conc_folder::String, channel::String)
    channel_path = joinpath(experiment_subfolder(protein, conc_folder), "Bliq VMS", channel)
    files = sort(filter(f -> !startswith(f, "."), readdir(channel_path)))
    return [joinpath(channel_path, f) for f in files]
end

total_counts(filepath::String) = sum(Float64.(load(filepath)))

# Ratio = C2/C1 (acceptor/donor) : ce sens croît avec [Cl⁻], comme la durée de
# vie ci-dessus, pour que les deux colonnes se lisent dans le même sens.
# Les fichiers de C1 et C2 sont triés puis appariés par position -- chaque
# paire (k-ième de C1, k-ième de C2) correspond à une des 3 acquisitions
# répétées, pas à un instant T identique entre les deux canaux.
function fret_ratios_of_folder(protein::String, conc_folder::String)
    c1_files = channel_image_files(protein, conc_folder, "C1")
    c2_files = channel_image_files(protein, conc_folder, "C2")
    length(c1_files) == length(c2_files) || error("Mismatched channel file counts for $protein/$conc_folder: $(length(c1_files)) vs $(length(c2_files))")
    return [total_counts(f2) / total_counts(f1) for (f1, f2) in zip(c1_files, c2_files)]
end

fret_results = Dict{String, Dict{String, Vector{Float64}}}()

for protein in proteins
    fret_results[protein] = Dict{String, Vector{Float64}}()
    for conc_folder in concentration_folders
        fret_results[protein][conc_folder] = fret_ratios_of_folder(protein, conc_folder)
        println("$protein / $conc_folder (FRET) : ", fret_results[protein][conc_folder])
    end
end

# -----------------------------------------------------------------------------
# Fit de Hill (partagé entre durée de vie et ratio FRET)
# -----------------------------------------------------------------------------

function gather_xy(data_dict::Dict{String, Dict{String, Vector{Float64}}}, protein::String)
    x_vals = Float64[]
    y_vals = Float64[]
    for (conc_value, conc_folder) in zip(concentration_values, concentration_folders)
        vals = data_dict[protein][conc_folder]
        append!(x_vals, fill(conc_value, length(vals)))
        append!(y_vals, vals)
    end
    return x_vals, y_vals
end

function fit_hill(x_vals::Vector{Float64}, y_vals::Vector{Float64})
    initial_params = [1.0, log(20.0), 1.0, minimum(y_vals)]
    fit = curve_fit(hill, x_vals, y_vals, initial_params; lower=[0.0, -16.0, 0.0, 0.0], upper=[Inf, 16.0, Inf, Inf])
    return fit.param
end

# `sub` distinguishes the derived lifetime columns (τₘ, τ₁, τ₂) from the plain
# τ/R columns without inventing a new "symbol" per caller: "" reproduces the
# original flat label ("τ", "τ₀", "τ_∞"), while e.g. sub="m" gives the compound
# label "τ_m", "τ_m₀", "τ_m_∞".
function hill_summary_lines(params::Vector{Float64}; base::String, sub::String="", unit::String)
    amplitude, ln_kd, hill_n, baseline = params
    label = isempty(sub) ? base : "$(base)_$(sub)"
    return [
        "Δ$(label) = $(round(amplitude, digits=3)) $(unit)",
        "Δ$(label) [%] = $(round(amplitude / (amplitude + baseline) * 100, digits=1)) %",
        "$(label)₀ = $(round(baseline, digits=3)) $(unit)",
        "$(label)_∞ = $(round(amplitude + baseline, digits=3)) $(unit)",
        "n = $(round(hill_n, digits=2))",
        "K_D = $(round(exp(ln_kd), digits=2)) mM",
    ]
end

# Version LaTeX (rendue par GR) du même résumé, pour le texte affiché dans les
# panneaux -- seule façon d'obtenir un vrai "∞" en indice (Unicode n'a pas de
# glyphe "∞ en indice"). `sub` (voir hill_summary_lines) devient un indice
# composé -- "m,0"/"m,∞" plutôt que d'empiler deux groupes d'indices distincts
# (non supporté proprement en LaTeX/mathtex). Les accolades autour de
# %$(symbol_display) évitent que GR fusionne p.ex. "\Delta" et "R" en un
# unique nom de macro invalide "\DeltaR".
function hill_box_latex(params::Vector{Float64}; base::String, sub::String="", unit::String)
    amplitude, ln_kd, hill_n, baseline = params
    tex_base = base == "τ" ? "\\tau" : base
    unit_tex = isempty(unit) ? "" : "\\ \\mathrm{$(unit)}"

    symbol_display = isempty(sub) ? tex_base : "$(tex_base)_{$(sub)}"
    zero_sub = isempty(sub) ? "0" : "$(sub),0"
    inf_sub  = isempty(sub) ? "\\infty" : "$(sub),\\infty"

    lines = [
        L"\Delta{%$(symbol_display)} = %$(round(amplitude, digits=3))%$(unit_tex)",
        L"\Delta{%$(symbol_display)}\ [\%] = %$(round(amplitude / (amplitude + baseline) * 100, digits=1))\%",
        L"%$(tex_base)_{%$(zero_sub)} = %$(round(baseline, digits=3))%$(unit_tex)",
        L"%$(tex_base)_{%$(inf_sub)} = %$(round(amplitude + baseline, digits=3))%$(unit_tex)",
        L"n = %$(round(hill_n, digits=2))",
        L"K_D = %$(round(exp(ln_kd), digits=2))\ \mathrm{mM}",
    ]
    return join(lines, "\n")
end

function compute_metric_stats(data_dict::Dict{String, Dict{String, Vector{Float64}}}; base::String, sub::String="", unit::String)
    stats = Dict{String, NamedTuple}()
    for protein in proteins
        x_vals, y_vals = gather_xy(data_dict, protein)
        params = fit_hill(x_vals, y_vals)

        println(protein)
        foreach(println, hill_summary_lines(params; base=base, sub=sub, unit=unit))

        stats[protein] = (x=x_vals, y=y_vals, yfit=hill(xfit, params), params=params)
    end
    return stats
end

# -----------------------------------------------------------------------------
# Graphique : panneau combiné (protéines superposées) + un panneau par protéine,
# pour chacune des deux mesures (durée de vie, ratio FRET)
# -----------------------------------------------------------------------------

marker_shapes = [:circle, :rect, :diamond]
xticks_log = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "∞"])
xfit = 10 .^ range(-2, stop=5, length=200)

#  Explicit margins here (rather than relying on Plots' auto-margin) since each
# panel is first built standalone and only later composited into the final 4x2
# grid -- without this, the y-axis label can end up clipped against the figure
# edge for panels that don't end up in the first column.
base_panel(xlabel, ylabel) = plot(xscale=:log10, framestyle=:box, grid=false,
    tickfontsize=11, guidefontsize=13, legendfontsize=10, titlefontsize=13,
    xlabel=xlabel, ylabel=ylabel, xticks=xticks_log, xlims=(0.005, 300000), palette=:Set1_3,
    left_margin=4Plots.mm, right_margin=3Plots.mm, top_margin=2Plots.mm, bottom_margin=2Plots.mm)

function make_combined_panel(stats::Dict{String, NamedTuple}; title::String, ylabel::String)
    panel = base_panel("", ylabel)
    plot!(panel, title=title, legend=:topleft)

    for (i, protein) in enumerate(proteins)
        d = stats[protein]
        scatter!(panel, d.x, d.y, label=protein, color=i, markershape=marker_shapes[i],
            markersize=5, markeralpha=0.5, markerstrokewidth=0.75, markerstrokecolor=:black)
        plot!(panel, xfit, d.yfit, label=nothing, lw=2, color=i)
    end

    # Un peu plus de marge que le padding automatique de Plots au-dessus/en
    # dessous des données, pour que les courbes ne touchent pas le cadre.
    all_y = reduce(vcat, (stats[p].y for p in proteins))
    ymin, ymax = extrema(all_y)
    yspan = ymax - ymin
    ylims!(panel, (ymin - 0.12 * yspan, ymax + 0.12 * yspan))

    return panel
end

function make_individual_panels(stats::Dict{String, NamedTuple}; ylabel::String, base::String, sub::String="", unit::String, extra_lines::Union{Nothing, Dict{String, String}}=nothing)
    panels = Plots.Plot[]

    for (i, protein) in enumerate(proteins)
        d = stats[protein]

        panel = base_panel(i == length(proteins) ? "[Cl⁻] [mM]" : "", ylabel)
        plot!(panel, title=protein, legend=false)

        scatter!(panel, d.x, d.y, color=i, markershape=marker_shapes[i],
            markersize=5, markeralpha=0.5, markerstrokewidth=0.75, markerstrokecolor=:black)
        plot!(panel, xfit, d.yfit, lw=2, color=i)

        # Cadrage vertical avec de la place au-dessus des données pour le texte
        # (un peu plus haut quand une ligne supplémentaire -- extra_lines -- s'ajoute).
        ymin, ymax = extrema(d.y)
        yspan = ymax - ymin
        top_pad = extra_lines === nothing ? 0.15 : 0.24
        ylims!(panel, (ymin - 0.15 * yspan, ymax + top_pad * yspan))

        box_text = hill_box_latex(d.params; base=base, sub=sub, unit=unit)
        if extra_lines !== nothing
            box_text = box_text * "\n" * extra_lines[protein]
        end

        xl, yl = xlims(panel), ylims(panel)
        text_x = xl[1] * (xl[2] / xl[1])^0.05
        text_y = yl[1] + 0.90 * (yl[2] - yl[1])
        annotate!(panel, text_x, text_y, text(box_text, 8, :left, :top))

        push!(panels, panel)
    end

    return panels
end

tau_stats      = compute_metric_stats(results_tau; base="τ", unit="ns")
tau_mean_stats = compute_metric_stats(results_tau_mean; base="τ", sub="m", unit="ns")
tau1_stats     = compute_metric_stats(results_tau1; base="τ", sub="1", unit="ns")
tau2_stats     = compute_metric_stats(results_tau2; base="τ", sub="2", unit="ns")
a1_stats       = compute_metric_stats(results_a1; base="a", sub="1", unit="")
fret_stats     = compute_metric_stats(fret_results; base="R", unit="")

tau_combined      = make_combined_panel(tau_stats; title="τ (1-exponential)", ylabel="τ [ns]")
tau_mean_combined = make_combined_panel(tau_mean_stats; title="τᵢ (2-exp mean)", ylabel="τᵢ [ns]")
tau1_combined     = make_combined_panel(tau1_stats; title="τ₁ (2-exponential)", ylabel="τ₁ [ns]")
tau2_combined     = make_combined_panel(tau2_stats; title="τ₂ (2-exponential)", ylabel="τ₂ [ns]")
a1_combined       = make_combined_panel(a1_stats; title="a₁ (2-exponential)", ylabel="a₁")
fret_combined     = make_combined_panel(fret_stats; title="FRET ratio", ylabel="FRET ratio (C2/C1)")

tau_panels      = make_individual_panels(tau_stats; ylabel="τ [ns]", base="τ", unit="ns", extra_lines=decay_model_line)
tau_mean_panels = make_individual_panels(tau_mean_stats; ylabel="τᵢ [ns]", base="τ", sub="i", unit="ns", extra_lines=decay_model_line)
tau1_panels     = make_individual_panels(tau1_stats; ylabel="τ₁ [ns]", base="τ", sub="1", unit="ns", extra_lines=decay_model_line)
tau2_panels     = make_individual_panels(tau2_stats; ylabel="τ₂ [ns]", base="τ", sub="2", unit="ns", extra_lines=decay_model_line)
a1_panels       = make_individual_panels(a1_stats; ylabel="a₁", base="a", sub="1", unit="", extra_lines=decay_model_line)
fret_panels     = make_individual_panels(fret_stats; ylabel="FRET ratio (C2/C1)", base="R", unit="")

# -----------------------------------------------------------------------------
# Export CSV : données brutes par réplicat, résumé des fits de Hill, et
# sélection du modèle de décroissance -- à l'exclusion des données brutes des
# histogrammes et des fits sur les histogrammes sommés des 3 réplicats (qui
# restent uniquement dans les PNG de histogram_fits/).
# -----------------------------------------------------------------------------

csv_folder = joinpath(@__DIR__, "csv_exports")
mkpath(csv_folder)

function write_csv(filepath::String, header::Vector{String}, rows::Vector{Vector{String}})
    open(filepath, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
end

# raw_data.csv : une ligne par protéine x concentration x réplicat. Le
# réplicat FRET associé est apparié par position (1er réplicat .sdt <-> 1er
# dossier TIFF trié, etc.), pas par identité d'acquisition -- même convention
# que l'appariement C1/C2 plus haut.
raw_data_rows = Vector{String}[]
for protein in proteins
    for (conc_value, conc_folder) in zip(concentration_values, concentration_folders)
        n_rep = length(results_tau[protein][conc_folder])
        for i in 1:n_rep
            push!(raw_data_rows, [
                protein, string(conc_value), string(i),
                string(results_tau[protein][conc_folder][i]),
                string(results_tau_mean[protein][conc_folder][i]),
                string(results_tau1[protein][conc_folder][i]),
                string(results_tau2[protein][conc_folder][i]),
                string(results_a1[protein][conc_folder][i]),
                string(results_chi2_1[protein][conc_folder][i]),
                string(results_chi2_2[protein][conc_folder][i]),
                string(fret_results[protein][conc_folder][i]),
            ])
        end
    end
end

write_csv(joinpath(csv_folder, "raw_data.csv"),
    ["protein", "concentration_mM", "replicate", "tau_1exp", "tau_mean", "tau_1", "tau_2", "a1", "chi2_1", "chi2_2", "fret_ratio"],
    raw_data_rows)

# hill_fit_summary.csv : une ligne par protéine x métrique
function hill_summary_row(protein::String, metric_name::String, params::Vector{Float64})
    amplitude, ln_kd, hill_n, baseline = params
    delta_pct = amplitude / (amplitude + baseline) * 100
    return [protein, metric_name, string(amplitude), string(delta_pct), string(baseline), string(amplitude + baseline), string(hill_n), string(exp(ln_kd))]
end

metric_stats_pairs = [
    ("tau", tau_stats), ("tau_mean", tau_mean_stats), ("tau_1", tau1_stats),
    ("tau_2", tau2_stats), ("a1", a1_stats), ("fret", fret_stats),
]
hill_summary_rows = Vector{String}[]
for (metric_name, stats) in metric_stats_pairs
    for protein in proteins
        push!(hill_summary_rows, hill_summary_row(protein, metric_name, stats[protein].params))
    end
end

write_csv(joinpath(csv_folder, "hill_fit_summary.csv"),
    ["protein", "metric", "delta", "delta_pct", "baseline", "asymptote", "n", "kd_mM"],
    hill_summary_rows)

# decay_model_selection.csv : une ligne par protéine
write_csv(joinpath(csv_folder, "decay_model_selection.csv"),
    ["protein", "mean_chi2_1exp", "mean_chi2_2exp", "selected_model"],
    decay_model_rows)

# -----------------------------------------------------------------------------
# Figure finale : 4 rangées (combiné + 3 protéines) x 6 colonnes
# (τ 1-exponentielle, τₘ, τ₁, τ₂, a₁, FRET)
# -----------------------------------------------------------------------------

plot(
    tau_combined, tau_mean_combined, tau1_combined, tau2_combined, a1_combined, fret_combined,
    tau_panels[1], tau_mean_panels[1], tau1_panels[1], tau2_panels[1], a1_panels[1], fret_panels[1],
    tau_panels[2], tau_mean_panels[2], tau1_panels[2], tau2_panels[2], a1_panels[2], fret_panels[2],
    tau_panels[3], tau_mean_panels[3], tau1_panels[3], tau2_panels[3], a1_panels[3], fret_panels[3],
    layout=(4, 6), size=(3360, 1700), dpi=100,
)
