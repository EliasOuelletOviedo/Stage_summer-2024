using Plots
using LsqFit
using FileIO
using Images
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
# proteins = ["mCer3_Topaz"]

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

# -----------------------------------------------------------------------------
# Sauvegarde d'un PNG (histogramme + ajustement) par fichier et par modèle
# -----------------------------------------------------------------------------

histogram_plots_folder = joinpath(@__DIR__, "histogram_fits")
mkpath(histogram_plots_folder)

sanitize_for_filename(s::String) = replace(s, " " => "")

# Les fichiers .sdt se terminent par "..._cN.sdt", N étant le réplicat (1-3) --
# on le réutilise tel quel plutôt que d'inventer un autre indexage.
function replicate_label(filename::String)
    m = match(r"_c(\d+)\.sdt$", filename)
    return m === nothing ? "c?" : "c$(m.captures[1])"
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

    savefig(p, filepath)
end

# Fits both a 1- and a 2-exponential decay model for one histogram, saving a
# PNG of each fit against the data. `tau` (the 1-exponential lifetime) drives
# the Hill-fit/plot below; chi2_1/chi2_2 (see `fit_reduced_chi_squared` in
# lifetime_analysis.jl) let the per-protein model-selection step further down
# report which model actually describes each protein's decay.
function lifetime_of_file(filepath::String)
    counts, _, histogram_resolution, _ = read_sdt_frame(filepath)
    x = Float64.(counts)

    path_parts = splitpath(filepath)
    protein, conc_folder, filename = path_parts[end-2], path_parts[end-1], path_parts[end]
    conc_label = sanitize_for_filename(conc_folder)
    rep_label = replicate_label(filename)

    params_1exp, data_xy_1exp = fit_decay(x, histogram_resolution, [3.0, 0.0, 1e-6])
    params_2exp, data_xy_2exp = fit_decay(x, histogram_resolution, [3.0, 0.5, 0.5, 0.0, 1e-6])

    chi2_1 = fit_reduced_chi_squared(params_1exp, data_xy_1exp)
    chi2_2 = fit_reduced_chi_squared(params_2exp, data_xy_2exp)

    save_decay_fit_plot(data_xy_1exp, params_1exp, chi2_1, "$protein $conc_folder $rep_label -- 1-exponential",
        joinpath(histogram_plots_folder, "$(protein)_$(conc_label)_$(rep_label)_1exp.png"))
    save_decay_fit_plot(data_xy_2exp, params_2exp, chi2_2, "$protein $conc_folder $rep_label -- 2-exponential",
        joinpath(histogram_plots_folder, "$(protein)_$(conc_label)_$(rep_label)_2exp.png"))

    return (tau=params_1exp[1], chi2_1=chi2_1, chi2_2=chi2_2)
end

function lifetimes_of_folder(folderpath::String)
    files = sort(filter(f -> endswith(f, ".sdt"), readdir(folderpath)))
    return [lifetime_of_file(joinpath(folderpath, f)) for f in files]
end

# -----------------------------------------------------------------------------
# Calcul des durées de vie pour chaque protéine / concentration
# -----------------------------------------------------------------------------

results = Dict{String, Dict{String, Vector{Float64}}}()
chi2_1_by_protein = Dict{String, Vector{Float64}}()
chi2_2_by_protein = Dict{String, Vector{Float64}}()

for protein in proteins
    results[protein] = Dict{String, Vector{Float64}}()
    chi2_1_by_protein[protein] = Float64[]
    chi2_2_by_protein[protein] = Float64[]

    for conc_folder in concentration_folders
        folderpath = joinpath(base_path, protein, conc_folder)
        file_results = lifetimes_of_folder(folderpath)

        results[protein][conc_folder] = [r.tau for r in file_results]
        append!(chi2_1_by_protein[protein], (r.chi2_1 for r in file_results))
        append!(chi2_2_by_protein[protein], (r.chi2_2 for r in file_results))

        println("$protein / $conc_folder : ", results[protein][conc_folder])
    end
end

# -----------------------------------------------------------------------------
# Choix du modèle de décroissance (1 vs 2 exponentielles) par protéine, à
# partir du χ² moyen de chaque modèle sur tous ses histogrammes (toutes
# concentrations et réplicats confondus) -- le modèle au χ² moyen le plus bas
# est retenu. Purement diagnostique (imprimé en console) ; n'affecte pas la
# durée de vie (1-exponentielle) qui alimente le graphique ci-dessous.
# -----------------------------------------------------------------------------

for protein in proteins
    mean_chi2_1 = mean(chi2_1_by_protein[protein])
    mean_chi2_2 = mean(chi2_2_by_protein[protein])
    n_exp = mean_chi2_2 < mean_chi2_1 ? 2 : 1

    println("$protein : ⟨χ²⟩ 1-exp = $(round(mean_chi2_1, digits=2)), ⟨χ²⟩ 2-exp = $(round(mean_chi2_2, digits=2)) -> $n_exp-exponential")
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

function hill_summary_lines(params::Vector{Float64}; symbol::String, unit::String)
    amplitude, ln_kd, hill_n, baseline = params
    return [
        "Δ$(symbol) = $(round(amplitude, digits=3)) $(unit)",
        "Δ$(symbol) [%] = $(round(amplitude / (amplitude + baseline) * 100, digits=1)) %",
        "$(symbol)₀ = $(round(baseline, digits=3)) $(unit)",
        "$(symbol)_∞ = $(round(amplitude + baseline, digits=3)) $(unit)",
        "n = $(round(hill_n, digits=2))",
        "K_D = $(round(exp(ln_kd), digits=2)) mM",
    ]
end

function compute_metric_stats(data_dict::Dict{String, Dict{String, Vector{Float64}}}; symbol::String, unit::String)
    stats = Dict{String, NamedTuple}()
    for protein in proteins
        x_vals, y_vals = gather_xy(data_dict, protein)
        params = fit_hill(x_vals, y_vals)

        println(protein)
        foreach(println, hill_summary_lines(params; symbol=symbol, unit=unit))

        stats[protein] = (x=x_vals, y=y_vals, yfit=hill(xfit, params), params=params)
    end
    return stats
end

# -----------------------------------------------------------------------------
# Graphique : un panneau par mesure (durée de vie, ratio FRET), les 3 protéines
# superposées dans chacun
# -----------------------------------------------------------------------------

marker_shapes = [:circle, :rect, :diamond]
xticks_log = ([0.01, 0.1, 1, 10, 100, 1000, 10000, 100000], ["0", "10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴", "∞"])
xfit = 10 .^ range(-2, stop=5, length=200)

base_panel(xlabel, ylabel) = plot(xscale=:log10, framestyle=:box, grid=false,
    tickfontsize=11, guidefontsize=13, legendfontsize=10, titlefontsize=13,
    xlabel=xlabel, ylabel=ylabel, xticks=xticks_log, xlims=(0.005, 300000), palette=:Set1_3,
    left_margin=7Plots.mm, right_margin=3Plots.mm, top_margin=2Plots.mm, bottom_margin=2Plots.mm)

function make_combined_panel(stats::Dict{String, NamedTuple}; title::String, ylabel::String)
    panel = base_panel("[Cl⁻] [mM]", ylabel)
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

lifetime_stats = compute_metric_stats(results; symbol="τ", unit="ns")
fret_stats     = compute_metric_stats(fret_results; symbol="R", unit="")

lifetime_combined = make_combined_panel(lifetime_stats; title="Lifetime variation of SCLM variants vs [Cl⁻]", ylabel="Lifetime [ns]")
fret_combined     = make_combined_panel(fret_stats; title="FRET ratio variation of SCLM variants vs [Cl⁻]", ylabel="FRET ratio (C2/C1)")

# -----------------------------------------------------------------------------
# Tableaux récapitulatifs (mêmes valeurs que hill_summary_lines, une ligne par
# protéine), un pour la durée de vie et un pour le ratio FRET, affichés sous
# leur graphique respectif
# -----------------------------------------------------------------------------

function table_row(protein::String, params::Vector{Float64})
    amplitude, ln_kd, hill_n, baseline = params
    return [
        protein,
        string(round(amplitude, digits=2)),
        # string(round(amplitude / (amplitude + baseline) * 100, digits=1)) * "%",
        string(round(baseline, digits=2)),
        string(round(amplitude + baseline, digits=2)),
        string(round(hill_n, digits=2)),
        string(round(exp(ln_kd), digits=1)),
    ]
end

# Pas de "τ₀"/"τ∞"/"R₀"/"R∞" ici (contrairement aux axes des graphiques, qui
# utilisent la police par défaut de GR) : la police "Courier" ci-dessous n'a
# pas ces glyphes indice Unicode ("GKS: glyph missing from current font") --
# d'où "Baseline"/"Asymptote" en toutes lettres, comme le tableau combiné
# original.
# lifetime_table_headers = ["Protein", "Δτ", "Δτ [%]", "τ_0", "τ_∞", "n", "K_D"]
# fret_table_headers      = ["Protein", "ΔR", "ΔR [%]", "R_0", "R_∞", "n", "K_D"]
lifetime_table_headers = ["Protein", "Δτ", "τ_0", "τ_∞", "n", "K_D"]
fret_table_headers      = ["Protein", "ΔR", "R_0", "R_∞", "n", "K_D"]

lifetime_table_rows = [table_row(protein, lifetime_stats[protein].params) for protein in proteins]
fret_table_rows     = [table_row(protein, fret_stats[protein].params) for protein in proteins]

# Un panneau sans axes contenant le texte formaté en colonnes alignées --
# Plots.jl n'a pas de type de série "tableau" ; "Courier" est la seule police
# à la fois monospace ET reconnue par GR sans avertissement (testé : "monospace"
# et "DejaVu Sans Mono" déclenchent tous deux un "GKS: could not find font").
# Police à 9pt (plutôt que 11pt comme le tableau combiné précédent) : chaque
# tableau n'occupe plus que la moitié de la largeur de la figure -- annotate!
# n'est pas rogné aux limites de son panneau, donc un texte trop large déborde
# carrément dans le panneau voisin plutôt que d'être coupé.
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

lifetime_table_panel = render_table_panel(lifetime_table_headers, lifetime_table_rows)
fret_table_panel     = render_table_panel(fret_table_headers, fret_table_rows)

# -----------------------------------------------------------------------------
# Figure finale : graphiques (durée de vie, ratio FRET) au-dessus, un tableau
# récapitulatif sous chacun
# -----------------------------------------------------------------------------

layout = @layout [a b; c{0.2h} d{0.2h}]
plot(lifetime_combined, fret_combined, lifetime_table_panel, fret_table_panel, layout=layout, size=(1400, 900), dpi=200)
