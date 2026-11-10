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

base_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-08-19_FLIM"

# Toujours aucun fichier IRF fourni avec ce jeu de données -- l'IRF caractérise
# l'électronique TCSPC du détecteur, pas l'échantillon, et la fenêtre TCSPC est
# identique ici (tac_r/tac_g = 50.03 ns / 4 = 12.5 ns sur 256 canaux) ; on
# réutilise donc celle de l'expérience SCLM, comme analyse_GECO.jl.
irf_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-07-20_SCLM/irf_20260710_16x_sclm.sdt"

irf_data = load_irf_from_sdt(irf_path)
RUNTIME[].irf            = irf_data
RUNTIME[].irf_bin_size   = compute_irf_bin_size(irf_data)
RUNTIME[].tcspc_window_size = round(irf_data[end, 1] + irf_data[2, 1], sigdigits=4)

# Compile la chaîne vec_to_lifetime -> Optim -> convolve une première fois :
# sans cet appel, le premier vrai fit ci-dessous a le temps du JIT (plusieurs
# secondes) qui dépasse le time_limit de l'optimiseur et renvoie NaN.
warmup_lifetime_fitting!()

# -----------------------------------------------------------------------------
# Plan d'expérience
#
# base_path/{GECO,Control}/{0..10 mM}/{750,800,850}_nm_c{1,2,3}.sdt
#
# Deux échantillons (FR-GECO et son contrôle), 11 concentrations de CaEGTA,
# 3 longueurs d'onde d'excitation x 3 réplicats = 9 fichiers par dossier.
#
# Le Control n'est pas une série affichée : il sert uniquement de bruit de fond
# à soustraire (voir "Soustraction du Control" plus bas). Les deux séries
# tracées sont donc GECO et GECO - Control.
#
# Comme dans analyse_GECO.jl, les dossiers sont nommés d'après la concentration
# de CaEGTA (0-10 mM) mais c'est le [Ca²⁺] libre résultant (μM, table de
# calibration CaEGTA) qui a un sens biologique pour FR-GECO et qui alimente le
# fit de Hill. "0 mM" (un vrai zéro) est placé à 1e-5 μM pour rester affichable
# sur l'axe log sans fausser le fit : vérifié empiriquement dans analyse_GECO.jl
# que le résultat est stable pour tout placeholder <= 1e-4.
# -----------------------------------------------------------------------------

geco_sample    = "GECO"
control_sample = "Control"

geco_label       = "GECO"
difference_label = "GECO - Control"
series_labels    = [geco_label, difference_label]

flim_wavelengths = ["750", "800", "850"]

concentration_folders  = ["0 mM", "1 mM", "2 mM", "3 mM", "4 mM", "5 mM", "6 mM", "7 mM", "8 mM", "9 mM", "10 mM"]
free_ca_by_caegta_conc = Dict("0 mM" => 0.0, "1 mM" => 0.017, "2 mM" => 0.038, "3 mM" => 0.065, "4 mM" => 0.100,
    "5 mM" => 0.150, "6 mM" => 0.225, "7 mM" => 0.351, "8 mM" => 0.602, "9 mM" => 1.35, "10 mM" => 39.0)
concentration_values = [c == 0.0 ? 1e-5 : c for c in getindex.(Ref(free_ca_by_caegta_conc), concentration_folders)]

hill(x, p) = (p[1] * x .^ p[3]) ./ (exp(p[2])^p[3] .+ x .^ p[3]) .+ p[4]
# p = [Δy, ln(K_D), n, y_0]

# -----------------------------------------------------------------------------
# Lecture des .sdt : somme spatiale
#
# Contrairement aux jeux de données précédents (histogrammes déjà réduits à un
# seul pixel), les fichiers de cette acquisition sont des piles FLIM complètes
# 1024 x 512 x 256. `read_sdt_frame` / `extract_channel_counts` ne conviennent
# donc PAS ici : ils font `reshape(vec(bloc), (256, :))`, ce qui suppose que
# l'axe temporel varie le plus vite, alors que SdtFile.read_sdt renvoie déjà un
# tableau Julia (image_y, image_x, adc_re) = (512, 1024, 256) où le temps est
# l'axe le plus LENT en mémoire column-major -- ce reshape mélangerait donc les
# pixels. On somme explicitement sur les deux axes spatiaux (dims=(1, 2)) pour
# obtenir le vecteur de 256 bins demandé.
#
# La somme est accumulée en Float64 : `extract_channel_counts` renvoie des
# UInt16, or additionner 512 x 1024 pixels dépasse largement 65535 sur le canal
# le plus brillant (pic ~7.2e4 comptes déjà à 750 nm / 10 mM).
#
# Les fichiers contiennent 2 canaux TCSPC ; on garde le canal 1, celui utilisé
# par toutes les analyses précédentes (read_sdt_frame renvoie `sdt.data[1]`).
# -----------------------------------------------------------------------------

# Canal TCSPC analysé. Les deux blocs sont de toute façon décompressés par
# read_sdt, donc analyser le canal 2 ne coûte qu'un second passage de fits.
# Surchargeable sans éditer le fichier :  SDT_CHANNEL=2 julia analyse_GECO_2.jl
# const SDT_CHANNEL = parse(Int, get(ENV, "SDT_CHANNEL", "1"))
const SDT_CHANNEL = 2

# Le canal 1 conserve les noms de sortie historiques (sans suffixe) pour que
# ses résultats déjà produits restent reproductibles à l'identique ; tout autre
# canal écrit à côté, dans des fichiers suffixés "_chN", sans rien écraser.
const CHANNEL_SUFFIX = SDT_CHANNEL == 1 ? "" : "_ch$(SDT_CHANNEL)"
const CHANNEL_TAG    = SDT_CHANNEL == 1 ? "" : " (ch $(SDT_CHANNEL))"

function spatially_summed_histogram(filepath::String)
    sdt = SdtFile.read_sdt(read(filepath), basename(filepath))

    isempty(sdt.data) && error("SDT file contains no data blocks: $filepath")
    isempty(sdt.measure_info) && error("SDT file contains no measurement info: $filepath")
    length(sdt.data) >= SDT_CHANNEL || error("SDT file has no channel $SDT_CHANNEL: $filepath")

    histogram_resolution = Int(sdt.measure_info[1].adc_re)
    block = sdt.data[SDT_CHANNEL]

    ndims(block) == 3 || error("Expected a 3-D (y, x, time) FLIM block, got size $(size(block)): $filepath")
    size(block, 3) == histogram_resolution ||
        error("Last axis ($(size(block, 3))) is not the $(histogram_resolution)-bin time axis: $filepath")

    counts = dropdims(sum(Float64, block; dims=(1, 2)), dims=(1, 2))
    return counts, histogram_resolution
end

# Les 3 réplicats d'une (échantillon, concentration, longueur d'onde) donnée,
# triés : "750_nm_c1.sdt", "750_nm_c2.sdt", "750_nm_c3.sdt".
function replicate_histograms(sample::String, conc_folder::String, wavelength::String)
    folderpath = joinpath(base_path, sample, conc_folder)
    files = sort(filter(f -> startswith(f, "$(wavelength)_nm") && endswith(f, ".sdt"), readdir(folderpath)))
    isempty(files) && error("No $(wavelength) nm .sdt file in $folderpath")

    histograms = Vector{Vector{Float64}}(undef, length(files))
    histogram_resolution = 0
    for (i, f) in enumerate(files)
        histograms[i], histogram_resolution = spatially_summed_histogram(joinpath(folderpath, f))
    end
    return histograms, histogram_resolution
end

# -----------------------------------------------------------------------------
# Soustraction du Control
#
# La soustraction se fait sur les HISTOGRAMMES, bin par bin, avant tout fit :
# la série "GECO - Control" est donc la décroissance du signal propre à GECO,
# le fond (diffusion, autofluorescence du tampon, dark counts) étant retiré du
# décompte de photons AVANT l'ajustement -- et non une différence de deux τ
# ajustés séparément.
#
# Le k-ième réplicat de GECO est apparié au k-ième du Control ; ce sont des
# acquisitions distinctes, l'appariement est positionnel et non simultané.
#
# Les bins négatifs (le fond soustrait dépasse le signal, essentiellement dans
# la queue et avant l'impulsion) sont ramenés à 0 pour le fit : l'objectif MLE
# de vec_to_lifetime est poissonien et n'a pas de sens pour un compte négatif.
# L'intensité, elle, reste la vraie différence non tronquée sum(GECO) -
# sum(Control), c.-à-d. le nombre net de photons attribuables à GECO.
# -----------------------------------------------------------------------------

subtract_background(geco::Vector{Float64}, control::Vector{Float64}) = max.(geco .- control, 0.0)

# -----------------------------------------------------------------------------
# Fits de décroissance
# -----------------------------------------------------------------------------

# vec_to_lifetime's default (first_fit=false) uses the tight ~20ms time_limit
# tuned for live acquisition display; under real system load a single fit can
# occasionally miss that budget and return NaN even after warmup_lifetime_fitting!.
# first_fit=true gives it the more generous ~60ms budget, and retrying on a NaN
# result (cheap here -- a few hundred fits total, one-off batch analysis) makes
# this robust to the rare case where even that isn't enough.
function fit_decay(counts::Vector{Float64}, histogram_resolution::Int, guess::Vector{Float64}; max_attempts::Int=10)
    for _ in 1:max_attempts
        params, data_xy = vec_to_lifetime(counts; guess=copy(guess), histogram_resolution=histogram_resolution, first_fit=true)
        isnan(params[1]) || return params, data_xy
    end
    error("Fit did not converge after $max_attempts attempts (guess length $(length(guess)))")
end

# Une seule durée de vie est tracée : le τ du modèle à 1 exponentielle.
fit_lifetime(counts::Vector{Float64}, histogram_resolution::Int) =
    fit_decay(counts, histogram_resolution, [3.0, 0.0, 1e-6])[1][1]

# -----------------------------------------------------------------------------
# Sauvegarde des histogrammes sommés : pour chaque série x concentration x
# longueur d'onde (2 x 11 x 3 = 66), les 3 réplicats sont sommés en un seul
# "gros" histogramme, ajusté avec les deux modèles (1 et 2 exponentielles) et
# sauvegardé -- 66 x 2 = 132 fichiers au total (même principe que
# analyse_GECO.jl, où c'était 22 x 2 = 44).
#
# Le fit à 2 exponentielles n'apparaît QUE dans ces PNG et dans le diagnostic
# de χ² plus bas : la figure de calibration n'affiche qu'une seule durée de vie.
# -----------------------------------------------------------------------------

histogram_plots_folder = joinpath(@__DIR__, "histogram_fits_GECO_2$(CHANNEL_SUFFIX)")
mkpath(histogram_plots_folder)

sanitize_for_filename(s::String) = replace(s, " " => "")

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

    # Cadrage sur la partie exploitable de la décroissance ; le seuil de 100
    # comptes peut n'être atteint par aucun bin sur les acquisitions les plus
    # faibles (850 nm), d'où le repli sur le plancher d'affichage.
    above_floor = filter(>(100), y_data_raw)
    ylims!(p, (isempty(above_floor) ? 0.5 : minimum(above_floor), ylims(p)[2]))

    savefig(p, filepath)
end

# Renvoie les deux χ² réduits (1 et 2 exponentielles) de l'histogramme sommé,
# consommés par le diagnostic de choix de modèle plus bas.
function save_summed_histogram_fits(series::String, wavelength::String, conc_folder::String,
                                    summed_counts::Vector{Float64}, histogram_resolution::Int)
    series_label = sanitize_for_filename(series)
    conc_label   = sanitize_for_filename(conc_folder)

    params_1exp, data_xy_1exp = fit_decay(summed_counts, histogram_resolution, [3.0, 0.0, 1e-6])
    params_2exp, data_xy_2exp = fit_decay(summed_counts, histogram_resolution, [3.0, 0.5, 0.5, 0.0, 1e-6])

    chi2_1 = fit_reduced_chi_squared(params_1exp, data_xy_1exp)
    chi2_2 = fit_reduced_chi_squared(params_2exp, data_xy_2exp)

    save_decay_fit_plot(data_xy_1exp, params_1exp, chi2_1,
        "$series$(CHANNEL_TAG) $(wavelength) nm $conc_folder (sum of 3) -- 1-exponential",
        joinpath(histogram_plots_folder, "$(series_label)_$(wavelength)nm_$(conc_label)_sum_1exp.png"))
    save_decay_fit_plot(data_xy_2exp, params_2exp, chi2_2,
        "$series$(CHANNEL_TAG) $(wavelength) nm $conc_folder (sum of 3) -- 2-exponential",
        joinpath(histogram_plots_folder, "$(series_label)_$(wavelength)nm_$(conc_label)_sum_2exp.png"))

    return chi2_1, chi2_2
end

# -----------------------------------------------------------------------------
# Analyse : durée de vie et intensité pour chaque fichier
#
# Chaque fichier n'est lu qu'une seule fois (la décompression d'une pile de
# 134 M de bins prend ~3 s, et il y en a 198) : la somme spatiale sert à la
# fois au fit par réplicat, à la soustraction du fond et à l'histogramme sommé
# des 3 réplicats.
#
# L'intensité est le nombre total de photons du fichier, exprimé en 10⁶ photons
# pour que les axes et les tableaux restent lisibles.
# -----------------------------------------------------------------------------

const PHOTONS_PER_UNIT = 1e6

# série -> longueur d'onde -> concentration -> valeur par réplicat
tau_results       = Dict(s => Dict(w => Dict{String, Vector{Float64}}() for w in flim_wavelengths) for s in series_labels)
intensity_results = Dict(s => Dict(w => Dict{String, Vector{Float64}}() for w in flim_wavelengths) for s in series_labels)

chi2_1_by_wavelength = Dict(w => Float64[] for w in flim_wavelengths)
chi2_2_by_wavelength = Dict(w => Float64[] for w in flim_wavelengths)

for wavelength in flim_wavelengths, conc_folder in concentration_folders
    geco_histograms,    histogram_resolution = replicate_histograms(geco_sample,    conc_folder, wavelength)
    control_histograms, _                    = replicate_histograms(control_sample, conc_folder, wavelength)

    length(geco_histograms) == length(control_histograms) ||
        error("Mismatched replicate counts for $(wavelength)nm / $conc_folder: " *
              "$(length(geco_histograms)) GECO vs $(length(control_histograms)) Control")

    difference_histograms = [subtract_background(g, c) for (g, c) in zip(geco_histograms, control_histograms)]

    histograms_by_series = Dict(geco_label => geco_histograms, difference_label => difference_histograms)

    # L'intensité de la différence est la vraie différence de photons (non
    # tronquée), pas la somme des bins écrêtés à 0 utilisés pour le fit.
    intensities_by_series = Dict(
        geco_label       => [sum(g) / PHOTONS_PER_UNIT for g in geco_histograms],
        difference_label => [(sum(g) - sum(c)) / PHOTONS_PER_UNIT for (g, c) in zip(geco_histograms, control_histograms)],
    )

    for series in series_labels
        histograms = histograms_by_series[series]

        tau_results[series][wavelength][conc_folder]       = [fit_lifetime(h, histogram_resolution) for h in histograms]
        intensity_results[series][wavelength][conc_folder] = intensities_by_series[series]

        chi2_1, chi2_2 = save_summed_histogram_fits(series, wavelength, conc_folder, sum(histograms), histogram_resolution)
        push!(chi2_1_by_wavelength[wavelength], chi2_1)
        push!(chi2_2_by_wavelength[wavelength], chi2_2)

        println("$series / $(wavelength)nm / $conc_folder : τ=", round.(tau_results[series][wavelength][conc_folder], digits=3),
                " I=", round.(intensity_results[series][wavelength][conc_folder], digits=3), " ×10⁶ ph")
    end
end

# -----------------------------------------------------------------------------
# Choix du modèle de décroissance (1 vs 2 exponentielles) par longueur d'onde,
# à partir du χ² moyen de chaque modèle sur tous les histogrammes sommés (les
# deux séries et toutes les concentrations confondues) -- le modèle au χ² moyen
# le plus bas est retenu. Purement diagnostique (imprimé en console) ; n'affecte
# pas la durée de vie (1-exponentielle) qui alimente les graphiques.
# -----------------------------------------------------------------------------

for wavelength in flim_wavelengths
    mean_chi2_1 = mean(chi2_1_by_wavelength[wavelength])
    mean_chi2_2 = mean(chi2_2_by_wavelength[wavelength])
    n_exp = mean_chi2_2 < mean_chi2_1 ? 2 : 1

    println("$(wavelength)nm : ⟨χ²⟩ 1-exp = $(round(mean_chi2_1, digits=2)), ⟨χ²⟩ 2-exp = $(round(mean_chi2_2, digits=2)) -> $n_exp-exponential")
end

# -----------------------------------------------------------------------------
# Fit de Hill (partagé entre durée de vie et intensité)
# -----------------------------------------------------------------------------

function gather_xy(data_dict::Dict{String, Dict{String, Vector{Float64}}}, wavelength::String)
    x_vals = Float64[]
    y_vals = Float64[]
    for (conc_value, conc_folder) in zip(concentration_values, concentration_folders)
        vals = data_dict[wavelength][conc_folder]
        append!(x_vals, fill(conc_value, length(vals)))
        append!(y_vals, vals)
    end
    return x_vals, y_vals
end

# Bornes volontairement plus larges que dans analyse_GECO.jl (qui imposait
# Δy >= 0 et y_0 >= 0) : après soustraction du fond, une série peut
# parfaitement décroître avec [Ca²⁺] ou passer sous zéro, et une borne à 0 y
# écraserait le fit contre la contrainte au lieu de suivre les données.
# L'amplitude et la ligne de base de départ sont dérivées des données
# elles-mêmes (moyennes aux deux extrémités de la gamme) plutôt que codées en
# dur, pour que le même code serve à τ (~ns) et à I (~10⁶ photons).
function fit_hill(x_vals::Vector{Float64}, y_vals::Vector{Float64})
    x_lo, x_hi = minimum(x_vals), maximum(x_vals)
    y_lo = mean(y_vals[x_vals .== x_lo])
    y_hi = mean(y_vals[x_vals .== x_hi])

    initial_params = [y_hi - y_lo, log(1.0), 1.0, y_lo]
    fit = curve_fit(hill, x_vals, y_vals, initial_params; lower=[-Inf, -16.0, 0.0, -Inf], upper=[Inf, 16.0, Inf, Inf])
    return fit.param
end

function hill_summary_lines(params::Vector{Float64}; base::String, unit::String)
    amplitude, ln_kd, hill_n, baseline = params
    return [
        "Δ$(base) = $(round(amplitude, digits=3)) $(unit)",
        "Δ$(base) [%] = $(round(amplitude / (amplitude + baseline) * 100, digits=1)) %",
        "$(base)₀ = $(round(baseline, digits=3)) $(unit)",
        "$(base)_∞ = $(round(amplitude + baseline, digits=3)) $(unit)",
        "n = $(round(hill_n, digits=2))",
        "K_D = $(round(exp(ln_kd), digits=3)) μM",
    ]
end

xfit = 10 .^ range(-5.3, stop=2, length=200)

function compute_metric_stats(data_dict::Dict{String, Dict{String, Vector{Float64}}}; row_label::String, base::String, unit::String)
    stats = Dict{String, NamedTuple}()
    for wavelength in flim_wavelengths
        x_vals, y_vals = gather_xy(data_dict, wavelength)
        params = fit_hill(x_vals, y_vals)

        println("$row_label / $(wavelength) nm")
        foreach(println, hill_summary_lines(params; base=base, unit=unit))

        stats[wavelength] = (x=x_vals, y=y_vals, yfit=hill(xfit, params), params=params)
    end
    return stats
end

tau_stats       = Dict(s => compute_metric_stats(tau_results[s];       row_label="τ -- $s", base="τ", unit="ns") for s in series_labels)
intensity_stats = Dict(s => compute_metric_stats(intensity_results[s]; row_label="I -- $s", base="I", unit="10⁶ ph") for s in series_labels)

# -----------------------------------------------------------------------------
# Graphiques : un panneau par (série, mesure), les 3 longueurs d'onde
# superposées dans chacun. Axe log comme dans analyse_GECO.jl -- le [Ca²⁺]
# libre testé couvre 0.017-39 μM, presque 4 ordres de grandeur.
# -----------------------------------------------------------------------------

marker_shapes = [:circle, :rect, :diamond]
xticks_log = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ["0", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹", "10²"])

base_panel(xlabel, ylabel) = plot(xscale=:log10, framestyle=:box, grid=false,
    tickfontsize=11, guidefontsize=13, legendfontsize=10, titlefontsize=13,
    xlabel=xlabel, ylabel=ylabel, xticks=xticks_log, xlims=(6e-6, 150), palette=:Set1_3,
    left_margin=10Plots.mm, right_margin=3Plots.mm, top_margin=2Plots.mm, bottom_margin=2Plots.mm)

function make_combined_panel(stats::Dict{String, NamedTuple}; title::String, ylabel::String, xlabel::String)
    panel = base_panel(xlabel, ylabel)
    plot!(panel, title=title, legend=:topleft)

    for (i, wavelength) in enumerate(flim_wavelengths)
        d = stats[wavelength]
        scatter!(panel, d.x, d.y, label="$(wavelength) nm", color=i, markershape=marker_shapes[i],
            markersize=5, markeralpha=0.5, markerstrokewidth=0.75, markerstrokecolor=:black)
        plot!(panel, xfit, d.yfit, label=nothing, lw=2, color=i)
    end

    # Un peu plus de marge que le padding automatique de Plots au-dessus/en
    # dessous des données, pour que les courbes ne touchent pas le cadre --
    # dissymétrique parce que la légende (:topleft) occupe le haut du panneau.
    all_y = reduce(vcat, (stats[w].y for w in flim_wavelengths))
    ymin, ymax = extrema(all_y)
    yspan = max(ymax - ymin, eps())
    ylims!(panel, (ymin - 0.10 * yspan, ymax + 0.24 * yspan))

    return panel
end

# Seule la dernière rangée porte l'étiquette de l'axe x (les deux rangées
# partagent la même gamme de concentrations).
tau_panels = [make_combined_panel(tau_stats[s];
    title="τ -- $s$(CHANNEL_TAG)", ylabel="τ [ns]",
    xlabel=(s == series_labels[end] ? "[Ca²⁺] [μM]" : "")) for s in series_labels]

intensity_panels = [make_combined_panel(intensity_stats[s];
    title="I -- $s$(CHANNEL_TAG)", ylabel="I [10⁶ photons]",
    xlabel=(s == series_labels[end] ? "[Ca²⁺] [μM]" : "")) for s in series_labels]

# -----------------------------------------------------------------------------
# Tableaux récapitulatifs (mêmes valeurs que hill_summary_lines, une ligne par
# longueur d'onde), un par panneau, affiché sous celui-ci -- "les descriptions
# sous les graphiques", comme dans analyse_GECO.jl.
# -----------------------------------------------------------------------------

function table_row(wavelength::String, params::Vector{Float64})
    amplitude, ln_kd, hill_n, baseline = params
    return [
        "$(wavelength) nm",
        string(round(amplitude, digits=3)),
        string(round(baseline, digits=3)),
        string(round(amplitude + baseline, digits=3)),
        string(round(hill_n, digits=2)),
        string(round(exp(ln_kd), digits=3)),
    ]
end

# Pas de "τ₀"/"τ∞" ici (contrairement aux axes des graphiques, qui utilisent la
# police par défaut de GR) : la police "Courier" ci-dessous n'a pas ces glyphes
# indice Unicode ("GKS: glyph missing from current font") -- d'où "τ_0"/"τ_inf".
# Les unités de la grandeur ajustée ne sont PAS répétées dans les en-têtes :
# elles figurent déjà sur l'axe y du panneau au-dessus, et les répéter rendait
# le tableau d'intensité plus large que son panneau -- annotate! n'étant pas
# rogné à son cadre, la colonne K_D débordait alors hors de la figure.
tau_table_headers       = ["Wavelength", "Δτ", "τ_0", "τ_inf", "n", "K_D [uM]"]
intensity_table_headers = ["Wavelength", "ΔI", "I_0", "I_inf", "n", "K_D [uM]"]

# Un panneau sans axes contenant le texte formaté en colonnes alignées --
# Plots.jl n'a pas de type de série "tableau" ; "Courier" est la seule police à
# la fois monospace ET reconnue par GR sans avertissement (testé : "monospace"
# et "DejaVu Sans Mono" déclenchent tous deux un "GKS: could not find font").
function render_table_panel(headers::Vector{String}, rows::Vector{Vector{String}})
    all_rows = vcat([headers], rows)
    ncols = length(headers)
    widths = [maximum(length(r[j]) for r in all_rows) for j in 1:ncols]
    fmt_row(r) = join([rpad(r[j], widths[j]) for j in 1:ncols], "  ")

    lines = [fmt_row(headers), fmt_row(["-"^widths[j] for j in 1:ncols])]
    append!(lines, [fmt_row(r) for r in rows])
    table_text = join(lines, "\n")

    panel = plot(axis=false, grid=false, framestyle=:none, xlims=(0, 1), ylims=(0, 1))
    annotate!(panel, 0.0, 1.0, text(table_text, 10, :left, :top, "Courier"))
    return panel
end

tau_table_panels = [render_table_panel(tau_table_headers,
    [table_row(w, tau_stats[s][w].params) for w in flim_wavelengths]) for s in series_labels]

intensity_table_panels = [render_table_panel(intensity_table_headers,
    [table_row(w, intensity_stats[s][w].params) for w in flim_wavelengths]) for s in series_labels]

# -----------------------------------------------------------------------------
# Export CSV : valeurs par réplicat et résumé des fits de Hill (les
# histogrammes bruts et les fits sur les sommes des 3 réplicats restent
# uniquement dans les PNG de histogram_fits_GECO_2/).
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

raw_data_rows = Vector{String}[]
for series in series_labels, wavelength in flim_wavelengths
    for (conc_value, conc_folder) in zip(concentration_values, concentration_folders)
        taus        = tau_results[series][wavelength][conc_folder]
        intensities = intensity_results[series][wavelength][conc_folder]
        for i in eachindex(taus)
            push!(raw_data_rows, [series, wavelength, conc_folder, string(conc_value), string(i),
                string(taus[i]), string(intensities[i])])
        end
    end
end

write_csv(joinpath(csv_folder, "geco_2_raw_data$(CHANNEL_SUFFIX).csv"),
    ["series", "wavelength_nm", "caegta_folder", "free_ca_uM", "replicate", "tau_ns", "intensity_Mphotons"],
    raw_data_rows)

function hill_summary_row(series::String, wavelength::String, metric_name::String, params::Vector{Float64})
    amplitude, ln_kd, hill_n, baseline = params
    delta_pct = amplitude / (amplitude + baseline) * 100
    return [series, wavelength, metric_name, string(amplitude), string(delta_pct), string(baseline),
        string(amplitude + baseline), string(hill_n), string(exp(ln_kd))]
end

hill_summary_rows = Vector{String}[]
for series in series_labels, wavelength in flim_wavelengths
    push!(hill_summary_rows, hill_summary_row(series, wavelength, "tau_ns", tau_stats[series][wavelength].params))
    push!(hill_summary_rows, hill_summary_row(series, wavelength, "intensity_Mphotons", intensity_stats[series][wavelength].params))
end

write_csv(joinpath(csv_folder, "geco_2_hill_fit_summary$(CHANNEL_SUFFIX).csv"),
    ["series", "wavelength_nm", "metric", "delta", "delta_pct", "baseline", "asymptote", "n", "kd_uM"],
    hill_summary_rows)

write_csv(joinpath(csv_folder, "geco_2_decay_model_selection$(CHANNEL_SUFFIX).csv"),
    ["wavelength_nm", "mean_chi2_1exp", "mean_chi2_2exp", "selected_model"],
    [[w, string(mean(chi2_1_by_wavelength[w])), string(mean(chi2_2_by_wavelength[w])),
      string(mean(chi2_2_by_wavelength[w]) < mean(chi2_1_by_wavelength[w]) ? 2 : 1)] for w in flim_wavelengths])

# -----------------------------------------------------------------------------
# Figure finale : 2 rangées (GECO, GECO - Control) x 2 colonnes (durée de vie,
# intensité), chaque panneau suivi de son tableau récapitulatif
# -----------------------------------------------------------------------------

layout = @layout [a b; c{0.11h} d{0.11h}; e f; g{0.11h} h{0.11h}]

final_plot = plot(
    tau_panels[1], intensity_panels[1], tau_table_panels[1], intensity_table_panels[1],
    tau_panels[2], intensity_panels[2], tau_table_panels[2], intensity_table_panels[2],
    layout=layout, size=(1500, 1350), dpi=150,
)

savefig(final_plot, joinpath(@__DIR__, "calibration_GECO_2$(CHANNEL_SUFFIX).png"))
final_plot
