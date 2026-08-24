using Plots
using LsqFit
using Statistics

# -----------------------------------------------------------------------------
# Figure comparative des deux canaux de détection, série GECO - Control.
#
# Ce script ne relit AUCUN fichier .sdt et ne refait AUCUN fit de décroissance :
# il repart des CSV produits par analyse_GECO_2.jl (un par canal), qui
# contiennent déjà, pour chaque longueur d'onde x concentration x réplicat, le
# τ à 1 exponentielle et le nombre total de photons. Seuls les fits de Hill
# sont refaits ici -- une régression à 4 paramètres sur 33 points, instantanée.
#
# Prérequis (déjà présents) :
#   csv_exports/geco_2_raw_data.csv       <- canal 1, produit par défaut
#   csv_exports/geco_2_raw_data_ch2.csv   <- canal 2, SDT_CHANNEL=2
# -----------------------------------------------------------------------------

csv_folder = joinpath(@__DIR__, "csv_exports")

# Les canaux sont désignés par leur filtre d'émission plutôt que par leur index
# TCSPC : c'est le filtre qui a un sens physique. Canal 1 = 675/65 (l'émission
# rouge lointaine de FR-GECO), canal 2 = 530/30.
channels = [
    (label="675/65", file="geco_2_raw_data.csv"),
    (label="530/30", file="geco_2_raw_data_ch2.csv"),
]

const PLOTTED_SERIES = "GECO - Control"

flim_wavelengths = ["750", "800", "850"]

hill(x, p) = (p[1] * x .^ p[3]) ./ (exp(p[2])^p[3] .+ x .^ p[3]) .+ p[4]
# p = [Δy, ln(K_D), n, y_0]

# -----------------------------------------------------------------------------
# Lecture des CSV
#
# Les colonnes n'ont ni guillemets ni virgules internes (write_csv de
# analyse_GECO_2.jl écrit des nombres et des étiquettes simples), un split sur
# la virgule suffit donc -- pas de dépendance à CSV.jl.
# -----------------------------------------------------------------------------

function read_raw_data(filepath::String)
    isfile(filepath) || error("Missing CSV: $filepath -- run analyse_GECO_2.jl first (SDT_CHANNEL=2 for channel 2)")

    lines  = readlines(filepath)
    header = split(lines[1], ',')
    col    = Dict(name => i for (i, name) in enumerate(header))

    rows = [split(l, ',') for l in lines[2:end] if !isempty(l)]
    filter!(r -> r[col["series"]] == PLOTTED_SERIES, rows)
    isempty(rows) && error("No \"$PLOTTED_SERIES\" rows in $filepath")

    return rows, col
end

# La gamme de concentrations est reconstruite depuis le CSV (colonnes
# caegta_folder / free_ca_uM) plutôt que redéclarée ici : la table de
# calibration CaEGTA -> [Ca²⁺] libre reste ainsi définie à un seul endroit,
# dans analyse_GECO_2.jl. Le "0 mM" y porte déjà son placeholder de 1e-5 μM,
# nécessaire pour l'axe log.
function concentration_axis(rows, col)
    pairs = unique([(r[col["caegta_folder"]], parse(Float64, r[col["free_ca_uM"]])) for r in rows])
    sort!(pairs, by=last)
    return first.(pairs), last.(pairs)
end

# métrique -> longueur d'onde -> concentration -> valeurs par réplicat
function group_values(rows, col, value_column::String, conc_folders::Vector{<:AbstractString})
    grouped = Dict(w => Dict{String, Vector{Float64}}(c => Float64[] for c in conc_folders) for w in flim_wavelengths)
    for r in rows
        wavelength = String(r[col["wavelength_nm"]])
        haskey(grouped, wavelength) || continue
        push!(grouped[wavelength][String(r[col["caegta_folder"]])], parse(Float64, r[col[value_column]]))
    end
    return grouped
end

# F/F_max : l'intensité est normalisée par sa valeur maximale, séparément pour
# chaque canal et chaque longueur d'onde -- c'est ce qui rend les trois
# longueurs d'onde superposables sur un même axe malgré des comptes de photons
# séparés par plus d'un ordre de grandeur.
#
# F_max est pris sur les MOYENNES par concentration, pas sur les réplicats
# individuels, pour qu'un seul point bruité ne fixe pas la normalisation ; un
# réplicat peut donc dépasser 1 de peu, ce qui est voulu.
#
# Le maximum ne tombe pas au même bout de la gamme selon le canal : à 675/65
# l'intensité croît avec [Ca²⁺] (F_max à 39 μM), à 530/30 elle décroît (F_max
# du côté des faibles [Ca²⁺]). Prendre le maximum plutôt que la valeur à une
# concentration fixée traite donc les deux canaux sans convention ad hoc.
function normalize_to_max(intensity_by_wavelength, conc_folders::Vector{<:AbstractString})
    normalized = Dict{String, Dict{String, Vector{Float64}}}()
    for wavelength in flim_wavelengths
        by_conc = intensity_by_wavelength[wavelength]
        f_max = maximum(mean(by_conc[c]) for c in conc_folders)
        normalized[wavelength] = Dict(c => by_conc[c] ./ f_max for c in conc_folders)
        println("F_max ($(wavelength) nm) = $(round(f_max, digits=4)) ×10⁶ photons")
    end
    return normalized
end

# -----------------------------------------------------------------------------
# Fit de Hill -- identique à celui de analyse_GECO_2.jl (bornes larges : après
# soustraction du fond une série peut décroître, ce que fait justement
# l'intensité du canal 530/30, et une borne à 0 sur l'amplitude écraserait le
# fit contre la contrainte).
# -----------------------------------------------------------------------------

function gather_xy(data_dict, wavelength::String, conc_folders, conc_values)
    x_vals = Float64[]
    y_vals = Float64[]
    for (conc_value, conc_folder) in zip(conc_values, conc_folders)
        vals = data_dict[wavelength][conc_folder]
        append!(x_vals, fill(conc_value, length(vals)))
        append!(y_vals, vals)
    end
    return x_vals, y_vals
end

function fit_hill(x_vals::Vector{Float64}, y_vals::Vector{Float64})
    x_lo, x_hi = minimum(x_vals), maximum(x_vals)
    y_lo = mean(y_vals[x_vals .== x_lo])
    y_hi = mean(y_vals[x_vals .== x_hi])

    initial_params = [y_hi - y_lo, log(1.0), 1.0, y_lo]
    fit = curve_fit(hill, x_vals, y_vals, initial_params; lower=[-Inf, -16.0, 0.0, -Inf], upper=[Inf, 16.0, Inf, Inf])
    return fit.param
end

xfit = 10 .^ range(-5.3, stop=2, length=200)

function compute_metric_stats(data_dict, conc_folders, conc_values; row_label::String, base::String, unit::String)
    stats = Dict{String, NamedTuple}()
    for wavelength in flim_wavelengths
        x_vals, y_vals = gather_xy(data_dict, wavelength, conc_folders, conc_values)
        params = fit_hill(x_vals, y_vals)
        amplitude, ln_kd, hill_n, baseline = params

        println("$row_label / $(wavelength) nm : Δ$(base) = $(round(amplitude, digits=3)) $unit, " *
                "$(base)₀ = $(round(baseline, digits=3)) $unit, $(base)_∞ = $(round(amplitude + baseline, digits=3)) $unit, " *
                "n = $(round(hill_n, digits=2)), K_D = $(round(exp(ln_kd), digits=3)) μM")

        stats[wavelength] = (x=x_vals, y=y_vals, yfit=hill(xfit, params), params=params)
    end
    return stats
end

# -----------------------------------------------------------------------------
# Panneaux -- mêmes conventions graphiques que analyse_GECO_2.jl
# -----------------------------------------------------------------------------

marker_shapes = [:circle, :rect, :diamond]
xticks_log = ([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100], ["0", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰", "10¹", "10²"])

base_panel(xlabel, ylabel) = plot(xscale=:log10, framestyle=:box, grid=false,
    tickfontsize=11, guidefontsize=13, legendfontsize=10, titlefontsize=13,
    xlabel=xlabel, ylabel=ylabel, xticks=xticks_log, xlims=(6e-6, 150), palette=:Set1_3,
    left_margin=10Plots.mm, right_margin=3Plots.mm, top_margin=2Plots.mm, bottom_margin=2Plots.mm)

function make_combined_panel(stats::Dict{String, NamedTuple}; title::String, ylabel::String, xlabel::String, legend_position::Symbol)
    panel = base_panel(xlabel, ylabel)
    plot!(panel, title=title, legend=legend_position)

    for (i, wavelength) in enumerate(flim_wavelengths)
        d = stats[wavelength]
        scatter!(panel, d.x, d.y, label="$(wavelength) nm", color=i, markershape=marker_shapes[i],
            markersize=5, markeralpha=0.5, markerstrokewidth=0.75, markerstrokecolor=:black)
        plot!(panel, xfit, d.yfit, label=nothing, lw=2, color=i)
    end

    # Marge dissymétrique : la légende occupe un coin du panneau, et le côté
    # qu'elle occupe dépend du sens de variation de la série (croissante à
    # 675/65, décroissante à 530/30) -- d'où legend_position en argument.
    all_y = reduce(vcat, (stats[w].y for w in flim_wavelengths))
    ymin, ymax = extrema(all_y)
    yspan = max(ymax - ymin, eps())
    pad_low, pad_high = legend_position == :bottomleft ? (0.24, 0.10) : (0.10, 0.24)
    ylims!(panel, (ymin - pad_low * yspan, ymax + pad_high * yspan))

    return panel
end

# Pas de glyphes indice Unicode ici : la police "Courier" (la seule monospace
# que GR accepte sans avertissement) n'a pas "τ₀"/"τ∞" -- d'où "τ_0"/"τ_inf".
# Les unités ne sont pas répétées dans les en-têtes, elles figurent déjà sur
# l'axe y du panneau au-dessus ; annotate! n'étant pas rogné à son cadre, un
# tableau trop large déborderait hors de la figure.
function render_table_panel(headers::Vector{String}, rows::Vector{Vector{String}})
    all_rows = vcat([headers], rows)
    ncols = length(headers)
    widths = [maximum(length(r[j]) for r in all_rows) for j in 1:ncols]
    fmt_row(r) = join([rpad(r[j], widths[j]) for j in 1:ncols], "  ")

    lines = [fmt_row(headers), fmt_row(["-"^widths[j] for j in 1:ncols])]
    append!(lines, [fmt_row(r) for r in rows])

    panel = plot(axis=false, grid=false, framestyle=:none, xlims=(0, 1), ylims=(0, 1))
    annotate!(panel, 0.0, 1.0, text(join(lines, "\n"), 10, :left, :top, "Courier"))
    return panel
end

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

make_table(stats::Dict{String, NamedTuple}, symbol::String) =
    render_table_panel(["Wavelength", "Δ$symbol", "$(symbol)_0", "$(symbol)_inf", "n", "K_D [uM]"],
        [table_row(w, stats[w].params) for w in flim_wavelengths])

# -----------------------------------------------------------------------------
# Construction de la figure : 2 rangées (canaux) x 3 colonnes (τ, photons,
# F/F_max), chaque panneau suivi de son tableau récapitulatif
# -----------------------------------------------------------------------------

plot_panels  = Plots.Plot[]
table_panels = Plots.Plot[]

for (row_index, channel) in enumerate(channels)
    is_last_row = row_index == length(channels)
    xlabel = is_last_row ? "[Ca²⁺] [μM]" : ""

    rows, col = read_raw_data(joinpath(csv_folder, channel.file))
    conc_folders, conc_values = concentration_axis(rows, col)

    println("\n=== $(channel.label) ($(PLOTTED_SERIES)) ===")

    tau_by_wavelength       = group_values(rows, col, "tau_ns", conc_folders)
    intensity_by_wavelength = group_values(rows, col, "intensity_Mphotons", conc_folders)
    normalized_by_wavelength = normalize_to_max(intensity_by_wavelength, conc_folders)

    tau_stats        = compute_metric_stats(tau_by_wavelength,        conc_folders, conc_values; row_label="$(channel.label) τ",       base="τ", unit="ns")
    intensity_stats  = compute_metric_stats(intensity_by_wavelength,  conc_folders, conc_values; row_label="$(channel.label) I",       base="I", unit="10⁶ ph")
    normalized_stats = compute_metric_stats(normalized_by_wavelength, conc_folders, conc_values; row_label="$(channel.label) F/F_max", base="F", unit="")

    # L'intensité et F/F_max décroissent avec [Ca²⁺] sur le canal 530/30 : la
    # légende y va en bas à gauche pour ne pas recouvrir les points.
    intensity_legend = mean(intensity_stats["850"].params[1]) < 0 ? :bottomleft : :topleft

    append!(plot_panels, [
        make_combined_panel(tau_stats;        title="$(channel.label) -- τ",       ylabel="τ [ns]",            xlabel=xlabel, legend_position=:topleft),
        make_combined_panel(intensity_stats;  title="$(channel.label) -- photons", ylabel="I [10⁶ photons]",   xlabel=xlabel, legend_position=intensity_legend),
        make_combined_panel(normalized_stats; title="$(channel.label) -- F/F_max", ylabel="F/F_max",           xlabel=xlabel, legend_position=intensity_legend),
    ])
    append!(table_panels, [make_table(tau_stats, "τ"), make_table(intensity_stats, "I"), make_table(normalized_stats, "F")])
end

layout = @layout [a b c; d{0.11h} e{0.11h} f{0.11h}; g h i; j{0.11h} k{0.11h} l{0.11h}]

final_plot = plot(
    plot_panels[1], plot_panels[2], plot_panels[3], table_panels[1], table_panels[2], table_panels[3],
    plot_panels[4], plot_panels[5], plot_panels[6], table_panels[4], table_panels[5], table_panels[6],
    layout=layout, size=(2100, 1350), dpi=150,
)

savefig(final_plot, joinpath(@__DIR__, "calibration_GECO_2_channels.png"))
final_plot
