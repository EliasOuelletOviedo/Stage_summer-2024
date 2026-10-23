using Plots

filepath = "/Users/eliasouellet-oviedo/Desktop/FPbase_Spectra2.csv.csv"

# ============================== RÉGLAGES ==============================

# Limites de l'axe des longueurs d'onde (nm)
xlims_nm = (400, 800)

# Seuil et transparence pour l'affichage des filtres (bande de fond)
filter_threshold = 0.5
filter_alpha = 0.25

# Noms des colonnes (protéines -> courbes) et leur couleur
protein_columns = [
    "mKelly2 EM"     => :red,
    "mCerulean3 EM"  => :dodgerblue,
    "Topaz EM"       => :gold,
]

# Noms des colonnes (filtres -> bandes de fond) et leur couleur
filter_columns = [
    "Semrock FF02-655/40" => :red,
    "Chroma D470/40m"     => :dodgerblue,
    "Chroma ET530/30m"    => :gold,
]

# Couleurs du fond et des axes
background_color = :black
axis_color = :white

# Highlights manuels (plage nm, couleur, alpha)
highlight_regions = [
    (570, 600, :orange, 0.25),
]

# ========================================================================

# ---- Lecture et parsing du fichier (";" comme séparateur, "," comme décimale) ----
raw = read(filepath, String)
raw = replace(raw, "﻿" => "")  # retire le BOM éventuel
lines = filter(!isempty, split(raw, r"\r?\n"))

header = [strip(h, '"') for h in split(lines[1], ';')]

parse_field(f) = isempty(f) ? missing : parse(Float64, replace(f, ',' => '.'))

data = [parse_field.(split(line, ';')) for line in lines[2:end]]

wavelengths = Float64[row[1] for row in data]

column_values(name) = begin
    idx = findfirst(==(name), header)
    idx === nothing && error("Colonne \"$name\" introuvable dans le fichier.")
    [coalesce(row[idx], NaN) for row in data]
end

# ---- Détection des plages contiguës où filtre > seuil ----
function contiguous_ranges(mask, x)
    ranges = Tuple{Float64,Float64}[]
    i, n = 1, length(mask)
    while i <= n
        if mask[i]
            j = i
            while j + 1 <= n && mask[j + 1]
                j += 1
            end
            push!(ranges, (x[i], x[j]))
            i = j + 1
        else
            i += 1
        end
    end
    return ranges
end

# ---- Construction du graphique ----
plt = plot(
    xlims = xlims_nm,
    ylims = (0, 1.05),
    framestyle = :axes,
    grid = false,
    legend = false,
    ticks = false,
    xlabel = "",
    ylabel = "",
    title = "",
    background_color = background_color,
    background_color_outside = background_color,
    foreground_color = axis_color,
    size = (900, 500),
)

# Bandes de fond pour les filtres (dessinées en premier, derrière les courbes)
for (name, color) in filter_columns
    values = column_values(name)
    mask = values .> filter_threshold
    ranges = contiguous_ranges(mask, wavelengths)
    for (x1, x2) in ranges
        band = Shape([x1, x2, x2, x1], [0, 0, 1.05, 1.05])
        plot!(plt, band, fillcolor = color, fillalpha = filter_alpha, linealpha = 0)
    end
end

# Highlights manuels
for (x1, x2, color, alpha) in highlight_regions
    band = Shape([x1, x2, x2, x1], [0, 0, 1.05, 1.05])
    plot!(plt, band, fillcolor = color, fillalpha = alpha, linealpha = 0)
end

# Courbes pour les protéines
for (name, color) in protein_columns
    values = column_values(name)
    plot!(plt, wavelengths, values, color = color, linewidth = 2)
end

display(plt)
