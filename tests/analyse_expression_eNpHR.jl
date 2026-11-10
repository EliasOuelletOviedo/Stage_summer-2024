using Plots
using FileIO
using Images

# -----------------------------------------------------------------------------
# Dossier de données : 6 champs, chacun imagé en fluorescence (YFP) et en
# transmission. Les fichiers d'une même paire partagent le numéro de champ :
#   2026-09-18_eNpHR_expression_10k_<i>.jpg        -> fluorescence
#   2026-09-18_eNpHR_expression_10k_trans_<i>.jpg  -> transmission
# -----------------------------------------------------------------------------

base_path = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-09-18"
prefix    = "2026-09-18_eNpHR_expression_10k"
n_fields  = 6

fluo_path(i)  = joinpath(base_path, "$(prefix)_$(i).jpg")
trans_path(i) = joinpath(base_path, "$(prefix)_trans_$(i).jpg")

image_panel(img, title_str) = plot(img, aspect_ratio=:equal, axis=nothing, border=:none,
    title=title_str, titlefontsize=13)

# Colonne vide entre les deux paires d'une même rangée, pour que chaque paire
# se lise comme un bloc.
spacer() = plot(framestyle=:none, grid=false, axis=nothing)

# La caméra enregistre la fluo YFP dans le canal vert : on réaffiche cette
# intensité en jaune (rouge + vert) pour correspondre au fluorophore.
to_yellow(img) = map(c -> RGB(green(c), green(c), zero(green(c))), img)

# -----------------------------------------------------------------------------
# Figure : 3 rangées × 2 paires (fluorescence | transmission) côte à côte
# -----------------------------------------------------------------------------

panels = Plots.Plot[]
for row in 1:3
    for (k, i) in enumerate((2row - 1, 2row))
        k == 2 && push!(panels, spacer())
        push!(panels, image_panel(to_yellow(load(fluo_path(i))), "Champ $i — YFP (fluorescence)"))
        push!(panels, image_panel(load(trans_path(i)), "Champ $i — transmission"))
    end
end

final_plot = plot(panels...,
    layout=grid(3, 5, widths=[0.24, 0.24, 0.04, 0.24, 0.24]),
    size=(2100, 1250), dpi=100, margin=2Plots.mm)

savefig(final_plot, joinpath(@__DIR__, "expression_eNpHR_2026-09-18.png"))
display(final_plot)
