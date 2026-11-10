# Étape 18 — histogrammes bruts sommés et IRF de la titration KIRIN du 2026-10-16.
#
# Même format que SCLM : un dossier par concentration, `<λ>_cNN.sdt` (10 acquisitions de
# 1 s), un seul bloc à deux voies utiles (1 bleu, 2 jaune), IRF `IRF_<λ>_16X.sdt` dans le
# dossier parent. Lecture par le chargeur de l'étape 15 (Charger.jl). `unknown` est écarté.
#
# Sorties (sorties/declins_purs_2026-10-16_KIRIN/) :
#   irf_bruts_KIRIN_<λ>nm.csv   une ligne par canal TCSPC : IRF brute et IRF utilisée (fond
#                               retiré, somme 1) de chaque voie, somme brute de chaque concentration
#   irf_bruts_KIRIN_<λ>nm.png   les deux voies côte à côte, sommes normalisées à 1 au maximum
#
#   julia --project=. step18_bruts_KIRIN.jl [--lambda 800]

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
include(joinpath(@__DIR__, "src", "deux_canaux", "DeuxCanaux.jl"))
using .SimFLIMFRET, .DeuxCanaux
const S = SimFLIMFRET
const D = DeuxCanaux
using Printf, Statistics, CSV, DataFrames

args = copy(ARGS)
valeur(nom, defaut) = (i = findfirst(==(nom), args); i === nothing ? defaut : args[i+1])
λ = parse(Int, valeur("--lambda", "800"))
donnees = joinpath(@__DIR__, "Data", "2026-10-16", "KIRIN")
dossier = joinpath(@__DIR__, "sorties", "declins_purs_2026-10-16_KIRIN")
mkpath(dossier)

# -----------------------------------------------------------------------------
# 1. Sommes par concentration et IRF
# -----------------------------------------------------------------------------
cellules, _ = D.charger(donnees; λs = (λ,))
etiquette(c) = replace(@sprintf("%g", c), "." => "p")
chemin_irf = joinpath(dirname(donnees), "IRF_$(λ)_16X.sdt")
irf_brute, pas, duree_irf, _ = S.lire_sdt_voies(chemin_irf)
n = size(irf_brute, 2)
tt = (0:n) .* pas
voies = ((D.BLEU, "bleu", "Canal bleu (1)"), (D.JAUNE, "jaune", "Canal jaune (2)"))
tab = DataFrame(canal = 0:n-1, t_debut_ns = tt[1:n], t_fin_ns = tt[2:n+1], dans_fenetre = [k in D.FENETRE for k in 1:n])
for (v, nom, _) in voies
    tab[!, "irf_$(nom)_brute"] = irf_brute[v, :]
    tab[!, "irf_$(nom)_utilisee"] = cellules[1].irf[v, :]
end
for (v, nom, _) in voies, c in cellules
    tab[!, "$(nom)_brut_$(etiquette(c.cl))mM"] = D.somme(c, v)
end
CSV.write(joinpath(dossier, "irf_bruts_KIRIN_$(λ)nm.csv"), tab)
@printf("KIRIN, %d nm : %d concentrations, IRF %s (%.0f s), pas %.4f ns\n", λ, length(cellules), basename(chemin_irf), duree_irf, pas)
for c in cellules
    @printf("  %6g mM : %2d acquisitions ; photons dans la fenêtre : bleu %.3g, jaune %.3g\n",
        c.cl, D.n_acq(c), D.photons(D.somme(c, D.BLEU)), D.photons(D.somme(c, D.JAUNE)))
end

# -----------------------------------------------------------------------------
# 2. Figure : les deux voies côte à côte, sommes normalisées à 1 au maximum
# -----------------------------------------------------------------------------
include(joinpath(@__DIR__, "src", "Figures.jl"))
using CairoMakie
const F = Figures
set_theme!(F.theme_sim())
couleurs = F.RAMPE[round.(Int, range(1, length(F.RAMPE), length(cellules)))]
tc = (tt[1:n] .+ tt[2:n+1]) ./ 2
W = collect(D.FENETRE)
fig = Figure(size = (1150, 540))
Label(fig[0, 1:2], "KIRIN · $(λ) nm · histogrammes bruts (somme des 10 acquisitions), normalisés à 1 au maximum";
    halign = :left, fontsize = 15, color = F.ENCRE)
# échelle des y partagée : du plus petit compte normalisé des deux voies jusqu'à 1
ymin = minimum(minimum(h[W]) / maximum(h[W]) for (v, _, _) in voies for c in cellules for h in (D.somme(c, v),))
g = filter(x -> ymin <= x <= 1, [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
for (j, (v, nom, titre)) in enumerate(voies)
    ax = Axis(fig[1, j]; yscale = log10, xlabel = "temps du TAC (ns)", ylabel = j == 1 ? "comptes normalisés (max = 1)" : "",
        title = titre, subtitle = F.virg(@sprintf("fenêtre utile %.2f–%.2f ns", tt[first(W)], tt[last(W)+1])))
    for (i, c) in enumerate(cellules)
        h = D.somme(c, v)
        ok = (h .> 0) .& tab.dans_fenetre
        lines!(ax, tc[ok], h[ok] ./ maximum(h[ok]); color = couleurs[i], linewidth = 1.5, label = F.virg(@sprintf("%g mM", c.cl)))
    end
    xlims!(ax, 0, 12.5)
    ylims!(ax, ymin, 1.0)
    ax.yticks = (g, [F.virg(@sprintf("%g", x)) for x in g])
    j == 2 && Legend(fig[1, 3], ax, "[Cl⁻]"; framevisible = false)
end
save(joinpath(dossier, "irf_bruts_KIRIN_$(λ)nm.png"), fig; px_per_unit = 2)
println("Sorties : $(dossier)")
