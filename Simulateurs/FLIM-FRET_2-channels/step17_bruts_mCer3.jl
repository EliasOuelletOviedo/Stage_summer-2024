# Étape 17 — histogrammes bruts sommés et IRF de la titration SCLM_mCer3 du 2026-10-16.
#
# Format propre à ce dossier : fichiers à plat `<[Cl⁻]>_cNNN.sdt` (virgule décimale :
# `2,5`), 100 acquisitions de 1 s par concentration, une IRF `irf.sdt` (100 s) dans le
# dossier même. Deux voies, mais rangées autrement que dans SCLM : chaque fichier a DEUX
# blocs de données d'une voie chacun, un par carte SPC-150N (bloc 1 = carte 3N0317 =
# voie 1, bloc 2 = carte 3N0318 = voie 2), là où SCLM a un seul bloc à plusieurs voies
# (carte SPC-QC-104). Pas du TAC 0,0489 ns : les 256 canaux couvrent une période laser
# (12,5 ns). L'échantillon `unknown` est écarté.
#
# Sorties (sorties/declins_purs_2026-10-16_SCLM_mCer3/) :
#   irf_bruts_mCer3.csv   une ligne par canal TCSPC, pour chaque voie : IRF brute, IRF
#                         utilisée (fond retiré, somme 1), somme brute de chaque concentration
#   irf_bruts_mCer3.png   les deux voies côte à côte, sommes normalisées à 1 au maximum
#   irf_bruts_mCer3_comptes.png   la même chose sans normalisation (comptes bruts)
#
#   julia --project=. step17_bruts_mCer3.jl

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, CSV, DataFrames
const SF = S.SdtFile

"""Un fichier à un bloc par carte : (voies × 256 canaux, pas_ns, durée_s, arrêt sur temps, n° de série des cartes)."""
function lire_blocs(chemin)
    s = SF.read_sdt(read(chemin), basename(chemin))
    mi = s.measure_info
    pas = [m.tac_r / (m.tac_g * m.adc_re) * 1e9 for m in mi]
    all(≈(pas[1]), pas) || error("$(chemin) : pas différents entre blocs $(pas)")
    d = reduce(vcat, [reshape(Float64.(vec(b)), 1, :) for b in s.data])
    return d, pas[1], Float64(mi[1].col_t), all(m -> m.stopt != 0, mi), [m.mod_ser_no for m in mi]
end

donnees = joinpath(@__DIR__, "Data", "2026-10-16", "SCLM_mCer3")
dossier = joinpath(@__DIR__, "sorties", "declins_purs_2026-10-16_SCLM_mCer3")
mkpath(dossier)

# -----------------------------------------------------------------------------
# 1. Sommes par concentration, deux voies
# -----------------------------------------------------------------------------
fichiers = filter(f -> endswith(f, ".sdt") && f != "irf.sdt", readdir(donnees))
concentrations = sort(unique(filter(!isnothing, [tryparse(Float64, replace(split(f, "_")[1], "," => ".")) for f in fichiers])))
nom_fichier(c) = replace(@sprintf("%g", c), "." => ",")
etiquette(c) = replace(@sprintf("%g", c), "." => "p")
irf_d, pas, duree_irf, _, cartes = lire_blocs(joinpath(donnees, "irf.sdt"))
nv, n = size(irf_d)
nv == 2 || @warn "IRF à $(nv) blocs (2 attendus)"
sommes = Dict{Float64,Matrix{Float64}}()
nacq = Dict{Float64,Int}()
for c in concentrations
    fs = filter(f -> startswith(f, nom_fichier(c) * "_"), fichiers)
    h = zeros(nv, n)
    for f in fs
        d, p, dur, arret, cf = lire_blocs(joinpath(donnees, f))
        size(d) == (nv, n) || error("$(f) : $(size(d)) au lieu de $((nv, n))")
        cf == cartes || error("$(f) : cartes $(cf), IRF $(cartes)")
        abs(p - pas) < 1e-9 || error("$(f) : pas du TAC $(p) ns, IRF à $(pas) ns")
        dur == 1.0 || @warn "$(f) : durée $(dur) s"
        arret || @warn "$(f) : acquisition sans arrêt sur temps"
        h .+= d
    end
    sommes[c] = h
    nacq[c] = length(fs)
end

tt = (0:n) .* pas
tab = DataFrame(canal = 0:n-1, t_debut_ns = tt[1:n], t_fin_ns = tt[2:n+1])
fenetres = UnitRange{Int}[]
@printf("SCLM_mCer3 : %d concentrations, %d voies (cartes %s), pas %.4f ns\n", length(concentrations), nv, join(cartes, ", "), pas)
for v in 1:nv
    # fenêtre : canaux que le TAC ne coupe pas (non nuls dans la somme de tout, IRF comprise)
    tot = sum(sommes[c][v, :] for c in concentrations) .+ irf_d[v, :]
    W = findfirst(>(0), tot):findlast(>(0), tot)
    push!(fenetres, W)
    # IRF utilisée : comme irf_banc (fond = moyenne des canaux valides avant le pic, négatifs à 0, somme 1)
    ib = irf_d[v, :]
    k = argmax(ib)
    fo = mean(ib[first(W):k-4])
    iu = max.(ib .- fo, 0.0)
    iu[setdiff(1:n, W)] .= 0.0
    iu ./= sum(iu)
    tab[!, "dans_fenetre_voie$(v)"] = [i in W for i in 1:n]
    tab[!, "irf_voie$(v)_brute"] = ib
    tab[!, "irf_voie$(v)_utilisee"] = iu
    @printf("  voie %d (%s) : fenêtre canaux %d–%d (base 0) = %.2f–%.2f ns ; IRF pic au canal %d, fond %.1f comptes/canal retiré\n",
        v, cartes[v], first(W) - 1, last(W) - 1, tt[first(W)], tt[last(W)+1], k - 1, fo)
end
for v in 1:nv, c in concentrations
    tab[!, "voie$(v)_brut_$(etiquette(c))mM"] = sommes[c][v, :]
end
CSV.write(joinpath(dossier, "irf_bruts_mCer3.csv"), tab)

# Échantillon inconnu : somme brute de ses acquisitions seulement (deux voies, sans IRF)
fu = filter(f -> startswith(f, "unknown_"), fichiers)
if !isempty(fu)
    hu = zeros(nv, n)
    for f in fu
        d, p, dur, arret, cf = lire_blocs(joinpath(donnees, f))
        (size(d) == (nv, n) && cf == cartes && abs(p - pas) < 1e-9) || error("$(f) : format différent de l'IRF")
        dur == 1.0 || @warn "$(f) : durée $(dur) s"
        hu .+= d
    end
    du = DataFrame(canal = 0:n-1, t_debut_ns = tt[1:n], t_fin_ns = tt[2:n+1])
    for v in 1:nv
        du[!, "dans_fenetre_voie$(v)"] = tab[!, "dans_fenetre_voie$(v)"]
        du[!, "voie$(v)_brut_inconnu"] = hu[v, :]
    end
    CSV.write(joinpath(dossier, "inconnu_bruts_mCer3.csv"), du)
    @printf("  inconnu : %d acquisitions ; photons dans la fenêtre : %s\n", length(fu),
        join([@sprintf("voie %d %.3g", v, sum(hu[v, fenetres[v]])) for v in 1:nv], ", "))
end
for c in concentrations
    @printf("  %6g mM : %3d acquisitions ; photons dans la fenêtre : %s\n", c, nacq[c],
        join([@sprintf("voie %d %.3g", v, sum(sommes[c][v, fenetres[v]])) for v in 1:nv], ", "))
end

# -----------------------------------------------------------------------------
# 2. Figure : les deux voies côte à côte, sommes normalisées à 1 au maximum
# -----------------------------------------------------------------------------
include(joinpath(@__DIR__, "src", "Figures.jl"))
using CairoMakie
const F = Figures
set_theme!(F.theme_sim())
couleurs = F.RAMPE[round.(Int, range(1, length(F.RAMPE), length(concentrations)))]
tc = (tt[1:n] .+ tt[2:n+1]) ./ 2
fig = Figure(size = (1150, 540))
Label(fig[0, 1:nv], "SCLM_mCer3 · histogrammes bruts (somme des acquisitions), normalisés à 1 au maximum";
    halign = :left, fontsize = 15, color = F.ENCRE)
# échelle des y partagée : du plus petit compte normalisé des voies jusqu'à 1
ymin = minimum(minimum(sommes[c][v, fenetres[v]]) / maximum(sommes[c][v, fenetres[v]]) for v in 1:nv for c in concentrations)
g = filter(x -> ymin <= x <= 1, [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
for v in 1:nv
    W = fenetres[v]
    ax = Axis(fig[1, v]; yscale = log10, xlabel = "temps du TAC (ns)", ylabel = v == 1 ? "comptes normalisés (max = 1)" : "",
        title = "Voie $(v) (carte $(cartes[v]))", subtitle = F.virg(@sprintf("fenêtre utile %.2f–%.2f ns", tt[first(W)], tt[last(W)+1])))
    for (i, c) in enumerate(concentrations)
        h = sommes[c][v, :]
        ok = (h .> 0) .& tab[!, "dans_fenetre_voie$(v)"]
        lines!(ax, tc[ok], h[ok] ./ maximum(h[ok]); color = couleurs[i], linewidth = 1.5, label = F.virg(@sprintf("%g mM", c)))
    end
    xlims!(ax, 0, 12.5)
    ylims!(ax, ymin, 1.0)
    ax.yticks = (g, [F.virg(@sprintf("%g", x)) for x in g])
    v == nv && Legend(fig[1, nv+1], ax, "[Cl⁻]"; framevisible = false)
end
save(joinpath(dossier, "irf_bruts_mCer3.png"), fig; px_per_unit = 2)

# Même figure sans normalisation : comptes bruts, échelle des y partagée du min au max
fig3 = Figure(size = (1150, 540))
Label(fig3[0, 1:nv], "SCLM_mCer3 · histogrammes bruts (somme des acquisitions), comptes non normalisés";
    halign = :left, fontsize = 15, color = F.ENCRE)
vals = [sommes[c][v, i] for v in 1:nv for c in concentrations for i in fenetres[v] if sommes[c][v, i] > 0]
ymin3, ymax3 = extrema(vals)
g3 = filter(x -> ymin3 <= x <= ymax3, [m * 10.0^e for e in 0:9 for m in (1, 2, 5)])
for v in 1:nv
    W = fenetres[v]
    ax = Axis(fig3[1, v]; yscale = log10, xlabel = "temps du TAC (ns)", ylabel = v == 1 ? "photons par canal TCSPC" : "",
        title = "Voie $(v) (carte $(cartes[v]))", subtitle = F.virg(@sprintf("fenêtre utile %.2f–%.2f ns", tt[first(W)], tt[last(W)+1])))
    for (i, c) in enumerate(concentrations)
        h = sommes[c][v, :]
        ok = (h .> 0) .& tab[!, "dans_fenetre_voie$(v)"]
        lines!(ax, tc[ok], h[ok]; color = couleurs[i], linewidth = 1.5, label = F.virg(@sprintf("%g mM", c)))
    end
    xlims!(ax, 0, 12.5)
    ylims!(ax, ymin3, ymax3)
    ax.yticks = (g3, [x >= 1e4 ? F.virg(@sprintf("%g", x)) : F.virg(@sprintf("%g", x)) for x in g3])
    v == nv && Legend(fig3[1, nv+1], ax, "[Cl⁻]"; framevisible = false)
end
save(joinpath(dossier, "irf_bruts_mCer3_comptes.png"), fig3; px_per_unit = 2)
println("Sorties : $(dossier)")
