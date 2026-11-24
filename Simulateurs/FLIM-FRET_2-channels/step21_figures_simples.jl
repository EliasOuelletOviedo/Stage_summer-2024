# Étape 21 — versions minimalistes des figures des étapes 19 (validation du simulateur) et
# 20 (biais sous autofluorescence), pour un lecteur qui ne connaît pas le projet.
#
# Aucun calcul : la figure relit les CSV des étapes 19 et 20 (les chiffres sont donc ceux des
# figures détaillées). Une idée par panneau, titres qui énoncent le résultat, pas de jargon
# interne ; les détails vont dans la légende (legende_simple.md, à côté de chaque figure).
#
#   julia --project=. step21_figures_simples.jl
#
# Sorties : sorties/validation_simulateur_2026-10-16_SCLM_850nm/validation_simulateur_simple.{png,pdf}
#           sorties/autofluorescence_2026-10-16_SCLM_850nm/autofluorescence_simple.{png,pdf}

include(joinpath(@__DIR__, "src", "Figures.jl"))
using CairoMakie, CSV, DataFrames, Printf, Statistics
const F = Figures
set_theme!(F.theme_sim())
virg1(x) = F.virg(@sprintf("%.1f", x))
signe(x) = (abs(x) < 0.05 ? "0" : (x > 0 ? "+" : "−") * virg1(abs(x)))
virgs(v) = [F.virg(@sprintf("%g", x)) for x in v]
milliers(n) = replace(string(round(Int, n)), r"(?<=\d)(?=(\d{3})+$)" => " ")
ORANGE = F.COULEUR["2"]

# =============================================================================
# Figure A.1 — validation du simulateur
# =============================================================================
dA1 = joinpath(@__DIR__, "sorties", "validation_simulateur_2026-10-16_SCLM_850nm")
ar = CSV.read(joinpath(dA1, "aller_retour.csv"), DataFrame)
rs = CSV.read(joinpath(dA1, "reel_simule.csv"), DataFrame)
N = round(Int, sum(rs.reel))
# origine du temps : mi-hauteur du front montant du modèle (≈ arrivée de l'impulsion laser)
im = argmax(rs.modele)
i50 = findfirst(>=(0.5rs.modele[im]), rs.modele)
t50 = rs.t_ns[i50-1] + (0.5rs.modele[im] - rs.modele[i50-1]) / (rs.modele[i50] - rs.modele[i50-1]) * (rs.t_ns[i50] - rs.t_ns[i50-1])
t = rs.t_ns .- t50
χ2r(r) = sum(abs2, r) / (length(r) - 1)

fig = Figure(size = (1400, 620), fontsize = 17)
axa = Axis(fig[1, 1]; title = "a   Le simulateur retrouve la fraction FRET injectée", titlesize = 19,
    subtitle = "moyenne ± écart-type de 500 simulations par point", subtitlesize = 15,
    xlabel = "fraction FRET injectée", ylabel = "fraction FRET retrouvée", aspect = 1,
    xticks = 0:0.25:1, yticks = 0:0.25:1, xtickformat = virgs, ytickformat = virgs)
lines!(axa, [0, 1], [0, 1]; color = F.MUET, linestyle = :dash, linewidth = 1.5)
coulN = Dict(1_000 => "#a9cdf2", 10_000 => "#2a78d6", 100_000 => "#0d3d7a")
dx = Dict(1_000 => -0.015, 10_000 => 0.0, 100_000 => 0.015)
for n in (1_000, 10_000, 100_000)
    l = ar[ar.N .== n, :]
    errorbars!(axa, l.alpha_injecte .+ dx[n], l.moyenne, l.ecart_type; color = coulN[n], linewidth = 2.2, whiskerwidth = 0)
    scatter!(axa, l.alpha_injecte .+ dx[n], l.moyenne; color = coulN[n], markersize = 11, strokewidth = 0)
end
limits!(axa, -0.03, 1.03, -0.08, 1.08)
Legend(fig[1, 1], [LineElement(color = F.MUET, linestyle = :dash), [MarkerElement(color = coulN[n], marker = :circle, markersize = 11) for n in (1_000, 10_000, 100_000)]...],
    ["retrouvée = injectée", "1 000 photons", "10 000 photons", "100 000 photons"];
    tellwidth = false, tellheight = false, halign = :left, valign = :top, margin = (14, 0, 0, 52), labelsize = 15, rowgap = 2)

gb = fig[1, 2] = GridLayout()
axb = Axis(gb[1, 1]; title = "b   Une décroissance simulée se superpose à la mesure", titlesize = 19,
    subtitle = "même échantillon (20 mM de Cl⁻), mêmes durées de vie, même instrument, $(milliers(N)) photons", subtitlesize = 15,
    yscale = log10, ylabel = "photons par canal")
lines!(axb, t, rs.simule; color = ORANGE, linewidth = 1.8)
scatter!(axb, t, rs.reel; color = F.ENCRE, markersize = 5, strokewidth = 0)
yt = [500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000]
axb.yticks = (yt, milliers.(yt))
Legend(gb[1, 1], [MarkerElement(color = F.ENCRE, marker = :circle, markersize = 8), LineElement(color = ORANGE, linewidth = 2.5),
        PolyElement(color = F.ZONE, strokecolor = F.AXE, strokewidth = 0.5)],
    ["mesure", "simulation", "bruit attendu (écart ±2)"]; tellwidth = false, tellheight = false, halign = :right, valign = :top,
    margin = (0, 12, 0, 12), labelsize = 16)
axr = Axis(gb[2, 1]; xlabel = "temps après l'impulsion laser (ns)", ylabel = "écart au modèle\n(unités de bruit)",
    xtickformat = virgs, yticks = [-6, -3, 0, 3, 6], ylabelsize = 15)
band!(axr, [t[1], t[end]], [-2.0, -2.0], [2.0, 2.0]; color = F.ZONE)
hlines!(axr, [0.0]; color = F.MUET, linewidth = 1)
scatter!(axr, t, rs.pearson_simule; color = (ORANGE, 0.85), markersize = 5, strokewidth = 0)
scatter!(axr, t, rs.pearson_reel; color = F.ENCRE, markersize = 5, strokewidth = 0)
text!(axr, t[end], 5.2; text = F.virg(@sprintf("χ²ᵣ (1 = bruit seul) : mesure %.1f · simulation %.1f", χ2r(rs.pearson_reel), χ2r(rs.pearson_simule))),
    align = (:right, :center), fontsize = 14, color = F.ENCRE2)
ylims!(axr, -7, 7)
linkxaxes!(axb, axr)
hidexdecorations!(axb; grid = false)
rowsize!(gb, 2, Relative(0.3))
rowgap!(gb, 8)
colsize!(fig.layout, 1, Relative(0.42))
colgap!(fig.layout, 40)
save(joinpath(dA1, "validation_simulateur_simple.png"), fig; px_per_unit = 2)
save(joinpath(dA1, "validation_simulateur_simple.pdf"), fig)
println("-> ", joinpath(dA1, "validation_simulateur_simple.png"))

# =============================================================================
# Figure A.2 — biais sous autofluorescence
# =============================================================================
dA2 = joinpath(@__DIR__, "sorties", "autofluorescence_2026-10-16_SCLM_850nm")
bi = CSV.read(joinpath(dA2, "biais_selon_cl.csv"), DataFrame)
pr = CSV.read(joinpath(dA2, "precision.csv"), DataFrame)
c_aff = 10.0
biais(sc, lect; connue = true) = only(bi[(bi.scenario .== sc) .& (bi.lecture .== lect) .&
    (bi.formes .== (connue ? "celles du modèle" : "différentes du modèle")) .& (bi.cl_vraie .== c_aff), :biais])
sigma(lect) = only(pr[(pr.lecture .== lect) .& (pr.cl_vraie .== c_aff), :ecart_type_robuste])
N_prec = only(unique(pr.N_bleu))
# (lecture, forme de l'autofluorescence connue du modèle ?, étiquette)
LIGNES = [("ratio", true, "Rapport\njaune / bleu"), ("mono", true, "Durée de vie\ndu donneur"),
          ("deux_auto", true, "Durées de vie, 2 canaux\n+ autofluorescence\n(forme connue)"),
          ("deux_auto", false, "Durées de vie, 2 canaux\n+ autofluorescence\n(forme mal connue)")]
ys = [4.0, 3.0, 2.0, 1.0]
BLEU, AMBRE = F.COULEUR["1"], F.COULEUR["4"]
h = 0.30                                   # épaisseur d'une barre

fig2 = Figure(size = (1400, 700), fontsize = 17)
ax1 = Axis(fig2[1, 1]; title = "a   Erreur sur la [Cl⁻] causée par l'autofluorescence", titlesize = 19,
    subtitle = "autofluorescence = 5 % des photons du canal · [Cl⁻] vraie : 10 mM", subtitlesize = 15,
    xlabel = "erreur sur la [Cl⁻] lue (mM)", yticks = (ys, last.(LIGNES)), yticklabelsize = 16,
    ygridvisible = false, xtickformat = virgs, xticks = -3:1:4)
vlines!(ax1, [0.0]; color = F.ENCRE2, linewidth = 1.2)
function barre!(ax, y, v; couleur, plein = true)
    plein ? poly!(ax, Rect(min(0, v), y - h / 2, abs(v), h); color = couleur) :
            poly!(ax, Rect(min(0, v), y - h / 2, abs(v), h); color = (couleur, 0.12), strokecolor = couleur, strokewidth = 2.2)
end
function etiquette!(ax, y, v; suffixe = "")
    text!(ax, v, y; text = signe(v) * suffixe, align = (v >= 0 ? :left : :right, :center), offset = (v >= 0 ? 7 : -7, 0),
        fontsize = 15, color = F.ENCRE)
end
for ((cle, connue, _), y) in zip(LIGNES, ys)
    vN, vF = biais("NADH", cle; connue), biais("FAD", cle; connue)
    barre!(ax1, y + 0.17, vN; couleur = BLEU, plein = connue)
    barre!(ax1, y - 0.17, vF; couleur = AMBRE, plein = connue)
    etiquette!(ax1, y + 0.17, vN)
    etiquette!(ax1, y - 0.17, vF)
end
limits!(ax1, -3.6, 4.2, 0.45, 4.55)
Legend(fig2[2, 1], [PolyElement(color = BLEU), PolyElement(color = AMBRE)], ["NADH (canal bleu)", "flavines (canal jaune)"];
    orientation = :horizontal, labelsize = 15, tellwidth = false)

ax2 = Axis(fig2[1, 2]; title = "b   Le prix : bruit de la mesure", titlesize = 19,
    subtitle = "$(milliers(N_prec)) photons bleus · [Cl⁻] vraie : 10 mM", subtitlesize = 15,
    xlabel = "écart-type de la [Cl⁻] lue (mM)", yticks = (ys, fill("", length(ys))), ygridvisible = false,
    xtickformat = virgs, xticks = 0:1:4)
for ((cle, _, _), y) in zip(LIGNES, ys)
    σ = sigma(cle)
    poly!(ax2, Rect(0, y - 0.22, σ, 0.44); color = F.MUET)
    text!(ax2, σ, y; text = virg1(σ), align = (:left, :center), offset = (7, 0), fontsize = 15, color = F.ENCRE)
end
limits!(ax2, 0, 4.4, 0.45, 4.55)
linkyaxes!(ax1, ax2)
colsize!(fig2.layout, 1, Relative(0.64))
colgap!(fig2.layout, 30)
save(joinpath(dA2, "autofluorescence_simple.png"), fig2; px_per_unit = 2)
save(joinpath(dA2, "autofluorescence_simple.pdf"), fig2)
println("-> ", joinpath(dA2, "autofluorescence_simple.png"))
