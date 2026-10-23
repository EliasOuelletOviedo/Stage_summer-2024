# Étape 10 — le « facteur 8 » sur les données réelles.
#
# La note de changements le place en tête : remplacer mono + Hill par un
# bi-exponentiel à durées de vie figées est à la fois un gain de précision et un
# test du modèle. Ici, sur les 27 histogrammes réels amincis par tirage binomial à
# N ∈ {10³, 3·10³, 10⁴, 3·10⁴} photons (200 tirages par histogramme et par N) :
#   2H  le pipeline actuel : mono-exponentiel (code de l'utilisateur) + Hill ;
#   4R  bi-exponentiel à τ_D0, τ_DA, t₀ figés (ajustement global), fond libre,
#       calibré par l'isotherme sur les histogrammes complets.
# Trois chiffres par cellule : le σ mesuré sur les données réelles, le σ du jumeau
# (mêmes opérations sur des tirages du modèle ajusté, fond compris) et la
# prédiction du simulateur (étape 4, normalisation C : photons bleus fixés).
#
#   julia --project=. -t 8 step10_facteur8.jl
#
# Sorties : sorties/reel/facteur8.csv (résumé), sorties/reel/amincis_tirages.csv
# (toutes les lectures, ignoré par git), sorties/figures/fig10_facteur8.png

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, CSV, DataFrames, TOML

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
prot = cfg["capteur"]["proteine"]
dossier = S.dossier_sorties("reel")
mkpath(dossier)
W = m.inst.fenetre

hs = S.lire_titration(prot)
C = [h.C for h in hs]
g = S.ajuster_global(hs, m.inst; τ_D0_init = m.p.τ_D0, τ_DA_init = m.p.τ_DA)
@printf("ajustement global : τ_D0 = %.4f, τ_DA = %.4f, t₀ = %.4f ns, χ²ᵣ = %.3f\n", g.τ_D0, g.τ_DA, g.t₀, g.chi2_global)
ctxs = S.creer_contextes(m, irf_mat)
ctx = ctxs[1]
e4 = S.BiExpFigeFond(g.τ_D0, g.τ_DA, g.t₀, m.inst)

"""Calibrations empiriques sur des histogrammes complets (réels ou modèle)."""
function calibrations(histos)
    f = [S.observable(e4, h, h, ctx)[1] for h in histos]
    τ = [S.ajuster_mono_utilisateur(ctx.utilisateur, h)[1][1] for h in histos]
    iso = S.ajuster_isotherme(C, f)
    hill = S.ajuster_hill([isinf(c) ? 1e5 : c for c in C], τ)     # « F » à 100 000 mM, comme analyse_SCLM_2.jl
    return iso, hill, f, τ
end
plein_reel = [h.n for h in hs]
plein_jumeau = map(1:length(hs)) do h
    μ = zeros(m.inst.n_bins)
    μ[W] .= g.μ[h]
    μ
end
iso_r, hill_r, f_r, τ_r = calibrations(plein_reel)
iso_j, hill_j, f_j, τ_j = calibrations(plein_jumeau)
@printf("calibration 4R : isotherme K = %.1f mM, n = %.2f (réel) ; %.1f mM, n = %.2f (jumeau)\n", iso_r.K_D, iso_r.n, iso_j.K_D, iso_j.n)
@printf("calibration 2H : Hill K = %.1f mM, n = %.2f (réel) ; %.1f mM, n = %.2f (jumeau)\n", hill_r[3], hill_r[4], hill_j[3], hill_j[4])

Ns = (1e3, 3e3, 1e4, 3e4)
R = 200
t = @elapsed begin
    dr = S.lectures_amincies(:reel, hs, nothing, Ns, R, m.inst, ctxs, e4, iso_r, hill_r; C_max, graine0 = 10_000_000)
    dj = S.lectures_amincies(:jumeau, hs, [g.μ[h] for h in eachindex(hs)], Ns, R, m.inst, ctxs, e4, iso_j, hill_j; C_max, graine0 = 20_000_000)
end
@printf("%d lectures en %.0f s\n", 2 * nrow(dr), t)
tirages = vcat(dr, dj)
CSV.write(joinpath(dossier, "amincis_tirages.csv"), tirages)
res = S.resumer_amincis(tirages)

# prédiction du simulateur (étape 4, normalisation C) aux concentrations communes
sim = CSV.read(S.dossier_sorties("resultats", "resultats_principaux.csv"), DataFrame; types = Dict(:estimateur => String, :normalisation => String))
cellules_sim = unique(sim[(sim.normalisation.=="C").&in.(sim.N_photons, Ref(collect(Ns))).&in.(sim.C_vrai, Ref(unique(C))), [:C_vrai, :N_photons]])
for r in eachrow(cellules_sim)
    push!(res, (source = "simulateur", C_vrai = r.C_vrai, N = r.N_photons); cols = :union)
end
sim_v(C_, N, k, col) = (d = sim[(sim.normalisation.=="C").&(sim.estimateur.==k).&(sim.N_photons.==N).&(sim.C_vrai.==C_), :]; nrow(d) == 1 ? d[1, col] : NaN)
for r in eachrow(res)
    r.source == "simulateur" || continue
    r.sigma_4R = sim_v(r.C_vrai, r.N, "4", :sigma_mM)
    r.sigma_2H = sim_v(r.C_vrai, r.N, "2H", :sigma_mM)
    r.sigma_robuste_4R = sim_v(r.C_vrai, r.N, "4", :sigma_robuste_mM)
    r.sigma_robuste_2H = sim_v(r.C_vrai, r.N, "2H", :sigma_robuste_mM)
    r.hors_plage_2H = sim_v(r.C_vrai, r.N, "2H", :taux_hors_plage)
    r.hors_plage_4R = sim_v(r.C_vrai, r.N, "4", :taux_hors_plage)
    r.biais_4R = sim_v(r.C_vrai, r.N, "4", :biais_mM)
    r.biais_2H = sim_v(r.C_vrai, r.N, "2H", :biais_mM)
    r.facteur = r.sigma_2H / r.sigma_4R
    r.facteur_robuste = r.sigma_robuste_2H / r.sigma_robuste_4R
end
sort!(res, [:C_vrai, :N, :source])
CSV.write(joinpath(dossier, "facteur8.csv"), res)

println("\nσ[Cl⁻] (mM, bruit de photons ; σ robuste entre parenthèses) et facteur 2H/4R :")
for c in (5.0, 20.0, 40.0), N in Ns
    @printf("C = %4.0f mM, N = %6.0f :", c, N)
    for src in ("reel", "jumeau", "simulateur")
        d = res[(res.source.==src).&(res.C_vrai.==c).&(res.N.==N), :]
        nrow(d) == 1 || continue
        @printf("  %-10s 4R %6.2f (%6.2f)  2H %6.2f (%6.2f) hp %3.0f %%  facteur %5.2f (%5.2f)", src, d.sigma_4R[1], d.sigma_robuste_4R[1],
            d.sigma_2H[1], d.sigma_robuste_2H[1], 100 * d.hors_plage_2H[1], d.facteur[1], d.facteur_robuste[1])
    end
    println()
end
resume = Dict{String,Any}("tau_D0" => g.τ_D0, "tau_DA" => g.τ_DA, "t0" => g.t₀, "R" => R,
    "calibration_4R_reel" => Dict("x0" => iso_r.x₀, "x_sat" => iso_r.x_sat, "K_D" => iso_r.K_D, "n" => iso_r.n),
    "calibration_2H_reel" => Dict("bas" => hill_r[1], "haut" => hill_r[2], "K" => hill_r[3], "n" => hill_r[4]),
    "calibration_4R_jumeau" => Dict("x0" => iso_j.x₀, "x_sat" => iso_j.x_sat, "K_D" => iso_j.K_D, "n" => iso_j.n),
    "calibration_2H_jumeau" => Dict("bas" => hill_j[1], "haut" => hill_j[2], "K" => hill_j[3], "n" => hill_j[4]))
open(io -> TOML.print(io, resume; sorted = true), joinpath(dossier, "facteur8_calibrations.toml"), "w")

include(joinpath(@__DIR__, "src", "Figures.jl"))
using .Figures, CairoMakie
set_theme!(Figures.theme_sim())
if isdefined(Figures, :figure_facteur8)
    Figures.figure_facteur8(res, S.dossier_sorties("figures", "fig10_facteur8.png"))
end
