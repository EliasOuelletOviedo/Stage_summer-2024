# Étape 6 — figures, à partir des fichiers plats des étapes 1, 3, 4, 5 et, s'ils
# existent, 7 et 9 à 13 (révision du 2026-10-15).
#
#   julia --project=. step6_figures.jl
#
# Sorties : sorties/figures/*.png et sorties/resultats/resume.md (tableaux chiffrés)

include(joinpath(@__DIR__, "src", "Figures.jl"))
include(joinpath(@__DIR__, "src", "Resume.jl"))
using .Figures, .Resume
const F = Figures
using CairoMakie, CSV, DataFrames, TOML

set_theme!(F.theme_sim())
racine = joinpath(@__DIR__, "sorties")
fig(n) = joinpath(racine, "figures", n)
mkpath(joinpath(racine, "figures"))
lire(f) = CSV.read(f, DataFrame; types = Dict(:estimateur => String, :normalisation => String, :mode_calibration => String))

cfg = TOML.parsefile(joinpath(@__DIR__, "parametres.toml"))
prot = cfg["capteur"]["proteine"]
tit = TOML.parsefile(joinpath(racine, "titration", "parametres_$(prot).toml"))

for p in ("mCer3_Topaz", "SCLM", "M2-SCLM")
    isfile(joinpath(racine, "titration", "titration_$(p).csv")) || continue
    F.figure_titration(p, joinpath(racine, "titration"), fig("titration_$(p).png"))
end

kapp = CSV.read(joinpath(racine, "calibration", "K_apparents.csv"), DataFrame; types = Dict(:estimateur => String))
courbes = CSV.read(joinpath(racine, "calibration", "courbes_calibration.csv"), DataFrame)
fmK = joinpath(racine, "incertitude", "mondes_K_apparents.csv")
mondes_K = isfile(fmK) ? CSV.read(fmK, DataFrame; types = Dict(:estimateur => String, :monde => String)) : nothing
F.figure5(kapp, courbes, fig("fig5_K_apparents.png"); K_vrai = tit["K_D"], mondes_K)

principal = joinpath(racine, "resultats", "resultats_principaux.csv")
fmg = joinpath(racine, "incertitude", "mondes_grille.csv")
mondes = isfile(fmg) ? lire(fmg) : nothing
if isfile(principal)
    df = lire(principal)
    F.figure1(df, fig("fig1_sigma_selon_C.png"); mondes)
    F.figure2(df, fig("fig2_sigma_selon_N.png"); mondes)
    F.figure3(df, fig("fig3_biais_selon_N.png"))
    F.figure6(df, fig("fig6_efficacite.png"))
    F.carte(df, fig("carte_meilleur_estimateur.png"))
end

nuis = joinpath(racine, "resultats", "nuisances.csv")
if isfile(nuis)
    nu = lire(nuis)
    eta = TOML.parsefile(joinpath(racine, "spectres", "eta.toml"))
    F.figure4(nu, fig("fig4_biais_selon_fuite.png"); η_nominal = eta["eta_YD"])
    F.figure4b(nu, fig("fig4b_autres_nuisances.png"))
    if isfile(principal) && any(nu.parametre_nuisance .== "tau_DA_capteur_ns")
        F.figure8(nu, df, fig("fig8_capteurs_hypothetiques.png"); τ_DA_nominal = tit["tau_DA"])
    end
    any(nu.parametre_nuisance .== "decalage_jaune_ps") && F.figure_decalage(nu, fig("fig13_tolerance_decalage.png"))
    any(nu.parametre_nuisance .== "blanchiment_accepteur") && F.figure_blanchiment(nu, fig("fig14_blanchiment.png"))
    if isfile(principal) && "variance_parametrique_mM2" in names(df) && any(nu.parametre_nuisance .== "Phi_A_relatif")
        F.figure11(df, nu, fig("fig11_incertitude_classement.png"))
    end
end

# révision du 2026-10-15 : batterie, facteur 8, canal jaune, fusion
isfile(joinpath(racine, "batterie", "batterie.csv")) && F.figure_batterie(joinpath(racine, "batterie"), fig("fig9_batterie_ingestion.png"))
f8 = joinpath(racine, "reel", "facteur8.csv")
isfile(f8) && F.figure_facteur8(CSV.read(f8, DataFrame; types = Dict(:source => String)), fig("fig10_facteur8.png"))
if isfile(joinpath(racine, "canal_jaune", "eta_YD_inversion.csv"))
    mes = CSV.read(normpath(joinpath(@__DIR__, "..", "..", "tests", "csv_exports", "hill_fit_summary.csv")), DataFrame)
    r = mes[(mes.protein.==prot).&(mes.metric.=="fret"), :]
    F.figure_canal_jaune(joinpath(racine, "canal_jaune"), fig("fig12_canal_jaune.png"); K_mes = r.kd_mM[1], plage_mes = r.asymptote[1] / r.baseline[1])
end
fr = joinpath(racine, "resultats", "fusion_resume.csv")
if isfile(fr)
    F.figure_fusion(CSV.read(fr, DataFrame; types = Dict(:capteur => String, :filtre => String)),
        CSV.read(joinpath(racine, "resultats", "fusion_trajectoire.csv"), DataFrame),
        CSV.read(joinpath(racine, "resultats", "fusion_saut.csv"), DataFrame; types = Dict(:capteur => String)), fig("fig7b_fusion.png"))
end
if isfile(principal)
    Resume.ecrire(lire(principal), isfile(nuis) ? lire(nuis) : nothing, joinpath(racine, "resultats", "resume.md"))
end
cm = joinpath(racine, "resultats", "cible_mobile.csv")
isfile(cm) && F.figure_cible_mobile(CSV.read(cm, DataFrame), fig("fig7_cible_mobile.png"))
println("figures écrites dans ", joinpath(racine, "figures"))
