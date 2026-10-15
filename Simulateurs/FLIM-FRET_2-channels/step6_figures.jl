# Étape 6 — figures, à partir des fichiers plats des étapes 1, 3, 4 et 5.
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
F.figure5(kapp, courbes, fig("fig5_K_apparents.png"); K_vrai = tit["K_D"])

principal = joinpath(racine, "resultats", "resultats_principaux.csv")
if isfile(principal)
    df = lire(principal)
    F.figure1(df, fig("fig1_sigma_selon_C.png"))
    F.figure2(df, fig("fig2_sigma_selon_N.png"))
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
end
if isfile(principal)
    Resume.ecrire(lire(principal), isfile(nuis) ? lire(nuis) : nothing, joinpath(racine, "resultats", "resume.md"))
end
cm = joinpath(racine, "resultats", "cible_mobile.csv")
isfile(cm) && F.figure_cible_mobile(CSV.read(cm, DataFrame), fig("fig7_cible_mobile.png"))
println("figures écrites dans ", joinpath(racine, "figures"))
