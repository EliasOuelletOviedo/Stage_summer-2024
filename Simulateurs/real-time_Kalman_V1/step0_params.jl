# step0_params.jl — le pont unique entre la titration et le simulateur.
#
# Lit, une fois, l'ajustement global de la titration SCLM du 2026-10-16 à 800 nm (modèle
# « libre » de l'étape 15 du simulateur FLIM-FRET à deux canaux), les IRF mesurées et le
# fond avant le pic, ajuste une Langmuir sur chaque quantité qui dépend du Cl⁻, et écrit le
# tout dans le bloc « step0 » de params.toml. Ensuite, le simulateur ne relit jamais la
# titration.
#
# Lancer depuis ce dossier, dans l'environnement du simulateur précédent (lecteur .sdt) :
#   julia --project=../FLIM-FRET_2-channels step0_params.jl

const ANCIEN = normpath(joinpath(@__DIR__, "..", "FLIM-FRET_2-channels"))
include(joinpath(ANCIEN, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET, CSV, DataFrames, Printf
const S = SimFLIMFRET

const LAMBDA = 800
const DATE = "2026-10-16"
const SORTIES = joinpath(ANCIEN, "sorties", "deux_canaux_2026-10-16_SCLM")
const IRF = joinpath(ANCIEN, "Data", DATE, "IRF_$(LAMBDA)_16X.sdt")

"""
Langmuir y = y₀ + (y∞ − y₀)·C/(K + C), moindres carrés non pondérés sur les 8 points (les
erreurs statistiques des sommes sont négligeables devant l'écart au modèle). K par recherche
sur une grille logarithmique, y₀ et y∞ par moindres carrés linéaires à K fixé.
"""
function langmuir(c, y)
    meilleur = (Inf, 0.0, 0.0, 0.0)
    for K in exp10.(range(-1, 4; length = 4001))
        X = hcat(ones(length(c)), c ./ (K .+ c))
        β = X \ y
        sse = sum(abs2, X * β .- y)
        sse < meilleur[1] && (meilleur = (sse, β[1], β[1] + β[2], K))
    end
    sse, y0, yinf, K = meilleur
    res = [y0 + (yinf - y0) * ci / (K + ci) for ci in c] .- y
    return (y0 = y0, yinf = yinf, K = K, rmse = sqrt(sse / length(c)), ecart_max = maximum(abs, res))
end

g = filter(r -> r.lambda_nm == LAMBDA && r.modele == "libre", CSV.read(joinpath(SORTIES, "global_parametres.csv"), DataFrame))[1, :]
t = sort(filter(r -> r.lambda_nm == LAMBDA, CSV.read(joinpath(SORTIES, "global_cellules.csv"), DataFrame)), :cl_mM)
cl = Float64.(t.cl_mM)
# J/B : photons des trois composantes (diffusion comprise) du jaune sur ceux du bleu, fond exclu.
JB = (t.photons_diffusion_jaune .+ t.photons_tauJ_jaune .+ t.photons_tauA_jaune) ./
     (t.photons_diffusion_bleu .+ t.photons_tau2_bleu .+ t.photons_tau3_bleu)
points = [("JB", JB), ("d_bleu", t.intensite_diffusion_bleu), ("p_bleu", t.amplitude_tau2_bleu),
          ("d_jaune", t.intensite_diffusion_jaune), ("q_jaune", t.amplitude_tauJ_jaune)]

_, pas, duree_irf, _ = S.lire_sdt_voies(IRF)
irf_b, _, niveau_b = S.irf_banc(IRF, 1)
irf_j, _, niveau_j = S.irf_banc(IRF, 2)

vec_toml(v) = "[" * join((@sprintf("%.10g", x) for x in v), ", ") * "]"
io = IOBuffer()
println(io, "# >>> step0 : écrit par step0_params.jl, ne pas éditer")
println(io, "\n[titration]")
println(io, "date = \"$DATE\"\nlongueur_onde_nm = $LAMBDA\nmodele_global = \"libre\"")
println(io, "source = \"FLIM-FRET_2-channels/sorties/deux_canaux_2026-10-16_SCLM (global_parametres.csv, global_cellules.csv)\"")
println(io, "irf = \"FLIM-FRET_2-channels/Data/$DATE/IRF_$(LAMBDA)_16X.sdt\"")
println(io, "\n[instrument]")
println(io, "pas_ns = $(pas)\nperiode_ns = 12.5\nn_canaux = $(length(irf_b))\nfenetre = [14, 195]   # le TAC coupe les canaux 2–13 et 196–256")
println(io, "t0_bleu_ns = $(g.t0_bleu_ps / 1000)\nt0_jaune_ns = $(g.t0_jaune_ps / 1000)")
println(io, "fond_bleu_cps_par_canal = $(niveau_b / duree_irf)    # niveau plat avant le pic de l'IRF, par seconde")
println(io, "fond_jaune_cps_par_canal = $(niveau_j / duree_irf)")
println(io, "irf_bleu = ", vec_toml(irf_b))
println(io, "irf_jaune = ", vec_toml(irf_j))
println(io, "\n[durees_de_vie]")
println(io, "tau_diffusion_ns = 0.001   # copie de l'IRF (lumière diffusée), comme dans l'ajustement global")
println(io, "tau2_ns = $(g.tau2_ns)\ntau3_ns = $(g.tau3_ns)\ntauJ_ns = $(g.tauJ_ns)\ntauA_ns = $(g.tauA_ns)")
println(io, "\n[points]   # les 8 concentrations de l'ajustement global, pour mémoire")
println(io, "cl_mM = ", vec_toml(cl))
println(io, "# JB : rapport des photons jaune/bleu (fond exclu) ; d : fraction d'intensité de la diffusion ;")
println(io, "# p : fraction d'amplitude de τ₂ dans le bleu ; q : fraction d'amplitude de τ_J dans le jaune")
for (nom, y) in points
    println(io, "$nom = ", vec_toml(y))
end
for (nom, y) in points
    l = langmuir(cl, Float64.(y))
    println(io, "\n[langmuir.$nom]")
    println(io, "y0 = $(l.y0)\nyinf = $(l.yinf)\nK_mM = $(l.K)\nrmse = $(l.rmse)\necart_max = $(l.ecart_max)")
    @printf("%-8s y0 = %8.4f  y∞ = %8.4f  K = %8.2f mM  rmse = %.4f  écart max = %.4f\n", nom, l.y0, l.yinf, l.K, l.rmse, l.ecart_max)
end
println(io, "# <<< step0")

chemin = joinpath(@__DIR__, "params.toml")
texte = read(chemin, String)
i, j = findfirst("# >>> step0", texte), findfirst("# <<< step0", texte)
(i === nothing || j === nothing) && error("marqueurs step0 introuvables dans params.toml")
write(chemin, texte[1:first(i)-1] * String(take!(io)) * texte[last(j)+1:end])
println("params.toml mis à jour (pas = $(round(pas, digits = 5)) ns, fond = $(round(niveau_b / duree_irf, digits = 2)) et $(round(niveau_j / duree_irf, digits = 2)) cps par canal)")
