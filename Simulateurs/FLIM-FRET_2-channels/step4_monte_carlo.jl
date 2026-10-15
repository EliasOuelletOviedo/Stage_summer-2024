# Étape 4 — grille principale : 13 concentrations × 7 budgets × 3 normalisations,
# 2 000 réalisations par cellule sous 10³ photons, 1 000 au-dessus.
#
# Tous les estimateurs sont calibrés sur le modèle vrai (sans erreur de
# calibration) : cette grille mesure la variance, le biais d'estimation à bas
# photons et la stabilité. Les erreurs de calibration sont à l'étape 5.
#
#   julia --project=. -t 8 step4_monte_carlo.jl [normalisations...]
#
# Sortie : sorties/resultats/resultats_principaux.csv (une ligne par cellule et
# par estimateur ; colonnes du plan + versions robustes et bornes propres).

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, CSV, DataFrames

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
println(S.decrire(m))
println("fils d'exécution : ", Threads.nthreads())

ctxs = S.creer_contextes(m, irf_mat)
ests = S.estimateurs_standard(m, ctxs[1].ws, irf_mat; C_max, ctx = ctxs[1])
cals = S.calibrer_tous(ests, m, ctxs[1]; C_max)
jeu = S.Jeu("exacte", [(k, e, cals[k]) for (k, e) in ests])

norms = isempty(ARGS) ? [:A, :B, :C] : Symbol.(ARGS)
cellules = S.Cellule[]
for (inorm, norm) in enumerate(norms), (iN, N) in enumerate(S.GRILLE_N), (iC, C) in enumerate(S.GRILLE_C)
    id = 100_000 * findfirst(==(norm), [:A, :B, :C]) + 100 * iN + iC
    push!(cellules, S.Cellule(id = id, C = C, N = N, norm = norm, R = S.repetitions(N)))
end
println(length(cellules), " cellules, ", sum(c.R for c in cellules), " réalisations")

fichier = S.dossier_sorties("resultats", "resultats_principaux.csv")
df = S.executer(cellules, _ -> m, _ -> [jeu], ctxs; fichier)
println("écrit : ", fichier, " (", nrow(df), " lignes)")
