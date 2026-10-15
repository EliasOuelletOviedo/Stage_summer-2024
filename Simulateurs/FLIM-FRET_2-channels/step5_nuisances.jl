# Étape 5 — balayages de nuisance, normalisation A, 10³ et 10⁴ photons.
#
# Pour chaque valeur d'un paramètre de nuisance, les données viennent du monde
# perturbé (« vrai »), et deux jeux d'estimation sont évalués sur les mêmes
# données :
#   « nominale » : estimateurs et calibrations construits sur les valeurs
#                  nominales (ce que tu crois) → biais d'une calibration
#                  imparfaite ou d'un changement entre la calibration et la tranche ;
#   « exacte »   : construits sur le monde perturbé lui-même → effet de la
#                  nuisance sur la précision seule.
# Balayages : fuite η_YD (0 à 1,5×), excitation directe ε (0 à 0,10), τ_D0 et
# τ_DA (±5, ±10 %), τ_A (±5, ±10, ±15 %), ε_dark·Φ_dark (0 à 0,05), et τ_DA
# du capteur (capteurs hypothétiques à FRET moins efficace, mode exact seul),
# et le gain du canal jaune (dérive de détecteur : change l'échelle du jaune,
# pas sa forme ; mode nominal seul).
# Dans le balayage des τ, l'estimateur conjoint en mode nominal est le 5b du
# plan (τ figés à des valeurs fausses).
#
#   julia --project=. -t 8 step5_nuisances.jl [balayages...]
#
# Sortie : sorties/resultats/nuisances.csv

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, CSV, DataFrames, StaticArrays

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
println(S.decrire(m))

avec_η(m, f) = S.remplacer(m; o = S.remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] m.o.η[1, 2]; f*m.o.η[2, 1] m.o.η[2, 2]])))
# gain du canal jaune (détecteur, ou tout ce qui change l'échelle du jaune sans changer sa forme)
avec_gain(m, g) = S.remplacer(m; o = S.remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] m.o.η[1, 2]; g*m.o.η[2, 1] g*m.o.η[2, 2]])))
avec_p(m; kw...) = S.remplacer(m; p = S.remplacer(m.p; kw...))
avec_o(m; kw...) = S.remplacer(m; o = S.remplacer(m.o; kw...))

# (nom, valeurs, monde perturbé, modes)
balayages = [
    ("eta_YD_relatif", [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5], v -> avec_η(m, v), ("nominale", "exacte")),
    ("epsilon", [0.0, 0.02, 0.04, 0.06, 0.08, 0.10], v -> avec_o(m; ε = v), ("nominale", "exacte")),
    ("tau_D0_relatif", [0.90, 0.95, 1.05, 1.10], v -> avec_p(m; τ_D0 = v * m.p.τ_D0), ("nominale", "exacte")),
    ("tau_DA_relatif", [0.90, 0.95, 1.05, 1.10], v -> avec_p(m; τ_DA = v * m.p.τ_DA), ("nominale", "exacte")),
    ("tau_A_relatif", [0.85, 0.90, 0.95, 1.05, 1.10, 1.15], v -> avec_p(m; τ_A = v * m.p.τ_A), ("nominale", "exacte")),
    ("epsilonPhi_dark", [0.0, 0.01, 0.02, 0.03, 0.04, 0.05], v -> avec_o(m; εΦ_dark = v), ("nominale", "exacte")),
    ("tau_DA_capteur_ns", [0.5, 1.0, 1.5, 2.1], v -> avec_p(m; τ_DA = v), ("exacte",)),
    ("gain_jaune_relatif", [0.8, 0.9, 0.95, 1.05, 1.1, 1.2], v -> avec_gain(m, v), ("nominale",)),
]
isempty(ARGS) || filter!(b -> b[1] in ARGS, balayages)

ctxs = S.creer_contextes(m, irf_mat)
ws = ctxs[1].ws
ests_nom = S.estimateurs_standard(m, ws, irf_mat; C_max)
cals_nom = S.calibrer_tous(ests_nom, m, ctxs[1]; C_max)
jeu_nominal(renomme5b) = S.Jeu("nominale", [(renomme5b && k == "5a" ? "5b" : k, e, cals_nom[k]) for (k, e) in ests_nom])

fichier = S.dossier_sorties("resultats", "nuisances.csv")
df = isfile(fichier) && !isempty(ARGS) ? filter(r -> !(r.parametre_nuisance in ARGS), CSV.read(fichier, DataFrame)) : DataFrame()
for (nomb, valeurs, monde, modes) in balayages
    println("\n=== balayage $nomb")
    for v in valeurs
        mv = monde(v)
        jeux = S.Jeu[]
        "nominale" in modes && push!(jeux, jeu_nominal(startswith(nomb, "tau_") && !endswith(nomb, "capteur_ns")))
        if "exacte" in modes
            ests_v = S.estimateurs_standard(mv, ws, irf_mat; C_max)
            cals_v = S.calibrer_tous(ests_v, mv, ctxs[1]; C_max)
            push!(jeux, S.Jeu("exacte", [(k, e, cals_v[k]) for (k, e) in ests_v]))
        end
        # Mêmes germes d'une valeur à l'autre (nombres aléatoires communs) : les
        # différences entre valeurs d'un balayage sont moins bruitées.
        cellules = [S.Cellule(id = 7_000_000 + 100iN + iC, C = C, N = N, norm = :A, R = 1000, nuisance = nomb, valeur = v)
                    for (iN, N) in enumerate((1e3, 1e4)) for (iC, C) in enumerate(S.GRILLE_C)]
        global df = S.executer(cellules, _ -> mv, _ -> jeux, ctxs; fichier, deja = df)
    end
end
println("écrit : ", fichier, " (", nrow(df), " lignes)")
