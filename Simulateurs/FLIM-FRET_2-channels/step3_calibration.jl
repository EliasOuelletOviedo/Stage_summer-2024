# Étape 3 — courbes de calibration sans bruit, K apparents, identifiabilité,
# bornes de Cramér-Rao.
#
#   julia --project=. step3_calibration.jl
#
# Sorties (sorties/calibration/) :
#   courbes_calibration.csv  observable sans bruit de chaque estimateur, 0–150 mM
#   K_apparents.csv          K_vrai → K_apparent (Hill ajustée à chaque courbe), et tes valeurs mesurées
#   identifiabilite.csv      corrélations à la borne, θ = (x, τ_D0, τ_DA, τ_A, t₀, échelles)
#   bornes.csv               σ[Cl⁻] minimal (mM) selon la donnée consommée, par C et N (normalisation A)

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, LinearAlgebra, CSV, DataFrames, TOML

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
println(S.decrire(m))
dossier = S.dossier_sorties("calibration")
mkpath(dossier)

ctx = S.Contexte(S.EspaceTravail(m.inst), S.charger_instance_utilisateur(irf_mat))
ests = S.estimateurs_standard(m, ctx.ws, irf_mat; C_max, ctx = ctx)
cals = S.calibrer_tous(ests, m, ctx; C_max)

# --- courbes de calibration ------------------------------------------------
grille = cals["1"].C
df = DataFrame(C_mM = grille, x = [S.fraction_active(C, m.p) for C in grille])
for (k, _) in ests
    df[!, k] = cals[k].O
end
CSV.write(joinpath(dossier, "courbes_calibration.csv"), df)

# --- K apparents -------------------------------------------------------------
"""B/J (l'inverse du ratio) : c'est ce que trace ton analyse TIFF (C2/C1, étiquettes probablement inversées)."""
struct RatioInverse <: S.Estimateur end
S.observable(::RatioInverse, nB, nY, ctx) = (sum(nB) / sum(nY), :ok)

lignes = []
for (k, e) in vcat(ests, [("1inv", RatioInverse())])
    k in ("2H", "5c") && continue     # 2H : même observable que 2 ; 5c : même courbe que 5a
    r = S.K_apparent(e, m, ctx)
    push!(lignes, (estimateur = k, nom = k == "1inv" ? "1 ratio inversé B/J" : S.NOMS[k], K_vrai = m.p.K_D, n_vrai = m.p.n,
        K_app = r.K_app, n_app = r.n_app, rapport_K = r.K_app / m.p.K_D, obs_0 = r.bas, obs_inf = r.haut))
    @printf("%-28s K_app = %7.2f mM (K_app/K_vrai = %.3f), n_app = %.3f\n", lignes[end].nom, r.K_app, r.K_app / m.p.K_D, r.n_app)
end
h = cals["2H"].hill
push!(lignes, (estimateur = "2H", nom = "2H Hill sur ta grille de titration", K_vrai = m.p.K_D, n_vrai = m.p.n,
    K_app = h[3], n_app = h[4], rapport_K = h[3] / m.p.K_D, obs_0 = h[1], obs_inf = h[2]))
@printf("%-28s K_app = %7.2f mM, n_app = %.3f (τ de %.3f à %.3f ns)\n", "2H Hill (grille titration)", h[3], h[4], h[1], h[2])
# Prédiction analytique pour une lecture pondérée par l'intensité : K_app = K·(D(0)/D(∞))^{1/n}
D(x) = x * m.p.τ_DA + (1 - x) * m.p.τ_D0
K_int = m.p.K_D * (D(m.p.x₀) / D(m.p.x_sat))^(1 / m.p.n)
@printf("prédiction ⟨τ⟩ pondérée en intensité : K_app = %.2f mM\n", K_int)
push!(lignes, (estimateur = "tau_intensite", nom = "⟨τ⟩ pondéré en intensité (analytique)", K_vrai = m.p.K_D, n_vrai = m.p.n,
    K_app = K_int, n_app = m.p.n, rapport_K = K_int / m.p.K_D, obs_0 = NaN, obs_inf = NaN))
# Tes ajustements de Hill expérimentaux (tests/csv_exports/hill_fit_summary.csv)
mes = CSV.read(normpath(joinpath(@__DIR__, "..", "..", "tests", "csv_exports", "hill_fit_summary.csv")), DataFrame)
prot = cfg["capteur"]["proteine"]
for (metric, lab) in (("tau", "mesuré : ton τ mono (Hill)"), ("fret", "mesuré : ton ratio TIFF C2/C1 (Hill)"))
    r = mes[(mes.protein.==prot).&(mes.metric.==metric), :]
    nrow(r) == 1 || continue
    push!(lignes, (estimateur = "mesure_" * metric, nom = lab, K_vrai = m.p.K_D, n_vrai = m.p.n, K_app = r.kd_mM[1], n_app = r.n[1],
        rapport_K = r.kd_mM[1] / m.p.K_D, obs_0 = r.baseline[1], obs_inf = r.asymptote[1]))
    @printf("%-28s K_app = %7.2f mM, n_app = %.3f\n", lab, r.kd_mM[1], r.n[1])
end
CSV.write(joinpath(dossier, "K_apparents.csv"), DataFrame(lignes))

# --- identifiabilité -----------------------------------------------------------
noms_θ = ["x", "τ_D0", "τ_DA", "τ_A", "t₀", "log S_B", "log S_J"]
id = DataFrame()
for (C, N) in ((5.0, 1e4), (20.0, 1e3), (20.0, 1e4), (40.0, 1e4))
    I, _ = S.information_fisher(C, m, N, :A, ctx.ws; params = (:x, :τ_D0, :τ_DA, :τ_A, :t0), canaux = (:B, :Y), echelles = :separees)
    R = S.correlations_crlb(I)
    for i in 1:7, j in i+1:7
        push!(id, (C_mM = C, N = N, parametre_1 = noms_θ[i], parametre_2 = noms_θ[j], correlation = R[i, j]); promote = true)
    end
    @printf("C=%5.1f N=%6.0f : corr(x,τ_A)=%+.3f corr(τ_DA,τ_A)=%+.3f corr(τ_D0,τ_A)=%+.3f corr(x,τ_D0)=%+.3f corr(x,τ_DA)=%+.3f corr(x,t₀)=%+.3f\n",
        C, N, R[1, 4], R[3, 4], R[2, 4], R[1, 2], R[1, 3], R[1, 5])
end
CSV.write(joinpath(dossier, "identifiabilite.csv"), id)

# --- bornes ----------------------------------------------------------------------
defs = (:absolue, :totaux, :conjoint, :conjoint_tauA_libre, :bleu, :bleu_tau_libres)
bornes = DataFrame()
for norm in (:A, :B, :C), N in S.GRILLE_N, C in sort(unique(vcat(S.GRILLE_C, 0.5:0.5:100.0)))
    ligne = (normalisation = String(norm), N_photons = N, C_mM = C)
    vals = [S.crlb_mM(d, C, m, N, norm, ctx.ws) for d in defs]
    push!(bornes, merge(ligne, NamedTuple{defs}(vals)))
end
CSV.write(joinpath(dossier, "bornes.csv"), bornes)
b = bornes[(bornes.normalisation.=="A").&(bornes.C_mM.==20.0), :]
println("\nBornes à 20 mM, normalisation A (mM) :")
show(stdout, b; allrows = true, allcols = true)
println()
