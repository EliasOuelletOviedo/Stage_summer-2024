# Étape 11 — propager l'incertitude des paramètres (changement 2).
#
# 1. Mondes IRF et t₀ : la titration réajustée (étape 9) avec l'IRF dilatée ×0,9 et
#    ×1,1, avec t₀ décalé de ± la gigue observée entre échantillons (3,5 ps), et avec
#    t₀ libre par histogramme. Chaque monde est cohérent : ses paramètres ET son IRF
#    (le code mono de l'utilisateur reçoit la même IRF). On y refait, en calibration
#    exacte, les cellules des figures 1 (N = 10⁴, normalisation A) et 2 (20 mM, trois
#    normalisations) et les K apparents de la figure 5 → bandes sur ces figures.
#    Le décalage de t₀ d'un canal entier (48,8 ps) demandé par la note est exclu :
#    les données le rejettent (χ²ᵣ ≈ 20 au lieu de 1,96, étape 9).
# 2. Covariance : 200 vecteurs (τ_D0, τ_DA, x₀, x_sat, K, n) tirés des covariances
#    de l'ajustement global (×χ²ᵣ) et de l'isotherme ; chaque estimateur, calibré sur
#    le point estimé, lit sans bruit chaque monde tiré → variance paramétrique, ajoutée
#    à côté de la variance photonique dans resultats_principaux.csv. Vérification par
#    un vrai Monte-Carlo sur 40 vecteurs à 20 mM et 10⁴ photons.
# 3. Classement : à partir des balayages τ_A (±20 %), Φ_A (±30 %), η_YD (±50 %) de
#    l'étape 5, le classement des estimateurs survit-il ?
#
#   julia --project=. -t 8 step11_incertitude.jl [mondes|covariance|classement]
#
# Sorties : sorties/incertitude/{mondes_grille,mondes_K_apparents,variance_parametrique,
#           decomposition_variance,classement}.csv ; colonnes ajoutées à resultats_principaux.csv

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, LinearAlgebra, CSV, DataFrames, StableRNGs, TOML

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
prot = cfg["capteur"]["proteine"]
dossier = S.dossier_sorties("incertitude")
mkpath(dossier)
parties = isempty(ARGS) ? ["mondes", "covariance", "classement"] : ARGS
lire(f) = CSV.read(f, DataFrame; types = Dict(:estimateur => String, :normalisation => String, :mode_calibration => String,
    :parametre_nuisance => String))

"""B/J : l'observable de l'analyse TIFF (C2/C1)."""
struct RatioInverse <: S.Estimateur end
S.observable(::RatioInverse, nB, nY, ctx) = (sum(nB) / sum(nY), :ok)

# -----------------------------------------------------------------------------
# 1. Mondes IRF et t₀
# -----------------------------------------------------------------------------
if "mondes" in parties
    reel = CSV.read(S.dossier_sorties("batterie", "reel_perturbations.csv"), DataFrame)
    choix = [("IRF ×0.90", "IRF × 0,9"), ("IRF ×1.10", "IRF × 1,1"), ("t₀ -3.5 ps", "t₀ − 3,5 ps"), ("t₀ +3.5 ps", "t₀ + 3,5 ps"),
        ("t₀ libre par histogramme", "t₀ par histogramme")]
    grille = DataFrame()
    kapp = DataFrame()
    for (prefixe, nom_monde) in choix
        i = findfirst(r -> startswith(r, prefixe), reel.perturbation)
        i === nothing && (println("monde absent : ", prefixe); continue)
        r = reel[i, :]
        mw, irf_w = S.monde_reajuste(m, irf_mat, r)
        println("\n=== monde « $nom_monde » : ", S.decrire(mw))
        ctxs = S.creer_contextes(mw, irf_w)
        ests = S.estimateurs_standard(mw, ctxs[1].ws, irf_w; C_max, ctx = ctxs[1])
        cals = S.calibrer_tous(ests, mw, ctxs[1]; C_max)
        jeu = S.Jeu("exacte", [(k, e, cals[k]) for (k, e) in ests])
        # mêmes germes que la grille principale (étape 4) : différences moins bruitées
        cellules = S.Cellule[]
        for (iC, C) in enumerate(S.GRILLE_C)
            push!(cellules, S.Cellule(id = 100_000 + 100 * 5 + iC, C = C, N = 1e4, norm = :A, R = S.repetitions(1e4)))
        end
        for (inorm, norm) in enumerate((:A, :B, :C)), (iN, N) in enumerate(S.GRILLE_N)
            (norm === :A && N == 1e4) && continue
            iC = findfirst(==(20.0), S.GRILLE_C)
            push!(cellules, S.Cellule(id = 100_000 * inorm + 100 * iN + iC, C = 20.0, N = N, norm = norm, R = S.repetitions(N)))
        end
        tmp = joinpath(dossier, "tmp_monde.csv")
        d = S.executer(cellules, _ -> mw, _ -> [jeu], ctxs; fichier = tmp)
        rm(tmp; force = true)
        d.monde .= nom_monde
        global grille = vcat(grille, d; cols = :union)
        CSV.write(joinpath(dossier, "mondes_grille.csv"), grille)
        # K apparents dans ce monde
        for (k, e) in vcat(ests, [("1inv", RatioInverse())])
            k in ("2H", "5c", "2F") && continue
            rk = S.K_apparent(e, mw, ctxs[1])
            push!(kapp, (monde = nom_monde, estimateur = k, K_vrai = mw.p.K_D, K_app = rk.K_app, rapport_K = rk.K_app / mw.p.K_D))
        end
        h = cals["2H"].hill
        push!(kapp, (monde = nom_monde, estimateur = "2H", K_vrai = mw.p.K_D, K_app = h[3], rapport_K = h[3] / mw.p.K_D))
        CSV.write(joinpath(dossier, "mondes_K_apparents.csv"), kapp)
        @printf("K vrai %.1f mM ; K apparents : τ mono (Hill, grille titration) %.1f mM, B/J %.1f mM, J/B %.1f mM\n", mw.p.K_D, h[3],
            kapp.K_app[(kapp.monde.==nom_monde).&(kapp.estimateur.=="1inv")][1], kapp.K_app[(kapp.monde.==nom_monde).&(kapp.estimateur.=="1")][1])
    end
    # résumé : σ à 20 mM et 10⁴ photons (normalisation A), par monde
    println("\nσ[Cl⁻] à 20 mM, 10⁴ photons, normalisation A (mM) :")
    princ = lire(S.dossier_sorties("resultats", "resultats_principaux.csv"))
    nomi = princ[(princ.C_vrai.==20.0).&(princ.N_photons.==1e4).&(princ.normalisation.=="A"), :]
    for k in ("1", "5d", "5a", "4", "2H")
        @printf("  %-3s nominal %6.2f |", k, nomi.sigma_mM[nomi.estimateur.==k][1])
        for (_, nom_monde) in choix
            d = grille[(grille.monde.==nom_monde).&(grille.C_vrai.==20.0).&(grille.N_photons.==1e4).&(grille.normalisation.=="A").&(grille.estimateur.==k), :]
            nrow(d) == 1 && @printf("  %s %6.2f", nom_monde, d.sigma_mM[1])
        end
        println()
    end
end

# -----------------------------------------------------------------------------
# 2. Covariance des ajustements → variance paramétrique
# -----------------------------------------------------------------------------
if "covariance" in parties
    hs = S.lire_titration(prot)
    C = [h.C for h in hs]
    g = S.ajuster_global(hs, m.inst; τ_D0_init = m.p.τ_D0, τ_DA_init = m.p.τ_DA)
    iso = S.ajuster_isotherme(C, g.x)
    cov_τ = g.cov[1:2, 1:2] .* g.chi2_global          # surdispersion : χ²ᵣ = 1,96
    @printf("σ(τ_D0) = %.4f ns, σ(τ_DA) = %.4f ns (×√χ²ᵣ) ; isotherme : σ(x₀) = %.4f, σ(x_sat) = %.4f, σ(K) = %.2f mM, σ(n) = %.3f\n",
        sqrt(cov_τ[1, 1]), sqrt(cov_τ[2, 2]), sqrt.(diag(iso.cov))...)
    J = 200
    tirages = S.tirer_parametres(m, cov_τ, iso, J, StableRNG(2026))
    ctxs = S.creer_contextes(m, irf_mat)
    ests = S.estimateurs_standard(m, ctxs[1].ws, irf_mat; C_max, ctx = ctxs[1])
    cals = S.calibrer_tous(ests, m, ctxs[1]; C_max)
    biais = Array{Float64}(undef, J, length(S.GRILLE_C), length(ests))
    suivant = Threads.Atomic{Int}(1)
    @sync for ctx in ctxs
        Threads.@spawn while (j = Threads.atomic_add!(suivant, 1)) <= J
            mj = S.remplacer(m; p = tirages[j])
            for (iC, Cv) in enumerate(S.GRILLE_C)
                b = S.biais_sans_bruit(ests, cals, mj, ctx, Cv)
                for (ie, (k, _)) in enumerate(ests)
                    biais[j, iC, ie] = b[k]
                end
            end
        end
    end
    vp = DataFrame()
    for (ie, (k, _)) in enumerate(ests), (iC, Cv) in enumerate(S.GRILLE_C)
        b = filter(isfinite, biais[:, iC, ie])
        push!(vp, (estimateur = k, C_vrai = Cv, n_tirages = length(b), biais_parametrique_moyen_mM = mean(b),
            variance_parametrique_mM2 = var(b), sigma_parametrique_mM = std(b), q05_mM = quantile(b, 0.05), q95_mM = quantile(b, 0.95)))
    end
    CSV.write(joinpath(dossier, "variance_parametrique.csv"), vp)
    println("σ paramétrique (mM) selon [Cl⁻] :")
    for k in ("1", "5d", "5a", "4", "2H", "6")
        d = vp[vp.estimateur.==k, :]
        @printf("  %-3s %s\n", k, join([@sprintf("%g:%.2f", c, s) for (c, s) in zip(d.C_vrai, d.sigma_parametrique_mM)], "  "))
    end
    # colonnes ajoutées à la grille principale : variance photonique et paramétrique séparées
    fp = S.dossier_sorties("resultats", "resultats_principaux.csv")
    princ = lire(fp)
    select!(princ, Not(intersect(names(princ), ["variance_photonique_mM2", "variance_parametrique_mM2", "biais_parametrique_moyen_mM", "sigma_total_mM"])))
    princ = leftjoin(princ, select(vp, :estimateur, :C_vrai, :variance_parametrique_mM2, :biais_parametrique_moyen_mM); on = [:estimateur, :C_vrai])
    princ.variance_photonique_mM2 = princ.sigma_mM .^ 2
    princ.sigma_total_mM = sqrt.(princ.variance_photonique_mM2 .+ coalesce.(princ.variance_parametrique_mM2, NaN))
    sort!(princ, [:normalisation, :N_photons, :C_vrai])
    CSV.write(fp, princ)
    println("colonnes variance_photonique_mM2, variance_parametrique_mM2, sigma_total_mM ajoutées à ", fp)

    # vérification : un vrai Monte-Carlo sur 40 vecteurs tirés, 20 mM, 10⁴ photons, normalisation A
    Jv, Rv = 40, 250
    jeu = S.Jeu("nominale", [(k, e, cals[k]) for (k, e) in ests])
    dec = DataFrame()
    par = Dict{String,Vector{Tuple{Float64,Float64}}}()
    for j in 1:Jv
        mj = S.remplacer(m; p = tirages[j])
        cell = S.Cellule(id = 11_000_000 + j, C = 20.0, N = 1e4, norm = :A, R = Rv)
        for l in S.simuler_cellule(cell, mj, [jeu], ctxs; crlb = false)
            push!(get!(par, l.estimateur, Tuple{Float64,Float64}[]), (l.biais_mM, l.sigma_mM))
        end
    end
    nomi = princ[(princ.C_vrai.==20.0).&(princ.N_photons.==1e4).&(princ.normalisation.=="A"), :]
    println("\nDécomposition de la variance à 20 mM, 10⁴ photons (mM²) :")
    for (ie, (k, _)) in enumerate(ests)
        v = par[k]
        moy_var = mean(last.(v) .^ 2)
        var_moy = var(first.(v))
        vph = nomi.variance_photonique_mM2[nomi.estimateur.==k][1]
        vpa = vp.variance_parametrique_mM2[(vp.estimateur.==k).&(vp.C_vrai.==20.0)][1]
        push!(dec, (estimateur = k, variance_photonique_MC_mM2 = moy_var, variance_parametrique_MC_mM2 = var_moy,
            variance_totale_MC_mM2 = moy_var + var_moy, variance_photonique_grille_mM2 = vph, variance_parametrique_sans_bruit_mM2 = vpa,
            variance_totale_predite_mM2 = vph + vpa))
        @printf("  %-3s MC : photonique %8.2f + paramétrique %7.2f = %8.2f   prédit : %8.2f + %7.2f = %8.2f\n", k, moy_var, var_moy,
            moy_var + var_moy, vph, vpa, vph + vpa)
    end
    CSV.write(joinpath(dossier, "decomposition_variance.csv"), dec)
end

# -----------------------------------------------------------------------------
# 3. Le classement survit-il aux paramètres supposés ?
# -----------------------------------------------------------------------------
if "classement" in parties
    nu = lire(S.dossier_sorties("resultats", "nuisances.csv"))
    princ = lire(S.dossier_sorties("resultats", "resultats_principaux.csv"))
    cles = ["1", "5d", "5a", "5c", "4", "3", "6", "2", "2H"]
    cl = DataFrame()
    for (param, valeurs) in (("tau_A_relatif", [0.8, 0.9, 1.1, 1.2]), ("Phi_A_relatif", [0.7, 0.85, 1.15, 1.3]),
                             ("eta_YD_relatif", [0.5, 0.75, 1.25, 1.5]))
        for v in valeurs, mode in ("exacte", "nominale"), N in (1e3, 1e4), Cv in (5.0, 20.0, 40.0)
            d = nu[(nu.parametre_nuisance.==param).&isapprox.(nu.valeur_nuisance, v).&(nu.mode_calibration.==mode).&(nu.N_photons.==N).&(nu.C_vrai.==Cv), :]
            r0 = princ[(princ.normalisation.=="A").&(princ.N_photons.==N).&(princ.C_vrai.==Cv), :]
            ks = [k for k in cles if any(d.estimateur .== (mode == "nominale" && k == "5a" && startswith(param, "tau") ? "5b" : k)) && any(r0.estimateur .== k)]
            length(ks) < 4 && continue
            kd(k) = mode == "nominale" && k == "5a" && startswith(param, "tau") ? "5b" : k
            a = [r0.RMSE_mM[r0.estimateur.==k][1] for k in ks]
            b = [d.RMSE_mM[d.estimateur.==kd(k)][1] for k in ks]
            ga, gb = ks[argmin(a)], ks[argmin(b)]
            # ex æquo à 5 % près avec le meilleur
            ex_a = Set(ks[a.<=1.05minimum(a)])
            ex_b = Set(ks[b.<=1.05minimum(b)])
            push!(cl, (parametre = param, valeur = v, mode = mode, N_photons = N, C_vrai = Cv, kendall = S.kendall(a, b),
                meilleur_nominal = join(sort(collect(ex_a)), "≈"), meilleur = join(sort(collect(ex_b)), "≈"),
                meilleur_inchange = !isempty(intersect(ex_a, ex_b)),
                ordre_nominal = join(ks[sortperm(a)], " < "), ordre = join(ks[sortperm(b)], " < "),
                RMSE_ratio = b[findfirst(==("1"), ks)], RMSE_5a = b[findfirst(==("5a"), ks)], RMSE_4 = b[findfirst(==("4"), ks)]))
        end
    end
    CSV.write(joinpath(dossier, "classement.csv"), cl)
    println("\nClassement par RMSE (20 mM, 10⁴ photons) :")
    for r in eachrow(cl[(cl.C_vrai.==20.0).&(cl.N_photons.==1e4), :])
        @printf("  %-15s ×%.2f %-9s  τ de Kendall %+.2f  meilleur %-6s (nominal %-6s)  %s\n", r.parametre, r.valeur, r.mode, r.kendall,
            r.meilleur, r.meilleur_nominal, r.ordre)
    end
    @printf("Kendall moyen (mode exact) : %.2f ; meilleur inchangé dans %d cas sur %d (exact) et %d sur %d (nominal)\n",
        mean(cl.kendall[cl.mode.=="exacte"]), count(cl.meilleur_inchange[cl.mode.=="exacte"]), count(cl.mode .== "exacte"),
        count(cl.meilleur_inchange[cl.mode.=="nominale"]), count(cl.mode .== "nominale"))
end
