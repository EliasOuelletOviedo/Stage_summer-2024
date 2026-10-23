# Étape 15 — analyse FLIM à deux canaux d'une titration de banc (2026-10-16 : SCLM, KIRIN).
#
# D'après « Analyse FLIM à deux canaux : modifications du code » (2026-10-16). Données :
# Data/2026-10-16/SCLM (défaut) ou KIRIN (--donnees), 8 concentrations × 3
# longueurs d'onde (750, 800, 850 nm) = 24 cellules, 10 acquisitions de 1 s chacune, deux
# canaux. Le dossier `unknown` est ignoré. Deux réponses, chacune avec sa règle écrite
# d'avance : la meilleure lecture (RMSE en mM) et la meilleure longueur d'onde.
#
# Étapes (une sortie persistée après chacune, dans sorties/deux_canaux_<date>/) :
#   charger   cellules, contrôles de Poisson et de tendance          (Charger.jl)
#   global    ajustement global par longueur d'onde, lien τ_J = τ₂,
#             diagnostics par point, signature du blanchiment         (Global.jl, Blanchir.jl)
#   calibrer  neuf lectures à pleins photons, courbes, K apparents,
#             borne de Fisher                                         (Lectures.jl, Calibrer.jl)
#   amincir   Monte-Carlo : sommes amincies et jumeau, 7 budgets      (Amincir.jl)
#   blanchir  scénarios A, B, C injectés à 1, 2, 5, 10 %             (Blanchir.jl)
#   agreger   fichier plat, règles, carte finale, figures             (Agreger.jl)
#   inconnue  concentration de l'échantillon inconnu (calib/unknown) par chaque lecture,
#             avec ses incertitudes                                    (Calibrer.jl)
#
#   julia --project=. -t 8 step15_analyse_deux_canaux.jl                    # tout (≈ 1 h)
#   julia --project=. -t 8 step15_analyse_deux_canaux.jl --essai            # 800 nm, 3 concentrations, peu de tirages (minutes)
#   julia --project=. -t 8 step15_analyse_deux_canaux.jl --donnees Data/2026-10-16/KIRIN charger global calibrer amincir blanchir agreger
#   julia --project=. -t 8 step15_analyse_deux_canaux.jl agreger            # une étape (reprend les sorties précédentes)
#   julia --project=. step15_analyse_deux_canaux.jl rejouer 800 20 1000 17 [jumeau] [A 0.05]   # un tirage seul
#
# Puissance laser : renseigner `puissance_mW` dans parametres_deux_canaux.toml (wattmètre)
# rend citable la comparaison à temps égal ; sans elle, elle reste conditionnelle.

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
include(joinpath(@__DIR__, "src", "deux_canaux", "DeuxCanaux.jl"))
using .SimFLIMFRET, .DeuxCanaux
const S = SimFLIMFRET
const D = DeuxCanaux
using Printf, Statistics, CSV, DataFrames, TOML, Serialization, Dates

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
args = copy(ARGS)
essai = "--essai" in args
filter!(!=("--essai"), args)
i_d = findfirst(==("--donnees"), args)
dossier_donnees = i_d === nothing ? joinpath(@__DIR__, "Data", "2026-10-16", "SCLM") : args[i_d+1]
i_d === nothing || deleteat!(args, i_d:i_d+1)
ETAPES = ["charger", "global", "calibrer", "amincir", "blanchir", "agreger", "inconnue"]
rejeu = !isempty(args) && args[1] == "rejouer"
etapes = rejeu ? String[] : isempty(args) ? ETAPES : args
all(in(ETAPES), etapes) || error("étapes connues : $(join(ETAPES, ", "))")

date_donnees = basename(dirname(normpath(dossier_donnees)))
jeu = basename(normpath(dossier_donnees))                 # le capteur : SCLM, KIRIN…
dossier = S.dossier_sorties("deux_canaux_" * date_donnees * "_" * jeu * (essai ? "_essai" : ""))
mkpath(dossier)
fig(nom) = S.dossier_sorties("figures", (essai ? "essai_" : "") * replace(nom, "deux_canaux" => "deux_canaux_" * jeu))
parametres = let f = joinpath(@__DIR__, "parametres_deux_canaux.toml")
    isfile(f) ? TOML.parsefile(f) : Dict{String,Any}()
end
puissances = Dict{Int,Float64}(parse(Int, k) => Float64(v) for (k, v) in get(parametres, "puissance_mW", Dict()) if v isa Real)
budget_vise = Float64(get(parametres, "budget_vise", D.BUDGET_VISE))
λs = essai ? (800,) : D.LONGUEURS_ONDE
concentrations = essai ? (cs = first.(D.dossiers_concentration(dossier_donnees)); [cs[1], 20.0, cs[end]]) : nothing
n_tirages = essai ? (N -> N < 1e3 ? 40 : 20) : D.n_tirages_defaut
n_tirages_blanchiment = essai ? (_ -> 10) : (_ -> Int(get(parametres, "tirages_blanchiment", 100)))
nf = Threads.nthreads()
println("Analyse à deux canaux — $(jeu), données $(dossier_donnees)", essai ? " (ESSAI : 800 nm, $(join(concentrations, "/")) mM)" : "", " — $(nf) fils")

# -----------------------------------------------------------------------------
# 1. Charger et contrôler
# -----------------------------------------------------------------------------
cellules, sources_fond = D.charger(dossier_donnees; λs, concentrations)
par_λ = Dict(λ => sort([c for c in cellules if c.λ == λ]; by = c -> c.cl) for λ in λs)
println(length(cellules), " cellules lues ; fond : ", join(unique(values(sources_fond)), " ; "))
if "charger" in etapes
    controles = DataFrame([D.controler(c) for c in cellules])
    CSV.write(joinpath(dossier, "controles_cellules.csv"), controles)
    CSV.write(joinpath(dossier, "acquisitions.csv"), D.table_acquisitions(cellules))
    println("\n1. Contrôles (10 acquisitions de 1 s par cellule) :")
    for r in eachrow(controles)
        @printf("  %d nm %6.1f mM : N_B %7.0f/s  J/B %.3f  σ(J/B) = %.2f × Poisson%s  pentes B %+.2f %% (t %+.1f), J %+.2f %% (t %+.1f), J/B %+.2f %% (t %+.1f)%s\n",
            r.lambda_nm, r.cl_mM, r.N_bleu / r.duree_s, r.JB, r.sigma_JB_sur_poisson, r.drapeau_poisson ? " [> 2 × Poisson]" : "",
            r.pente_bleu_pct, r.t_bleu, r.pente_jaune_pct, r.t_jaune, r.pente_JB_pct, r.t_JB,
            r.drapeau_tendance ? "  [tendance : $(r.tendances)]" : "")
    end
end
controles = CSV.read(joinpath(dossier, "controles_cellules.csv"), DataFrame)

# -----------------------------------------------------------------------------
# 2. Ajustement global, lien des composantes rapides, diagnostics, signature du blanchiment
# -----------------------------------------------------------------------------
if "global" in etapes
    println("\n2. Ajustement global (8 concentrations × 2 canaux par longueur d'onde) :")
    modeles = Dict{Int,D.ModeleGlobal}()
    params = DataFrame()
    tab = DataFrame()
    liens = DataFrame()
    svd_tab = DataFrame()
    svd_poids = DataFrame()
    signatures = DataFrame()
    for λ in λs
        cs = par_λ[λ]
        lie = D.ajuster_global(cs; lien = true)
        libre = D.ajuster_global(cs; lien = false, τ_init = (lie.τ2, lie.τ3, lie.τA, 1.0), t0_init = (lie.t0B, lie.t0J))
        test = D.tester_lien(lie, libre)
        push!(liens, test)
        mg = test.lien_tient ? lie : libre
        modeles[λ] = mg
        for m in (lie, libre)
            push!(params, (lambda_nm = λ, modele = m.lien ? "lie" : "libre", retenu = m === mg, tau2_ns = m.τ2, tau3_ns = m.τ3,
                tauA_ns = m.τA, tauJ_ns = m.τJ, sigma_tau2 = m.σ_τ[1], sigma_tau3 = m.σ_τ[2], sigma_tauA = m.σ_τ[3], sigma_tauJ = m.σ_τ[4],
                t0_bleu_ps = 1000m.t0B, t0_jaune_ps = 1000m.t0J, chi2r = m.χ2, deviance = m.deviance, n_parametres = m.nparam, converge = m.converge))
        end
        t = D.table_global(mg)
        t0B, t0J = D.t0_par_cellule(cs, mg)
        t.t0_bleu_par_cellule_ps = 1000 .* t0B
        t.t0_jaune_par_cellule_ps = 1000 .* t0J
        mel = D.tester_melange(cs, mg)
        t.f_melange = [x.f for x in mel]
        t.chi2r_melange = [x.chi2r for x in mel]
        t.exces_melange = t.chi2r_melange .- [x.chi2r_global for x in mel]
        append!(tab, t)
        for (canal, s) in D.residus_svd(cs, mg)
            append!(svd_tab, DataFrame(lambda_nm = λ, canal = canal, canal_tcspc = collect(D.FENETRE), vecteur = s.vecteur, part = s.part))
            append!(svd_poids, DataFrame(lambda_nm = λ, canal = canal, cl_mM = mg.cl, poids = s.poids, part = s.part))
        end
        ctx = D.contexte(mg, cs[1]; mono = false)
        for c in cs
            push!(signatures, D.signature_blanchiment(c, ctx))
        end
        fu = D.tester_fuite(t)
        χ2 = t.chi2r_bleu .+ t.chi2r_jaune
        @printf("  %d nm : lié χ²ᵣ %.3f, libre χ²ᵣ %.3f ; Δdéviance/χ²ᵣ = %.0f (seuil 3,84) → le lien %s (τ₂ bleu %.3f ns, τ_J libre %.3f ns) ; modèle retenu : %s\n",
            λ, lie.χ2, libre.χ2, test.statistique, test.lien_tient ? "TIENT" : "FORCE", libre.τ2, libre.τJ, mg.lien ? "lié" : "libre")
        @printf("        τ₂ %.3f, τ₃ %.3f, τ_A %.3f, τ_J %.3f ns ; t₀ bleu %+.1f ps, jaune %+.1f ps ; diffusion %.0f photons/s (bleu, %.0f–%.0f)\n",
            mg.τ2, mg.τ3, mg.τA, mg.τJ, 1000mg.t0B, 1000mg.t0J, mean(t.diffusion_bleu_par_s), extrema(t.diffusion_bleu_par_s)...)
        @printf("        λ_fuite %.3f → %.3f (Spearman %+.2f, variation %+.0f %% sur la plage) : %s ; profil χ²ᵣ global « %s », excès du mélange à deux états « %s » (max %+.2f à %g mM)\n",
            t.lambda_fuite[1], t.lambda_fuite[end], fu.spearman, fu.variation_pct_sur_plage, fu.constant ? "constant" : "dérive avec le Cl⁻",
            D.forme_profil(t.cl_mM, χ2), D.forme_profil(t.cl_mM, t.exces_melange), maximum(t.exces_melange), t.cl_mM[argmax(t.exces_melange)])
        @printf("        t₀ par cellule : bleu %s ps (Spearman %+.2f), jaune %s ps (Spearman %+.2f)\n",
            join([@sprintf("%.1f", x) for x in t.t0_bleu_par_cellule_ps], " "), D.spearman(t.cl_mM, t.t0_bleu_par_cellule_ps),
            join([@sprintf("%.1f", x) for x in t.t0_jaune_par_cellule_ps], " "), D.spearman(t.cl_mM, t.t0_jaune_par_cellule_ps))
        for r in eachrow(t)
            @printf("        %6.1f mM : bleu I %.3f/%.3f/%.3f (diff/τ₂/τ₃) α %.3f/%.3f ; jaune I %.3f/%.3f/%.3f ; χ²ᵣ %.2f/%.2f ; mélange f %.3f χ²ᵣ %.2f\n",
                r.cl_mM, r.intensite_diffusion_bleu, r.intensite_tau2_bleu, r.intensite_tau3_bleu, r.amplitude_tau2_bleu, r.amplitude_tau3_bleu,
                r.intensite_diffusion_jaune, r.intensite_tauJ_jaune, r.intensite_tauA_jaune, r.chi2r_bleu, r.chi2r_jaune, r.f_melange, r.chi2r_melange)
        end
    end
    CSV.write(joinpath(dossier, "global_parametres.csv"), params)
    CSV.write(joinpath(dossier, "global_lien.csv"), liens)
    CSV.write(joinpath(dossier, "global_cellules.csv"), tab)
    CSV.write(joinpath(dossier, "global_residus_svd.csv"), svd_tab)
    CSV.write(joinpath(dossier, "global_residus_poids.csv"), svd_poids)
    CSV.write(joinpath(dossier, "tendances_blanchiment.csv"), signatures)
    serialize(joinpath(dossier, "modeles_globaux.jls"), modeles)
    println("  signature du blanchiment : ", join(["$(r.lambda_nm)/$(r.cl_mM) $(r.scenario)" for r in eachrow(signatures) if r.scenario != "aucun"], ", "),
        isempty(filter(!=("aucun"), signatures.scenario)) ? "aucune cellule ne dérive" : "")
end
modeles = deserialize(joinpath(dossier, "modeles_globaux.jls"))

besoin_lectures = any(in(etapes), ("calibrer", "amincir", "blanchir", "inconnue")) || rejeu
contextes = Dict{Int,D.Contexte}()
if besoin_lectures
    for λ in λs
        contextes[λ] = D.contexte(modeles[λ], par_λ[λ][1]; fils = nf)
    end
end

# -----------------------------------------------------------------------------
# 3. Lectures à pleins photons, calibration, K apparents, borne de Fisher
# -----------------------------------------------------------------------------
if "calibrer" in etapes
    println("\n3. Lectures à pleins photons et calibration :")
    refs = DataFrame()
    for λ in λs, c in par_λ[λ]
        ctx = contextes[λ]
        mg = modeles[λ]
        o, e = D.lire(ctx, D.somme(c, D.BLEU), D.somme(c, D.JAUNE), mg.duree_somme_s[D.indice(mg, c.cl)], 1)
        par_acq = [D.lire(ctx, D.acquisition(c, D.BLEU, i), D.acquisition(c, D.JAUNE, i), c.duree_s, 1)[1] for i in 1:D.n_acq(c)]
        for (q, k) in enumerate(D.LECTURES)
            v = filter(isfinite, [p[q] for p in par_acq])
            push!(refs, (lambda_nm = λ, cl_mM = c.cl, lecture = k, o = o[q], code = e[q], sigma_ref = length(v) >= 3 ? std(v) / sqrt(length(v)) : NaN,
                o_moyenne_acquisitions = isempty(v) ? NaN : mean(v), n_acquisitions_valides = length(v)))
        end
    end
    CSV.write(joinpath(dossier, "lectures_pleins_photons.csv"), refs)
    courbes = Dict{Tuple{Int,String},D.Courbe}()
    cal = DataFrame()
    for λ in λs, k in D.LECTURES
        r = sort(refs[(refs.lambda_nm.==λ).&(refs.lecture.==k), :], :cl_mM)
        C = D.courbe(k, λ, Vector{Float64}(r.cl_mM), Vector{Float64}(r.o), Vector{Float64}(r.sigma_ref))
        courbes[(λ, k)] = C
        L = D.langmuir(C.c, C.o, C.σ_ref)
        σc = D.sigma_calibration(C)
        zv, ouv = D.violation_monotonie(C)
        push!(cal, (lambda_nm = λ, lecture = k, nom = D.NOMS_LECTURES[k], monotone = C.monotone, violation_sigma = zv, violation_entre_mM = ouv,
            o_0mM = C.o[1], o_max = C.o[end],
            plage = C.o[1] / C.o[end], K_apparent_mM = L.K, K_bas_mM = L.K_bas, K_haut_mM = L.K_haut, chi2r_langmuir = L.chi2r,
            sigma_calibration_20mM = σc[findfirst(==(20.0), C.c) === nothing ? 1 : findfirst(==(20.0), C.c)]))
        @printf("  %d nm %-38s %s  %10.4f → %10.4f  K app. %6.1f mM [%5.1f ; %5.1f]  σ_calib(20 mM) %.3f mM\n", λ, D.NOMS_LECTURES[k],
            C.monotone ? "monotone   " : @sprintf("NON MONOTONE (%s mM, %.1f σ)", ouv, zv), C.o[1], C.o[end], L.K, L.K_bas, L.K_haut,
            cal.sigma_calibration_20mM[end])
    end
    CSV.write(joinpath(dossier, "calibration.csv"), cal)
    serialize(joinpath(dossier, "courbes.jls"), courbes)
    fisher = DataFrame()
    for λ in λs
        append!(fisher, D.bornes_fisher(modeles[λ], contextes[λ], D.BUDGETS))
    end
    # R du Kalman : borne de toute l'information au nombre de photons d'une acquisition de 1 s
    kal = DataFrame()
    for λ in λs, c in par_λ[λ]
        nb = D.photons(D.somme(c, D.BLEU)) / D.n_acq(c)
        b = D.bornes_fisher(modeles[λ], contextes[λ], [nb])
        r = b[b.cl_mM.==c.cl, :]
        push!(kal, (lambda_nm = λ, cl_mM = c.cl, photons_bleus_par_acquisition = nb, sigma_fisher_toute_mM = r.toute[1],
            R_kalman_mM2 = r.toute[1]^2, sigma_fisher_rapport_mM = r.rapport[1]))
    end
    CSV.write(joinpath(dossier, "fisher.csv"), fisher)
    CSV.write(joinpath(dossier, "kalman_R.csv"), kal)
end

# -----------------------------------------------------------------------------
# 4. Amincissement et Monte-Carlo
# -----------------------------------------------------------------------------
if "amincir" in etapes
    println("\n4. Monte-Carlo : sommes amincies (réel) et jumeau, budgets $(join([@sprintf("%g", N) for N in D.BUDGETS], ", ")) photons bleus")
    est = D.monte_carlo(cellules, modeles, contextes; n_tirages, titre = "intacts")
    CSV.write(joinpath(dossier, "estimations_tirages.csv"), est)
    println("  ", nrow(est), " tirages lus")
end

# -----------------------------------------------------------------------------
# 5. Blanchiment injecté
# -----------------------------------------------------------------------------
if "blanchir" in etapes
    taux = essai ? (0.02, 0.10) : D.TAUX_BLANCHIMENT
    println("\n5. Blanchiment : scénarios A, B, C à ", join([@sprintf("%g %%", 100t) for t in taux], ", "))
    scen = [(sc, β) for sc in D.SCENARIOS for β in taux]
    estb = D.monte_carlo(cellules, modeles, contextes; sources = (:reel,), n_tirages = n_tirages_blanchiment, scenarios = scen,
        titre = "blanchis")
    CSV.write(joinpath(dossier, "estimations_blanchiment.csv"), estb)
    println("  ", nrow(estb), " tirages blanchis lus")
end

# -----------------------------------------------------------------------------
# 6. Agrégation
# -----------------------------------------------------------------------------
if "agreger" in etapes
    println("\n6. Agrégation")
    courbes = deserialize(joinpath(dossier, "courbes.jls"))
    est = CSV.read(joinpath(dossier, "estimations_tirages.csv"), DataFrame)
    fb = joinpath(dossier, "estimations_blanchiment.csv")
    isfile(fb) && (est = vcat(est, CSV.read(fb, DataFrame); cols = :union))
    est.source = String.(est.source)
    est.scenario = String.(est.scenario)
    long = D.inverser_tirages(est, courbes)
    stats = D.statistiques(long, courbes)
    flat = D.fichier_plat(stats, controles)
    D.ecrire_avec_en_tete(joinpath(dossier, "resultats_deux_canaux.csv"), flat, D.EN_TETE_FICHIER_PLAT)
    CSV.write(joinpath(dossier, "jumeau_contre_reel.csv"), D.jumeau_contre_reel(stats))
    fisher = CSV.read(joinpath(dossier, "fisher.csv"), DataFrame)
    CSV.write(joinpath(dossier, "efficacites.csv"), D.efficacites(flat, fisher))
    regles = D.regle_lecture(flat)
    CSV.write(joinpath(dossier, "regle_lecture.csv"), regles)
    signatures = CSV.read(joinpath(dossier, "tendances_blanchiment.csv"), DataFrame)
    pentes = any(flat.scenario_blanchiment .!= "aucun") ? D.pentes_blanchiment(flat) : DataFrame()
    if nrow(pentes) > 0
        CSV.write(joinpath(dossier, "pentes_blanchiment.csv"), pentes)
        CSV.write(joinpath(dossier, "blanchiment_tolerable.csv"), D.blanchiment_tolerable(pentes; budget = budget_vise))
        CSV.write(joinpath(dossier, "regle_lecture_blanchie.csv"), D.regle_lecture_blanchie(flat, pentes, signatures, controles))
    end
    refs = CSV.read(joinpath(dossier, "lectures_pleins_photons.csv"), DataFrame)
    globaux = CSV.read(joinpath(dossier, "global_cellules.csv"), DataFrame)
    lo = D.regle_longueur_onde(flat, regles, controles, refs, globaux; budget = budget_vise, puissances_mW = puissances)
    CSV.write(joinpath(dossier, "regle_longueur_onde.csv"), lo)
    carte = D.carte_finale(flat)
    CSV.write(joinpath(dossier, "carte_finale.csv"), carte)
    println("  règle de lecture au budget visé ($(budget_vise) photons bleus) :")
    for r in eachrow(regles[(regles.budget_photons.==budget_vise), :])
        @printf("    %d nm %-38s échecs max %5.1f %%  |biais| max %7.2f mM  RMSE médiane 5–35 %7.2f mM  pire %7.2f mM (%g mM)  %s\n",
            r.lambda_nm, D.NOMS_LECTURES[r.lecture], 100r.echec_max, r.biais_max_mM, r.rmse_mediane_5_35_mM, r.pire_cas_mM, r.pire_cas_cl_mM,
            r.survit ? (r.gagnante ? "GAGNANTE" : "rang $(r.rang)") : "éliminée ($(r.eliminee_par))")
    end
    for r in eachrow(lo)
        @printf("  %d nm : plage J/B ×%.2f ; photons égaux : %s (%.2f mM) ; temps égal (%.0f ms) : %s (%.2f mM) ; %.0f photons bleus/s à 20 mM\n",
            r.lambda_nm, r.plage_JB, r.gagnante_photons_egaux, r.rmse_mediane_photons_egaux_mM, 1000r.duree_temps_egal_s,
            r.gagnante_temps_egal, r.rmse_mediane_temps_egal_mM, r.photons_bleus_par_s_20mM)
    end
    resume = Dict{String,Any}("date" => string(Dates.today()), "donnees" => dossier_donnees, "budget_vise" => budget_vise,
        "essai" => essai, "fond" => join(unique(values(sources_fond)), " ; "),
        "puissance_mW" => Dict(string(k) => v for (k, v) in puissances))
    for r in eachrow(regles[regles.gagnante, :])
        resume["gagnante_$(r.lambda_nm)nm_$(round(Int, r.budget_photons))"] = r.lecture
    end
    open(joinpath(dossier, "resume.toml"), "w") do io
        TOML.print(io, resume; sorted = true)
    end
end
if "agreger" in etapes
    isdefined(Main, :Figures) || include(joinpath(@__DIR__, "src", "Figures.jl"))
    if isdefined(Main.Figures, :figures_deux_canaux)
        Base.invokelatest(Main.Figures.figures_deux_canaux, dossier, fig; budget = budget_vise, essai, jeu)
    end
end

# -----------------------------------------------------------------------------
# 7. L'échantillon inconnu : sa concentration par chaque lecture, avec ses incertitudes
# -----------------------------------------------------------------------------
if "inconnue" in etapes
    println("\n7. Échantillon inconnu ($(jeu)/unknown)")
    # courbes reconstruites depuis leurs nœuds (pentes de Hermite à jour, même non monotones)
    courbes = Dict(k => D.courbe(C.lecture, C.λ, C.c, C.o, C.σ_ref) for (k, C) in deserialize(joinpath(dossier, "courbes.jls")))
    gc = CSV.read(joinpath(dossier, "global_cellules.csv"), DataFrame)
    inconnues = D.charger_inconnue(dossier_donnees; λs)
    # tension rapport/forme des cellules de la titration, pour comparaison
    tension_titration = Dict(λ => [x.tension_rapport_forme for x in D.tester_melange(par_λ[λ], modeles[λ])] for λ in λs)
    tab = DataFrame()
    diag = DataFrame()
    for cel in inconnues
        λ = cel.λ
        ctx, mg = contextes[λ], modeles[λ]
        T = D.n_acq(cel) * cel.duree_s
        hB, hJ = D.somme(cel, D.BLEU), D.somme(cel, D.JAUNE)
        o, e = D.lire(ctx, hB, hJ, T, 1)
        uniques = [D.lire(ctx, D.acquisition(cel, D.BLEU, i), D.acquisition(cel, D.JAUNE, i), cel.duree_s, 1)[1] for i in 1:D.n_acq(cel)]
        # saut de t₀ entre cellules, mesuré sur la titration (écart-type des t₀ par point), appliqué à l'inconnu
        g = gc[gc.lambda_nm.==λ, :]
        σtB, σtJ = std(g.t0_bleu_par_cellule_ps) / 1000, std(g.t0_jaune_par_cellule_ps) / 1000
        lus(b, j) = D.lire(ctx, b, j, T, 1)[1]
        oBp, oBm = lus(D.decaler(hB, σtB, cel.pas_ns), hJ), lus(D.decaler(hB, -σtB, cel.pas_ns), hJ)
        oJp, oJm = lus(hB, D.decaler(hJ, σtJ, cel.pas_ns)), lus(hB, D.decaler(hJ, -σtJ, cel.pas_ns))
        # diagnostics : contrôles, t₀ et mélange à deux états de l'inconnu
        ct = D.controler(cel)
        r1 = D.ajuster_global([cel]; lien = mg.lien, τ_init = (mg.τ2, mg.τ3, mg.τA, mg.τJ), t0_init = (mg.t0B, mg.t0J), figer_tau = true)
        mel = D.tester_melange([cel], mg)[1]
        push!(diag, (lambda_nm = λ, N_bleu_par_s = ct.N_bleu / ct.duree_s, JB = ct.JB, sigma_JB_sur_poisson = ct.sigma_JB_sur_poisson,
            pente_JB_pct = ct.pente_JB_pct, t_JB = ct.t_JB, tendances = ct.tendances,
            t0_bleu_ps = 1000r1.t0B, t0_bleu_titration_ps = mean(g.t0_bleu_par_cellule_ps), sigma_t0_bleu_ps = 1000σtB,
            t0_jaune_ps = 1000r1.t0J, t0_jaune_titration_ps = mean(g.t0_jaune_par_cellule_ps), sigma_t0_jaune_ps = 1000σtJ,
            f_melange = mel.f, f_melange_commune = mel.f_commune, exces_melange = mel.chi2r - (r1.χ2B[1] + r1.χ2J[1]) / 2,
            tension_rapport_forme = mel.tension_rapport_forme, tension_titration_mediane = median(tension_titration[λ]),
            tension_titration_max = maximum(tension_titration[λ]),
            alpha_tau2_bleu = r1.aB[1, 2] / (r1.aB[1, 2] + r1.aB[1, 3]), diffusion_bleu_par_s = r1.IB[1, 1] / T))
        for (q, k) in enumerate(D.LECTURES)
            C = courbes[(λ, k)]
            v = filter(isfinite, [u[q] for u in uniques])
            σo = length(v) >= 3 ? std(v) / sqrt(length(v)) : NaN
            p = D.predire(C, o[q], σo)
            inv1(x) = (s1 = D.inverser_toutes(C, x); length(s1) == 1 ? s1[1] : NaN)
            σt0 = sqrt(((inv1(oBp[q]) - inv1(oBm[q])) / 2)^2 + ((inv1(oJp[q]) - inv1(oJm[q])) / 2)^2)
            h = D.hill(C.c, C.o, C.σ_ref)
            push!(tab, (lambda_nm = λ, lecture = k, nom = D.NOMS_LECTURES[k], monotone = C.monotone, observable = o[q], code = e[q],
                sigma_observable = σo, c_mM = p.c, solutions = join([@sprintf("%.2f", x) for x in p.solutions], " ; "),
                mediane_mM = p.mediane, bas68_mM = p.bas68, haut68_mM = p.haut68, bas95_mM = p.bas95, haut95_mM = p.haut95,
                sigma_stat_mM = p.sigma, sigma_mesure_mM = p.sigma_mesure, sigma_calibration_mM = p.sigma_calibration,
                sigma_t0_mM = σt0, sigma_interpolation_mM = abs(p.c - p.c_lineaire) / 2,
                sigma_total_mM = sqrt(p.sigma^2 + σt0^2 + (abs(p.c - p.c_lineaire) / 2)^2), c_lineaire_mM = p.c_lineaire,
                c_hill_mM = D.inverser_hill(h, o[q]), chi2r_hill = h.chi2r,
                hors_plage_mc = p.hors_plage_mc))
        end
    end
    # une valeur par lecture : moyenne des trois longueurs d'onde (mesures et calibrations indépendantes),
    # pondérée par 1/σ_total² ; σ gonflé par √χ²ᵣ si les longueurs d'onde ne s'accordent pas
    comb = DataFrame()
    for k in D.LECTURES
        d = tab[(tab.lecture.==k).&tab.monotone.&isfinite.(tab.c_mM).&isfinite.(tab.sigma_total_mM), :]
        nrow(d) == 0 && (push!(comb, (lecture = k, nom = D.NOMS_LECTURES[k], n_lambda = 0, c_mM = NaN, sigma_mM = NaN, chi2r = NaN,
            sigma_gonfle_mM = NaN, lambdas = "")); continue)
        w = 1 ./ d.sigma_total_mM .^ 2
        m = sum(w .* d.c_mM) / sum(w)
        σ = 1 / sqrt(sum(w))
        χ2r = nrow(d) > 1 ? sum(w .* (d.c_mM .- m) .^ 2) / (nrow(d) - 1) : NaN
        push!(comb, (lecture = k, nom = D.NOMS_LECTURES[k], n_lambda = nrow(d), c_mM = m, sigma_mM = σ, chi2r = χ2r,
            sigma_gonfle_mM = σ * max(1.0, isfinite(χ2r) ? sqrt(χ2r) : 1.0), lambdas = join(string.(d.lambda_nm), " ")))
    end
    CSV.write(joinpath(dossier, "inconnue.csv"), tab)
    CSV.write(joinpath(dossier, "inconnue_combinee.csv"), comb)
    CSV.write(joinpath(dossier, "inconnue_diagnostics.csv"), diag)
    for r in eachrow(diag)
        @printf("  %d nm : J/B %.4f (%.2f × Poisson, pente %+.3f %%/acq, t %+.1f) ; t₀ bleu %.1f ps (titration %.1f ± %.1f), jaune %.1f ps (%.1f ± %.1f)\n",
            r.lambda_nm, r.JB, r.sigma_JB_sur_poisson, r.pente_JB_pct, r.t_JB, r.t0_bleu_ps, r.t0_bleu_titration_ps, r.sigma_t0_bleu_ps,
            r.t0_jaune_ps, r.t0_jaune_titration_ps, r.sigma_t0_jaune_ps)
        @printf("          mélange : f %.3f (forme seule) contre %.3f (J/B imposé) ; tension rapport/forme %.0f (titration : médiane %.0f, max %.0f ; seuil 3,84) ; diffusion %.0f photons/s\n",
            r.f_melange, r.f_melange_commune, r.tension_rapport_forme, r.tension_titration_mediane, r.tension_titration_max, r.diffusion_bleu_par_s)
    end
    for r in eachrow(tab)
        @printf("  %d nm %-38s %s %7.2f mM  [68 %% : %6.2f ; %6.2f]  σ mesure %.2f, calib. %.2f, t₀ %.2f, interp. %.2f → total %.2f  (linéaire %.2f, Hill %.2f)\n",
            r.lambda_nm, r.nom, r.monotone ? " " : "*", r.c_mM, r.bas68_mM, r.haut68_mM, r.sigma_mesure_mM,
            r.sigma_calibration_mM, r.sigma_t0_mM, r.sigma_interpolation_mM, r.sigma_total_mM, r.c_lineaire_mM, r.c_hill_mM)
    end
    println("  combinées sur les longueurs d'onde (lectures monotones) :")
    for r in eachrow(comb)
        @printf("    %-38s %7.2f ± %.2f mM (χ²ᵣ %.2f, ± %.2f gonflé ; %s nm)\n", r.nom, r.c_mM, r.sigma_mM, r.chi2r, r.sigma_gonfle_mM, r.lambdas)
    end
    isdefined(Main, :Figures) || include(joinpath(@__DIR__, "src", "Figures.jl"))
    if isdefined(Main.Figures, :figure_dc_inconnue)
        Base.invokelatest(Main.Figures.figure_dc_inconnue, tab, comb, fig("fig23_deux_canaux_inconnue.png"); essai, jeu)
    end
end

# -----------------------------------------------------------------------------
# Rejouer un tirage seul
# -----------------------------------------------------------------------------
if rejeu
    λ, cl, N, j = parse(Int, args[2]), parse(Float64, args[3]), parse(Float64, args[4]), parse(Int, args[5])
    src = length(args) >= 6 && args[6] in ("reel", "jumeau") ? Symbol(args[6]) : :reel
    k0 = length(args) >= 6 && args[6] in ("reel", "jumeau") ? 7 : 6
    sc = length(args) >= k0 ? Symbol(args[k0]) : :aucun
    β = length(args) >= k0 + 1 ? parse(Float64, args[k0+1]) : 0.0
    cel = only(filter(c -> c.cl == cl, par_λ[λ]))
    r = D.rejouer(cel, modeles[λ], contextes[λ], N, j; source = src, scenario = sc, taux = β)
    @printf("Tirage %d nm, %g mM, %g photons, n° %d (%s, %s %g) : p = %.4f, %d photons bleus, %d jaunes\n", λ, cl, N, j, src, sc, β, r.p,
        sum(r.bleu[D.FENETRE]), sum(r.jaune[D.FENETRE]))
    for k in D.LECTURES
        o, e = r.lectures[k]
        @printf("  %-38s %12.5f  code %d\n", D.NOMS_LECTURES[k], o, e)
    end
end
