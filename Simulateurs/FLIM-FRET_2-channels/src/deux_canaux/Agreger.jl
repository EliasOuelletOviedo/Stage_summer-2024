# Agreger.jl — étape 7 : des estimations aux décisions.
#
# Chaque tirage est inversé par la courbe de sa lecture (mM), puis résumé par cellule de
# la grille (λ, concentration, lecture, budget, scénario, taux) : biais, σ corrigé de
# l'amincissement, RMSE, taux d'échec. Viennent ensuite les deux règles de décision,
# écrites avant de voir les courbes (document du 2026-10-16, étape 5), la carte finale,
# le blanchiment tolérable, l'efficacité de chaque lecture et l'accord réel/jumeau.

const FENETRE_CLAMPING = (5.0, 35.0)
const INTERIEUR = (2.5, 80.0)
const BUDGET_VISE = 1e4
const SEUIL_ECHEC = 0.01
const SEUIL_BIAIS_MM = 1.0

dans(c, intervalle) = intervalle[1] <= c <= intervalle[2]

"""Une ligne par tirage et par lecture : concentration estimée (mM) et code final."""
function inverser_tirages(est::DataFrame, courbes::Dict{Tuple{Int,String},Courbe})
    n = nrow(est)
    parts = DataFrame[]
    for k in LECTURES
        ĉ = fill(NaN, n)
        code = fill(Int8(OK), n)
        o, e = est[!, Symbol("o_" * k)], est[!, Symbol("e_" * k)]
        for i in 1:n
            C = courbes[(est.lambda_nm[i], k)]
            if !C.monotone
                code[i] = NON_MONOTONE
            elseif e[i] != OK
                code[i] = e[i]
            else
                v, c = inverser(C, o[i])
                ĉ[i], code[i] = v, c
            end
        end
        push!(parts, DataFrame(source = est.source, lambda_nm = est.lambda_nm, cl_mM = est.cl_mM, budget_photons = est.budget_photons,
            p = est.p, scenario = est.scenario, taux = est.taux, tirage = est.tirage, lecture = fill(k, n), c_hat = ĉ, code = code))
    end
    return reduce(vcat, parts)
end

"""
    statistiques(long, courbes) -> DataFrame

Par cellule de la grille : biais = moyenne des estimations réussies − concentration (la
valeur à pleins photons, exacte aux nœuds de la calibration) ; σ_tirages ; σ = σ_tirages
/√(1 − p) pour les tirages réels (le jumeau, poissonien, n'est pas corrigé) ; RMSE ;
taux d'échec = non-convergence + paramètre sur une borne + hors de la plage calibrée.
`sigma_total` ajoute l'incertitude de la calibration (8 points à photons finis) : c'est
la vraie ; `sigma` isole le bruit du tirage.
"""
function statistiques(long::DataFrame, courbes::Dict{Tuple{Int,String},Courbe})
    σcal = Dict{Tuple{Int,String,Float64},Float64}()
    for ((λ, k), C) in courbes
        for (i, c) in enumerate(C.c)
            σcal[(λ, k, c)] = sigma_calibration(C)[i]
        end
    end
    out = DataFrame()
    for g in groupby(long, [:source, :lambda_nm, :cl_mM, :lecture, :budget_photons, :scenario, :taux]; sort = true)
        n = nrow(g)
        c = g.cl_mM[1]
        p = g.p[1]
        codes = g.code
        mono = !any(==(NON_MONOTONE), codes)
        n_aj = count(in((NON_CONVERGENCE, BORNE, AUTRE)), codes)
        n_hp = count(==(HORS_PLAGE), codes)
        ĉ = g.c_hat[codes.==OK]
        biais = isempty(ĉ) ? NaN : mean(ĉ) - c
        σt = length(ĉ) >= 2 ? std(ĉ) : NaN
        σ = g.source[1] == "reel" ? sigma_vraie(σt, p) : σt
        σc = get(σcal, (g.lambda_nm[1], g.lecture[1], c), NaN)
        σtot = sqrt(σ^2 + σc^2)
        mono || (biais = σt = σ = σtot = NaN)
        push!(out, (lambda_nm = g.lambda_nm[1], cl_mM = c, lecture = g.lecture[1], budget_photons = g.budget_photons[1],
            scenario_blanchiment = g.scenario[1], taux_blanchiment = g.taux[1], biais_mM = biais, sigma_mM = σ,
            rmse_mM = sqrt(biais^2 + σ^2), taux_echec = mono ? (n_aj + n_hp) / n : NaN, n_tirages = n,
            source = g.source[1], sigma_tirages_mM = σt, p_amincissement = p, sigma_calibration_mM = σc,
            sigma_total_mM = σtot, rmse_total_mM = sqrt(biais^2 + σtot^2), taux_echec_ajustement = n_aj / n,
            taux_hors_plage = n_hp / n, n_reussis = length(ĉ), monotone = mono))
    end
    return out
end

const EN_TETE_FICHIER_PLAT = [
    "Analyse FLIM à deux canaux — titration du 2026-10-16 (step15_analyse_deux_canaux.jl). Une ligne par cellule de la grille.",
    "Vérité de référence : la somme des 10 acquisitions de chaque cellule. Elle mesure la reproductibilité, pas la justesse : un biais commun à la référence et aux tirages est invisible ici.",
    "Les valeurs absolues (mM, K apparents) ne sont pas citables si le capteur est dégradé ; le classement des lectures et les rapports de σ le sont.",
    "sigma_mM : écart-type des estimations, corrigé de l'amincissement (σ_tirages/√(1−p)) ; sigma_total_mM et rmse_total_mM ajoutent l'incertitude de la calibration (la vraie).",
    "taux_echec = non-convergence + paramètre sur une borne + observable hors de la plage calibrée (aux deux bornes de la titration, la moitié des tirages sort par construction).",
    "scenario_blanchiment : aucun = tirages intacts ; A, B, C = blanchiment injecté à taux_blanchiment (fraction blanchie pendant l'acquisition), voir le README.",
    "drapeau_poisson : J/B disperse à plus de 2 × Poisson d'une acquisition à l'autre ; drapeau_tendance : tendance > 2σ sur les 10 acquisitions (intensités ou J/B)."]

"""Écrit un CSV précédé de lignes de commentaire (« # »)."""
function ecrire_avec_en_tete(chemin::AbstractString, df::DataFrame, lignes)
    open(chemin, "w") do io
        for l in lignes
            println(io, "# ", l)
        end
        CSV.write(io, df; append = true, writeheader = true)
    end
    return chemin
end

"""Le fichier plat : tirages réels seulement, drapeaux des contrôles joints."""
function fichier_plat(stats::DataFrame, controles::DataFrame)
    f = stats[stats.source.=="reel", Not(:source)]
    drap = select(controles, :lambda_nm, :cl_mM, :drapeau_poisson, :drapeau_tendance)
    out = leftjoin(f, drap; on = [:lambda_nm, :cl_mM])
    ordre_sc = Dict("aucun" => 0, "A" => 1, "B" => 2, "C" => 3)
    out.ordre_lecture = [findfirst(==(k), LECTURES) for k in out.lecture]
    out.ordre_scenario = [ordre_sc[s] for s in out.scenario_blanchiment]
    sort!(out, [:lambda_nm, :cl_mM, :ordre_lecture, :budget_photons, :ordre_scenario, :taux_blanchiment])
    return select!(out, Not([:ordre_lecture, :ordre_scenario]))
end

"""Réel contre jumeau : rapport des σ (tirages intacts)."""
function jumeau_contre_reel(stats::DataFrame)
    s = stats[stats.scenario_blanchiment.=="aucun", :]
    r = s[s.source.=="reel", [:lambda_nm, :cl_mM, :lecture, :budget_photons, :sigma_mM, :biais_mM, :taux_echec]]
    j = s[s.source.=="jumeau", [:lambda_nm, :cl_mM, :lecture, :budget_photons, :sigma_mM, :biais_mM, :taux_echec]]
    rename!(r, :sigma_mM => :sigma_reel_mM, :biais_mM => :biais_reel_mM, :taux_echec => :echec_reel)
    rename!(j, :sigma_mM => :sigma_jumeau_mM, :biais_mM => :biais_jumeau_mM, :taux_echec => :echec_jumeau)
    out = innerjoin(r, j; on = [:lambda_nm, :cl_mM, :lecture, :budget_photons])
    out.rapport_sigma = out.sigma_reel_mM ./ out.sigma_jumeau_mM
    return out
end

"""Efficacité de chaque lecture : (borne de Fisher de l'information qu'elle exploite / σ)²."""
function efficacites(flat::DataFrame, fisher::DataFrame)
    f = flat[flat.scenario_blanchiment.=="aucun", [:lambda_nm, :cl_mM, :lecture, :budget_photons, :sigma_mM, :rmse_mM]]
    out = innerjoin(f, fisher; on = [:lambda_nm, :cl_mM, :budget_photons])
    out.information = [INFORMATION_LECTURE[k] for k in out.lecture]
    out.borne_mM = [r[Symbol(r.information)] for r in eachrow(out)]
    out.efficacite = (out.borne_mM ./ out.sigma_mM) .^ 2
    return select(out, :lambda_nm, :cl_mM, :lecture, :budget_photons, :information, :borne_mM, :sigma_mM, :efficacite,
        :rapport, :toute, :forme, :bleu, :jaune)
end

"""
Biais en mM par point de pourcentage blanchi : pente des moindres carrés du biais contre
le taux (0, 1, 2, 5, 10 %), par lecture, cellule, budget et scénario.
"""
function pentes_blanchiment(flat::DataFrame)
    base = flat[flat.scenario_blanchiment.=="aucun", :]
    bl = flat[flat.scenario_blanchiment.!="aucun", :]
    out = DataFrame()
    for g in groupby(bl, [:lambda_nm, :cl_mM, :lecture, :budget_photons, :scenario_blanchiment]; sort = true)
        b0 = base[(base.lambda_nm.==g.lambda_nm[1]).&(base.cl_mM.==g.cl_mM[1]).&(base.lecture.==g.lecture[1]).&(base.budget_photons.==g.budget_photons[1]), :]
        nrow(b0) == 1 || continue
        x = vcat(0.0, 100 .* g.taux_blanchiment)
        y = vcat(b0.biais_mM[1], g.biais_mM)
        ok = isfinite.(y)
        pente = count(ok) >= 3 ? cov(x[ok], y[ok]) / var(x[ok]) : NaN
        push!(out, (lambda_nm = g.lambda_nm[1], cl_mM = g.cl_mM[1], lecture = g.lecture[1], budget_photons = g.budget_photons[1],
            scenario_blanchiment = g.scenario_blanchiment[1], pente_mM_par_pct = pente,
            biais_a_10pct_mM = g.biais_mM[findfirst(==(0.10), g.taux_blanchiment)] - b0.biais_mM[1]))
    end
    return out
end

"""
Taux de blanchiment tolérable pour un biais sous 1 mM, par lecture, longueur d'onde et
scénario : 1 mM / (plus forte pente sur 5–35 mM) au budget visé ; aussi sur 2,5–80 mM.
"""
function blanchiment_tolerable(pentes::DataFrame; budget = BUDGET_VISE)
    out = DataFrame()
    p = pentes[pentes.budget_photons.==budget, :]
    for g in groupby(p, [:lambda_nm, :lecture, :scenario_blanchiment]; sort = true)
        f = g[dans.(g.cl_mM, Ref(FENETRE_CLAMPING)), :]
        i = g[dans.(g.cl_mM, Ref(INTERIEUR)), :]
        pmax(v) = (v = filter(isfinite, abs.(v)); isempty(v) ? NaN : maximum(v))
        push!(out, (lambda_nm = g.lambda_nm[1], lecture = g.lecture[1], scenario_blanchiment = g.scenario_blanchiment[1],
            budget_photons = budget, pente_max_5_35_mM_par_pct = pmax(f.pente_mM_par_pct),
            taux_tolerable_5_35_pct = SEUIL_BIAIS_MM / pmax(f.pente_mM_par_pct),
            pente_max_interieur_mM_par_pct = pmax(i.pente_mM_par_pct), taux_tolerable_interieur_pct = SEUIL_BIAIS_MM / pmax(i.pente_mM_par_pct)))
    end
    return out
end

"""
    regle_lecture(flat; colonne_biais, colonne_rmse) -> DataFrame

La règle de choix de la lecture, écrite avant de voir les courbes. À une longueur d'onde
et un budget, trois critères dans l'ordre, chacun éliminatoire :
  1. taux d'échec < 1 % à chaque concentration jugée ;
  2. |biais| < 1 mM à chaque concentration jugée ;
  3. la plus petite RMSE médiane (vraie, calibration comprise) sur 5–35 mM, la fenêtre de
     clamping (nœuds 5, 10 et 20 mM), le pire cas rapporté à côté.
Le gagnant est la meilleure RMSE parmi les survivants. Concentrations jugées : 2,5 à
80 mM. Aux deux bornes de la calibration (0 et 131 mM), la moitié des tirages sort de
l'intervalle calibré par construction, quelle que soit la lecture ; elles sont
rapportées dans le fichier plat, pas jugées.
"""
function regle_lecture(flat::DataFrame; colonne_biais::Symbol = :biais_mM, colonne_rmse::Symbol = :rmse_total_mM)
    f = flat[flat.scenario_blanchiment.=="aucun", :]
    out = DataFrame()
    for g in groupby(f, [:lambda_nm, :budget_photons]; sort = true)
        lignes = NamedTuple[]
        for h in groupby(g, :lecture; sort = false)
            di = h[dans.(h.cl_mM, Ref(INTERIEUR)), :]
            d5 = h[dans.(h.cl_mM, Ref(FENETRE_CLAMPING)), :]
            nrow(di) == 0 && continue
            echec = maximum(di.taux_echec)
            ie = argmax(replace(di.taux_echec, NaN => -1.0))
            biais = maximum(abs.(di[!, colonne_biais]))
            med = median(d5[!, colonne_rmse])
            k = argmax(replace(di[!, colonne_rmse], NaN => Inf))
            pire = di[k, colonne_rmse]
            elim = !all(h.monotone) ? "non monotone" : !(echec < SEUIL_ECHEC) ? "échecs" : !(biais < SEUIL_BIAIS_MM) ? "biais" :
                   !isfinite(med) ? "RMSE indéfinie" : ""
            push!(lignes, (lambda_nm = g.lambda_nm[1], budget_photons = g.budget_photons[1], lecture = h.lecture[1],
                echec_max = echec, echec_max_cl_mM = di.cl_mM[ie], echec_max_hors_plage = di.taux_hors_plage[ie],
                biais_max_mM = biais, rmse_mediane_5_35_mM = med, pire_cas_mM = pire,
                pire_cas_cl_mM = di.cl_mM[k], survit = isempty(elim), eliminee_par = elim))
        end
        d = DataFrame(lignes)
        d.rang = zeros(Int, nrow(d))
        s = findall(d.survit)
        for (r, i) in enumerate(s[sortperm(d.rmse_mediane_5_35_mM[s])])
            d.rang[i] = r
        end
        d.gagnante = d.rang .== 1
        append!(out, d)
    end
    return out
end

"""
Second classement, sous blanchiment au taux mesuré. Le scénario d'une longueur d'onde est
le plus fréquent parmi ses cellules (« aucun » si aucune ne dérive ; le pire des trois,
lecture par lecture, s'il est indéterminé). Pendant une mesure de N photons, la fraction
blanchie vaut taux_mesuré × N / débit bleu de la cellule ; le biais ajouté est cette
fraction (en %) × la pente en mM par point de pourcentage.
"""
function regle_lecture_blanchie(flat::DataFrame, pentes::DataFrame, signatures::DataFrame, controles::DataFrame)
    f = copy(flat[flat.scenario_blanchiment.=="aucun", :])
    f.biais_blanchi_mM = copy(f.biais_mM)
    f.rmse_blanchie_mM = copy(f.rmse_total_mM)
    scen = Dict{Int,String}()
    for g in groupby(signatures, :lambda_nm)
        v = filter(!=("aucun"), g.scenario)
        scen[g.lambda_nm[1]] = isempty(v) ? "aucun" : (cm = Dict(s => count(==(s), v) for s in unique(v)); first(sort(collect(cm); by = x -> -x[2]))[1])
    end
    for (i, r) in enumerate(eachrow(f))
        sg = signatures[(signatures.lambda_nm.==r.lambda_nm).&(signatures.cl_mM.==r.cl_mM), :]
        ct = controles[(controles.lambda_nm.==r.lambda_nm).&(controles.cl_mM.==r.cl_mM), :]
        (nrow(sg) == 1 && nrow(ct) == 1) || continue
        k_bl = sg.taux_mesure_par_s[1]
        k_bl > 0 || continue
        β_pct = 100 * k_bl * r.budget_photons / (ct.N_bleu[1] / ct.duree_s[1])
        sc = scen[r.lambda_nm]
        cand = pentes[(pentes.lambda_nm.==r.lambda_nm).&(pentes.cl_mM.==r.cl_mM).&(pentes.lecture.==r.lecture).&(pentes.budget_photons.==r.budget_photons), :]
        nrow(cand) == 0 && continue
        pente = sc in ("A", "B", "C") ? cand.pente_mM_par_pct[findfirst(==(sc), cand.scenario_blanchiment)] :
                cand.pente_mM_par_pct[argmax(abs.(replace(cand.pente_mM_par_pct, NaN => 0.0)))]
        f.biais_blanchi_mM[i] = r.biais_mM + β_pct * pente
        f.rmse_blanchie_mM[i] = sqrt(f.biais_blanchi_mM[i]^2 + r.sigma_total_mM^2)
    end
    regle = regle_lecture(f; colonne_biais = :biais_blanchi_mM, colonne_rmse = :rmse_blanchie_mM)
    regle.scenario_blanchiment = [scen[λ] for λ in regle.lambda_nm]
    return regle
end

"""RMSE interpolée en log-log entre les budgets simulés (NaN hors de la grille)."""
function rmse_a(Ns::Vector{Float64}, r::Vector{Float64}, N::Real)
    ok = isfinite.(r) .& (r .> 0)
    Ns, r = Ns[ok], r[ok]
    length(Ns) >= 2 || return NaN
    (N < Ns[1] || N > Ns[end]) && return NaN
    i = clamp(searchsortedlast(Ns, N), 1, length(Ns) - 1)
    t = (log(N) - log(Ns[i])) / (log(Ns[i+1]) - log(Ns[i]))
    return exp(log(r[i]) + t * (log(r[i+1]) - log(r[i])))
end

"""
    regle_longueur_onde(flat, regles, controles, references, globaux; budget, puissances_mW) -> DataFrame

Trois quantités par longueur d'onde : la plage dynamique de J/B (0 → 131 mM), la RMSE
médiane sur 5–35 mM à photons bleus égaux (le budget visé) et à temps égal (la durée qui
donne le budget visé au débit bleu médian de la titration ; chaque cellule reçoit alors
débit × durée photons, RMSE interpolée entre les budgets), et les photons bleus par
seconde, par milliwatt si la puissance est mesurée. La comparaison à photons égaux ne
dépend d'aucune hypothèse ; à temps égal, elle dépend de la puissance laser.
"""
function regle_longueur_onde(flat::DataFrame, regles::DataFrame, controles::DataFrame, references::DataFrame, globaux::DataFrame;
                             budget = BUDGET_VISE, puissances_mW = Dict{Int,Float64}())
    f = flat[flat.scenario_blanchiment.=="aucun", :]
    debit = controles.N_bleu ./ controles.duree_s
    t_vise = budget / median(debit)
    out = DataFrame()
    for λ in sort(unique(f.lambda_nm))
        rg = regles[(regles.lambda_nm.==λ).&(regles.budget_photons.==budget), :]
        g = rg[rg.gagnante, :]
        gag = nrow(g) == 1 ? g.lecture[1] : "aucune"
        survivantes = rg.lecture[rg.survit]
        ref = references[(references.lambda_nm.==λ).&(references.lecture.=="1"), :]
        sort!(ref, :cl_mM)
        ct = controles[controles.lambda_nm.==λ, :]
        # à temps égal : RMSE médiane 5–35 mM de chaque lecture survivante
        temps = Dict{String,Float64}()
        for k in survivantes
            vals = Float64[]
            for c in sort(unique(f.cl_mM[(f.lambda_nm.==λ)]))
                dans(c, FENETRE_CLAMPING) || continue
                d = sort(f[(f.lambda_nm.==λ).&(f.cl_mM.==c).&(f.lecture.==k), :], :budget_photons)
                r_c = ct[ct.cl_mM.==c, :]
                nrow(r_c) == 1 || continue
                push!(vals, rmse_a(d.budget_photons, d.rmse_total_mM, r_c.N_bleu[1] / r_c.duree_s[1] * t_vise))
            end
            temps[k] = median(vals)
        end
        best_t = isempty(temps) ? "aucune" : first(sort(collect(temps); by = last))[1]
        gl = globaux[globaux.lambda_nm.==λ, :]
        P = get(puissances_mW, λ, NaN)
        d20 = ct[ct.cl_mM.==20.0, :]
        debit20 = nrow(d20) == 1 ? d20.N_bleu[1] / d20.duree_s[1] : NaN
        push!(out, (lambda_nm = λ, plage_JB = ref.o[1] / ref.o[end], gagnante_photons_egaux = gag,
            rmse_mediane_photons_egaux_mM = nrow(g) == 1 ? g.rmse_mediane_5_35_mM[1] : NaN,
            duree_temps_egal_s = t_vise, gagnante_temps_egal = best_t, rmse_mediane_temps_egal_mM = get(temps, best_t, NaN),
            photons_bleus_par_s_20mM = debit20, photons_bleus_par_s_median = median(ct.N_bleu ./ ct.duree_s),
            puissance_mW = P, photons_bleus_par_s_par_mW = debit20 / P, photons_bleus_par_s_par_mW2 = debit20 / P^2,
            jaune_diffusion_0mM = gl.intensite_diffusion_jaune[1], jaune_rapide_0mM = gl.intensite_tauJ_jaune[1],
            jaune_accepteur_0mM = gl.intensite_tauA_jaune[1], cmax_mM = gl.cl_mM[end], jaune_diffusion_cmax = gl.intensite_diffusion_jaune[end],
            jaune_rapide_cmax = gl.intensite_tauJ_jaune[end], jaune_accepteur_cmax = gl.intensite_tauA_jaune[end]))
    end
    # aucune gagnante à aucune longueur d'onde : aucune n'est « meilleure »
    out.meilleure_photons_egaux = out.rmse_mediane_photons_egaux_mM .== minimum(vcat(filter(isfinite, out.rmse_mediane_photons_egaux_mM), Inf))
    out.meilleure_temps_egal = out.rmse_mediane_temps_egal_mM .== minimum(vcat(filter(isfinite, out.rmse_mediane_temps_egal_mM), Inf))
    return out
end

"""
La carte finale : pour chaque (λ, concentration, budget), la lecture de plus petite RMSE
vraie parmi celles qui échouent moins d'une fois sur cent et dont le biais reste sous
1 mM à ce point. Aux bornes de la calibration, seuls les échecs d'ajustement comptent
(la sortie de plage y est structurelle).
"""
function carte_finale(flat::DataFrame)
    f = flat[flat.scenario_blanchiment.=="aucun", :]
    bornes = extrema(f.cl_mM)
    out = DataFrame()
    for g in groupby(f, [:lambda_nm, :cl_mM, :budget_photons]; sort = true)
        ech = g.cl_mM[1] in bornes ? g.taux_echec_ajustement : g.taux_echec
        ok = g.monotone .& (ech .< SEUIL_ECHEC) .& (abs.(g.biais_mM) .< SEUIL_BIAIS_MM) .& isfinite.(g.rmse_total_mM)
        if any(ok)
            i = findall(ok)[argmin(g.rmse_total_mM[ok])]
            push!(out, (lambda_nm = g.lambda_nm[1], cl_mM = g.cl_mM[1], budget_photons = g.budget_photons[1], gagnante = g.lecture[i],
                rmse_total_mM = g.rmse_total_mM[i], n_candidates = count(ok)))
        else
            push!(out, (lambda_nm = g.lambda_nm[1], cl_mM = g.cl_mM[1], budget_photons = g.budget_photons[1], gagnante = "aucune",
                rmse_total_mM = NaN, n_candidates = 0))
        end
    end
    return out
end
