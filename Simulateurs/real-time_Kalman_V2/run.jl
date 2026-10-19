# run.jl — étape 2, la boucle fermée : les sept courses du plan, dans l'ordre, puis les critères
# de sortie et les quatre figures. Chaque course prend ses réglages dans les verdicts des
# précédentes (β retenu à la course 2, calibration de P₁₁ à la course 3, coude à la course 5).
#
# Lancer depuis ce dossier :  julia --project=. -t auto run.jl  [--rapide]
# (--rapide : 2 répétitions partout, pour vérifier que tout tourne)

using Printf, Statistics, LinearAlgebra, StableRNGs, Base.Threads
for f in ("params", "generate", "estimate", "physio", "kalman", "regulateur", "boucle")
    include(joinpath(@__DIR__, "src", f * ".jl"))
end

const SORTIES = joinpath(@__DIR__, "sorties")

# ----------------------------------------------------------------------------- masques et métriques

moy(x) = (x = filter(isfinite, x); isempty(x) ? NaN : mean(x))
fenetres(t, ws) = reduce(.|, [(w[1] .<= t .< w[2]) for w in ws])
masque_phase(t, p, nom) = (ph = p.toml["scenario"]["phases"][findfirst(x -> x[1] == nom, p.toml["scenario"]["phases"])]; ph[2] .<= t .< ph[3])
noms_phases(p) = [x[1] for x in p.toml["scenario"]["phases"]]

"""Transitoires et plateaux d'après la pente du vrai C (seuils de params.toml)."""
function transitoires(com, p)
    s = p.toml["scenario"]
    d = abs.([diff(com[:C_vrai]) ./ p.dt; NaN])
    return d .> s["seuil_transitoire_mM_par_s"], d .< s["seuil_plateau_mM_par_s"]
end

"""Autocorrélation à un pas, sur les seules paires (k, k+1) toutes deux dans le masque."""
function autocorr(x, m)
    i = [k for k in 1:length(x)-1 if m[k] && m[k+1] && isfinite(x[k]) && isfinite(x[k+1])]
    length(i) < 5 && return NaN
    a, b = x[i], x[i.+1]
    return cor(a, b)
end

"""Métriques d'un filtre sur une course."""
function metriques(com, T, p)
    t, C, cons = com[:t], com[:C_vrai], com[:consigne]
    e = T[:C_hat] .- C
    tr, pl = transitoires(com, p)
    et = fenetres(t, p.toml["scenario"]["etabli"])
    fin = t .>= t[end] + p.dt - p.toml["course"]["fin_biais_s"]
    ph = Dict(n => masque_phase(t, p, n) for n in noms_phases(p))
    νR = T[:innov_ratio]
    return merge((rmse = sqrt(moy(e .^ 2)), biais_30s = moy(e[fin]), rmse_plateaux = sqrt(moy(e[pl] .^ 2)),
        P11_sur_err2 = moy(T[:P11]) / moy(e .^ 2)),
        (; (Symbol("biais_", n) => moy(e[ph[n]]) for n in noms_phases(p))...),
        (; (Symbol("P11_err2_", n) => moy(T[:P11][ph[n]]) / moy(e[ph[n]] .^ 2) for n in noms_phases(p))...),
        (innov_transitoires = moy(νR[tr]), innov_plateaux = moy(νR[pl]),
        autocorr_transitoires = autocorr(νR, tr), autocorr_plateaux = autocorr(νR, pl),
        innov_5a_transitoires = moy(T[:innov_5a][tr]),
        erreur_regulation = moy(cons[et] .- T[:C_hat][et]), erreur_reelle = moy(cons[et] .- C[et]),
        rho12_etabli = moy(T[:rho12][et]), biais_g_30s = moy(T[:g_hat][fin] .- com[:g_vrai][fin])))
end

"""Les quatre signatures du couplage filtre–régulateur, sur les fenêtres de régime établi."""
function signatures(com, T, p)
    t, u = com[:t], com[:u]
    ws = p.toml["scenario"]["etabli"]
    lag = round(Int, 1.0 / p.dt)
    ag, resp, cors = Float64[], Float64[], Float64[]
    du = [NaN; diff(u)]
    for w in ws
        i = findall(w[1] .<= t .< w[2])
        push!(ag, std(du[i[2:end]]))
        P = T[:P11][i]
        push!(resp, (maximum(P) - minimum(P)) / mean(P))
    end
    a = abs.(du)
    for L in -lag:lag
        pairs = [(T[:P11][k], a[k+L]) for w in ws for k in findall(w[1] .<= t .< w[2]) if 1 <= k + L <= length(t) && w[1] <= t[k+L] < w[2] && isfinite(a[k+L])]
        x, y = first.(pairs), last.(pairs)
        push!(cors, std(x) > 0 && std(y) > 0 ? cor(x, y) : 0.0)
    end
    return (agitation = mean(ag), correlation_P11_agitation = cors[argmax(abs.(cors))], respiration_P11 = mean(resp),
            rho12_etabli = moy(T[:rho12][fenetres(t, ws)]))
end

"""Moyenne de NamedTuples (champ par champ)."""
moyenne(M) = NamedTuple{keys(M[1])}(Tuple(moy([getfield(x, k) for x in M]) for k in keys(M[1])))
ecart_type(M, k) = std(filter(isfinite, [getfield(x, k) for x in M]))

function repeter(f, n)
    out = Vector{Any}(undef, n)
    @threads for r in 1:n
        out[r] = f(r)
    end
    return out
end

function ecrire_csv(chemin, lignes)
    isempty(lignes) && return
    open(chemin, "w") do io
        println(io, join(string.(keys(first(lignes))), ","))
        for l in lignes
            println(io, join((x isa Number ? @sprintf("%.6g", x) : string(x) for x in l), ","))
        end
    end
end

nom_derive(forme, α, p) = forme == "saut" ? @sprintf("saut %g %%", 100p.derive.saut) : α == 0 ? "aucune" : @sprintf("rampe %g %%", 100α)
derives(p) = [[("rampe", α) for α in p.toml["derive"]["balayage_alpha"]]; [("saut", 0.0)]]
regulateur(p, ph; Kp = Kp_nominal(p, ph), kw...) = PI(Kp, Kp * p.toml["regulateur"]["Ki_sur_Kp_par_s"]; kw...)

"""Fisher contre Monte-Carlo à C_reference (repris de l'étape 1), et échecs de 5a à 5 mM."""
function verifier_R(p::Params, lec)
    Cr, M = p.toml["lectures"]["C_reference_mM"], p.toml["course"]["verification_R_tirages"]
    tir(C, i) = histogrammes_c(C, 1.0, p; rng = StableRNG(p.toml["course"]["graine"] + 5_000_000 + i))
    L = filter(l -> l[1].ok && l[2].ok, map(i -> (h = tir(Cr, i); (estimer(h..., lec.JB.cfg), estimer(h..., lec.S.cfg))), 1:M))
    Z = hcat([l[1].z for l in L], [l[2].z for l in L])
    Ce, Rf = cov(Z), mean(R_jointe(l...) for l in L)
    e5 = count(i -> isnan(inverser(lec.S.courbe, estimer(tir(5.0, i)..., lec.S.cfg).z)), 1:M) / M
    return (N_photons_bleus = p.N, var_ratio_reelle = Ce[1, 1], var_ratio_fisher = Rf[1, 1], var_5a_reelle = Ce[2, 2],
            var_5a_fisher = Rf[2, 2], cov_reelle = Ce[1, 2], cov_fisher = Rf[1, 2], echecs_5a_5mM_pct = 100 * e5)
end

include(joinpath(@__DIR__, "src", "figures.jl"))

# ----------------------------------------------------------------------------- les sept courses

function main(; rapide = "--rapide" in ARGS)
    mkpath(SORTIES)
    p0 = charger_params()
    np = rapide ? 2 : p0.toml["course"]["repetitions_performance"]
    ns = rapide ? 2 : p0.toml["course"]["repetitions_statistiques"]
    K = p0.toml["kalman"]
    Ns = p0.toml["photons"]["balayage"]
    βs = K["balayage_beta"]
    σs = K["balayage_sigma_C0_diagnostic"]
    ph0 = physio(p0)
    lecs = Dict(N => lecteurs(avec(p0; N)) for N in Ns)
    verif = [verifier_R(avec(p0; N), lecs[N]) for N in Ns]
    nomβ(β, σ = K["sigma_C0_mM"]) = @sprintf("fusion β=%g σC0=%g", β, σ)

    # 1 et 2 — boucle ouverte (régulateur sur le vrai C), tous les filtres passifs sur les mêmes histogrammes
    passifs = [[spec(nomβ(β, σ), :fusion, β, σ) for σ in σs for β in βs]; spec("ratio", :ratio, 0.0); spec("ratio_libre", :ratio_libre, 0.0)]
    BO = []
    for N in Ns, (forme, α) in derives(p0)
        p = avec(p0; N, forme, alpha = α)
        res = repeter(r -> course(p, lecs[N], ph0, r; pilote = :oracle, passifs, reg = regulateur(p, ph0)), ns)
        for s in passifs
            M = [metriques(c, T[s.nom], p) for (c, T) in res]
            push!(BO, merge((N_photons_bleus = N, derive = nom_derive(forme, α, p), filtre = s.nom, beta = s.β,
                             sigma_C0 = get(s, :σC0, NaN)), moyenne(M), (biais_30s_sd = ecart_type(M, :biais_30s),)))
        end
        @printf("course 1–2  N = %g  %-10s  fait\n", N, nom_derive(forme, α, p))
    end
    ligne(L, N, d, f) = L[findfirst(x -> x.N_photons_bleus == N && x.derive == d && x.filtre == f, L)]
    cr = p0.toml["criteres"]
    passe(β) = all(N -> (l = ligne(BO, N, "aucune", nomβ(β)); abs(l.innov_transitoires) < cr["innovation_transitoire_max"] &&
                        abs(l.autocorr_transitoires) < cr["autocorr_par_phase_max"] && abs(l.autocorr_plateaux) < cr["autocorr_par_phase_max"]), Ns)
    # repli si aucun β ne passe tout le test 1 : le plus petit β qui corrige l'innovation moyenne de transitoire
    # (le défaut que β vise), l'autocorrélation étant alors rapportée comme indépendante de β
    innov_ok(β) = all(N -> abs(ligne(BO, N, "aucune", nomβ(β)).innov_transitoires) < cr["innovation_transitoire_max"], Ns)
    iβ = findfirst(passe, βs)
    iβ2 = findfirst(innov_ok, βs)
    β_ret = iβ !== nothing ? βs[iβ] : iβ2 !== nothing ? βs[iβ2] : βs[argmin([abs(ligne(BO, 1000.0, "aucune", nomβ(β)).innov_transitoires) for β in βs])]
    note_β = iβ !== nothing ? "plus petit β qui passe le test 1" : iβ2 !== nothing ? "aucun β ne passe tout le test 1 ; plus petit β qui corrige l'innovation moyenne de transitoire" :
             "aucun β ne corrige l'innovation : β à l'innovation la plus faible"
    @printf("β retenu : %g (%s)\n", β_ret, note_β)
    fus = spec("fusion", :fusion, β_ret)

    # 3 et 4 — boucle fermée, variante A : fusion contre ratio seul (g figé) et ratio à g libre, mêmes germes
    BF, traces4 = [], Dict()
    for N in Ns, (forme, α) in derives(p0), s in (fus, spec("ratio", :ratio, β_ret), spec("ratio_libre", :ratio_libre, β_ret))
        p = avec(p0; N, forme, alpha = α)
        res = repeter(r -> course(p, lecs[N], ph0, r; pilote = s, reg = regulateur(p, ph0)), np)
        M = [metriques(c, T[s.nom], p) for (c, T) in res]
        push!(BF, merge((N_photons_bleus = N, derive = nom_derive(forme, α, p), filtre = s.nom), moyenne(M),
                        (erreur_reelle_sd = ecart_type(M, :erreur_reelle), erreur_regulation_sd = ecart_type(M, :erreur_regulation))))
        (N == 1000.0 && α == 0.2) && (traces4[s.nom] = res)
        @printf("course 3–4  N = %g  %-10s  %-12s erreur affichée %+.2f  réelle %+.2f mM\n", N, nom_derive(forme, α, p), s.nom,
                BF[end].erreur_regulation, BF[end].erreur_reelle)
    end
    c_cal = Dict(N => ligne(BF, N, "aucune", "fusion").P11_sur_err2 for N in Ns)

    # 5 — balayage de Kp (N = 1000, sans dérive), quatre signatures ; contrôle à β = 0 ; puis β au coude − 30 %
    p = avec(p0; N = 1000.0, alpha = 0.0)
    rg = p0.toml["regulateur"]
    Kn = Kp_nominal(p, ph0)
    d0, d1 = rg["balayage_Kp_decades"]
    Kps = Kn .* exp10.(range(d0, d1; length = round(Int, (d1 - d0) * rg["points_par_decade"]) + 1))
    SIG = []
    sig(β, Kp; variante = :A, c = 1.0, σ0² = 1.0) = begin
        s = spec("fusion", :fusion, β)
        res = repeter(r -> course(p, lecs[1000.0], ph0, r; pilote = s, reg = regulateur(p, ph0; Kp, variante, c, σ0²)), ns)
        merge(moyenne([signatures(c, T["fusion"], p) for (c, T) in res]), moyenne([metriques(c, T["fusion"], p) for (c, T) in res]))
    end
    for β in unique([β_ret, 0.0, maximum(βs)]), Kp in Kps
        push!(SIG, merge((essai = "Kp", beta = β, Kp = Kp, Kp_sur_nominal = Kp / Kn), sig(β, Kp)))
        @printf("course 5  β = %g  Kp = %.4f (×%.2f)  agitation %.4f\n", β, Kp, Kp / Kn, SIG[end].agitation)
    end
    coude(β) = begin
        L = filter(x -> x.essai == "Kp" && x.beta == β, SIG)
        a = [x.agitation / x.Kp for x in L]
        i = findfirst(>(cr["coude_rapport"] * a[1]), a)
        i === nothing ? (NaN, Kps[end]) : (L[i].Kp, L[i].Kp)
    end
    Kc, Kc_borne = coude(β_ret)
    Kc0, _ = coude(0.0)
    for β in βs
        push!(SIG, merge((essai = "beta", beta = β, Kp = rg["facteur_coude_beta"] * Kc_borne, Kp_sur_nominal = rg["facteur_coude_beta"] * Kc_borne / Kn),
                         sig(β, rg["facteur_coude_beta"] * Kc_borne)))
    end

    # 6 — variantes B et C du régulateur (Kp nominal), P₁₁ corrigé par la calibration de la course 3
    σ0² = median(filter(isfinite, reduce(vcat, [T["fusion"][:P11][fenetres(c[:t], p.toml["scenario"]["etabli"])] for (c, T) in traces4["fusion"]]))) / c_cal[1000.0]
    VAR = []
    for (v, Kp) in ((:A, Kn), (:B, Kn), (:C, rg["Kp_variante_C_facteur"] * Kn))
        push!(VAR, merge((variante = string(v), Kp = Kp), sig(β_ret, Kp; variante = v, c = c_cal[1000.0], σ0²)))
        @printf("course 6  variante %s : agitation %.4f  erreur réelle %+.2f\n", v, VAR[end].agitation, VAR[end].erreur_reelle)
    end

    # 7 — filtre à quatre états : k_pompe estimable ? (boucle fermée sur ce filtre, sans dérive)
    K7, traces7 = [], Dict()
    for k in p0.toml["physiologie"]["balayage_k_pompe"]
        ph = physio(p0; k)
        s = spec("fusion4", :fusion4, β_ret)
        res = repeter(r -> course(p, lecs[1000.0], ph, r; pilote = s, reg = regulateur(p, ph0)), np)
        traces7[k] = res[1]
        for (c, T) in res
            t = c[:t]
            push!(K7, (k_vrai = k, k_final = T["fusion4"][:k_hat][end], sigma_k_final = sqrt(T["fusion4"][:P44][end]),
                       P44_repos_debut = T["fusion4"][:P44][findlast(<(10.0), t)], P44_fin_creneau = T["fusion4"][:P44][findlast(<(30.0), t)],
                       P44_fin_descente = T["fusion4"][:P44][findlast(<(65.0), t)], erreur_reelle = metriques(c, T["fusion4"], p).erreur_reelle))
        end
        @printf("course 7  k vrai %.2f : k̂ final %.3f ± %.3f\n", k, mean(x.k_final for x in K7 if x.k_vrai == k), mean(x.sigma_k_final for x in K7 if x.k_vrai == k))
    end
    f4 = filtre(p, lecs[1000.0], ph0, spec("f4", :fusion4, β_ret))
    rang = [(C_mM = C, rang = rank(observabilite(f4, C, 1.0); rtol = 1e-10), cond = cond(observabilite(f4, C, 1.0))) for C in (5.0, 5.5, 10.0, 25.0)]

    # robustesse — la conclusion de la course 4 sous k_pompe et J_max balayés (rampe 20 %, N = 1000)
    ROB = []
    pr = avec(p0; N = 1000.0, alpha = 0.2)
    for (k, J) in unique(vcat([(k, ph0.J_max) for k in p0.toml["physiologie"]["balayage_k_pompe"]], [(ph0.k, J) for J in p0.toml["physiologie"]["balayage_J_max"]])),
        s in (fus, spec("ratio", :ratio, β_ret))
        ph = physio(p0; k, J_max = J)
        res = repeter(r -> course(pr, lecs[1000.0], ph, r; pilote = s, reg = regulateur(pr, ph)), np)
        M = [metriques(c, T[s.nom], pr) for (c, T) in res]
        push!(ROB, merge((k_pompe = k, J_max = J, filtre = s.nom), moyenne(M), (erreur_reelle_sd = ecart_type(M, :erreur_reelle),)))
        @printf("robustesse k %.2f J %g %-7s : affichée %+.2f  réelle %+.2f\n", k, J, s.nom, ROB[end].erreur_regulation, ROB[end].erreur_reelle)
    end

    # ------------------------------------------------------------------------- critères de sortie
    l3 = ligne(BF, 1000.0, "aucune", "fusion")
    bo_ret = Dict(N => ligne(BO, N, "aucune", nomβ(β_ret)) for N in Ns)
    bo_0 = Dict(N => ligne(BO, N, "aucune", nomβ(0.0)) for N in Ns)
    phs = noms_phases(p0)
    Kp5 = (1 / 5 - ph0.k) / dflux(u_maintien(p0.toml["regulateur"]["C_ref_reglage_mM"], ph0), ph0)
    marge = Kc_borne / Kn
    o(b) = b ? "oui" : "non"
    f2(x) = @sprintf("%+.2f", x)
    vr(N) = verif[findfirst(v -> v.N_photons_bleus == N, verif)]
    crit = Any[
        (critere = "innovation moyenne, transitoires, N = 1000 (β retenu)", valeur = @sprintf("%+.3f", bo_ret[1000.0].innov_transitoires), seuil = "< 0,05", adosse = "référence : −0,274 à l'étape 1", tenu = o(abs(bo_ret[1000.0].innov_transitoires) < cr["innovation_transitoire_max"])),
        (critere = "autocorrélation d'innovation, transitoires / plateaux, N = 1000", valeur = @sprintf("%+.3f / %+.3f", bo_ret[1000.0].autocorr_transitoires, bo_ret[1000.0].autocorr_plateaux), seuil = "< 0,10", adosse = "blancheur théorique", tenu = o(abs(bo_ret[1000.0].autocorr_transitoires) < cr["autocorr_par_phase_max"] && abs(bo_ret[1000.0].autocorr_plateaux) < cr["autocorr_par_phase_max"])),
        (critere = "hausse de la RMSE sur les plateaux, β retenu contre β = 0", valeur = join((@sprintf("%+.1f %%", 100(bo_ret[N].rmse_plateaux / bo_0[N].rmse_plateaux - 1)) for N in Ns), " / "), seuil = "≤ 5 %", adosse = "référence : course 1, même budget", tenu = o(all(bo_ret[N].rmse_plateaux <= (1 + cr["hausse_rmse_plateaux_max"]) * bo_0[N].rmse_plateaux for N in Ns))),
        (critere = "P₁₁/err² par phase, boucle fermée, N = 1000", valeur = join((@sprintf("%.2f", getfield(l3, Symbol("P11_err2_", n))) for n in phs), " / "), seuil = "dans [0,5 ; 2]", adosse = "conservé de l'étape 1", tenu = o(all(cr["P11_sur_err2"][1] <= getfield(l3, Symbol("P11_err2_", n)) <= cr["P11_sur_err2"][2] for n in phs))),
        (critere = "biais du filtre, 5 phases, boucle fermée sans dérive, N = 1000 (mM)", valeur = join((f2(getfield(l3, Symbol("biais_", n))) for n in phs), " / "), seuil = "< ±0,30", adosse = "remesuré contre la bonne dynamique", tenu = o(all(abs(getfield(l3, Symbol("biais_", n))) < cr["biais_filtre_phases_mM"] for n in phs))),
        (critere = "erreur de régulation établie, sans dérive, fusion (mM)", valeur = join((f2(ligne(BF, N, "aucune", "fusion").erreur_regulation) for N in Ns), " / "), seuil = "< 0,5", adosse = "10 % du créneau de 20 mM", tenu = o(all(abs(ligne(BF, N, "aucune", "fusion").erreur_regulation) < cr["erreur_regulation_etablie_mM"] for N in Ns))),
        (critere = "erreur réelle, rampe 20 %, fusion (mM)", valeur = join((f2(ligne(BF, N, "rampe 20 %", "fusion").erreur_reelle) for N in Ns), " / "), seuil = "< 1,5", adosse = "la fusion doit tenir sa promesse", tenu = o(all(abs(ligne(BF, N, "rampe 20 %", "fusion").erreur_reelle) < cr["erreur_reelle_fusion_rampe20_mM"] for N in Ns))),
        (critere = "erreur réelle, rampe 20 %, ratio seul (mM)", valeur = join((f2(ligne(BF, N, "rampe 20 %", "ratio").erreur_reelle) for N in Ns), " / "), seuil = "> 4 attendu", adosse = "sinon la question 2 répond non", tenu = o(all(abs(ligne(BF, N, "rampe 20 %", "ratio").erreur_reelle) > cr["erreur_reelle_ratio_rampe20_mM"] for N in Ns))),
        (critere = "accord Fisher / Monte-Carlo, var J/B et var 5a (%)", valeur = join((@sprintf("%.1f / %.1f", 100abs(vr(N).var_ratio_fisher / vr(N).var_ratio_reelle - 1), 100abs(vr(N).var_5a_fisher / vr(N).var_5a_reelle - 1)) for N in Ns), " ; "), seuil = "< 5", adosse = "tenu à 0,6 et 0,4 % à l'étape 1", tenu = o(all(abs(vr(N).var_ratio_fisher / vr(N).var_ratio_reelle - 1) < cr["accord_fisher_pct"] / 100 && abs(vr(N).var_5a_fisher / vr(N).var_5a_reelle - 1) < cr["accord_fisher_pct"] / 100 for N in Ns))),
        (critere = "ρ₁₂ en boucle fermée contre boucle ouverte, régime établi, N = 1000", valeur = @sprintf("%.2f contre %.2f", l3.rho12_etabli, bo_ret[1000.0].rho12_etabli), seuil = "≤ BO + 0,1", adosse = "référence : boucle ouverte au même budget", tenu = o(l3.rho12_etabli <= bo_ret[1000.0].rho12_etabli + cr["rho12_hausse_max"])),
        (critere = "marge de gain avant le coude (Kp coude / Kp nominal)", valeur = isnan(Kc) ? @sprintf("> %.1f (pas de coude dans le balayage)", marge) : @sprintf("%.1f", marge), seuil = "≥ 3", adosse = "opinion, à réviser après la course 5", tenu = o(marge >= cr["marge_coude_min"])),
        (critere = "échecs de 5a, C = 5 mM, N = 1000 (%)", valeur = @sprintf("%.1f", vr(1000.0).echecs_5a_5mM_pct), seuil = "à rapporter", adosse = "32,7 % mesuré à l'étape 1", tenu = "—"),
    ]
    issue = (abs(ligne(BF, 1000.0, "rampe 20 %", "fusion").erreur_reelle) < cr["erreur_reelle_fusion_rampe20_mM"] &&
             abs(ligne(BF, 1000.0, "rampe 20 %", "ratio").erreur_reelle) > cr["erreur_reelle_ratio_rampe20_mM"]) ? "A — la fusion gagne" :
            (!isnan(Kc) && Kc < Kp5) ? "C — le couplage domine" : "B — la boucle fermée annule l'avantage"
    @printf("\nβ retenu %g ; Kp nominal %.4f ; coude %s (β = 0 : %s) ; Kp pour τ_BF = 5 s : %.4f\nissue : %s\n", β_ret, Kn,
            isnan(Kc) ? "aucun" : @sprintf("%.4f", Kc), isnan(Kc0) ? "aucun" : @sprintf("%.4f", Kc0), Kp5, issue)
    foreach(c -> @printf("  %-70s %-30s %-10s %s\n", c.critere, c.valeur, c.seuil, c.tenu), crit)
    resume = [(cle = "beta_retenu", valeur = string(β_ret), note = note_β), (cle = "Kp_nominal", valeur = @sprintf("%.5f", Kn), note = "τ_BF = 3 s au maintien de 25 mM"),
              (cle = "Kp_coude", valeur = isnan(Kc) ? "aucun" : @sprintf("%.5f", Kc), note = "β retenu"), (cle = "Kp_coude_beta0", valeur = isnan(Kc0) ? "aucun" : @sprintf("%.5f", Kc0), note = "contrôle, Q constant"),
              (cle = "Kp_tau_BF_5s", valeur = @sprintf("%.5f", Kp5), note = "gain nécessaire pour τ_BF = 5 s"), (cle = "calibration_P11_N1000", valeur = @sprintf("%.3f", c_cal[1000.0]), note = "c des variantes B et C"),
              (cle = "sigma0_carre_variante_C", valeur = @sprintf("%.4f", σ0²), note = "médiane de P₁₁/c en régime établi, course 3"), (cle = "issue", valeur = issue, note = "posée d'avance dans le plan")]
    for (nom, l) in (("boucle_ouverte", BO), ("boucle_fermee", BF), ("couplage", SIG), ("variantes", VAR), ("k_pompe", K7), ("observabilite", rang),
                     ("robustesse", ROB), ("verification_R", verif), ("criteres", crit), ("resume", resume))
        ecrire_csv(joinpath(SORTIES, nom * ".csv"), l)
    end
    figure_article(traces4, p0, joinpath(SORTIES, "fig1_erreur_affichee_reelle.png"))
    figure_blancheur(BO, p0, β_ret, joinpath(SORTIES, "fig2_blancheur.png"))
    figure_couplage(SIG, Kn, Kc, β_ret, p0, joinpath(SORTIES, "fig3_coude.png"))
    figure_k_pompe(traces7, ph0, p0, joinpath(SORTIES, "fig4_k_pompe.png"))
    figure_fisher(verif, joinpath(SORTIES, "fig_fisher.png"))
    rapide || video_boucle(traces4, p0, joinpath(SORTIES, "video_rampe20"))
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
