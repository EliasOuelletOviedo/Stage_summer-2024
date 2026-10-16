# run.jl — la course : 900 pas à 10 Hz, trois lectures (J/B, 5a, Kalman). Tous les réglages du
# filtre (marche aléatoire au balayage de σ_w,g, modèle (g, ġ), variante Ornstein-Uhlenbeck)
# tournent sur les mêmes lectures, aux trois budgets. Sorties : un fichier plat par pas, un
# résumé par scénario et par filtre, les vérifications, les critères de sortie, trois vidéos.
# Simule tout d'abord, rend ensuite.
#
# Lancer depuis ce dossier :  julia --project=. -t auto run.jl  [--rapide]
# (--rapide : 2 répétitions, sans vidéo, pour vérifier que tout tourne)

using Printf, Statistics, LinearAlgebra, StableRNGs, Base.Threads
for f in ("params", "generate", "estimate", "kalman")
    include(joinpath(@__DIR__, "src", f * ".jl"))
end

const SORTIES = joinpath(@__DIR__, "sorties")
const COLONNES = [:t, :C_vrai, :g_vrai, :C_ratio, :C_5a, :C_kalman, :g_chapeau, :P11, :P22, :rho12,
                  :innov_ratio, :innov_5a, :N_bleu, :N_jaune, :gp_chapeau, :z_ratio, :z_5a]

nom_ma(σ) = @sprintf("MA σg=%.0e", σ)
reference(p::Params) = nom_ma(p.toml["kalman"]["sigma_wg"])

"""Les réglages du filtre comparés sur les mêmes lectures."""
function configurations(p::Params)
    k, v = p.toml["kalman"], p.toml["kalman"]["vitesse"]
    c = [(nom = nom_ma(σ), modele = :marche, σwg = σ, σwgp = NaN, τgp = Inf) for σ in k["balayage_sigma_wg"]]
    for σ in v["balayage_sigma_wgp"]
        push!(c, (nom = @sprintf("vitesse σġ=%.0e", σ), modele = :vitesse, σwg = v["sigma_wg"], σwgp = σ, τgp = Inf))
    end
    push!(c, (nom = @sprintf("OU σġ=%.0e τ=%gs", v["sigma_wgp_rappel"], v["tau_rappel_s"]), modele = :vitesse,
              σwg = v["sigma_wg"], σwgp = v["sigma_wgp_rappel"], τgp = v["tau_rappel_s"]))
    return c
end

function filtre(p::Params, lec, c)
    k, r = p.toml["kalman"], p.toml["trajectoire"]
    cR, cS = lec.JB.courbe, lec.S.courbe
    v = c.modele === :vitesse
    return EKF(modele = c.modele, x0 = v ? [k["C0_mM"], k["g0"], 0.0] : [k["C0_mM"], k["g0"]],
        σ0 = v ? [k["sigma_C0_mM"], k["sigma_g0"], k["vitesse"]["sigma_gp0"]] : [k["sigma_C0_mM"], k["sigma_g0"]],
        σw = v ? [k["sigma_wC_mM"], c.σwg, c.σwgp] : [k["sigma_wC_mM"], c.σwg], dt = p.dt, C_repos = r["C_repos_mM"],
        τ = r["tau_KCC2_s"], τ_gp = c.τgp, h = Function[c -> valeur(cR, c), c -> valeur(cS, c)],
        dh = Function[c -> derivee(cR, c), c -> derivee(cS, c)], gain = [true, false])
end

"""Une course (répétition r) : une table par réglage du filtre, et les histogrammes si demandés (vidéo)."""
function course(p::Params, lec, r; histos = false, configs = configurations(p))
    ts = collect(temps(p))
    n, W = length(ts), p.fenetre
    f = Dict(c.nom => filtre(p, lec, c) for c in configs)
    T = Dict(c.nom => Dict(col => fill(NaN, n) for col in COLONNES) for c in configs)
    H = Vector{NTuple{2,Vector{Float64}}}(undef, histos ? n : 0)
    for (k, t) in enumerate(ts)
        hb, hj = histogrammes(t, p; rng = rng_pas(k, r, p))
        histos && (H[k] = (hb, hj))
        lR, lS = estimer(hb, hj, lec.JB.cfg), estimer(hb, hj, lec.S.cfg)
        Cr = lR.ok ? inverser(lec.JB.courbe, lR.z) : NaN
        Cs = lS.ok ? inverser(lec.S.courbe, lS.z) : NaN
        for (nom, fi) in f
            predire!(fi)
            R = R_attendue(lec, fi.x[1], fi.x[2])          # Fisher à l'estimation courante du filtre
            ν = mettre_a_jour!(fi, [lR.z, lS.z], R, [lR.ok, lS.ok])
            d = diagnostics(fi)
            v = (t, concentration(t, p), derive(t, p), Cr, Cs, d.C, d.g, d.P11, d.P22, d.rho12, ν[1], ν[2],
                 sum(hb[W]), sum(hj[W]), d.gp, lR.z, lS.z)
            foreach(((c, x),) -> T[nom][c][k] = x, zip(COLONNES, v))
        end
    end
    return T, H
end

autocorr1(x) = (x = filter(!isnan, x); y = x .- mean(x); sum(y[1:end-1] .* y[2:end]) / sum(abs2, y))
moy_ok(x) = (x = filter(!isnan, x); isempty(x) ? NaN : mean(x))

"""
Métriques d'une course pour une estimation de C (colonne `Cest`). Pour une lecture, aussi le
biais « linéaire » de son observable brute, (z − h(C))/h′(C), sans inversion ni troncature ;
pour le Kalman, les diagnostics du filtre.
"""
function metriques(T, Cest, p::Params; courbe = nothing, z = :z_ratio)
    t, C = T[:t], T[:C_vrai]
    tr, d = p.toml["trajectoire"], p.derive
    fin = t .>= t[end] + p.dt - p.toml["course"]["fin_biais_s"]
    montee = tr["t_debut_montee_s"] .<= t .< tr["t_fin_montee_s"]
    e = T[Cest] .- C
    v = .!isnan.(e)
    lin = courbe === nothing ? NaN : mean(((T[z] .- valeur.(Ref(courbe), C)) ./ derivee.(Ref(courbe), C))[fin])
    out = (rmse = sqrt(mean(e[v] .^ 2)), biais_30s = mean(e[v.&fin]), err_max_montee = maximum(abs.(e[v.&montee]); init = -Inf),
           echecs_pct = 100 * count(!, v) / length(v), biais_lineaire_30s = lin)
    courbe === nothing || return merge(out, (t_P22_moitie = NaN, P11_sur_err2 = NaN, biais_g_30s = NaN, t_recup_g = NaN,
        rho12_moyen = NaN, var_innov_ratio = NaN, var_innov_5a = NaN, autocorr_innov_ratio = NaN, innov_ratio_montee = NaN,
        innov_ratio_hors_montee = NaN))
    i = findfirst(<=(p.toml["kalman"]["sigma_g0"]^2 / 2), T[:P22])
    eg = T[:g_chapeau] .- T[:g_vrai]
    j = d.forme == "saut" ? findfirst(k -> t[k] >= d.t_saut && abs(eg[k]) <= 0.5d.saut, eachindex(t)) : nothing
    νR, νS = T[:innov_ratio], T[:innov_5a]
    return merge(out, (t_P22_moitie = i === nothing ? NaN : t[i], P11_sur_err2 = mean(T[:P11]) / mean(e .^ 2),
        biais_g_30s = mean(eg[fin]), t_recup_g = j === nothing ? NaN : t[j] - d.t_saut, rho12_moyen = mean(abs.(T[:rho12])),
        var_innov_ratio = var(filter(!isnan, νR)), var_innov_5a = var(filter(!isnan, νS)), autocorr_innov_ratio = autocorr1(νR),
        innov_ratio_montee = moy_ok(νR[montee]), innov_ratio_hors_montee = moy_ok(νR[.!montee])))
end

function ecrire_csv(chemin, lignes)
    open(chemin, "w") do io
        println(io, join(string.(keys(first(lignes))), ","))
        for l in lignes
            println(io, join((x isa Number ? @sprintf("%.6g", x) : string(x) for x in l), ","))
        end
    end
end

nom_derive(forme, α, p) = forme == "saut" ? @sprintf("saut %g %%", 100p.derive.saut) : α == 0 ? "aucune" : @sprintf("rampe %g %%", 100α)

tirage(p, C, i) = histogrammes_c(C, 1.0, p; rng = StableRNG(p.toml["course"]["graine"] + 5_000_000 + i))

"""Dispersion réelle de (z_R, z_S) sur des tirages à C_reference, contre le R de Fisher moyen."""
function verifier_R(p::Params, lec)
    Cr, M = reglage(p, "lectures", "C_reference_mM"), p.toml["course"]["verification_R_tirages"]
    L = filter(l -> l[1].ok && l[2].ok, map(i -> (h = tirage(p, Cr, i); (estimer(h..., lec.JB.cfg), estimer(h..., lec.S.cfg))), 1:M))
    Z = hcat([l[1].z for l in L], [l[2].z for l in L])
    Ce, Rf = cov(Z), mean(R_jointe(l...) for l in L)
    return (N_photons_bleus = p.N, tirages = length(L), var_ratio_reelle = Ce[1, 1], var_ratio_fisher = Rf[1, 1],
            var_5a_reelle = Ce[2, 2], var_5a_fisher = Rf[2, 2], cov_reelle = Ce[1, 2], cov_fisher = Rf[1, 2],
            rho_reel = Ce[1, 2] / sqrt(Ce[1, 1] * Ce[2, 2]), rho_fisher = Rf[1, 2] / sqrt(Rf[1, 1] * Rf[2, 2]))
end

"""
D'où vient le biais de 5a ? Sur des tirages à C fixé (g = 1) : le biais de l'observable brute
(z − h(C))/h′(C), que reçoit le filtre, contre le biais de la concentration inversée, calculé
sur les seuls tirages dans la plage (les autres sont des échecs).
"""
function biais_lectures(p::Params, lec)
    M = p.toml["course"]["verification_R_tirages"]
    lignes = []
    for C in p.toml["course"]["concentrations_biais_mM"], (nom, l) in (("J/B", lec.JB), ("5a", lec.S))
        z = zeros(M)
        @threads for i in 1:M
            z[i] = estimer(tirage(p, C, i)..., l.cfg).z
        end
        Ĉ = inverser.(Ref(l.courbe), z)
        ok = .!isnan.(Ĉ)
        push!(lignes, (N_photons_bleus = p.N, C_mM = C, lecture = nom, biais_observable_mM = (mean(z) - valeur(l.courbe, C)) / derivee(l.courbe, C),
                       sigma_mM = std(z) / abs(derivee(l.courbe, C)), biais_inverse_mM = mean(Ĉ[ok]) - C, echecs_pct = 100 * count(!, ok) / M))
    end
    return lignes
end

"""Copie de `p` sans diffusion dans aucun canal."""
sans_diffusion(p::Params) = (L = copy(p.lang); L["d_bleu"] = L["d_jaune"] = Langmuir(0.0, 0.0, 1.0);
                             Params((f === :lang ? L : getfield(p, f) for f in fieldnames(Params))...))

"""
Test d'amplitudes sans bruit sur toute la plage : p et q que 5a retrouve contre p et q injectés,
avec la diffusion du générateur et sans elle. Écart constant : normalisation ; écart qui suit
la diffusion et disparaît sans elle : la diffusion, que 5a ne modélise pas.
"""
function amplitudes(p::Params)
    a, b, _ = reglage(p, "lectures", "grille_calibration_mM")
    ps = (avec = p, sans = sans_diffusion(p))
    cp, cq = Cfg5a(p, (1.0, 0.0)), Cfg5a(p, (0.0, 1.0))
    map(collect(a:p.toml["course"]["pas_amplitudes_mM"]:b)) do C
        r = Dict(k => (estimer(histogrammes_c(C, 1.0, q)..., cp).z, estimer(histogrammes_c(C, 1.0, q)..., cq).z) for (k, q) in pairs(ps))
        (C_mM = C, p_injecte = p.lang["p_bleu"](C), p_avec_diffusion = r[:avec][1], p_sans_diffusion = r[:sans][1],
         q_injecte = p.lang["q_jaune"](C), q_avec_diffusion = r[:avec][2], q_sans_diffusion = r[:sans][2],
         d_bleu = p.lang["d_bleu"](C), d_jaune = p.lang["d_jaune"](C))
    end
end

ligne(resume, N, der, filtre) = resume[findfirst(r -> r.N_photons_bleus == N && r.derive == der && r.filtre == filtre, resume)]

"""Critères de sortie vers l'étape 2 (chiffrés d'avance dans params.toml), pour chaque réglage du filtre."""
function criteres(resume, p::Params)
    c, ref = p.toml["criteres"], reference(p)
    ders = unique(r.derive for r in resume)
    out = []
    for N in c["budgets"], cf in configurations(p)
        K = [ligne(resume, N, d, cf.nom) for d in ders]
        J = [ligne(resume, N, d, "J/B") for d in ders]
        Rf = [ligne(resume, N, d, ref) for d in ders]
        b10 = ligne(resume, N, "rampe 10 %", cf.nom).biais_30s
        P = [k.P11_sur_err2 for k in K]
        q = maximum(k.rmse / j.rmse for (k, j) in zip(K, J))
        bat = cf.modele === :marche ? "—" : all(begin
            r = Rf[i]; k = K[i]
            startswith(d, "rampe") ? k.rmse <= r.rmse :
            startswith(d, "saut") ? k.rmse <= (1 + c["tolerance_saut"]) * r.rmse && abs(k.biais_30s) <= abs(r.biais_30s) + c["tolerance_biais_saut_mM"] :
            k.rmse <= (1 + c["tolerance_saut"]) * r.rmse end for (i, d) in enumerate(ders)) ? "oui" : "non"
        push!(out, (N_photons_bleus = N, filtre = cf.nom, biais_rampe10_mM = b10, P11_sur_err2_min = minimum(P), P11_sur_err2_max = maximum(P),
            rmse_kalman_sur_JB_max = q, c1_biais = abs(b10) < c["biais_rampe10_max_mM"] ? "oui" : "non",
            c2_calibration = all(c["P11_sur_err2"][1] .<= P .<= c["P11_sur_err2"][2]) ? "oui" : "non",
            c3_rmse = q < 1 ? "oui" : "non", bat_la_reference = bat))
    end
    return out
end

"""Correctif 4 : le σ_w,g dont P₁₁/err² sous la dérive la plus forte est ≥ 1 et le plus proche de 1."""
function sigma_recommande(resume, p::Params, N)
    pire = nom_derive("rampe", maximum(p.toml["derive"]["balayage_alpha"]), p)
    c = [(σ, ligne(resume, N, pire, nom_ma(σ)).P11_sur_err2) for σ in p.toml["kalman"]["balayage_sigma_wg"]]
    c = filter(x -> x[2] >= 1, c)
    return isempty(c) ? NaN : c[argmin(last.(c))][1]
end

include(joinpath(@__DIR__, "src", "figures.jl"))

function main(; rapide = "--rapide" in ARGS)
    mkpath(SORTIES)
    p0 = charger_params()
    nrep = rapide ? 2 : p0.toml["course"]["repetitions"]
    vid, ref = p0.toml["video"], reference(p0)
    configs = configurations(p0)
    plat = [ref; [c.nom for c in configs if c.modele === :vitesse && c.τgp == Inf && c.σwgp == p0.toml["kalman"]["vitesse"]["sigma_wgp_rappel"]]]
    lignes_pas, resume, coherence, verif, biais5a, saut, videos = [], [], [], [], [], Dict(), Dict()
    for N in p0.toml["photons"]["balayage"]
        pN = avec(p0; N)
        lec = lecteurs(pN)
        Cr = reglage(pN, "lectures", "C_reference_mM")
        μ = histogrammes_c(Cr, 1.0, pN)
        lR, lS, lp, lq = estimer(μ..., lec.JB.cfg), estimer(μ..., lec.S.cfg), estimer(μ..., Cfg5a(pN, (1.0, 0.0))), estimer(μ..., Cfg5a(pN, (0.0, 1.0)))
        σJB = sqrt(lR.σ²) / abs(derivee(lec.JB.courbe, Cr))
        push!(coherence, (N_photons_bleus = N, C_JB_mM = inverser(lec.JB.courbe, lR.z), C_5a_mM = inverser(lec.S.courbe, lS.z),
            p_5a = lp.z, p_titration = pN.lang["p_bleu"](Cr), q_5a = lq.z, q_titration = pN.lang["q_jaune"](Cr), w_p = lec.w[1], w_q = lec.w[2],
            sigmaC_JB = σJB, sigmaC_p = lec.σC_reference.p, sigmaC_q = lec.σC_reference.q, sigmaC_5a = lec.σC_reference.pq,
            sigma_5a_sur_JB = lec.σC_reference.pq / σJB))
        @printf("N = %5.0f  σ_C à %g mM : J/B %.2f, 5a p %.2f, q %.2f, p+q %.2f mM\n", N, Cr, σJB, lec.σC_reference...)
        push!(verif, verifier_R(pN, lec))
        append!(biais5a, biais_lectures(pN, lec))
        for (forme, α) in [[("rampe", α) for α in p0.toml["derive"]["balayage_alpha"]]; [("saut", 0.0)]]
            p = avec(pN; forme, alpha = α)
            nom = nom_derive(forme, α, p)
            nv = N != vid["N_bleu"] ? nothing : forme == "saut" ? "video_saut" : α == 0 ? "video_sans_derive" : α == vid["alpha"] ? "video_rampe" : nothing
            res = Vector{Any}(undef, nrep)
            @threads for r in 1:nrep
                res[r] = course(p, lec, r; histos = !rapide && nv !== nothing && r == 1, configs)
            end
            T1 = res[1][1]
            nv === nothing || rapide || (videos[nv] = (T1[ref], res[1][2], p))
            forme == "saut" && N == vid["N_bleu"] && (saut = Dict(c.nom => (t = T1[ref][:t], g = T1[ref][:g_vrai],
                gm = mean(res[r][1][c.nom][:g_chapeau] for r in 1:nrep)) for c in configs))
            for f in plat, k in eachindex(T1[f][:t])
                push!(lignes_pas, merge((N_photons_bleus = N, derive = nom, filtre = f), NamedTuple{Tuple(COLONNES)}(Tuple(T1[f][c][k] for c in COLONNES))))
            end
            moy(M) = merge(NamedTuple{keys(M[1])}(Tuple(mean(getfield(x, k) for x in M) for k in keys(M[1]))),
                           (rmse_sd = std(getfield.(M, :rmse)), biais_30s_sd = std(getfield.(M, :biais_30s))))
            tete(l, f) = (N_photons_bleus = N, derive = nom, lecture = l, filtre = f)
            push!(resume, merge(tete("J/B", "J/B"), moy([metriques(x[1][ref], :C_ratio, p; courbe = lec.JB.courbe, z = :z_ratio) for x in res])))
            push!(resume, merge(tete("5a", "5a"), moy([metriques(x[1][ref], :C_5a, p; courbe = lec.S.courbe, z = :z_5a) for x in res])))
            for c in configs
                push!(resume, merge(tete("Kalman", c.nom), moy([metriques(x[1][c.nom], :C_kalman, p) for x in res])))
            end
            k = ligne(resume, N, nom, ref)
            @printf("  %-12s biais 30 s : J/B %+6.2f  Kalman %+6.2f mM ; RMSE %.2f ; P11/err² %.2f\n", nom,
                    ligne(resume, N, nom, "J/B").biais_30s, k.biais_30s, k.rmse, k.P11_sur_err2)
        end
    end
    amp = amplitudes(avec(p0; N = vid["N_bleu"]))
    crit = criteres(resume, p0)
    for N in p0.toml["photons"]["balayage"]
        @printf("σ_w,g recommandé (P11/err² ≥ 1 et le plus proche de 1 sous la dérive la plus forte), N = %g : %.0e\n",
                N, sigma_recommande(resume, p0, N))
    end
    @printf("amplitudes sans diffusion : max |p̂ − p| = %.1e, max |q̂ − q| = %.1e\n",
            maximum(abs(a.p_sans_diffusion - a.p_injecte) for a in amp), maximum(abs(a.q_sans_diffusion - a.q_injecte) for a in amp))
    for (nom, l) in (("courses", lignes_pas), ("resume", resume), ("coherence", coherence), ("verification_R", verif),
                     ("biais_lectures", biais5a), ("amplitudes", amp), ("criteres", crit))
        ecrire_csv(joinpath(SORTIES, nom * ".csv"), l)
    end
    figure_resume(resume, saut, p0, joinpath(SORTIES, "fig_resume.png"))
    figure_reglage(resume, p0, joinpath(SORTIES, "fig_reglage.png"))
    figure_fisher(verif, joinpath(SORTIES, "fig_fisher.png"))
    figure_amplitudes(amp, joinpath(SORTIES, "fig_amplitudes.png"))
    for (nom, (T, H, p)) in videos
        video(T, H, p, joinpath(SORTIES, nom))
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
