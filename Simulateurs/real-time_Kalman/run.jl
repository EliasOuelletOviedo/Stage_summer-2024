# run.jl — la course : 900 pas à 10 Hz, trois lectures (J/B, 5a, Kalman), le Kalman avec R
# complet et avec R diagonal sur les mêmes lectures ; un fichier plat par pas, un résumé par
# scénario, la vérification de R, deux vidéos. Simule tout d'abord, rend ensuite.
#
# Lancer depuis ce dossier :  julia --project=. -t auto run.jl  [--rapide]
# (--rapide : 2 répétitions, sans vidéo, pour vérifier que tout tourne)

using Printf, Statistics, LinearAlgebra, StableRNGs, Base.Threads
for f in ("params", "generate", "estimate", "kalman")
    include(joinpath(@__DIR__, "src", f * ".jl"))
end

const SORTIES = joinpath(@__DIR__, "sorties")
const COLONNES = [:t, :C_vrai, :g_vrai, :C_ratio, :C_5a, :C_kalman, :g_chapeau, :P11, :P22, :rho12,
                  :innov_ratio, :innov_5a, :N_bleu, :N_jaune]
const MODES = (:complet, :diagonal)

function filtre(p::Params, lec; σwg = p.toml["kalman"]["sigma_wg"])
    k, r = p.toml["kalman"], p.toml["trajectoire"]
    cR, cS = lec.JB.courbe, lec.S.courbe
    return EKF(x0 = [k["C0_mM"], k["g0"]], P0 = Diagonal([k["sigma_C0_mM"]^2, k["sigma_g0"]^2]),
        σw = [k["sigma_wC_mM"], σwg], dt = p.dt, C_repos = r["C_repos_mM"], τ = r["tau_KCC2_s"],
        h = Function[c -> valeur(cR, c), c -> valeur(cS, c)], dh = Function[c -> derivee(cR, c), c -> derivee(cS, c)],
        gain = [true, false])
end

"""Une course (répétition r) : une table par mode de R, et les histogrammes si demandés (vidéo)."""
function course(p::Params, lec, r; histos = false, σwg = p.toml["kalman"]["sigma_wg"])
    ts = collect(temps(p))
    n, W = length(ts), p.fenetre
    f = Dict(m => filtre(p, lec; σwg) for m in MODES)
    T = Dict(m => Dict(c => fill(NaN, n) for c in COLONNES) for m in MODES)
    H = Vector{NTuple{2,Vector{Float64}}}(undef, histos ? n : 0)
    for (k, t) in enumerate(ts)
        hb, hj = histogrammes(t, p; rng = rng_pas(k, r, p))
        histos && (H[k] = (hb, hj))
        lR, lS = estimer(hb, hj, lec.JB.cfg), estimer(hb, hj, lec.S.cfg)
        Cr = lR.ok ? inverser(lec.JB.courbe, lR.z) : NaN
        Cs = lS.ok ? inverser(lec.S.courbe, lS.z) : NaN
        for m in MODES
            predire!(f[m])
            R = R_attendue(lec, f[m].x...)          # Fisher à l'estimation courante du filtre
            ν = mettre_a_jour!(f[m], [lR.z, lS.z], m == :complet ? R : Matrix(Diagonal(R)), [lR.ok, lS.ok])
            d = diagnostics(f[m])
            v = (t, concentration(t, p), derive(t, p), Cr, Cs, d.C, d.g, d.P11, d.P22, d.rho12, ν[1], ν[2], sum(hb[W]), sum(hj[W]))
            foreach(((c, x),) -> T[m][c][k] = x, zip(COLONNES, v))
        end
    end
    return T, H
end

autocorr1(x) = (x = filter(!isnan, x); y = x .- mean(x); sum(y[1:end-1] .* y[2:end]) / sum(abs2, y))

"""Métriques d'une course pour une estimation de C (colonne `Cest`) et, pour le Kalman, les diagnostics du filtre."""
function metriques(T, Cest, p::Params; kalman = false)
    t, C = T[:t], T[:C_vrai]
    tr, d = p.toml["trajectoire"], p.derive
    fin = t .>= t[end] + p.dt - p.toml["course"]["fin_biais_s"]
    montee = tr["t_debut_montee_s"] .<= t .< tr["t_fin_montee_s"]
    e = Cest .- C
    v = .!isnan.(e)
    out = (rmse = sqrt(mean(e[v] .^ 2)), biais_fin = mean(e[v.&fin]), err_max_montee = maximum(abs.(e[v.&montee]); init = -Inf),
           echecs_pct = 100 * count(!, v) / length(v))
    kalman || return merge(out, (t_P22_moitie = NaN, P11_sur_err2 = NaN, biais_g_fin = NaN, t_recup_g = NaN, rho12_moyen = NaN,
                                 var_innov_ratio = NaN, var_innov_5a = NaN, autocorr_innov_ratio = NaN, innov_ratio_montee = NaN))
    i = findfirst(<=(p.toml["kalman"]["sigma_g0"]^2 / 2), T[:P22])
    eg = T[:g_chapeau] .- T[:g_vrai]
    j = d.forme == "saut" ? findfirst(k -> t[k] >= d.t_saut && abs(eg[k]) <= 0.5d.saut, eachindex(t)) : nothing
    νR, νS = T[:innov_ratio], T[:innov_5a]
    return merge(out, (t_P22_moitie = i === nothing ? NaN : t[i], P11_sur_err2 = mean(T[:P11]) / mean(e .^ 2),
        biais_g_fin = mean(eg[fin]), t_recup_g = j === nothing ? NaN : t[j] - d.t_saut, rho12_moyen = mean(abs.(T[:rho12])),
        var_innov_ratio = var(filter(!isnan, νR)), var_innov_5a = var(filter(!isnan, νS)),
        autocorr_innov_ratio = autocorr1(νR), innov_ratio_montee = mean(filter(!isnan, νR[montee]))))
end

function ecrire_csv(chemin, entetes, lignes)
    open(chemin, "w") do io
        println(io, join(entetes, ","))
        for l in lignes
            println(io, join((x isa AbstractString ? x : @sprintf("%.6g", x) for x in l), ","))
        end
    end
end

nom_derive(forme, α, p) = forme == "saut" ? @sprintf("saut %g %%", 100p.derive.saut) : α == 0 ? "aucune" : @sprintf("rampe %g %%", 100α)

"""Dispersion réelle de (z_R, z_S) sur des tirages à C_reference, contre le R de Fisher moyen."""
function verifier_R(p::Params, lec)
    Cr, M = reglage(p, "lectures", "C_reference_mM"), p.toml["course"]["verification_R_tirages"]
    L = map(1:M) do i
        hb, hj = histogrammes_c(Cr, 1.0, p; rng = StableRNG(p.toml["course"]["graine"] + 5_000_000 + i))
        estimer(hb, hj, lec.JB.cfg), estimer(hb, hj, lec.S.cfg)
    end
    L = filter(l -> l[1].ok && l[2].ok, L)
    Z = hcat([l[1].z for l in L], [l[2].z for l in L])
    Ce, Rf = cov(Z), mean(R_jointe(l...) for l in L)
    return (p.N, Ce[1, 1], Rf[1, 1], Ce[2, 2], Rf[2, 2], Ce[1, 2], Rf[1, 2], Ce[1, 2] / sqrt(Ce[1, 1] * Ce[2, 2]),
            Rf[1, 2] / sqrt(Rf[1, 1] * Rf[2, 2]), 100 * (1 - length(L) / M))
end

include(joinpath(@__DIR__, "src", "figures.jl"))

function main(; rapide = "--rapide" in ARGS)
    mkpath(SORTIES)
    p0 = charger_params()
    nrep = rapide ? 2 : p0.toml["course"]["repetitions"]
    vid = p0.toml["video"]
    σ0 = p0.toml["kalman"]["sigma_wg"]
    lignes_pas, resume, coherence, verif, saut = [], [], [], [], Dict()
    videos = Dict()
    for N in p0.toml["photons"]["balayage"]
        pN = avec(p0; N)
        lec = lecteurs(pN)
        Cr = reglage(pN, "lectures", "C_reference_mM")
        μ = histogrammes_c(Cr, 1.0, pN)
        lR, lS, lp, lq = estimer(μ..., lec.JB.cfg), estimer(μ..., lec.S.cfg), estimer(μ..., Cfg5a(pN, (1.0, 0.0))), estimer(μ..., Cfg5a(pN, (0.0, 1.0)))
        σJB = sqrt(lR.σ²) / abs(derivee(lec.JB.courbe, Cr))
        push!(coherence, (N, inverser(lec.JB.courbe, lR.z), inverser(lec.S.courbe, lS.z), lp.z, pN.lang["p_bleu"](Cr), lq.z,
                          pN.lang["q_jaune"](Cr), lec.w..., σJB, lec.σC_reference..., lec.σC_reference.pq / σJB))
        @printf("N = %5.0f  σ_C à %g mM : J/B %.2f, 5a p %.2f, q %.2f, p+q %.2f mM ; w = (%.3g, %.3g)\n", N, Cr, σJB, lec.σC_reference..., lec.w...)
        push!(verif, verifier_R(pN, lec))
        # le réglage de Q (σ_w,g) ne se balaye qu'au budget de la vidéo
        σs = N == vid["N_bleu"] ? p0.toml["kalman"]["balayage_sigma_wg"] : [σ0]
        scen = [[(f, α, σ) for (f, α) in [[("rampe", α) for α in p0.toml["derive"]["balayage_alpha"]]; [("saut", 0.0)]]] for σ in σs]
        for (forme, α, σwg) in reduce(vcat, scen)
            p = avec(pN; forme, alpha = α)
            nom = nom_derive(forme, α, p)
            nom_video = forme != "rampe" || N != vid["N_bleu"] ? nothing :
                        α == 0 && σwg == σ0 ? "video_sans_derive" : α == vid["alpha"] && σwg == σ0 ? "video_avec_derive" :
                        α == vid["alpha"] && σwg == vid["sigma_wg_regle"] ? "video_avec_derive_sigma_wg_regle" : nothing
            pour_video = !rapide && nom_video !== nothing
            res = Vector{Any}(undef, nrep)
            @threads for r in 1:nrep
                res[r] = course(p, lec, r; histos = pour_video && r == 1, σwg)
            end
            T1 = res[1][1]
            pour_video && (videos[nom_video] = (T1[:complet], res[1][2], p))
            forme == "saut" && (saut[(N, σwg)] = T1[:complet])
            σwg == σ0 && for m in MODES, k in eachindex(T1[m][:t])
                push!(lignes_pas, (N, nom, "R $(m)", (T1[m][c][k] for c in COLONNES)...))
            end
            for (lect, Cest, m, kal) in (("J/B", :C_ratio, :complet, false), ("5a", :C_5a, :complet, false),
                                      ("Kalman R complet", :C_kalman, :complet, true), ("Kalman R diagonal", :C_kalman, :diagonal, true))
                M = [metriques(res[r][1][m], res[r][1][m][Cest], p; kalman = kal) for r in 1:nrep]
                moy(f) = mean(getfield(x, f) for x in M)
                push!(resume, (N, nom, σwg, lect, moy(:rmse), std(getfield.(M, :rmse)), moy(:biais_fin), std(getfield.(M, :biais_fin)),
                               (moy(f) for f in keys(M[1])[3:end])...))
            end
            k = findlast(r -> r[4] == "Kalman R complet", resume)
            @printf("  %-12s σ_w,g = %.0e  biais fin : J/B %+6.2f  5a %+6.2f  Kalman %+6.2f mM ; P11/err² complet %.2f, diagonal %.2f\n", nom,
                    σwg, resume[k-2][7], resume[k-1][7], resume[k][7], resume[k][12], resume[k+1][12])
        end
    end
    ecrire_csv(joinpath(SORTIES, "courses.csv"), vcat(["N_photons_bleus", "derive", "R"], string.(COLONNES)), lignes_pas)
    ecrire_csv(joinpath(SORTIES, "resume.csv"), ["N_photons_bleus", "derive", "sigma_wg", "lecture", "rmse_mM", "rmse_sd", "biais_30s_mM", "biais_30s_sd",
        "err_max_montee_mM", "echecs_pct", "t_P22_moitie_s", "P11_sur_err2", "biais_g_30s", "t_recup_g_s", "rho12_moyen",
        "var_innov_ratio", "var_innov_5a", "autocorr_innov_ratio", "innov_ratio_montee"], resume)
    ecrire_csv(joinpath(SORTIES, "coherence.csv"), ["N_photons_bleus", "C_JB_mM", "C_5a_mM", "p_5a", "p_titration", "q_5a", "q_titration",
        "w_p", "w_q", "sigmaC_JB", "sigmaC_p", "sigmaC_q", "sigmaC_5a", "sigma_5a_sur_JB"], coherence)
    ecrire_csv(joinpath(SORTIES, "verification_R.csv"), ["N_photons_bleus", "var_ratio_reelle", "var_ratio_fisher", "var_5a_reelle",
        "var_5a_fisher", "cov_reelle", "cov_fisher", "rho_reel", "rho_fisher", "echecs_pct"], verif)
    figure_resume(resume, saut, p0, joinpath(SORTIES, "fig_resume.png"))
    for (nom, (T, H, p)) in videos
        titre = Dict("video_sans_derive" => "sans dérive", "video_avec_derive" => "rampe de 10 %, σ_w,g = 1e-4",
                     "video_avec_derive_sigma_wg_regle" => "rampe de 10 %, σ_w,g = 1e-3")[nom]
        figure_ecarts(T, p, joinpath(SORTIES, "ecarts_" * replace(nom, "video_" => "") * ".png"),
                      @sprintf("Écart mesuré − attendu sur les 900 pas de la vidéo (%s, N = %g)", titre, p.N))
        video(T, H, p, joinpath(SORTIES, nom))
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
