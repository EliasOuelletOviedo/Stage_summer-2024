# Étape 14 — les analyses du simulateur sur le banc à deux canaux du 2026-10-16.
#
# Données : Data/2026-10-16 — deux échantillons (test_0 : 0 mM ; test_131 : 131 mM de Cl⁻,
# d'après les noms de dossiers), excités à 750, 800 et 850 nm, 10 acquisitions de 1 s
# (temps fixé) chacune, les deux canaux enregistrés ensemble (voie 1 bleu, voie 2 jaune) ;
# une IRF par longueur d'onde, vue par les deux détecteurs (diffuseur, 30 s).
#   1. acquisitions : photons, J/B, stabilité sur les 10 répétitions contre Poisson ;
#   2. IRF des deux détecteurs : largeur, queue, retard du jaune sur le bleu ;
#   3. canal bleu : l'ajustement global de la titration (M2, τ liés, x par histogramme) ;
#   4. ajustement conjoint bleu + jaune (le module du changement 3), IRF propre à chaque
#      détecteur : ρ_D, ρ_A, τ_A, ε, retard résiduel, montée ;
#   5. les tests de la note : accepteur émetteur (H2), plage de J/B, fuite, montée, ε(λ) ;
#   6. précision des lectures sur les histogrammes réels amincis (discriminabilité
#      0 contre 131 mM), le jumeau simulé et la prédiction de Fisher ;
#   7. le monde du simulateur reconstruit avec les paramètres mesurés : bornes à 20 mM.
#
#   julia --project=. -t 8 step14_banc.jl [dossier]
#
# Sorties : sorties/banc_2026-10-16/*.csv, *.toml ; figures fig15–fig17.

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, Random, CSV, DataFrames, StableRNGs, TOML, Dates, LinearAlgebra

dossier_donnees = isempty(ARGS) ? joinpath(@__DIR__, "Data", "2026-10-16") : ARGS[1]
nom_banc = "banc_" * basename(normpath(dossier_donnees))
dossier = S.dossier_sorties(nom_banc)
mkpath(dossier)
m, irf_mat, cfg = S.charger_configuration()
W = S.FENETRE_BANC
λs = (750, 800, 850)

hs = S.lire_banc(dossier_donnees)
println(length(hs), " acquisitions lues ; durées ", unique([h.duree_s for h in hs]), " s")

# -----------------------------------------------------------------------------
# 1. Acquisitions et stabilité
# -----------------------------------------------------------------------------
acq = DataFrame()
for ech in unique([h.echantillon for h in hs]), λ in λs
    g = [h for h in hs if h.echantillon == ech && h.λ == λ]
    NB = [sum(h.bleu[W]) for h in g]
    NJ = [sum(h.jaune[W]) for h in g]
    JB = NJ ./ NB
    σp = mean(JB) * sqrt(mean(1 ./ NB .+ 1 ./ NJ))
    pente = cov(collect(1.0:length(JB)), JB) / var(collect(1.0:length(JB)))
    push!(acq, (echantillon = ech, C_mM = g[1].C, lambda_nm = λ, n = length(g), N_bleu = mean(NB), sigma_N_bleu = std(NB),
        sigma_poisson_bleu = sqrt(mean(NB)), N_jaune = mean(NJ), sigma_N_jaune = std(NJ), sigma_poisson_jaune = sqrt(mean(NJ)),
        JB = mean(JB), sigma_JB = std(JB), sigma_poisson_JB = σp, pente_JB_pct_par_acq = 100 * pente / mean(JB)))
end
CSV.write(joinpath(dossier, "acquisitions.csv"), acq)
println("\nAcquisitions (fenêtre $(first(W))–$(last(W))) :")
for r in eachrow(acq)
    @printf("  %-8s %d nm : N_B %7.0f (σ %5.0f, Poisson %4.0f)  N_J %7.0f (σ %5.0f, Poisson %4.0f)  J/B %.4f ± %.4f (Poisson %.4f)  pente %+.3f %%/acq\n",
        r.echantillon, r.lambda_nm, r.N_bleu, r.sigma_N_bleu, r.sigma_poisson_bleu, r.N_jaune, r.sigma_N_jaune, r.sigma_poisson_jaune,
        r.JB, r.sigma_JB, r.sigma_poisson_JB, r.pente_JB_pct_par_acq)
end
for λ in λs
    a0 = acq[(acq.C_mM.==0).&(acq.lambda_nm.==λ), :]
    a1 = acq[(acq.C_mM.!=0).&(acq.lambda_nm.==λ), :]
    @printf("  plage de J/B à %d nm : %.3f / %.3f = ×%.2f\n", λ, a0.JB[1], a1.JB[1], a0.JB[1] / a1.JB[1])
end

# -----------------------------------------------------------------------------
# 2. IRF des deux détecteurs
# -----------------------------------------------------------------------------
irfs = Dict{Tuple{Int,Int},Vector{Float64}}()
pas = NaN
irf_tab = DataFrame()
for λ in λs, voie in 1:2
    irf, p, fond = S.irf_banc(joinpath(dossier_donnees, "IRF_$(λ)_16X.sdt"), voie)
    global pas = p
    irfs[(λ, voie)] = irf
    k = argmax(irf)
    push!(irf_tab, (lambda_nm = λ, voie = voie == 1 ? "bleu" : "jaune", fond_par_canal = fond, canal_pic = k,
        largeur_mi_hauteur_ps = 1000 * S.largeur_mi_hauteur(irf, pas), queue_apres_1ns = sum(irf[k+16:end]),
        retard_sur_bleu_ps = voie == 1 ? 0.0 : 1000 * S.retard_irf(irfs[(λ, 1)], irf, pas)))
end
CSV.write(joinpath(dossier, "irf.csv"), irf_tab)
@printf("\nPas du TAC : %.5f ns ; période = %.2f canaux\n", pas, 12.5 / pas)
for r in eachrow(irf_tab)
    @printf("  IRF %d nm %-5s : fond %.0f/canal, pic au canal %d, largeur à mi-hauteur %.0f ps, queue au-delà de 1 ns %.1f %%, retard sur le bleu %+.1f ps\n",
        r.lambda_nm, r.voie, r.fond_par_canal, r.canal_pic, r.largeur_mi_hauteur_ps, 100 * r.queue_apres_1ns, r.retard_sur_bleu_ps)
end

# -----------------------------------------------------------------------------
# 3. Canal bleu : ajustement global (comme la titration)
# -----------------------------------------------------------------------------
bleu = DataFrame()
glob = Dict{Int,Any}()
for λ in λs
    inst = S.instrument_banc(irfs[(λ, 1)], pas)
    g_hs = [h for h in hs if h.λ == λ]
    histos = [S.Histogramme("banc", h.echantillon, h.C, h.cycle, h.bleu) for h in g_hs]
    g = S.ajuster_global(histos, inst; τ_D0_init = m.p.τ_D0, τ_DA_init = m.p.τ_DA)
    gt = S.ajuster_global(histos, inst; τ_D0_init = g.τ_D0, τ_DA_init = g.τ_DA, t0_par_histo = true)
    glob[λ] = (g = g, gt = gt, histos = histos, inst = inst)
    C = [h.C for h in histos]
    E = 1 - g.τ_DA / g.τ_D0
    x0, x1 = mean(g.x[C.==0]), mean(g.x[C.!=0])
    Nn = [sum(h.bleu[W]) - g.b[i] * length(W) for (i, h) in enumerate(g_hs)]
    I_mes = mean(Nn[C.!=0]) / mean(Nn[C.==0])
    I_pred = (1 - E * x1) / (1 - E * x0)
    A_rel = mean(g.A[C.!=0]) / mean(g.A[C.==0])
    push!(bleu, (lambda_nm = λ, tau_D0 = g.τ_D0, sigma_tau_D0 = g.σ_glob[1], tau_DA = g.τ_DA, sigma_tau_DA = g.σ_glob[2], E = E, t0 = g.t₀,
        x_0mM = x0, sd_x_0mM = std(g.x[C.==0]), x_131mM = x1, sd_x_131mM = std(g.x[C.!=0]), chi2_r = g.chi2_global,
        chi2_r_t0_par_histo = gt.chi2_global, tau_DA_t0_par_histo = gt.τ_DA, gigue_t0_ps = 1000 * std(gt.t0s),
        intensite_mesuree_131_sur_0 = I_mes, intensite_predite_131_sur_0 = I_pred, amplitude_131_sur_0 = A_rel))
end
CSV.write(joinpath(dossier, "ajustement_bleu.csv"), bleu)
tit = TOML.parsefile(S.dossier_sorties("titration", "parametres_mCer3_Topaz.toml"))
x131_pred = S.fraction_active(131.0, m.p)
println("\nCanal bleu, ajustement global M2 (20 histogrammes par longueur d'onde) — titration du 2026-07-20 : τ_D0 = $(round(m.p.τ_D0, digits = 3)), τ_DA = $(round(m.p.τ_DA, digits = 4)), x₀ = $(round(m.p.x₀, digits = 3)), x(131 mM) prédit = $(round(x131_pred, digits = 3))")
for r in eachrow(bleu)
    @printf("  %d nm : τ_D0 = %.3f ± %.3f, τ_DA = %.4f ± %.4f ns (E = %.3f), x(0) = %.3f ± %.3f, x(131) = %.3f ± %.3f, χ²ᵣ = %.3f (t₀ par histo : %.3f, τ_DA %.4f, gigue %.1f ps) ; intensité 131/0 mesurée ×%.3f, prédite ×%.3f ; amplitude ×%.3f\n",
        r.lambda_nm, r.tau_D0, r.sigma_tau_D0, r.tau_DA, r.sigma_tau_DA, r.E, r.x_0mM, r.sd_x_0mM, r.x_131mM, r.sd_x_131mM, r.chi2_r,
        r.chi2_r_t0_par_histo, r.tau_DA_t0_par_histo, r.gigue_t0_ps, r.intensite_mesuree_131_sur_0, r.intensite_predite_131_sur_0, r.amplitude_131_sur_0)
end

# -----------------------------------------------------------------------------
# 4. Ajustement conjoint bleu + jaune (histogrammes cumulés sur les 10 répétitions)
# -----------------------------------------------------------------------------
conj = DataFrame()
conjoints = Dict{Tuple{Int,String},Any}()
for λ in λs
    instB = S.instrument_banc(irfs[(λ, 1)], pas)
    instJ = S.instrument_banc(irfs[(λ, 2)], pas)
    g_hs = [h for h in hs if h.λ == λ]
    Cs = sort(unique([h.C for h in g_hs]))
    hB = [sum(h.bleu for h in g_hs if h.C == c) for c in Cs]
    hY = [sum(h.jaune for h in g_hs if h.C == c) for c in Cs]
    b = glob[λ].g
    for (variante, kw) in (("τ du bleu, montée = τ_DA", (τ_fixes = true, τ_init = (b.τ_D0, b.τ_DA))),
                           ("τ du bleu, montée libre", (τ_fixes = true, τ_init = (b.τ_D0, b.τ_DA), montee_libre = true)),
                           ("tout libre", (τ_fixes = false, τ_init = (b.τ_D0, b.τ_DA))))
        t = @elapsed r = S.ajuster_bleu_jaune(hB, hY, instB, m; inst_jaune = instJ, Δt_libre = true, kw...)
        conjoints[(λ, variante)] = (r = r, hB = hB, hY = hY, Cs = Cs, instB = instB, instJ = instJ)
        # ⟨τ⟩ du jaune (modèle ajusté, sans fond) aux deux concentrations : ∑ t·μ / ∑ μ au-delà du pic
        push!(conj, (lambda_nm = λ, variante = variante, rho_D = r.ρ_D, sigma_rho_D = r.σ[:ρ_D], rho_A = r.ρ_A, sigma_rho_A = r.σ[:ρ_A],
            rho_A_sur_rho_D = r.ρ_A / r.ρ_D, tau_A = r.τ_A, sigma_tau_A = r.σ[:τ_A], epsilon = r.ε, sigma_epsilon = r.σ[:ε],
            retard_ps = 1000 * r.Δt, sigma_retard_ps = 1000 * r.σ[:Δt], tau_montee = r.τ_r, sigma_tau_montee = r.σ[:τ_r],
            tau_D0 = r.τ_D0, tau_DA = r.τ_DA, x_0mM = r.x[1], x_131mM = r.x[end], chi2_r_bleu = r.chi2_bleu, chi2_r_jaune = r.chi2_jaune,
            converge = r.converge, secondes = t))
        @printf("  %d nm, %-26s ρ_D' = %.4f ± %.4f, ρ_A' = %.4f ± %.4f (ρ_A/ρ_D = %.3f), τ_A = %.3f ± %.3f ns, ε = %.3f ± %.3f, retard %+.1f ± %.1f ps, montée %.4f ± %.4f ns, x = %.3f / %.3f, χ²ᵣ bleu %.2f jaune %.2f (%.0f s)\n",
            λ, variante, r.ρ_D, r.σ[:ρ_D], r.ρ_A, r.σ[:ρ_A], r.ρ_A / r.ρ_D, r.τ_A, r.σ[:τ_A], r.ε, r.σ[:ε], 1000 * r.Δt, 1000 * r.σ[:Δt],
            r.τ_r, r.σ[:τ_r], r.x[1], r.x[end], r.chi2_bleu, r.chi2_jaune, t)
    end
end
CSV.write(joinpath(dossier, "ajustement_conjoint.csv"), conj)

# -----------------------------------------------------------------------------
# 5. Les tests de la note
# -----------------------------------------------------------------------------
ws_nom = S.EspaceTravail(m.inst)
JB_nom(C) = (μ = S.formes(C, m, ws_nom); sum(μ[2]) / sum(μ[1]))
plage_nom = JB_nom(0.0) / JB_nom(131.0)
rho_spectres = m.o.η[2, 2] * S.phi_A(m.p) / (m.o.η[2, 1] * m.p.Φ_D)
tests = Dict{String,Any}("plage_JB_predite_nominal_0_131" => plage_nom, "rho_A_sur_rho_D_spectres" => rho_spectres,
    "x_131_predit_isotherme_titration" => x131_pred)
println(@sprintf("\nTests — plage de J/B prédite par le modèle nominal (0 contre 131 mM) : ×%.2f ; sous H2 (accepteur non émissif) : ×1,00 ; ρ_A/ρ_D des spectres et de Φ_A = 0,60 : %.2f", plage_nom, rho_spectres))
for λ in λs
    a0 = acq[(acq.C_mM.==0).&(acq.lambda_nm.==λ), :JB][1]
    a1 = acq[(acq.C_mM.!=0).&(acq.lambda_nm.==λ), :JB][1]
    r = conjoints[(λ, "τ du bleu, montée = τ_DA")].r
    rl = conjoints[(λ, "τ du bleu, montée libre")].r
    q = r.ρ_A / r.ρ_D
    # si les spectres de fuite sont justes : Φ_A effectif ; si Φ_A = 0,60 : η_YD effectif
    ΦA_eff = q / rho_spectres * S.phi_A(m.p)
    ηYD_eff = rho_spectres / q * m.o.η[2, 1]
    tests["$(λ)nm"] = Dict("plage_JB_mesuree" => a0 / a1, "rho_A_sur_rho_D" => q, "Phi_A_eff_si_fuite_spectres" => ΦA_eff,
        "eta_YD_si_Phi_A_catalogue" => ηYD_eff, "facteur_fuite_si_Phi_A_catalogue" => ηYD_eff / m.o.η[2, 1],
        "tau_A" => r.τ_A, "epsilon" => r.ε, "retard_residuel_ps" => 1000r.Δt, "tau_montee_libre" => rl.τ_r,
        "sigma_tau_montee" => rl.σ[:τ_r], "gain_jaune_relatif_si_fuite_spectres" => r.ρ_D / (m.o.η[2, 1] / m.o.η[1, 1]))
    @printf("  %d nm : plage mesurée ×%.2f ; ρ_A/ρ_D = %.2f → Φ_A effectif %.2f (fuite des spectres) ou fuite η_YD × %.2f (Φ_A = 0,60) ; τ_A = %.2f ns ; ε = %.3f ; montée libre %.3f ± %.3f ns (τ_DA %.3f) ; gain relatif du jaune ×%.2f (fuite des spectres)\n",
        λ, a0 / a1, q, ΦA_eff, ηYD_eff / m.o.η[2, 1], r.τ_A, r.ε, rl.τ_r, rl.σ[:τ_r], r.τ_DA, r.ρ_D / (m.o.η[2, 1] / m.o.η[1, 1]))
end
open(joinpath(dossier, "tests.toml"), "w") do io
    TOML.print(io, tests; sorted = true)
end

# -----------------------------------------------------------------------------
# 6. Précision des lectures : histogrammes réels amincis, jumeau, Fisher
# -----------------------------------------------------------------------------
println("\nPrécision : discriminabilité 0 contre 131 mM sur les histogrammes amincis…")
ctxs = S.creer_contextes(m, irf_mat; avec_utilisateur = false)
nthr = length(ctxs)
"""Instances du code de l'utilisateur avec l'IRF du banc (bleu) tronquée à la période (195 canaux)."""
function instances_banc(λ)
    irf = irfs[(λ, 1)][1:195]
    mat = hcat(collect(0:194) .* pas, irf .* 1e5)
    return [S.charger_instance_utilisateur(mat) for _ in 1:nthr]
end
Ns = (1e3, 3e3, 1e4, 3e4)
Rr = 40
prec = DataFrame()
for λ in λs
    cj = conjoints[(λ, "τ du bleu, montée = τ_DA")]
    r = cj.r
    F = S.formes_banc(r, cj.instB, cj.instJ)
    insts_u = instances_banc(λ)
    g_hs = [h for h in hs if h.λ == λ]
    # jumeau : le modèle ajusté de chaque concentration, remis aux photons de chaque répétition
    μ_mod = Dict(c => (zeros(256), zeros(256)) for c in cj.Cs)
    for (i, c) in enumerate(cj.Cs)
        μ_mod[c][1][W] .= r.μ_bleu[i]
        μ_mod[c][2][W] .= r.μ_jaune[i]
    end
    taches = [(src, N, k) for src in (:reel, :jumeau) for N in Ns for k in eachindex(g_hs)]
    lectures = Vector{Any}(undef, length(taches))
    suivant = Threads.Atomic{Int}(1)
    @sync for t in 1:nthr
        Threads.@spawn begin
            iu = insts_u[t]
            bB, bJ = zeros(256), zeros(256)
            while (i = Threads.atomic_add!(suivant, 1)) <= length(taches)
                src, N, k = taches[i]
                h = g_hs[k]
                ic = findfirst(==(h.C), cj.Cs)
                NB = sum(h.bleu[W])
                p = N / NB
                out = NamedTuple[]
                for j in 1:Rr
                    rng = StableRNG(S.graine(16_000_000 + 100k + round(Int, 10log10(N)) + (src === :jumeau ? 50_000 : 0), j + 1000λ))
                    if src === :reel
                        S.amincir_banc!(bB, bJ, h, p, rng)
                    else
                        μB = μ_mod[h.C][1] .* (N / sum(μ_mod[h.C][1][W]))
                        μJ = μ_mod[h.C][2] .* (N / sum(μ_mod[h.C][1][W]))
                        S.tirer_deux!(bB, bJ, μB, μJ, rng)
                    end
                    sB, sJ = sum(bB[W]), sum(bJ[W])
                    fB, fJ = r.b_bleu[ic] * N / sum(r.μ_bleu[ic]), r.b_jaune[ic] * N / sum(r.μ_bleu[ic])
                    x5d = try S.lire_forme(F, bB, bJ, fB, fJ; mode = :commune) catch; NaN end
                    x5a = try S.lire_forme(F, bB, bJ, fB, fJ; mode = :separees) catch; NaN end
                    x4 = try S.lire_forme(F, bB, bJ, fB, fJ; mode = :bleu) catch; NaN end
                    τm = try S.ajuster_mono_banc(iu, copy(bB); pas) catch; NaN end
                    push!(out, (source = String(src), lambda_nm = λ, N_bleu = N, C_mM = h.C, cycle = h.cycle, tirage = j,
                        ratio = sJ / sB, x_5d = x5d, x_5a = x5a, x_4 = x4, tau_mono = τm))
                end
                lectures[i] = out
            end
        end
    end
    d = DataFrame(reduce(vcat, lectures))
    CSV.write(joinpath(dossier, "lectures_amincies_$(λ)nm.csv"), d)
    # Fisher : σ de chaque lecture au point ajusté, remis à l'échelle de N
    for N in Ns, src in ("reel", "jumeau")
        dd = d[(d.source.==src).&(d.N_bleu.==N), :]
        a, b2 = dd[dd.C_mM.==0, :], dd[dd.C_mM.!=0, :]
        push!(prec, (source = src, lambda_nm = λ, N_bleu = N,
            d_ratio = S.discriminabilite(a.ratio, b2.ratio), d_5d = S.discriminabilite(a.x_5d, b2.x_5d),
            d_5a = S.discriminabilite(a.x_5a, b2.x_5a), d_4 = S.discriminabilite(a.x_4, b2.x_4),
            d_mono = S.discriminabilite(a.tau_mono, b2.tau_mono), echec_mono = mean(isnan.(dd.tau_mono))))
    end
    # prédiction de Fisher (modèle ajusté, fonds connus)
    for N in Ns
        σx = Dict{Symbol,Vector{Float64}}()
        Δ = Dict{Symbol,Float64}()
        for (mode, params) in ((:commune, 2), (:separees, 3), (:bleu, 2))
            v = map(eachindex(cj.Cs)) do ic
                x = r.x[ic]
                Sc = N / sum(r.μ_bleu[ic])
                fB, fJ = r.b_bleu[ic] * Sc, r.b_jaune[ic] * Sc
                A = r.A[ic] * Sc
                μf(θ) = mode === :bleu ? θ[2] .* (F.uB .+ θ[1] .* F.vB) .+ fB :
                    vcat(θ[2] .* (F.uB .+ θ[1] .* F.vB) .+ fB, (mode === :commune ? θ[2] : θ[3]) .* (F.uJ .+ θ[1] .* F.vJ) .+ fJ)
                θ = mode === :separees ? [x, A, A] : [x, A]
                I = S.fisher_numerique(μf, θ)
                sqrt(inv(Symmetric(I))[1, 1])
            end
            σx[mode] = v
            Δ[mode] = abs(r.x[end] - r.x[1])
        end
        Rv = [sum(r.μ_jaune[ic]) / sum(r.μ_bleu[ic]) for ic in eachindex(cj.Cs)]
        σR = [Rv[ic] * sqrt(1 / N + 1 / (N * Rv[ic])) for ic in eachindex(cj.Cs)]
        push!(prec, (source = "Fisher", lambda_nm = λ, N_bleu = N,
            d_ratio = abs(Rv[end] - Rv[1]) / sqrt(mean(σR .^ 2)), d_5d = Δ[:commune] / sqrt(mean(σx[:commune] .^ 2)),
            d_5a = Δ[:separees] / sqrt(mean(σx[:separees] .^ 2)), d_4 = Δ[:bleu] / sqrt(mean(σx[:bleu] .^ 2)), d_mono = NaN,
            echec_mono = NaN))
    end
end
CSV.write(joinpath(dossier, "precision.csv"), prec)
println("Discriminabilité d' (0 contre 131 mM) :")
for r in eachrow(sort(prec, [:lambda_nm, :N_bleu, :source]))
    @printf("  %d nm N_B = %6.0f %-8s ratio %7.1f  5d %7.1f  5a %6.1f  4 bleu %6.1f  mono %6.1f   (ratio/5a %.1f, ratio/4 %.1f, ratio/mono %.1f)\n",
        r.lambda_nm, r.N_bleu, r.source, r.d_ratio, r.d_5d, r.d_5a, r.d_4, r.d_mono, r.d_ratio / r.d_5a, r.d_ratio / r.d_4, r.d_ratio / r.d_mono)
end

# -----------------------------------------------------------------------------
# 7. Le monde du simulateur reconstruit avec les paramètres mesurés (800 nm)
# -----------------------------------------------------------------------------
r8 = conjoints[(800, "τ du bleu, montée = τ_DA")].r
b8 = glob[800].g
Kn = m.p.K_D^m.p.n
rr = Kn / (Kn + 131.0^m.p.n)
x0m, x1m = r8.x[1], r8.x[end]
xsat_m = (x1m - x0m * rr) / (1 - rr)
ηBD = m.o.η[1, 1]
o_m = S.remplacer(m.o; η = S.SMatrix{2,2}([ηBD m.o.η[1, 2]; r8.ρ_D*ηBD r8.ρ_A*ηBD*m.p.Φ_D/S.phi_A(m.p)]), ε = r8.ε)
p_m = S.remplacer(m.p; τ_D0 = r8.τ_D0, τ_DA = r8.τ_DA, τ_A = r8.τ_A, x₀ = x0m, x_sat = max(xsat_m, 0.001))
inst_m = S.InstrumentParams(n_bins = 256, T_rep = 12.5, pas = pas, irf = irfs[(800, 1)], t₀ = r8.t₀, fenetre = W, Δt_jaune = r8.Δt)
mm = S.Modele(p_m, o_m, inst_m)
println("\nMonde mesuré (800 nm) : ", S.decrire(mm))
wsm = S.EspaceTravail(inst_m)
bornes = DataFrame()
for (nomm, mod, wsx) in (("nominal (titration, spectres)", m, ws_nom), ("mesuré au banc", mm, wsm))
    for N in (1e3, 1e4), C in (5.0, 20.0, 40.0)
        push!(bornes, (monde = nomm, N_photons = N, C_mM = C,
            absolue = S.crlb_mM(:absolue, C, mod, N, :A, wsx), totaux = S.crlb_mM(:totaux, C, mod, N, :A, wsx),
            conjoint = S.crlb_mM(:conjoint, C, mod, N, :A, wsx), bleu = S.crlb_mM(:bleu, C, mod, N, :A, wsx),
            JB = (μ = S.attendus(C, mod, N, :A, wsx); sum(μ[2]) / sum(μ[1]))))
    end
end
CSV.write(joinpath(dossier, "bornes_monde_mesure.csv"), bornes)
println("Bornes de Cramér-Rao (mM, normalisation A) :")
for r in eachrow(bornes)
    @printf("  %-30s N = %-6g C = %4.0f mM : tout %.2f  totaux J/B %.2f  forme 2 canaux %.2f  bleu %.2f   (J/B = %.2f)\n",
        r.monde, r.N_photons, r.C_mM, r.absolue, r.totaux, r.conjoint, r.bleu, r.JB)
end
open(joinpath(dossier, "monde_mesure.toml"), "w") do io
    TOML.print(io, Dict("date" => string(Dates.today()), "source" => dossier_donnees, "lambda_nm" => 800,
        "tau_D0" => r8.τ_D0, "tau_DA" => r8.τ_DA, "tau_A" => r8.τ_A, "epsilon" => r8.ε, "dt_jaune_residuel_ns" => r8.Δt,
        "rho_D_avec_gain" => r8.ρ_D, "rho_A_avec_gain" => r8.ρ_A, "x_0mM" => x0m, "x_131mM" => x1m,
        "x_sat_avec_K_n_titration" => xsat_m, "K_D_titration" => m.p.K_D, "n_titration" => m.p.n,
        "eta_YD_effectif_avec_gain" => r8.ρ_D * ηBD, "eta_YA_effectif_avec_gain" => r8.ρ_A * ηBD * m.p.Φ_D / S.phi_A(m.p)); sorted = true)
end
println("\nécrit dans ", dossier)
