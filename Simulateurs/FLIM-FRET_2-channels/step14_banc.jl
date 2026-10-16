# Étape 14 — banc à deux canaux du 2026-10-16 : une autre version de SCLM, traitée seule.
#
# Données : Data/2026-10-16 — deux échantillons (test_0 et test_131 : 0 et 131 mM de Cl⁻
# d'après les noms de dossiers), excités à 750, 800 et 850 nm, 10 acquisitions de 1 s
# (temps fixé) chacune, les deux canaux ensemble (voie 1 bleu, voie 2 jaune) ; une IRF par
# longueur d'onde, vue par les deux détecteurs (30 s). Cette version de SCLM n'est pas
# comparée aux autres acquisitions : tout ce qui suit vient de ces fichiers seuls.
#   1. les 10 acquisitions : photons, J/B, stabilité d'une acquisition à l'autre ;
#   2. les IRF des deux détecteurs : largeur, queue lente, retard du jaune sur le bleu ;
#   3. la somme des 10 acquisitions de chaque échantillon (10 fois plus de photons) :
#      ce que contiennent les déclins (trois exponentielles libres par canal), phasors,
#      rapport J/B, et les valeurs de référence de chaque technique à pleins photons ;
#   4. la stabilité des techniques d'analyse quand les photons manquent : chaque somme
#      est amincie par tirage binomial de 10² à 10⁵ photons bleus (le jaune suit), et
#      chaque technique lit chaque tirage. Lectures normalisées sur l'écart entre les deux
#      échantillons (0 à 0 mM, 1 à 131 mM, valeurs à pleins photons) : bruit, biais,
#      échecs, discriminabilité ; jumeau simulé et borne de Fisher pour comparaison ;
#   5. la longueur d'onde d'excitation.
# Techniques (les estimateurs du simulateur) : 1 rapport J/B ; 2 ton mono sur le bleu ;
# 2J ton mono sur le jaune ; 3 bi-exponentiel à τ libres (bleu) ; 4 bi-exponentiel à τ
# figés (bleu) ; 5a forme des deux canaux, échelles séparées ; 5d forme des deux canaux,
# échelle commune ; 6 phasor (bleu).
#
#   julia --project=. -t 8 step14_banc.jl [dossier]
#
# Sorties : sorties/banc_2026-10-16/*.csv, *.toml ; sorties/figures/fig15_banc_2026-10-16_*.png

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, Random, CSV, DataFrames, StableRNGs, TOML, Dates, LinearAlgebra

dossier_donnees = isempty(ARGS) ? joinpath(@__DIR__, "Data", "2026-10-16") : ARGS[1]
nom_banc = "banc_" * basename(normpath(dossier_donnees))
dossier = S.dossier_sorties(nom_banc)
mkpath(dossier)
# fichiers des premiers essais (comparaisons à d'autres acquisitions, abandonnées : cette version se traite seule)
foreach(f -> rm(joinpath(dossier, f)), filter(f -> startswith(f, "lectures_amincies_") || f in ("ajustement_conjoint.csv",
    "ajustement_bleu.csv", "ajustement_M2_bleu.csv", "tests.toml", "bornes_monde_mesure.csv", "monde_mesure.toml", "phasors.csv",
    "precision.csv", "trois_exponentielles.csv"), readdir(dossier)))
W = S.FENETRE_BANC
λs = (750, 800, 850)
etats = ["test_0" => 0.0, "test_131" => 131.0]
hs = S.lire_banc(dossier_donnees)
println(length(hs), " acquisitions lues ; durées ", unique([h.duree_s for h in hs]), " s")
resume = Dict{String,Any}("date" => string(Dates.today()), "donnees" => dossier_donnees)

# -----------------------------------------------------------------------------
# 1. Les 10 acquisitions : stabilité d'une acquisition à l'autre
# -----------------------------------------------------------------------------
acq = DataFrame()
for (ech, C) in etats, λ in λs
    g = [h for h in hs if h.echantillon == ech && h.λ == λ]
    NB = [sum(h.bleu[W]) for h in g]
    NJ = [sum(h.jaune[W]) for h in g]
    JB = NJ ./ NB
    σp = mean(JB) * sqrt(mean(1 ./ NB .+ 1 ./ NJ))
    pente = cov(collect(1.0:length(JB)), JB) / var(collect(1.0:length(JB)))
    push!(acq, (echantillon = ech, C_mM = C, lambda_nm = λ, n = length(g), N_bleu = mean(NB), sigma_N_bleu = std(NB),
        sigma_poisson_bleu = sqrt(mean(NB)), N_jaune = mean(NJ), sigma_N_jaune = std(NJ), sigma_poisson_jaune = sqrt(mean(NJ)),
        JB = mean(JB), sigma_JB = std(JB), sigma_poisson_JB = σp, pente_JB_pct_par_acq = 100 * pente / mean(JB)))
end
CSV.write(joinpath(dossier, "acquisitions.csv"), acq)
println("\n1. Les 10 acquisitions de 1 s (fenêtre $(first(W))–$(last(W))) :")
for r in eachrow(acq)
    @printf("  %-8s %d nm : N_B %7.0f (σ %5.0f = %.1f × Poisson)  N_J %7.0f (σ %5.0f = %.1f × Poisson)  J/B %.4f ± %.4f (%.2f × Poisson)  pente %+.3f %%/acq\n",
        r.echantillon, r.lambda_nm, r.N_bleu, r.sigma_N_bleu, r.sigma_N_bleu / r.sigma_poisson_bleu, r.N_jaune, r.sigma_N_jaune,
        r.sigma_N_jaune / r.sigma_poisson_jaune, r.JB, r.sigma_JB, r.sigma_JB / r.sigma_poisson_JB, r.pente_JB_pct_par_acq)
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
        queue_apres_4ns = sum(irf[k+63:end]), retard_sur_bleu_ps = voie == 1 ? 0.0 : 1000 * S.retard_irf(irfs[(λ, 1)], irf, pas)))
end
CSV.write(joinpath(dossier, "irf.csv"), irf_tab)
@printf("\n2. IRF — pas du TAC %.5f ns, période = %.2f canaux\n", pas, 12.5 / pas)
for r in eachrow(irf_tab)
    @printf("  %d nm %-5s : largeur à mi-hauteur %.0f ps, queue au-delà de 1 ns %.1f %%, de 4 ns %.1f %%, retard sur le bleu %+.1f ps\n",
        r.lambda_nm, r.voie, r.largeur_mi_hauteur_ps, 100 * r.queue_apres_1ns, 100 * r.queue_apres_4ns, r.retard_sur_bleu_ps)
end

# -----------------------------------------------------------------------------
# 3. La somme des 10 acquisitions de chaque échantillon
# -----------------------------------------------------------------------------
somme(ech, λ, voie) = sum(voie == 1 ? h.bleu : h.jaune for h in hs if h.echantillon == ech && h.λ == λ)
sommes = Dict((ech, λ) => (somme(ech, λ, 1), somme(ech, λ, 2)) for (ech, _) in etats for λ in λs)
mexp = DataFrame()
lisses = Dict{Tuple{Int,Int},Any}()
deux = Dict{Int,Any}()
for λ in λs, (voie, nomv) in ((1, "bleu"), (2, "jaune"))
    inst = S.instrument_banc(irfs[(λ, voie)], pas)
    H = [sommes[(e, λ)][voie] for (e, _) in etats]
    r = S.ajuster_multiexp(H, inst; τ_init = (0.1, 0.9, 3.0), montee = voie == 2)
    lisses[(λ, voie)] = r
    if voie == 1
        deux[λ] = S.ajuster_multiexp(H, inst; τ_init = (0.6, 2.5))     # base de l'estimateur 4
    end
    for (s, (_, C)) in enumerate(etats)
        push!(mexp, (lambda_nm = λ, canal = nomv, C_mM = C, N = sum(H[s][W]), tau_1 = r.τ[1], tau_2 = r.τ[2], tau_3 = r.τ[3],
            fraction_1 = r.fractions_intensite[s][1], fraction_2 = r.fractions_intensite[s][2], fraction_3 = r.fractions_intensite[s][3],
            chi2_r = r.chi2[s], t0 = r.t₀))
    end
end
CSV.write(joinpath(dossier, "sommes_trois_exponentielles.csv"), mexp)
println("\n3. Sommes des 10 acquisitions — trois exponentielles libres (τ communs aux deux échantillons ; fractions d'intensité) :")
for g in groupby(mexp, [:lambda_nm, :canal]; sort = false)
    @printf("  %d nm %-5s (N = %.2g / %.2g) : τ = %.3f, %.3f, %.3f ns ; 0 mM : %.2f %.2f %.2f ; 131 mM : %.2f %.2f %.2f ; χ²ᵣ %.1f / %.1f\n",
        g.lambda_nm[1], g.canal[1], g.N[1], g.N[2], g.tau_1[1], g.tau_2[1], g.tau_3[1], g.fraction_1[1], g.fraction_2[1], g.fraction_3[1],
        g.fraction_1[2], g.fraction_2[2], g.fraction_3[2], g.chi2_r[1], g.chi2_r[2])
end
println("  bi-exponentiel du bleu (base de l'estimateur 4) :")
for λ in λs
    d = deux[λ]
    @printf("    %d nm : τ = %.3f et %.3f ns ; fraction d'amplitude courte 0 mM %.3f, 131 mM %.3f ; χ²ᵣ %.1f / %.1f\n", λ, d.τ[1], d.τ[2],
        d.amplitudes[1][1] / sum(d.amplitudes[1]), d.amplitudes[2][1] / sum(d.amplitudes[2]), d.chi2[1], d.chi2[2])
end
ph = DataFrame()
for λ in λs, (ech, C) in etats, (voie, nomv) in ((1, "bleu"), (2, "jaune"))
    z = S.phasor_brut(sommes[(ech, λ)][voie], pas) / S.phasor_brut(irfs[(λ, voie)], pas; fenetre = 1:256)
    push!(ph, (lambda_nm = λ, C_mM = C, canal = nomv, g = real(z), s = imag(z), tau_phase_ns = imag(z) / real(z) / (2π / 12.5)))
end
CSV.write(joinpath(dossier, "sommes_phasors.csv"), ph)
println("  durée de vie de phase (IRF de chaque détecteur corrigée) :")
for g in groupby(ph, [:lambda_nm, :canal]; sort = false)
    @printf("    %d nm %-5s : 0 mM τ_φ = %.2f ns, 131 mM τ_φ = %.2f ns\n", g.lambda_nm[1], g.canal[1], g.tau_phase_ns[1], g.tau_phase_ns[2])
end
for λ in λs
    JB = [sum(sommes[(e, λ)][2][W]) / sum(sommes[(e, λ)][1][W]) for (e, _) in etats]
    @printf("  J/B des sommes à %d nm : %.4f à 0 mM, %.4f à 131 mM (plage ×%.2f)\n", λ, JB[1], JB[2], JB[1] / JB[2])
    resume["JB_$(λ)nm"] = JB
end

# -----------------------------------------------------------------------------
# 4. Stabilité des techniques quand les photons manquent
# -----------------------------------------------------------------------------
const TECH = ["1", "2", "2J", "3", "4", "5a", "5d", "6"]
const NOM_TECH = Dict("1" => "1 rapport J/B", "2" => "2 mono, bleu (ton code)", "2J" => "2J mono, jaune (ton code)",
    "3" => "3 bi-exp. τ libres (bleu)", "4" => "4 bi-exp. τ figés (bleu)", "5a" => "5a forme 2 canaux, échelles séparées",
    "5d" => "5d forme 2 canaux, échelle commune", "6" => "6 phasor (bleu)")
nthr = Threads.nthreads()
inst_u(λ, voie) = [S.charger_instance_utilisateur(hcat(collect(0:194) .* pas, irfs[(λ, voie)][1:195] .* 1e5)) for _ in 1:nthr]
Ns = (1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5)
R = 400
tirages = DataFrame()
refs = DataFrame()
for λ in λs
    println("\n4. ", λ, " nm : lectures des sommes, puis des tirages amincis…")
    rB, rJ = lisses[(λ, 1)], lisses[(λ, 2)]
    G = S.gabarits_banc(rB, rJ)
    d2 = deux[λ]
    instB = S.instrument_banc(irfs[(λ, 1)], pas)
    wsB = [S.EspaceTravail(instB) for _ in 1:nthr]
    b2 = zeros(256); b1 = zeros(256)
    S.base!(b1, wsB[1], d2.τ[1], d2.t₀); S.base!(b2, wsB[1], d2.τ[2], d2.t₀)
    F4 = S.FormesBanc(b2[W], b1[W] .- b2[W], zeros(length(W)), zeros(length(W)), W)
    P0 = S.phasor_brut(vcat(zeros(first(W) - 1), rB.μ[1] .- rB.fonds[1], zeros(256 - last(W))), pas)
    P1 = S.phasor_brut(vcat(zeros(first(W) - 1), rB.μ[2] .- rB.fonds[2], zeros(256 - last(W))), pas)
    iuB, iuJ = inst_u(λ, 1), inst_u(λ, 2)
    # lecture de toutes les techniques sur un couple d'histogrammes ; fonds connus remis à l'échelle
    function lire_tout(bB, bJ, s, N, t)
        Sc = N / sum(rB.μ[s])
        fB, fJ = rB.fonds[s] * Sc, rJ.fonds[s] * Sc
        f4 = d2.fonds[s] * Sc
        sB, sJ = sum(bB[W]), sum(bJ[W])
        lit(f) = try f() catch; NaN end
        return Dict("1" => sB > 0 ? sJ / sB : NaN,
            "2" => lit(() -> S.ajuster_mono_banc(iuB[t], copy(bB); pas)),
            "2J" => lit(() -> S.ajuster_mono_banc(iuJ[t], copy(bJ); pas)),
            "3" => lit(() -> S.lire_biexp_libre(bB, wsB[t], W, (d2.τ[1], d2.τ[2]), d2.t₀, f4)),
            "4" => lit(() -> S.lire_forme(F4, bB, bJ, f4, 0.0; mode = :bleu)),
            "5a" => lit(() -> S.lire_forme(G, bB, bJ, fB, fJ; mode = :separees)),
            "5d" => lit(() -> S.lire_forme(G, bB, bJ, fB, fJ; mode = :commune)),
            "6" => sB > 0 ? S.lecture_phasor(bB, P0, P1, pas) : NaN)
    end
    # références à pleins photons (les sommes)
    ref = Dict{String,Vector{Float64}}(k => Float64[] for k in TECH)
    for (s, (e, _)) in enumerate(etats)
        hB, hJ = sommes[(e, λ)]
        v = lire_tout(hB, hJ, s, sum(hB[W]), 1)
        for k in TECH
            push!(ref[k], v[k])
        end
    end
    for k in TECH
        push!(refs, (lambda_nm = λ, technique = k, lecture_0mM = ref[k][1], lecture_131mM = ref[k][2]))
        @printf("  %-38s %10.4f à 0 mM, %10.4f à 131 mM\n", NOM_TECH[k], ref[k][1], ref[k][2])
    end
    # tirages : réels (sommes amincies) et jumeau (Poisson du lissage à trois exponentielles)
    lisse(s) = (vcat(zeros(first(W) - 1), rB.μ[s], zeros(256 - last(W))), vcat(zeros(first(W) - 1), rJ.μ[s], zeros(256 - last(W))))
    taches = [(src, N, s, lot) for src in (:reel, :jumeau) for N in Ns for s in 1:2 for lot in 1:8]
    sorties = Vector{Any}(undef, length(taches))
    suivant = Threads.Atomic{Int}(1)
    @sync for t in 1:nthr
        Threads.@spawn begin
            bB, bJ = zeros(256), zeros(256)
            while (i = Threads.atomic_add!(suivant, 1)) <= length(taches)
                src, N, s, lot = taches[i]
                e = etats[s][1]
                hB, hJ = sommes[(e, λ)]
                p = N / sum(hB[W])
                μB, μJ = lisse(s)
                cB = N / sum(μB[W])
                out = NamedTuple[]
                for j in (lot-1)*(R÷8)+1:lot*(R÷8)
                    rng = StableRNG(S.graine(17_000_000 + 1000s + round(Int, 10log10(N)) + (src === :jumeau ? 500 : 0) + 10λ, j))
                    if src === :reel
                        S.amincir!(bB, hB, p, rng)
                        S.amincir!(bJ, hJ, p, rng)
                    else
                        S.tirer_deux!(bB, bJ, μB .* cB, μJ .* cB, rng)
                    end
                    v = lire_tout(bB, bJ, s, N, t)
                    push!(out, (source = String(src), lambda_nm = λ, N_bleu = N, C_mM = etats[s][2], tirage = j,
                        N_jaune = sum(bJ[W]), (Symbol("o_" * k) => v[k] for k in TECH)...))
                end
                sorties[i] = out
            end
        end
    end
    global tirages = vcat(tirages, DataFrame(reduce(vcat, sorties)); cols = :union)
end
CSV.write(joinpath(dossier, "tirages_amincis.csv"), tirages)
CSV.write(joinpath(dossier, "references_pleins_photons.csv"), refs)

# statistiques : lectures normalisées f = (o − o₀)/(o₁₃₁ − o₀), o₀ et o₁₃₁ à pleins photons
stab = DataFrame()
for g in groupby(tirages, [:source, :lambda_nm, :N_bleu]; sort = true), k in TECH
    r = refs[(refs.lambda_nm.==g.lambda_nm[1]).&(refs.technique.==k), :]
    o0, o1 = r.lecture_0mM[1], r.lecture_131mM[1]
    col = Symbol("o_" * k)
    f0 = (g[g.C_mM.==0, col] .- o0) ./ (o1 - o0)
    f1 = (g[g.C_mM.!=0, col] .- o0) ./ (o1 - o0)
    ok0, ok1 = filter(isfinite, f0), filter(isfinite, f1)
    σ0, σ1 = std(ok0), std(ok1)
    rob(v) = 1.4826 * median(abs.(v .- median(v)))
    push!(stab, (source = g.source[1], lambda_nm = g.lambda_nm[1], N_bleu = g.N_bleu[1], technique = k, nom = NOM_TECH[k],
        sigma_0mM = σ0, sigma_131mM = σ1, sigma = sqrt((σ0^2 + σ1^2) / 2), sigma_robuste = sqrt((rob(ok0)^2 + rob(ok1)^2) / 2),
        biais_0mM = mean(ok0), biais_131mM = mean(ok1) - 1, echec = 1 - (length(ok0) + length(ok1)) / (length(f0) + length(f1)),
        d_prime = abs(mean(ok1) - mean(ok0)) / sqrt((σ0^2 + σ1^2) / 2)))
end
# borne de Fisher (bruit de Poisson seul, modèle lissé, fonds connus) pour 1, 4, 5a, 5d
for λ in λs, N in Ns
    rB, rJ = lisses[(λ, 1)], lisses[(λ, 2)]
    G = S.gabarits_banc(rB, rJ)
    R0 = [sum(rJ.μ[s]) / sum(rB.μ[s]) for s in 1:2]
    σR = [R0[s] * sqrt(1 / N + 1 / (N * R0[s])) for s in 1:2]
    σf = Dict("1" => sqrt(mean(σR .^ 2)) / abs(R0[2] - R0[1]))
    for (k, mode) in (("5a", :separees), ("5d", :commune))
        v = map(1:2) do s
            Sc = N / sum(rB.μ[s])
            bB0, bJ0 = rB.fonds[s] * Sc, rJ.fonds[s] * Sc
            μf(θ) = vcat(θ[2] .* (G.uB .+ θ[1] .* G.vB) .+ bB0, (mode === :commune ? θ[2] : θ[3]) .* (G.uJ .+ θ[1] .* G.vJ) .+ bJ0)
            θ = mode === :separees ? [s - 1.0, Float64(N), Float64(N)] : [s - 1.0, Float64(N)]
            sqrt(inv(Symmetric(S.fisher_numerique(μf, θ)))[1, 1])
        end
        σf[k] = sqrt(mean(v .^ 2))
    end
    for (k, v) in σf
        push!(stab, (source = "Fisher", lambda_nm = λ, N_bleu = N, technique = k, nom = NOM_TECH[k], sigma_0mM = NaN, sigma_131mM = NaN,
            sigma = v, sigma_robuste = NaN, biais_0mM = NaN, biais_131mM = NaN, echec = NaN, d_prime = 1 / v))
    end
end
CSV.write(joinpath(dossier, "stabilite.csv"), stab)
println("\nStabilité à 800 nm (lectures normalisées : 0 → 0 mM, 1 → 131 mM ; σ et biais en fraction de l'écart) :")
for N in Ns
    @printf("  N_B = %-6g", N)
    for k in TECH
        d = stab[(stab.source.=="reel").&(stab.lambda_nm.==800).&(stab.N_bleu.==N).&(stab.technique.==k), :]
        @printf("  %s σ %.3f b %+.3f e %2.0f%%", k, d.sigma[1], (d.biais_0mM[1] + d.biais_131mM[1]) / 2, 100 * d.echec[1])
    end
    println()
end
println("\nDiscriminabilité d' à 10⁴ photons bleus (réel / jumeau / Fisher) :")
for λ in λs
    @printf("  %d nm :", λ)
    for k in TECH
        v(src) = (d = stab[(stab.source.==src).&(stab.lambda_nm.==λ).&(stab.N_bleu.==1e4).&(stab.technique.==k), :]; nrow(d) == 1 ? d.d_prime[1] : NaN)
        @printf("  %s %.1f/%.1f/%.1f", k, v("reel"), v("jumeau"), v("Fisher"))
    end
    println()
end
open(joinpath(dossier, "resume.toml"), "w") do io
    TOML.print(io, resume; sorted = true)
end

include(joinpath(@__DIR__, "src", "Figures.jl"))
using .Figures, CairoMakie
set_theme!(Figures.theme_sim())
if isdefined(Figures, :figure_banc_stabilite)
    Figures.figure_banc_declins(dossier, dossier_donnees, S.dossier_sorties("figures", "fig15_banc_2026-10-16_declins.png"))
    Figures.figure_banc_stabilite(dossier, S.dossier_sorties("figures", "fig16_banc_2026-10-16_stabilite.png"))
end
