# Étape 2 — validation interne : les huit tests du plan, plus trois tests de
# fidélité (noyau contre intégrale exacte, FFT contre convolution directe,
# adaptateur contre tes propres résultats). Tant qu'ils ne passent pas tous,
# aucun résultat du simulateur ne vaut la peine d'être interprété.
#
#   julia --project=. -t 8 step2_validation.jl
#
# Sortie : sorties/validation/rapport_validation.md

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, LinearAlgebra, Random, StableRNGs, ForwardDiff, QuadGK, StaticArrays, CSV, DataFrames

m, irf_mat, cfg = S.charger_configuration()
ws = S.EspaceTravail(m.inst)
W = m.inst.fenetre
Δ, T, n = S.pas_temps(m.inst), m.inst.T_rep, m.inst.n_bins
println(S.decrire(m))

resultats = Tuple{String,String,String,Bool}[]
function verdict!(test, critere, valeur, ok)
    push!(resultats, (test, critere, valeur, ok))
    @printf("%-4s %-62s %s\n", ok ? "OK" : "ÉCHEC", test, valeur)
end
relerr(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), 1e-300)
x20 = S.fraction_active(20.0, m.p)

# -----------------------------------------------------------------------------
# 1. Aller-retour sans bruit
# -----------------------------------------------------------------------------
# (a) modèle complet, tous les paramètres libres : θ = (x, τ_D0, τ_DA, τ_A, t₀, log S_B, log S_Y)
let
    μB, μY = S.attendus(20.0, m, 1e6, :C, ws)
    nobs = vcat(μB[W], μY[W])
    mW = length(W)
    function modele!(μ, J, θ)
        mm = S.remplacer(m; p = S.remplacer(m.p; τ_D0 = θ[2], τ_DA = θ[3], τ_A = θ[4]))
        mB, mY, dB, dY = S.formes_et_derivees(θ[1], mm, ws; t₀ = θ[5])
        SB, SY = exp(θ[6]), exp(θ[7])
        for (i, k) in enumerate(W)
            μ[i] = SB * mB[k]
            μ[mW+i] = SY * mY[k]
            for j in 1:5
                J[i, j] = SB * dB[k, j]
                J[mW+i, j] = SY * dY[k, j]
            end
            J[i, 6], J[i, 7] = μ[i], 0.0
            J[mW+i, 6], J[mW+i, 7] = 0.0, μ[mW+i]
        end
    end
    mB, mY = S.formes(20.0, m, ws)
    vrai = [x20, m.p.τ_D0, m.p.τ_DA, m.p.τ_A, 0.0, log(sum(μB) / sum(mB)), log(sum(μY) / sum(mY))]
    depart = vrai .* [1.1, 1.03, 0.95, 1.04, 1, 1, 1] .+ [0, 0, 0, 0, 0.02, 0.1, -0.1]
    r = S.lm_poisson(modele!, depart, [-1, 0.5, 0.01, 0.5, -1, -Inf, -Inf], [2, 10, 3, 10, 1, Inf, Inf], nobs;
        maxiter = 1000, tol_dec = 1e-24)
    e = maximum(abs.(r.θ[1:5] .- vrai[1:5]) ./ max.(abs.(vrai[1:5]), [1, 1, 1, 1, Δ]))
    verdict!("1a aller-retour sans bruit, modèle complet (7 paramètres libres)", "écart relatif max < 1e-8",
        @sprintf("%.1e (%d itérations)", e, r.iterations), e < 1e-8)
end
# (b) estimateurs à τ figés sur des données sans bruit : x exactement
let
    ctx = S.Contexte(ws, nothing)
    for (k, e) in (("4 (canal bleu sans fuite d'accepteur)", nothing), ("5a", S.Conjoint(m, ws)), ("5d", S.Conjoint(m, ws; commune = true)))
        if e === nothing
            mm = S.remplacer(m; o = S.remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] 0.0; m.o.η[2, 1] m.o.η[2, 2]])))
            est = S.BiExpFige(mm, ws)
            μB, μY = S.attendus(20.0, mm, 1e4, :C, ws)
        else
            est = e
            μB, μY = S.attendus(20.0, m, 1e4, :C, ws)
        end
        o, _ = S.observable(est, μB, μY, ctx)
        verdict!("1b estimateur $k sans bruit rend x", "écart sur x < 1e-10", @sprintf("%.1e", abs(o - x20)), abs(o - x20) < 1e-10)
    end
end

# -----------------------------------------------------------------------------
# 2. Conservation des photons
# -----------------------------------------------------------------------------
let
    inst_plein = S.remplacer(m.inst; fenetre = 1:n)
    p = m.p
    k0 = 1 / p.τ_D0
    donneur = S.ExpSum([p.Φ_D * k0 * x20, p.Φ_D * k0 * (1 - x20)], [p.τ_DA, p.τ_D0])
    s = p.Φ_A / p.τ_A * x20 * (1 / p.τ_DA - k0) / (1 / p.τ_DA - 1 / p.τ_A)
    sens = S.ExpSum([-s, s], [p.τ_DA, p.τ_A])
    eD = abs(sum(S.forward(donneur, inst_plein, ws)) / S.integrale_donneur(x20, p) - 1)
    eA = abs(sum(S.forward(sens, inst_plein, ws)) / S.integrale_sensibilisee(x20, p) - 1)
    verdict!("2  conservation : donneur Φ_D[x τ_DA/τ_D0 + 1 − x]", "écart relatif < 1e-12", @sprintf("%.1e", eD), eD < 1e-12)
    verdict!("2  conservation : émission sensibilisée Φ_A·x·E", "écart relatif < 1e-12", @sprintf("%.1e", eA), eA < 1e-12)
end

# -----------------------------------------------------------------------------
# 3. Amplitudes à t = 0 égales aux fractions (canal bleu idéal)
# -----------------------------------------------------------------------------
let
    o = S.remplacer(m.o; η = SMatrix{2,2}([1.0 0.0; 0.0 1.0]))
    esB, _ = S.expsums(x20, m.p, o)
    f = esB.α[1] / (esB.α[1] + esB.α[2])
    verdict!("3  amplitudes du bleu idéal : x et 1 − x, rien sur τ_A", "écart < 1e-14",
        @sprintf("%.1e ; α(τ_A) = %.1e", abs(f - x20), esB.α[3]), abs(f - x20) < 1e-14 && esB.α[3] == 0)
end

# -----------------------------------------------------------------------------
# 4. Limite à Cl⁻ saturant (x_sat = 0)
# -----------------------------------------------------------------------------
let
    p0 = S.remplacer(m.p; x_sat = 0.0)
    esB, esY = S.expsums(S.fraction_active(Inf, p0), p0, m.o)
    mono = abs(esB.α[1]) < 1e-15 && abs(esB.α[3]) < 1e-15 && esB.α[2] > 0
    copie = maximum(abs.(esY.α .- (m.o.η[2, 1] / m.o.η[1, 1]) .* esB.α)) / esB.α[2]
    o0 = S.remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] m.o.η[1, 2]; 0.0 m.o.η[2, 2]]))
    _, esY0 = S.expsums(0.0, p0, o0)
    verdict!("4  Cl⁻ saturant : bleu mono-exponentiel à τ_D0", "seule l'amplitude τ_D0 non nulle", mono ? "oui" : "non", mono)
    verdict!("4  Cl⁻ saturant : jaune = (η_YD/η_BD)·bleu, nul si η_YD = 0", "écart < 1e-14",
        @sprintf("%.1e ; plus grande amplitude du jaune (η_YD=0) = %.1e", copie, maximum(abs.(esY0.α))), copie < 1e-14 && maximum(abs.(esY0.α)) == 0)
    iu = S.charger_instance_utilisateur(irf_mat)
    μB, _ = S.attendus(Inf, S.remplacer(m; p = p0), 1e6, :C, ws)
    τu = S.ajuster_mono_utilisateur(iu, μB)[1][1]
    verdict!("4  ton ajustement mono sur ce déclin rend τ_D0", "écart sur τ < 0,01 ns",
        @sprintf("%.4f contre %.4f ns (écart %.2f %% ≈ rapport des pas 0,048861/0,048828)", τu, m.p.τ_D0, 100 * (τu / m.p.τ_D0 - 1)),
        abs(τu - m.p.τ_D0) < 0.01)
end

# -----------------------------------------------------------------------------
# 5. Limite dégénérée : η_YD = 0, ε = 0 → jaune = émission sensibilisée pure
# -----------------------------------------------------------------------------
let
    o = S.remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] m.o.η[1, 2]; 0.0 m.o.η[2, 2]]), ε = 0.0, εΦ_dark = 0.0)
    md = S.remplacer(m; o = o)
    _, μY = S.attendus(20.0, md, 1e6, :A, ws)
    nobs = μY[W]
    B1, d1, B2, d2, poub = zeros(n), zeros(n), zeros(n), zeros(n), zeros(n)
    function modele!(μ, J, θ)
        S.bases_derivees!(B1, d1, poub, ws, θ[1], 0.0)
        S.bases_derivees!(B2, d2, poub, ws, θ[2], 0.0)
        for (i, k) in enumerate(W)
            μ[i] = θ[3] * B1[k] + θ[4] * B2[k]
            J[i, 1], J[i, 2], J[i, 3], J[i, 4] = θ[3] * d1[k], θ[4] * d2[k], B1[k], B2[k]
        end
    end
    S.base!(B2, ws, m.p.τ_A, 0.0)
    A = sum(nobs) / sum(B2[W])
    r = S.lm_poisson(modele!, [0.9 * m.p.τ_DA, 1.05 * m.p.τ_A, -0.7A, 1.2A], [0.01, 0.01, -Inf, -Inf], [20, 20, Inf, Inf], nobs;
        maxiter = 1000, tol_dec = 1e-24)
    τc, τl, Ac, Al = r.θ[1] < r.θ[2] ? (r.θ[1], r.θ[2], r.θ[3], r.θ[4]) : (r.θ[2], r.θ[1], r.θ[4], r.θ[3])
    e = max(abs(τc / m.p.τ_DA - 1), abs(τl / m.p.τ_A - 1), abs(Ac / Al + 1))
    verdict!("5  jaune pur sensibilisé : −A à τ_DA, +A à τ_A, rapport −1", "écart relatif < 1e-8",
        @sprintf("τ = (%.6f, %.6f) ns, A_court/A_long = %.10f ; écart %.1e", τc, τl, Ac / Al, e), e < 1e-8)
end

# -----------------------------------------------------------------------------
# 6. Efficacité asymptotique à 10⁶ photons
# -----------------------------------------------------------------------------
function efficacite_asymptotique(mm, est, Mrep, def; C = 20.0, N = 1e6)
    wss = [S.EspaceTravail(mm.inst) for _ in 1:Threads.nthreads()]
    μB, μY = S.attendus(C, mm, N, :C, wss[1])
    xs = zeros(Mrep)
    suivant = Threads.Atomic{Int}(1)
    @sync for w in wss
        Threads.@spawn begin
            ctx = S.Contexte(w, nothing)
            nB, nY = zeros(n), zeros(n)
            while (r = Threads.atomic_add!(suivant, 1)) <= Mrep
                S.tirer_deux!(nB, nY, μB, μY, StableRNG(S.graine(999, r)))
                xs[r] = S.observable(est, nB, nY, ctx)[1]
            end
        end
    end
    I = def === :conjoint ? S.information_fisher(C, mm, N, :C, wss[1]; params = (:x,), canaux = (:B, :Y), echelles = :separees)[1] :
        def === :absolue ? S.information_fisher(C, mm, N, :C, wss[1]; params = (:x,), canaux = (:B, :Y), echelles = :commune)[1] :
        def === :bleu ? S.information_fisher(C, mm, N, :C, wss[1]; params = (:x,), canaux = (:B,), echelles = :separees)[1] :
        S.information_fisher(C, mm, N, :C, wss[1]; params = (:x, :τ_D0, :τ_DA), canaux = (:B,), echelles = :separees)[1]
    eff = S.variance_crlb(I) / var(xs)
    ic = eff .* (1 .+ [-1.96, 1.96] .* sqrt(2 / (Mrep - 1)))   # IC à 95 % de var(x̂)/CRLB
    return eff, ic, mean(xs) - S.fraction_active(C, mm.p)
end
let
    Mrep = 20_000
    sansBA = S.remplacer(m; o = S.remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] 0.0; m.o.η[2, 1] m.o.η[2, 2]])))
    for (lab, mm, est, def) in (("5a conjoint", m, S.Conjoint(m, ws), :conjoint),
                                ("5d conjoint échelle commune", m, S.Conjoint(m, ws; commune = true), :absolue),
                                ("4 biexp τ figés (bleu sans fuite d'accepteur)", sansBA, S.BiExpFige(sansBA, ws), :bleu),
                                ("3 biexp τ libres (idem)", sansBA, S.BiExpLibre((m.p.τ_DA, m.p.τ_D0), 0.0, W, (0.02, 20.0)), :bleu_tau_libres))
        t = @elapsed eff, ic, biais = efficacite_asymptotique(mm, est, Mrep, def)
        ok = ic[1] <= 1.02 && ic[2] >= 0.98
        verdict!("6  efficacité à 10⁶ photons, $lab", "IC 95 % de CRLB/var recoupe [0,98 ; 1,02]",
            @sprintf("%.3f [%.3f ; %.3f], biais sur x %.1e (%d réal., %.0f s)", eff, ic..., biais, Mrep, t), ok)
    end
end

# -----------------------------------------------------------------------------
# 7. Dérivées analytiques contre différentiation automatique
# -----------------------------------------------------------------------------
let
    rng = Xoshiro(7)
    pire = 0.0
    for essai in 1:6
        p = S.remplacer(m.p; τ_D0 = 2.5 + 2rand(rng), τ_DA = 0.1 + 1.5rand(rng), τ_A = 2.5 + rand(rng))
        o = S.remplacer(m.o; ε = 0.1rand(rng), εΦ_dark = 0.05rand(rng))
        x = rand(rng)
        t₀ = 0.2 * (rand(rng) - 0.5)
        mm = S.remplacer(m; p = p, o = o)
        mB, mY, dB, dY = S.formes_et_derivees(x, mm, ws; t₀ = t₀)
        function gen(θ)
            α = S.amplitudes_M2_lent(θ[1], θ[2], θ[3], θ[4], p.Φ_D, p.Φ_A, o.η, o.ε, o.εΦ_dark)
            τ = [θ[3], θ[2], θ[4]]
            b = S.forward_direct(S.ExpSum(α[1, :], τ), m.inst.irf, θ[5], n, Δ, T)
            y = S.forward_direct(S.ExpSum(α[2, :], τ), m.inst.irf, θ[5], n, Δ, T)
            return vcat(b[W], y[W])
        end
        θ = [x, p.τ_D0, p.τ_DA, p.τ_A, t₀]
        Jad = ForwardDiff.jacobian(gen, θ)
        Jan = vcat(dB[W, :], dY[W, :])
        pire = max(pire, relerr(Jan, Jad), relerr(vcat(mB[W], mY[W]), gen(θ)))
    end
    verdict!("7  dérivées analytiques contre ForwardDiff (6 jeux aléatoires)", "écart relatif < 1e-10", @sprintf("%.1e", pire), pire < 1e-10)
end

# -----------------------------------------------------------------------------
# 8. Linéarité phasor
# -----------------------------------------------------------------------------
let
    Cs = vcat(0.0, exp.(range(log(0.5), log(500.0); length = 30)))
    dmax = Dict{Symbol,Float64}()
    for mec in (:M2_lent, :M2_rapide)
        mm = S.remplacer(m; p = S.remplacer(m.p; mecanisme = mec))
        for (c, nomc) in ((1, :B), (2, :Y))
            P = [S.phasor(S.formes(C, mm, ws)[c], m.inst) for C in Cs]
            dmax[Symbol(mec, :_, nomc)] = maximum(S.distance_droite(z, P[1], P[end]) for z in P)
        end
    end
    ok = dmax[:M2_lent_B] < 1e-10 && dmax[:M2_lent_Y] < 1e-10
    verdict!("8  phasor : la titration M2 lent est un segment (bleu et jaune)", "distance max < 1e-10 (en longueur de segment)",
        @sprintf("bleu %.1e, jaune %.1e ; M2 rapide : bleu %.1e, jaune %.1e", dmax[:M2_lent_B], dmax[:M2_lent_Y],
            dmax[:M2_rapide_B], dmax[:M2_rapide_Y]), ok)
end

# -----------------------------------------------------------------------------
# Fidélité : noyau, convolution, adaptateur
# -----------------------------------------------------------------------------
let
    function F_exact(u, τ)
        p(v) = (v2 = mod(v, T); exp(-v2 / τ) / (1 - exp(-T / τ)))
        g(a) = begin
            pts = [0.0, Δ]
            for k in -1:1
                b = u + a - k * T
                0 < b < Δ && push!(pts, b)
            end
            sort!(pts)
            sum(quadgk(b -> p(u + a - b), pts[i], pts[i+1]; rtol = 1e-14, atol = 1e-18)[1] for i in 1:length(pts)-1)
        end
        return quadgk(g, 0, Δ; rtol = 1e-13, atol = 1e-18)[1] / Δ
    end
    pire = maximum(abs(S.noyau(u, τ, Δ, T) / F_exact(u, τ) - 1) for τ in (0.15, 3.7), u in (-0.04, -0.01, 0.0, 0.02, Δ, 0.5, 12.0))
    verdict!("9  noyau contre l'intégrale double exacte (quadrature)", "écart relatif < 1e-11", @sprintf("%.1e", pire), pire < 1e-11)
    es = S.ExpSum([1.0, -0.4, 0.7], [m.p.τ_DA, m.p.τ_D0, m.p.τ_A])
    inst2 = S.remplacer(m.inst; t₀ = 0.0371)
    e = relerr(S.forward(es, inst2, ws), S.forward_direct(es, inst2))
    verdict!("10 convolution circulaire FFT contre convolution directe", "écart relatif < 1e-12", @sprintf("%.1e", e), e < 1e-12)
    ref = CSV.read(normpath(joinpath(@__DIR__, "..", "..", "tests", "csv_exports", "raw_data.csv")), DataFrame)
    ligne = ref[(ref.protein.=="SCLM").&(ref.concentration_mM.==0.01).&(ref.replicate.==1), :]
    h = S.lire_titration("SCLM")[1]
    iu = S.charger_instance_utilisateur(irf_mat)
    τu = S.ajuster_mono_utilisateur(iu, h.n)[1][1]
    verdict!("11 adaptateur : ton τ mono de SCLM 0 mM c1 (raw_data.csv)", "identique", @sprintf("%.16f contre %.16f", τu, ligne.tau_1exp[1]),
        τu == ligne.tau_1exp[1])
end

# -----------------------------------------------------------------------------
dossier = S.dossier_sorties("validation")
mkpath(dossier)
open(joinpath(dossier, "rapport_validation.md"), "w") do io
    println(io, "# Validation interne du simulateur\n")
    println(io, "Généré par `step2_validation.jl`. Modèle : ", S.decrire(m), "\n")
    println(io, "| Résultat | Test | Critère | Valeur |")
    println(io, "| --- | --- | --- | --- |")
    for (t, c, v, ok) in resultats
        cellule(x) = replace(x, "|" => "\\|")
        println(io, "| ", ok ? "OK" : "**ÉCHEC**", " | ", cellule(t), " | ", cellule(c), " | ", cellule(v), " |")
    end
end
nok = count(r -> r[4], resultats)
println("\n$nok / $(length(resultats)) tests réussis — rapport : $(joinpath(dossier, "rapport_validation.md"))")
nok == length(resultats) || exit(1)
