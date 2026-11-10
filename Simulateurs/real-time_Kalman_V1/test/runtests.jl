# Les quatre tests à faire passer avant le reste, plus deux tests du filtre seul sur des mesures
# fabriquées à la main (sans histogramme) et le test d'amplitudes avec et sans diffusion. Lancer : julia --project=. test/runtests.jl

using Test, Statistics, StableRNGs
include(joinpath(@__DIR__, "..", "run.jl"))

const p = charger_params()
const lec = lecteurs(p)

# un ensemble englobant : un groupe qui échoue n'empêche pas les suivants de tourner
@testset "simulateur Kalman temps réel" begin

@testset "1. aller-retour sans bruit : J/B et 5a retrouvent C" begin
    for C in (20.0, 23.7)        # 20 mM est un nœud de la grille de calibration, 23,7 mM non
        μ = histogrammes_c(C, 1.0, p)
        @test inverser(lec.JB.courbe, estimer(μ..., lec.JB.cfg).z) ≈ C atol = 0.01
        @test inverser(lec.S.courbe, estimer(μ..., lec.S.cfg).z) ≈ C atol = 0.01
    end
    μ = histogrammes_c(20.0, 1.0, p)
    W = p.fenetre
    @test sum(μ[1][W]) - sum(p.fond_b[W]) ≈ p.N rtol = 1e-9          # le bleu porte N photons
end

@testset "2. invariance de 5a : jaune × 1,2" begin
    # (a) l'hypothèse structurelle : sur un histogramme que le modèle de 5a décrit exactement
    # (les colonnes de 5a + fond), A_J absorbe g au bruit numérique près
    c = Cfg5a(p, (0.0, 1.0))
    W = p.fenetre
    poids(X, x) = size(X, 2) == 3 ? [5.0 / 300; x] : x       # 5 photons de diffusion si 5a la modélise
    μb = zeros(p.n); μb[W] .= 300 .* c.Xb * poids(c.Xb, [0.6, 0.4]); μb .+= p.fond_b
    jaune(g) = (μ = zeros(p.n); μ[W] .= g .* 300 .* c.Xj * poids(c.Xj, [0.12, 0.88]); μ .+ p.fond_j)
    @test abs(estimer(μb, jaune(1.0), c).z - estimer(μb, jaune(1.2), c).z) < 1e-8
    # (b) le générateur complet : avec la diffusion dans 5a, l'invariance y est aussi exacte ;
    # sans elle, le reste doit rester négligeable devant le bruit d'un pas (moins de 1/1000 de σ)
    μ1, μ2 = histogrammes_c(20.0, 1.0, p), histogrammes_c(20.0, 1.2, p)
    for w in ((1.0, 0.0), (0.0, 1.0), lec.w)                          # p, q et l'observable transmise
        c = Cfg5a(p, w)
        l1 = estimer(μ1..., c)
        @test abs(l1.z - estimer(μ2..., c).z) < 1e-3 * sqrt(l1.σ²)
    end
end

@testset "3. sensibilité de J/B : jaune × 1,2" begin
    z1, z2 = (estimer(histogrammes_c(20.0, g, p)..., lec.JB.cfg).z for g in (1.0, 1.2))
    @test z2 / z1 ≈ 1.2 rtol = 1e-9
    @test z1 ≈ jb(20.0, p) rtol = 1e-9
end

@testset "4. filtre à g constant : le Kalman converge vers J/B" begin
    # une seule course varie de ±5 % (ĝ erre autour de 1) : la convergence se juge sur 8 courses
    res = Vector{Any}(undef, 8)
    Threads.@threads for r in 1:8
        res[r] = course(avec(p; alpha = 0.0), lec, r)[1]
    end
    for m in (reference(p), ma_reglee(p, p.N), nom_ma(1e-4))   # la référence OU, la marche réglée, le premier défaut
        ecart(sel) = mean(begin
            T = x[m]; t, K, R = T[:t], T[:C_kalman], T[:C_ratio]
            v = .!isnan.(R) .& sel(t)
            (mean(K[v]) - mean(R[v])) / mean(R[v]) end for x in res)
        @test abs(ecart(t -> t .>= t[end] - 30)) < 0.05
        @test abs(ecart(t -> trues(length(t)))) < 0.05
        @test abs(mean(mean(x[m][:g_chapeau][x[m][:t] .>= 60]) for x in res) - 1) < 0.02
    end
end

@testset "5. filtre seul, mesures fabriquées : g s'apprend, et J/B seul ne le peut pas" begin
    hR, dhR = c -> 6 / (1 + c / 30), c -> -0.2 / (1 + c / 30)^2
    hS, dhS = c -> c / (c + 40), c -> 40 / (c + 40)^2
    σR, σS = 0.03, 0.02
    fab(ok; modele = :marche) = EKF(modele = modele, x0 = modele === :marche ? [20.0, 1.0] : [20.0, 1.0, 0.0],
        σ0 = modele === :marche ? [2.0, 0.05] : [2.0, 0.05, 2e-3], σw = modele === :marche ? [0.1, 1e-4] : [0.1, 1e-4, 1e-5],
        dt = 0.1, C_repos = 20.0, τ = 1e9, h = Function[hR, hS], dh = Function[dhR, dhS], gain = [true, false])
    # g constant à 1,05, inconnu (le filtre part de 1) : un filtre juste l'apprend et reste
    # cohérent avec son P. Suivre une rampe est une question de réglage de Q, pas de code :
    # c'est la course qui le mesure.
    rng = StableRNG(1)
    f2, f1 = fab(true), fab(true)
    R = Diagonal([σR^2, σS^2])
    g = 1.05
    for k in 1:900
        z = [g * hR(20.0) + σR * randn(rng), hS(20.0) + σS * randn(rng)]
        mettre_a_jour!(predire!(f2), z, R, [true, true])
        mettre_a_jour!(predire!(f1), z, R, [true, false])
    end
    @test abs(f2.x[2] - g) < 3 * sqrt(f2.P[2, 2]) + 1e-3
    @test abs(f2.x[1] - 20.0) < 3 * sqrt(f2.P[1, 1]) + 0.1
    # sans 5a, H n'a qu'une ligne : l'écart de gain passe dans Ĉ ; g étant dans l'état, le
    # filtre l'avoue par un P₁₁ qui ne se resserre pas
    @test abs(f1.x[1] - 20.0) > 3 * abs(f2.x[1] - 20.0)
    @test f1.P[1, 1] > 10 * f2.P[1, 1]
end

@testset "6. modèle (g, ġ), mesures fabriquées : une rampe se suit sans traînage" begin
    hR, dhR = c -> 6 / (1 + c / 30), c -> -0.2 / (1 + c / 30)^2
    hS, dhS = c -> c / (c + 40), c -> 40 / (c + 40)^2
    mk(m) = EKF(modele = m, x0 = m === :marche ? [20.0, 1.0] : [20.0, 1.0, 0.0], σ0 = m === :marche ? [2.0, 0.05] : [2.0, 0.05, 2e-3],
                σw = m === :marche ? [0.1, 1e-4] : [0.1, 1e-4, 1e-5], dt = 0.1, C_repos = 20.0, τ = 1e9,
                h = Function[hR, hS], dh = Function[dhR, dhS], gain = [true, false])
    rng = StableRNG(2)
    fv, fm = mk(:vitesse), mk(:marche)
    R = Diagonal([0.03^2, 0.02^2])
    g = 1.0
    for k in 1:900
        g = 1 + 0.10 * k / 900                        # rampe : ġ = 1,11e-3 s⁻¹
        z = [g * hR(20.0) + 0.03 * randn(rng), hS(20.0) + 0.02 * randn(rng)]
        mettre_a_jour!(predire!(fv), z, R, [true, true])
        mettre_a_jour!(predire!(fm), z, R, [true, true])
    end
    @test abs(fv.x[2] - g) < 3 * sqrt(fv.P[2, 2]) + 1e-3      # le modèle à vitesse suit, et le sait
    @test abs(fv.x[3] - 0.10 / 90) < 3 * sqrt(fv.P[3, 3])     # et retrouve la pente
    @test abs(fm.x[2] - g) > 5 * abs(fv.x[2] - g)             # la marche aléatoire raide traîne
end

@testset "7. amplitudes : 5a retrouve exactement p et q injectés, avec et sans diffusion" begin
    # avec la diffusion dans 5a (correctif A), l'écart doit tomber au niveau machine dans les
    # deux cas ; sans elle, seul le générateur sans diffusion est décrit exactement
    for q0 in (p, sans_diffusion(p)), C in (0.0, 20.0, 80.0, 131.0)
        μ = histogrammes_c(C, 1.0, q0)
        @test estimer(μ..., Cfg5a(p, (1.0, 0.0))).z ≈ p.lang["p_bleu"](C) atol = 1e-9
        @test estimer(μ..., Cfg5a(p, (0.0, 1.0))).z ≈ p.lang["q_jaune"](C) atol = 1e-9
    end
end

end
