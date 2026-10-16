# Les quatre tests à faire passer avant le reste, plus un test du filtre seul sur des mesures
# fabriquées à la main (sans histogramme). Lancer : julia --project=. test/runtests.jl

using Test, Statistics, StableRNGs
include(joinpath(@__DIR__, "..", "run.jl"))

const p = charger_params()
const lec = lecteurs(p)

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
    # (deux exponentielles + fond, sans diffusion), A_J absorbe g au bruit numérique près
    c = Cfg5a(p, (0.0, 1.0))
    W = p.fenetre
    μb = zeros(p.n); μb[W] .= 300 .* c.Xb * [0.6, 0.4]; μb .+= p.fond_b
    jaune(g) = (μ = zeros(p.n); μ[W] .= g .* 300 .* c.Xj * [0.12, 0.88]; μ .+ p.fond_j)
    @test abs(estimer(μb, jaune(1.0), c).z - estimer(μb, jaune(1.2), c).z) < 1e-8
    # (b) le générateur complet ajoute la diffusion, que 5a ne modélise pas : le reste est
    # mesuré et doit rester négligeable devant le bruit d'un pas (moins de 1/1000 de σ)
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
    T, _ = course(avec(p; alpha = 0.0), lec, 1)
    for m in MODES
        t, K, R = T[m][:t], T[m][:C_kalman], T[m][:C_ratio]
        v = .!isnan.(R)
        fin = v .& (t .>= t[end] - 30)
        @test abs(mean(K[fin]) - mean(R[fin])) / mean(R[fin]) < 0.05
        @test abs(mean(K[v]) - mean(R[v])) / mean(R[v]) < 0.05
        @test abs(mean(T[m][:g_chapeau][fin]) - 1) < 0.02
    end
end

@testset "5. filtre seul, mesures fabriquées : g s'apprend, et J/B seul ne le peut pas" begin
    hR, dhR = c -> 6 / (1 + c / 30), c -> -0.2 / (1 + c / 30)^2
    hS, dhS = c -> c / (c + 40), c -> 40 / (c + 40)^2
    σR, σS = 0.03, 0.02
    fab(ok) = EKF(x0 = [20.0, 1.0], P0 = Diagonal([4.0, 0.05^2]), σw = [0.1, 1e-4], dt = 0.1, C_repos = 20.0, τ = 1e9,
                  h = Function[hR, hS], dh = Function[dhR, dhS], gain = [true, false])
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
