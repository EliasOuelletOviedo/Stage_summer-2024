# Tests de l'étape 2. Lancer : julia --project=. -t auto test/runtests.jl

using Test, Statistics, StableRNGs, LinearAlgebra
include(joinpath(@__DIR__, "..", "run.jl"))

const p = charger_params()
const lec = lecteurs(p)
const ph = physio(p)

@testset "étape 2" begin

@testset "1. reprise de l'étape 1 : J/B et 5a retrouvent C sans bruit" begin
    for C in (5.0, 20.0, 23.7)
        μ = histogrammes_c(C, 1.0, p)
        @test inverser(lec.JB.courbe, estimer(μ..., lec.JB.cfg).z) ≈ C atol = 0.01
        @test inverser(lec.S.courbe, estimer(μ..., lec.S.cfg).z) ≈ C atol = 0.01
    end
end

@testset "2. physiologie : flux saturant, équilibre, unidirectionnalité" begin
    @test flux(0.0, ph) == 0.0
    @test flux(1.0, ph) < ph.J_max
    @test flux(ph.u_demi, ph) ≈ ph.J_max / 2
    @test dflux(0.0, ph) > dflux(1.0, ph)                       # le gain chute à forte commande
    for Cc in (10.0, 25.0)
        u = u_maintien(Cc, ph)
        @test flux(u, ph) ≈ ph.k * (Cc - ph.C_repos) rtol = 1e-12
        sans_bruit = Physio(ph.C_repos, ph.k, ph.J_max, ph.u_demi, 0.0, ph.dt)
        C = Cc
        for _ in 1:2000
            C = pas_cellule(C, u, sans_bruit, StableRNG(1))
        end
        @test C ≈ Cc rtol = 1e-9
    end
    @test u_maintien(3.0, ph) == 0.0                           # sous le repos : la commande ne peut rien
end

@testset "3. PI : anti-emballement et saturation" begin
    r = PI(0.02, 0.01)
    for _ in 1:100
        @test commande!(r, -5.0, 1.0, 0.1) == 0.0              # consigne sous le chlore : u saturé à 0
    end
    @test r.I == 0.0                                           # l'intégrateur n'a pas accumulé de dette
    u = commande!(r, 5.0, 1.0, 0.1)
    @test 0 < u < 1 && r.I > 0
    rB = PI(0.02, 0.01; variante = :B, c = 1.0)
    @test commande!(rB, 0.5, 1.0, 0.1) == 0.0                  # zone morte : |e| < √P₁₁
    rC = PI(0.02, 0.0; variante = :C, c = 1.0, σ0² = 1.0)
    @test commande!(rC, 5.0, 1.0, 0.1) ≈ 0.02 * 5.0 / 2        # gain divisé par 2 quand P₁₁ = σ₀²
end

@testset "4. observabilité de k_pompe : rang 4 hors équilibre, 3 à l'équilibre" begin
    f = filtre(p, lec, ph, spec("f4", :fusion4, 0.0))
    @test rank(observabilite(f, 25.0, 1.0); rtol = 1e-10) == 4
    @test rank(observabilite(f, ph.C_repos, 1.0); rtol = 1e-10) == 3
end

@testset "5. filtre, mesures fabriquées : la commande connue, g appris ; J/B seul à g figé biaisé" begin
    hR, dhR = c -> 6 / (1 + c / 30), c -> -0.2 / (1 + c / 30)^2
    hS, dhS = c -> c / (c + 40), c -> 40 / (c + 40)^2
    σR, σS = 0.03, 0.02
    mk(m) = EKF(m, m === :ratio ? [5.0, 1.0] : [5.0, 1.0, 0.0], Matrix(Diagonal(m === :ratio ? [4.0, 0.0] : [4.0, 0.05^2, 1e-6])),
                ph, 0.3, 0.0, m === :ratio ? [0.0] : [1e-4, 1e-5], 1 - 0.1 / 30, Function[hR, hS], Function[dhR, dhS],
                m === :ratio ? [true, false] : [true, true])
    fF, fR = mk(:fusion), mk(:ratio)
    rng = StableRNG(3)
    sans_bruit = Physio(ph.C_repos, ph.k, ph.J_max, ph.u_demi, 0.0, ph.dt)
    C, g, u = ph.C_repos, 1.10, 0.0
    for k in 1:900
        t = k * ph.dt
        u = 10 <= t < 50 ? u_maintien(25.0, ph) + (t < 13 ? 0.5 : 0.0) : 0.0
        z = [g * hR(C) + σR * randn(rng), hS(C) + σS * randn(rng)]
        R = Diagonal([σR^2, σS^2])
        mettre_a_jour!(fF, z, R, [true, true])
        mettre_a_jour!(fR, z, R, [true, true])
        predire!(fF, u); predire!(fR, u)
        C = pas_cellule(C, u, sans_bruit, rng)
    end
    @test abs(fF.x[2] - g) < 3 * sqrt(fF.P[2, 2]) + 1e-3        # la fusion apprend g = 1,10
    @test abs(fF.x[1] - C) < 3 * sqrt(fF.P[1, 1]) + 0.2
    @test abs(fR.x[1] - C) > 3 * sqrt(fR.P[1, 1])               # J/B seul à g figé : biaisé, et confiant
end

end
