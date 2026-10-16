# kalman.jl — filtre de Kalman étendu sur l'état (C, g) ou (C, g, ġ).
#
# Le filtre ne connaît rien de la photophysique : il reçoit des observables, leur
# covariance R et, pour chaque mesure, une fonction h(C) et sa dérivée, multipliée ou non
# par g. Ici, le rapport vaut g·h_R(C) et 5a vaut h_S(C) : H a un zéro dans la colonne de
# g pour 5a, c'est ce qui rend C et g identifiables.
#
# Propagation (la trajectoire est inconnue) : C relaxe vers son repos avec τ_KCC2. Pour g,
# deux modèles :
#   :marche   g_{k+1} = g_k + w_g                       (marche aléatoire)
#   :vitesse  g_{k+1} = g_k + Δt·ġ_k + w_g,  ġ_{k+1} = ρ·ġ_k + w_ġ
# avec ρ = 1 (vitesse constante : suit une rampe sans traînage) ou ρ = 1 − Δt/τ_ġ
# (Ornstein-Uhlenbeck : oublie une pente apprise après un saut). Q est diagonale.

using LinearAlgebra

mutable struct EKF
    x::Vector{Float64}
    P::Matrix{Float64}
    Q::Matrix{Float64}
    F::Matrix{Float64}          # x_{k+1} = F·x_k + b
    b::Vector{Float64}
    h::Vector{Function}
    dh::Vector{Function}
    gain::Vector{Bool}          # la mesure i vaut g·hᵢ(C) si gain[i], hᵢ(C) sinon
end

function EKF(; modele = :marche, x0::Vector{Float64}, σ0::Vector{Float64}, σw::Vector{Float64}, dt, C_repos, τ,
             τ_gp = Inf, h, dh, gain)
    a = 1 - dt / τ
    F = modele === :marche ? [a 0.0; 0.0 1.0] :
        modele === :vitesse ? [a 0.0 0.0; 0.0 1.0 dt; 0.0 0.0 (isinf(τ_gp) ? 1.0 : 1 - dt / τ_gp)] :
        error("modèle de g inconnu : $modele")
    n = size(F, 1)
    length(x0) == length(σ0) == length(σw) == n || error("état de dimension $n attendu")
    b = zeros(n)
    b[1] = (1 - a) * C_repos
    return EKF(copy(x0), Matrix(Diagonal(σ0 .^ 2)), Matrix(Diagonal(σw .^ 2)), F, b, h, dh, gain)
end

"""Prédiction : x ← F·x + b, P ← F·P·Fᵀ + Q."""
function predire!(f::EKF)
    f.x = f.F * f.x .+ f.b
    f.P = f.F * f.P * f.F' + f.Q
    return f
end

"""
Mise à jour avec les mesures valables (`ok`) ; R peut être complet ou diagonal. Forme de
Joseph pour garder P symétrique définie positive. Rend les innovations normalisées
νᵢ/√Sᵢᵢ (NaN pour une mesure absente) : elles doivent être blanches et de variance 1.
"""
function mettre_a_jour!(f::EKF, z::AbstractVector, R::AbstractMatrix, ok::AbstractVector{Bool})
    out = fill(NaN, length(z))
    idx = findall(ok)
    isempty(idx) && return out
    C, g = f.x[1], f.x[2]
    zr = zeros(1, length(f.x) - 2)
    hx = [f.gain[i] ? g * f.h[i](C) : f.h[i](C) for i in idx]
    H = reduce(vcat, [hcat(f.gain[i] ? [g * f.dh[i](C) f.h[i](C)] : [f.dh[i](C) 0.0], zr) for i in idx])
    Ri = R[idx, idx]
    ν = z[idx] .- hx
    S = Symmetric(H * f.P * H' + Ri)
    K = f.P * H' / S
    f.x .+= K * ν
    A = I - K * H
    f.P = A * f.P * A' + K * Ri * K'
    out[idx] .= ν ./ sqrt.(diag(S))
    return out
end

"""Les diagnostics d'un pas : Ĉ, ĝ, la pente ĝ̇ (NaN sans état de vitesse), P₁₁, P₂₂, ρ₁₂ = P₁₂/√(P₁₁P₂₂)."""
diagnostics(f::EKF) = (C = f.x[1], g = f.x[2], gp = length(f.x) > 2 ? f.x[3] : NaN, P11 = f.P[1, 1], P22 = f.P[2, 2],
                       rho12 = f.P[1, 2] / sqrt(f.P[1, 1] * f.P[2, 2]))
