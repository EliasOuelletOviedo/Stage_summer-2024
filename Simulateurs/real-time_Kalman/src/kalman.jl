# kalman.jl — filtre de Kalman étendu sur l'état x = (C, g).
#
# Le filtre ne connaît rien de la photophysique : il reçoit des observables, leur
# covariance R et, pour chaque mesure, une fonction h(C) et sa dérivée, multipliée ou non
# par g. Ici, le rapport vaut g·h_R(C) et 5a vaut h_S(C) : H a un zéro dans la colonne de
# g pour 5a, c'est ce qui rend C et g identifiables.
#
# Propagation minimale (la trajectoire est inconnue) : C relaxe vers son repos avec τ_KCC2,
# g suit une marche aléatoire. Q = diag(σ²_w,C, σ²_w,g) est le seul vrai réglage.

using LinearAlgebra

mutable struct EKF
    x::Vector{Float64}
    P::Matrix{Float64}
    Q::Matrix{Float64}
    dt::Float64
    C_repos::Float64
    τ::Float64
    h::Vector{Function}
    dh::Vector{Function}
    gain::Vector{Bool}          # la mesure i vaut g·hᵢ(C) si gain[i], hᵢ(C) sinon
end

EKF(; x0, P0, σw, dt, C_repos, τ, h, dh, gain) =
    EKF(float.(copy(x0)), Matrix{Float64}(P0), Matrix(Diagonal(float.(σw) .^ 2)), dt, C_repos, τ, h, dh, gain)

"""C_{k+1} = C_k + Δt·(C_repos − C_k)/τ + w_C ; g_{k+1} = g_k + w_g."""
function predire!(f::EKF)
    a = 1 - f.dt / f.τ
    f.x[1] = f.C_repos + a * (f.x[1] - f.C_repos)
    F = [a 0.0; 0.0 1.0]
    f.P = F * f.P * F' + f.Q
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
    C, g = f.x
    hx = [f.gain[i] ? g * f.h[i](C) : f.h[i](C) for i in idx]
    H = reduce(vcat, [(f.gain[i] ? [g * f.dh[i](C) f.h[i](C)] : [f.dh[i](C) 0.0]) for i in idx])
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

"""Les diagnostics d'un pas : Ĉ, ĝ, P₁₁, P₂₂ et la corrélation ρ₁₂ = P₁₂/√(P₁₁P₂₂)."""
diagnostics(f::EKF) = (C = f.x[1], g = f.x[2], P11 = f.P[1, 1], P22 = f.P[2, 2],
                       rho12 = f.P[1, 2] / sqrt(f.P[1, 1] * f.P[2, 2]))
