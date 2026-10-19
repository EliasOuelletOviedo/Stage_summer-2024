# kalman.jl — filtre de Kalman étendu, étape 2 : la prédiction connaît la commande.
#
# Trois modèles d'état :
#   :fusion   x = (C, g, ġ)           J/B et 5a ; C suit la physiologie, g et ġ l'OU de l'étape 1
#   :fusion4  x = (C, g, ġ, k_pompe)  idem, k_pompe estimé (course 7)
#   :ratio    x = (C, g)              J/B seul ; g en marche aléatoire, sans ancre (course 4)
#
# Prédiction : C' = C + Δt·[J(u) − k·(C − C_repos)] ; g' = g + Δt·ġ ; ġ' = (1 − Δt/τ_OU)·ġ ; k' = k.
# Q dépend de la commande : σ_C²(u) = σ_C,0² + β·(J(u)·Δt)², β l'erreur relative sur le flux au carré.
# Mesure (inchangée) : z_R = g·h_R(C), z_S = h_S(C). Le filtre ne connaît rien de la photophysique.

using LinearAlgebra

mutable struct EKF
    modele::Symbol
    x::Vector{Float64}
    P::Matrix{Float64}
    ph::Physio
    σC0::Float64
    β::Float64
    σw::Vector{Float64}         # écarts-types par pas des états autres que C
    ρ::Float64                  # 1 − Δt/τ_OU
    h::Vector{Function}
    dh::Vector{Function}
    mesures::Vector{Bool}       # lectures que ce filtre utilise (J/B, 5a)
end

"""Prédiction avec la commande u appliquée pendant le pas. Rend le filtre."""
function predire!(f::EKF, u)
    ph, n = f.ph, length(f.x)
    C = f.x[1]
    k = f.modele === :fusion4 ? f.x[4] : ph.k
    J = flux(u, ph)
    F = Matrix{Float64}(I, n, n)
    F[1, 1] = 1 - ph.dt * k
    xn = copy(f.x)
    xn[1] = C + ph.dt * (J - k * (C - ph.C_repos))
    if n >= 3
        xn[2] = f.x[2] + ph.dt * f.x[3]
        xn[3] = f.ρ * f.x[3]
        F[2, 3] = ph.dt
        F[3, 3] = f.ρ
    end
    f.modele === :fusion4 && (F[1, 4] = -ph.dt * (C - ph.C_repos))
    Q = Diagonal([f.σC0^2 + f.β * (J * ph.dt)^2; f.σw .^ 2])
    f.x = xn
    f.P = F * f.P * F' + Q
    return f
end

"""
Mise à jour avec les lectures valables et utilisées par ce filtre ; forme de Joseph. Rend les
innovations normalisées νᵢ/√Sᵢᵢ (NaN pour une lecture absente ou non utilisée).
"""
function mettre_a_jour!(f::EKF, z::AbstractVector, R::AbstractMatrix, ok::AbstractVector{Bool})
    out = fill(NaN, length(z))
    idx = findall(ok .& f.mesures)
    isempty(idx) && return out
    C, g = f.x[1], f.x[2]
    zr = zeros(1, length(f.x) - 2)
    gain = (true, false)            # J/B est multiplié par g, 5a non
    hx = [gain[i] ? g * f.h[i](C) : f.h[i](C) for i in idx]
    H = reduce(vcat, [hcat(gain[i] ? [g * f.dh[i](C) f.h[i](C)] : [f.dh[i](C) 0.0], zr) for i in idx])
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

"""Diagnostics d'un pas : Ĉ, ĝ, P₁₁, P₂₂, ρ₁₂, et k̂_pompe, P₄₄ pour le filtre à quatre états."""
diagnostics(f::EKF) = (C = f.x[1], g = f.x[2], P11 = f.P[1, 1], P22 = f.P[2, 2], rho12 = f.P[1, 2] / sqrt(f.P[1, 1] * f.P[2, 2]),
                       k = length(f.x) >= 4 ? f.x[4] : NaN, P44 = length(f.x) >= 4 ? f.P[4, 4] : NaN)

"""
Matrice d'observabilité O = [H; HF; HF²; HF³] du filtre à quatre états, linéarisée au point
(C, g, ġ = 0, k) : son rang vaut 4 si et seulement si f(C)·p′(C)·(C − C_repos) ≠ 0.
"""
function observabilite(f::EKF, C, g)
    ph = f.ph
    k = f.x[4]
    F = Matrix{Float64}(I, 4, 4)
    F[1, 1], F[1, 4], F[2, 3], F[3, 3] = 1 - ph.dt * k, -ph.dt * (C - ph.C_repos), ph.dt, f.ρ
    H = [g*f.dh[1](C) f.h[1](C) 0 0; f.dh[2](C) 0 0 0]
    return reduce(vcat, [H * F^i for i in 0:3])
end
