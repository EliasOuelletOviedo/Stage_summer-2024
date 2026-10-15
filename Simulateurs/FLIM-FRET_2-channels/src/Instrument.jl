# Instrument.jl — périodisation, intégration par canal, convolution circulaire.
#
# Le noyau F(u) est le nombre de photons attendus dans un canal de détection
# dont le début est décalé de u par rapport au début d'un canal de l'IRF, pour
# un déclin e^{-t/τ} :
#   - déclin périodisé exactement (série géométrique des impulsions
#     précédentes, aucune troncature) ;
#   - IRF supposée uniforme à l'intérieur de chacun de ses canaux (c'est un
#     histogramme : on ne sait rien de plus fin) ;
#   - intégré sur le canal de détection, pas échantillonné au centre.
# Ces trois choix donnent F(u) = (1/Δ)∫(Δ−|w|) p(u+w) dw, le déclin périodisé
# p lissé par un triangle de demi-largeur Δ : F est C¹ partout (y compris au
# saut de p à l'impulsion), donc dérivable en τ et en t₀, et Σₘ F(mΔ − t₀) = τ
# exactement pour tout t₀ (conservation des photons).
#
# Avec Q = 1/(e^{T/τ} − 1), c = τ²·4sinh²(Δ/2τ)/Δ et E = e^{-u/τ} :
#   F(u) = Q·c·E + G(u)        (Q·c·E : impulsions précédentes)
#   G(u) = c·E                                               si u ≥ Δ
#   G(u) = [E τ² φ₂(Δ/τ) + Δτ(1−E) − τ² φ₂(u/τ)]/Δ           si 0 ≤ u < Δ
#   G(u) = τ² φ₂((Δ+u)/τ)/Δ                                  si −Δ ≤ u < 0
# où φ₂(z) = e^{-z} − 1 + z. u est ramené dans [−Δ, T−Δ) par périodicité.

"""φ₂(z) = e^{-z} − 1 + z, calculé sans annulation catastrophique pour z petit."""
@inline function phi2(z)
    if abs(z) < 0.02
        z2 = z * z
        return z2 * (0.5 - z / 6 + z2 / 24 - z2 * z / 120 + z2 * z2 / 720 - z2 * z2 * z / 5040)
    end
    return expm1(-z) + z
end

"""Ramène u dans [−Δ, T−Δ) (floor plutôt que mod : dérivable par l'AD)."""
@inline reduire(u, Δ, T) = u - T * floor((u + Δ) / T)

"""
    noyau(u, τ, Δ, T)

F(u) pour une exponentielle de durée de vie τ (générique, pour l'AD).
"""
function noyau(u, τ, Δ, T)
    u = reduire(u, Δ, T)
    x = Δ / τ
    c = τ^2 * 4 * sinh(x / 2)^2 / Δ
    Q = 1 / expm1(T / τ)
    E = exp(-u / τ)
    if u >= Δ
        G = c * E
    elseif u >= 0
        G = (E * τ^2 * phi2(x) - Δ * τ * expm1(-u / τ) - τ^2 * phi2(u / τ)) / Δ
    else
        G = τ^2 * phi2((Δ + u) / τ) / Δ
    end
    return Q * c * E + G
end

"""
    noyau_derivees(u, τ, Δ, T) -> (F, ∂F/∂u, ∂F/∂τ)

Dérivées analytiques du noyau (le test 7 les compare à l'AD).
"""
function noyau_derivees(u::Float64, τ::Float64, Δ::Float64, T::Float64)
    u = reduire(u, Δ, T)
    x = Δ / τ
    s = 4 * sinh(x / 2)^2
    c = τ^2 * s / Δ
    dc = (τ / Δ) * (2s - 2x * sinh(x))
    Q = 1 / expm1(T / τ)
    dQ = Q * (1 + Q) * T / τ^2
    E = exp(-u / τ)
    dE_du = -E / τ
    dE_dτ = E * u / τ^2
    if u >= Δ
        G = c * E
        dG_du = c * dE_du
        dG_dτ = dc * E + c * dE_dτ
    elseif u >= 0
        a = u / τ
        p2x, p2a = phi2(x), phi2(a)
        unmE = -expm1(-a)                 # 1 − E
        G = (E * τ^2 * p2x + Δ * τ * unmE - τ^2 * p2a) / Δ
        dG_du = (-E * τ * p2x + Δ * E - τ * unmE) / Δ
        dG_dτ = (E * τ * (a * p2x + 2p2x + x * expm1(-x)) + Δ * unmE - Δ * E * a - 2τ * p2a + τ * a * unmE) / Δ
    else
        l = (Δ + u) / τ
        p2l = phi2(l)
        G = τ^2 * p2l / Δ
        dG_du = -τ * expm1(-l) / Δ
        dG_dτ = (2τ * p2l + τ * l * expm1(-l)) / Δ
    end
    F = Q * c * E + G
    dF_du = Q * c * dE_du + dG_du
    dF_dτ = (dQ * c + Q * dc) * E + Q * c * dE_dτ + dG_dτ
    return F, dF_du, dF_dτ
end

"""f[m+1] = F(mΔ − t₀), m = 0…n−1."""
function noyau!(f::AbstractVector, τ, t₀, n::Int, Δ, T)
    @inbounds for m in 0:n-1
        f[m+1] = noyau(m * Δ - t₀, τ, Δ, T)
    end
    return f
end

function noyau_derivees!(f, dfτ, dfu, τ::Float64, t₀::Float64, n::Int, Δ::Float64, T::Float64)
    @inbounds for m in 0:n-1
        f[m+1], dfu[m+1], dfτ[m+1] = noyau_derivees(m * Δ - t₀, τ, Δ, T)
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Convolution circulaire par FFT, planifiée une fois par espace de travail
# -----------------------------------------------------------------------------

"""
    EspaceTravail(inst)

Plans FFT et tampons d'un fil d'exécution. Les plans sont créés une fois (la
planification FFTW n'est pas sûre en parallèle : créer les espaces de travail
avant les boucles parallèles), puis réutilisés.
"""
struct EspaceTravail{P1,P2}
    n::Int
    Δ::Float64
    T::Float64
    irf::Vector{Float64}
    plan_r::P1
    plan_ir::P2
    H::Vector{ComplexF64}      # rfft de l'IRF
    f::Vector{Float64}
    dfτ::Vector{Float64}
    dfu::Vector{Float64}
    Fk::Vector{ComplexF64}
    tmp::Vector{ComplexF64}
end

function EspaceTravail(inst::InstrumentParams)
    n = inst.n_bins
    length(inst.irf) == n || error("IRF de $(length(inst.irf)) canaux pour une grille de $n")
    plan_r = plan_rfft(zeros(n))
    plan_ir = plan_irfft(zeros(ComplexF64, n ÷ 2 + 1), n)
    H = plan_r * inst.irf
    return EspaceTravail(n, pas_temps(inst), inst.T_rep, copy(inst.irf), plan_r, plan_ir, H,
        zeros(n), zeros(n), zeros(n), zeros(ComplexF64, n ÷ 2 + 1), zeros(ComplexF64, n ÷ 2 + 1))
end

"""out = IRF ⊛ g (circulaire, sur la période)."""
function convoluer!(out::AbstractVector{Float64}, ws::EspaceTravail, g::AbstractVector{Float64})
    mul!(ws.Fk, ws.plan_r, g)
    ws.tmp .= ws.H .* ws.Fk
    mul!(out, ws.plan_ir, ws.tmp)   # c2r détruit tmp : tampon jetable
    return out
end

"""out = IRF ⊛ F(·; τ, t₀) : la base convoluée d'une exponentielle d'amplitude 1."""
function base!(out, ws::EspaceTravail, τ::Float64, t₀::Float64)
    noyau!(ws.f, τ, t₀, ws.n, ws.Δ, ws.T)
    return convoluer!(out, ws, ws.f)
end

"""Base convoluée et ses dérivées ∂/∂τ et ∂/∂t₀ (= −IRF ⊛ ∂F/∂u)."""
function bases_derivees!(B, dBτ, dBt, ws::EspaceTravail, τ::Float64, t₀::Float64)
    noyau_derivees!(ws.f, ws.dfτ, ws.dfu, τ, t₀, ws.n, ws.Δ, ws.T)
    convoluer!(B, ws, ws.f)
    convoluer!(dBτ, ws, ws.dfτ)
    convoluer!(dBt, ws, ws.dfu)
    dBt .*= -1
    return nothing
end

"""
    forward(es, inst, ws) -> Vector{Float64}

Comptages attendus par canal sur toute la période (avant la fenêtre TAC), pour
une excitation du donneur : Σᵢ αᵢ·(IRF ⊛ Fᵢ). La somme vaut Σᵢ αᵢτᵢ.
"""
function forward(es::ExpSum, inst::InstrumentParams, ws::EspaceTravail)
    μ = zeros(inst.n_bins)
    B = zeros(inst.n_bins)
    for (α, τ) in zip(es.α, es.τ)
        base!(B, ws, Float64(τ), inst.t₀)
        μ .+= α .* B
    end
    return μ
end

"""Même modèle par convolution directe O(n²), générique (différentiation automatique)."""
function forward_direct(es::ExpSum, irf::AbstractVector, t₀, n::Int, Δ, T)
    Tt = promote_type(eltype(es.α), eltype(es.τ), typeof(t₀))
    g = zeros(Tt, n)
    for (α, τ) in zip(es.α, es.τ)
        for m in 0:n-1
            g[m+1] += α * noyau(m * Δ - t₀, τ, Δ, T)
        end
    end
    μ = zeros(Tt, n)
    for k in 0:n-1, j in 0:n-1
        μ[k+1] += irf[j+1] * g[mod(k - j, n)+1]
    end
    return μ
end

forward_direct(es::ExpSum, inst::InstrumentParams; t₀ = inst.t₀) =
    forward_direct(es, inst.irf, t₀, inst.n_bins, pas_temps(inst), inst.T_rep)

"""Masque la fenêtre TAC : les canaux hors fenêtre ne reçoivent aucun photon."""
function masquer!(μ::AbstractVector, fenetre::UnitRange{Int})
    for k in eachindex(μ)
        k in fenetre || (μ[k] = 0.0)
    end
    return μ
end
