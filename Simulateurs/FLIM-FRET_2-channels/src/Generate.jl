# Generate.jl — de [Cl⁻] aux deux histogrammes entiers.
#
# Les sept étapes du plan : (1) x(C) ; (2) amplitudes analytiques par canal ;
# (3) déclin périodisé ; (4) convolution circulaire avec l'IRF ; (5) intégration
# par canal [3–5 : Instrument.jl] ; (6) mise à l'échelle photonique ; (7) Poisson.
# Plus la fenêtre TAC : les canaux 1–13 et 243–256 de tes .sdt sont vides.

const C_REF_NORMALISATION_A = 20.0   # mM ; voir `echelle`

"""
    formes(C, m, ws; t₀) -> (mB, mY)

Comptages attendus par excitation du donneur dans les deux canaux, fenêtre TAC
appliquée (zéros hors fenêtre). Garde-fou : erreur si un μ_k est négatif, ce
qui signalerait des paramètres incohérents (τ_DA > τ_D0 par exemple).
"""
function formes(C::Real, m::Modele, ws::EspaceTravail; t₀::Float64 = m.inst.t₀)
    x = fraction_active(C, m.p)
    return formes_x(x, m, ws; t₀)
end

function formes_x(x::Real, m::Modele, ws::EspaceTravail; t₀::Float64 = m.inst.t₀)
    esB, esY = expsums(x, m.p, m.o)
    inst = t₀ == m.inst.t₀ ? m.inst : remplacer(m.inst; t₀ = t₀)
    mB = forward(esB, inst, ws)
    # le détecteur jaune peut être en retard de Δt_jaune sur le bleu
    mY = forward(esY, inst.Δt_jaune == 0 ? inst : remplacer(inst; t₀ = t₀ + inst.Δt_jaune), ws)
    for μ in (mB, mY)
        lim = 1e-10 * maximum(abs, μ)
        minimum(μ) < -lim && error("comptage attendu négatif (min $(minimum(μ))) : paramètres incohérents")
        μ .= max.(μ, 0.0)
        masquer!(μ, m.inst.fenetre)
    end
    return mB, mY
end

"""
    echelle(norm, N, mB, mY, mB_ref)

Facteur d'échelle photonique selon la normalisation du budget :
- `:A` temps d'acquisition fixe : N = photons bleus attendus à
  C_REF_NORMALISATION_A (20 mM) ; à toute autre concentration, les deux canaux
  suivent la physique. Les estimateurs à un canal y voient tes photons bleus
  actuels, le canal jaune est un bonus.
- `:B` photons totaux : N_bleu + N_jaune = N.
- `:C` photons du canal bleu : N_bleu = N.
"""
function echelle(norm::Symbol, N::Real, mB, mY, mB_ref)
    norm === :A && return N / sum(mB_ref)
    norm === :B && return N / (sum(mB) + sum(mY))
    norm === :C && return N / sum(mB)
    error("normalisation inconnue : $norm")
end

"""
    attendus(C, m, N, norm, ws; C_ref) -> (μB, μY)

Comptages attendus des deux canaux (le rapport bleu/jaune n'est pas libre :
il sort du modèle).
"""
function attendus(C::Real, m::Modele, N::Real, norm::Symbol, ws::EspaceTravail;
                  C_ref::Real = C_REF_NORMALISATION_A, t₀::Float64 = m.inst.t₀)
    mB, mY = formes(C, m, ws; t₀)
    # normalisation A : le temps d'acquisition est fixé par le capteur intact (un
    # échantillon blanchi émet moins pour le même temps)
    mB_ref = norm === :A ? formes(C_ref, sans_blanchiment(m), ws; t₀)[1] : mB
    S = echelle(norm, N, mB, mY, mB_ref)
    return S .* mB, S .* mY
end

"""
Germe d'une réalisation : même germe pour les deux canaux d'une réalisation,
germes indépendants entre réalisations et entre cellules de la grille.
StableRNG : flux identique d'une version de Julia à l'autre.
"""
graine(cellule::Integer, r::Integer) = UInt64(1_000_003) * UInt64(cellule) + UInt64(r) + UInt64(0x5EED)

"""Tirage de Poisson canal par canal (les canaux hors fenêtre ont μ = 0)."""
function tirer!(n::AbstractVector{Float64}, μ::AbstractVector{Float64}, rng::AbstractRNG)
    @inbounds for k in eachindex(n, μ)
        n[k] = μ[k] > 0 ? Float64(rand(rng, Poisson(μ[k]))) : 0.0
    end
    return n
end

function tirer_deux!(nB, nY, μB, μY, rng::AbstractRNG)
    tirer!(nB, μB, rng)
    tirer!(nY, μY, rng)
    return nothing
end

"""Le même monde sans blanchiment (référence de la normalisation A)."""
sans_blanchiment(m::Modele) = (m.p.β_A == 0 && m.p.β_D == 0) ? m : remplacer(m; p = remplacer(m.p; β_A = 0.0, β_D = 0.0))

"""
    derive_gain(rng, n, dt; τ_g, amplitude) -> Vector{Float64}

Gain du canal jaune relatif au bleu au cours du temps (changement 4) : marche
aléatoire sur ℓ = log g, ℓ_{k+1} = ℓ_k + σ√dt·ξ_k, avec σ = amplitude/√τ_g. Après
τ_g secondes, l'écart-type de ℓ vaut `amplitude` (5 % par défaut) : τ_g est
l'échelle de temps de la dérive des détecteurs (quelques minutes à l'heure).
g(0) = 1 (gain calibré au départ).
"""
function derive_gain(rng::AbstractRNG, n::Integer, dt::Real; τ_g::Real, amplitude::Real = 0.05)
    σ = amplitude / sqrt(τ_g)
    ℓ = zeros(n)
    for k in 2:n
        ℓ[k] = ℓ[k-1] + σ * sqrt(dt) * randn(rng)
    end
    return exp.(ℓ)
end
