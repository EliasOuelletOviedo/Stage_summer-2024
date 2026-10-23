# physio.jl — le modèle physiologique (étape 2) : le chlore devient un bilan entre KCC2 qui
# l'extrait et eNpHR qui le fait entrer.
#
#   dC/dt = J(u) − k_pompe·(C − C_repos) + ξ,   J(u) = J_max·u/(u + u_½),   u ∈ [0, 1]
#
# Discrétisé au pas d'acquisition : C_{k+1} = C_k + Δt·[J(u_k) − k_pompe·(C_k − C_repos)] + ξ_k.
# La commande ne peut que faire monter le chlore ; pour le faire descendre, il faut couper et
# attendre KCC2.

using StableRNGs

struct Physio
    C_repos::Float64
    k::Float64              # k_pompe (s⁻¹)
    J_max::Float64          # mM·s⁻¹
    u_demi::Float64
    σ_ξ::Float64            # bruit physiologique, par pas
    dt::Float64
end

"""La physiologie de params.toml ; `k` et `J_max` se remplacent pour les balayages de robustesse."""
function physio(p::Params; k = p.toml["physiologie"]["k_pompe_par_s"], J_max = p.toml["physiologie"]["J_max_mM_par_s"])
    ph = p.toml["physiologie"]
    return Physio(ph["C_repos_mM"], Float64(k), Float64(J_max), ph["u_demi"], ph["sigma_xi_mM"], p.dt)
end

flux(u, ph::Physio) = ph.J_max * u / (u + ph.u_demi)
dflux(u, ph::Physio) = ph.J_max * ph.u_demi / (u + ph.u_demi)^2

"""Un pas de la cellule : C_{k+1} à partir de C_k et de la commande u_k."""
pas_cellule(C, u, ph::Physio, rng) = C + ph.dt * (flux(u, ph) - ph.k * (C - ph.C_repos)) + ph.σ_ξ * randn(rng)

"""La commande qui maintient C à l'équilibre : J(u) = k·(C − C_repos)."""
function u_maintien(C, ph::Physio)
    J = ph.k * (C - ph.C_repos)
    J <= 0 && return 0.0
    J >= ph.J_max && return 1.0
    return ph.u_demi * J / (ph.J_max - J)
end

"""
Gain proportionnel qui place le pôle dominant de la boucle fermée à τ_BF (au premier ordre,
au maintien de C_ref) : 1/τ_BF ≈ k_pompe + Kp·∂J/∂u.
"""
function Kp_nominal(p::Params, ph::Physio)
    r = p.toml["regulateur"]
    return (1 / r["tau_BF_cible_s"] - ph.k) / dflux(u_maintien(r["C_ref_reglage_mM"], ph), ph)
end

"""La consigne à l'instant t (créneaux de params.toml)."""
function consigne(t, p::Params)
    c = p.toml["scenario"]["consigne"]
    i = findlast(x -> t >= x[1] - 1e-9, c)
    return c[i][2]
end

"""Le nom de la phase de consigne à l'instant t."""
phase(t, p::Params) = (ph = p.toml["scenario"]["phases"]; ph[findlast(x -> t >= x[2] - 1e-9, ph)][1])

"""Le bruit physiologique du pas k de la répétition r : un flux à part de celui des photons."""
rng_cellule(k, r, p::Params) = StableRNG(p.toml["course"]["graine"] + 7_000_000 + 10_000r + k)
