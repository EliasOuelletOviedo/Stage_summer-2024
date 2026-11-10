# Calibrate.jl — courbes de calibration sans bruit, inversion, K apparents.
#
# Chaque estimateur reçoit SA courbe : le générateur en mode sans bruit sur une
# grille de concentrations, l'observable de l'estimateur à chaque point. C'est
# ce que tu ferais au banc, et c'est ce qui produit les K apparents divergents.

"""Grille de calibration : dense à bas [Cl⁻], jusqu'à C_max."""
function grille_calibration(C_max::Real = 150.0)
    g = vcat(0.0:0.05:2.0, 2.1:0.1:10.0, 10.25:0.25:40.0, 40.5:0.5:C_max)
    return unique(sort(vcat(g, C_max)))
end

"""
    Calibration

`C` croissant, `O` observable sans bruit correspondante (monotone), et pour
l'estimateur 2H les paramètres de Hill (inversion analytique, comme ton
pipeline) au lieu de la courbe exacte.
"""
struct Calibration
    C::Vector{Float64}
    O::Vector{Float64}
    croissante::Bool
    C_max::Float64
    hill::Union{Nothing,NTuple{4,Float64}}   # (bas, haut, K, n)
end

"""
    calibrer(e, m, ctx; C_max, N=1e6, hill=false) -> Calibration

Observable de `e` sur les comptages attendus (sans Poisson) du modèle `m`,
normalisation C à N photons bleus (les ajustements sans bruit ne dépendent
pas de N, sauf par les tolérances de ton optimiseur, d'où N grand).
`hill = true` (estimateur 2H) : inversion par une Hill ajustée comme dans ton
pipeline plutôt que par la courbe exacte.
"""
function calibrer(e::Estimateur, m::Modele, ctx::Contexte; C_max::Real = 150.0, N::Real = 1e6,
                  grille = grille_calibration(C_max), hill::Bool = false)
    O = similar(grille)
    for (i, C) in enumerate(grille)
        μB, μY = attendus(C, m, N, :C, ctx.ws)
        o, st = observable(e, μB, μY, ctx)
        st === :echec && error("calibration de $(nom(e)) : échec sans bruit à C = $C mM")
        O[i] = o
    end
    d = diff(O)
    croissante = sum(d) > 0
    if !(all(>(0), d) || all(<(0), d))
        # bruit numérique de l'optimiseur sur un plateau : on impose la monotonie
        n_viol = count(v -> croissante ? v <= 0 : v >= 0, d)
        @warn "courbe de calibration non strictement monotone pour $(nom(e)) ($n_viol violations) : monotonie imposée"
        for i in 2:length(O)
            O[i] = croissante ? max(O[i], O[i-1] + 1e-12) : min(O[i], O[i-1] - 1e-12)
        end
    end
    h = hill ? calibration_hill_utilisateur(e, m, ctx) : nothing
    return Calibration(collect(grille), O, croissante, Float64(C_max), h)
end

"""
    inverser(cal, o) -> (Ĉ, hors_plage)

Interpolation linéaire de la courbe inverse ; hors de la plage de la courbe,
Ĉ est ramené à 0 ou à C_max et `hors_plage` est vrai.
"""
function inverser(cal::Calibration, o::Real)
    isnan(o) && return (NaN, false)
    cal.hill !== nothing && return inverser_hill(cal.hill, o, cal.C_max)
    O, C = cal.O, cal.C
    if cal.croissante
        o <= O[1] && return (0.0, o < O[1])
        o >= O[end] && return (cal.C_max, true)
        j = searchsortedlast(O, o)
    else
        o >= O[1] && return (0.0, o > O[1])
        o <= O[end] && return (cal.C_max, true)
        j = searchsortedlast(-O, -o)
    end
    w = (o - O[j]) / (O[j+1] - O[j])
    return (C[j] + w * (C[j+1] - C[j]), false)
end

# -----------------------------------------------------------------------------
# Hill : K apparents, et calibration du pipeline actuel (estimateur 2H)
# -----------------------------------------------------------------------------

"""Hill à 4 paramètres, forme de analyse_SCLM_2.jl : bas + (haut − bas)·Cⁿ/(Kⁿ + Cⁿ)."""
hill_modele(C, θ) = [θ[1] + (θ[2] - θ[1]) * (isinf(c) ? 1.0 : c^θ[4] / (θ[3]^θ[4] + c^θ[4])) for c in C]

"""
    ajuster_hill(C, O) -> (bas, haut, K, n)

Moindres carrés non pondérés (comme ton `fit_hill`). Pour les K apparents, on
l'applique à la courbe sans bruit sur une grille large (plateaux bien définis).
"""
function ajuster_hill(C::AbstractVector, O::AbstractVector; n_libre::Bool = true)
    i0, i1 = argmin(C), argmax(C)
    θ0 = [O[i0], O[i1], 20.0, 1.0]
    lo = [-Inf, -Inf, 1e-3, n_libre ? 0.1 : 1.0]
    hi = [Inf, Inf, 1e5, n_libre ? 10.0 : 1.0]
    θ, _, _ = lm_moindres_carres(hill_modele, θ0, lo, hi, collect(Float64, C), collect(Float64, O))
    return (θ[1], θ[2], θ[3], θ[4])
end

"""
Inversion analytique d'une Hill : C = K·((O − bas)/(haut − O))^{1/n}, ramenée à
[0, C_max] ; au-delà d'un plateau, l'estimation colle à ce plateau (0 ou C_max).
"""
function inverser_hill(h::NTuple{4,Float64}, o::Real, C_max::Real)
    bas, haut, K, n = h
    if haut > bas
        o <= bas && return (0.0, o < bas)
        o >= haut && return (Float64(C_max), true)
    else
        o >= bas && return (0.0, o > bas)
        o <= haut && return (Float64(C_max), true)
    end
    C = K * ((o - bas) / (haut - o))^(1 / n)
    C > C_max && return (Float64(C_max), true)
    return (C, false)
end

"""
Ce que fait ton pipeline : τ mono-exponentiel sans bruit (équivalent de tes
millions de photons) aux concentrations de ta titration, F placé à
100 000 mM comme dans analyse_SCLM_2.jl, puis Hill à 4 paramètres.
"""
function calibration_hill_utilisateur(e::Estimateur, m::Modele, ctx::Contexte)
    Cs = [0.01, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 146.0, 100000.0]
    τs = map(Cs) do C
        μB, μY = attendus(C, m, 1e6, :C, ctx.ws)
        o, st = observable(e, μB, μY, ctx)
        st === :echec && error("calibration de Hill : échec sans bruit à C = $C")
        o
    end
    return ajuster_hill(Cs, τs)
end

"""
    K_apparent(e, m, ctx) -> (K_app, n_app, bas, haut)

Hill ajustée à la courbe sans bruit de l'estimateur sur une grille large
(0 à 50·K_D, log), pour la cartographie K_vrai → K_apparent.
"""
function K_apparent(e::Estimateur, m::Modele, ctx::Contexte; facteur = 50.0, npts = 120)
    K = m.p.K_D
    grille = vcat(0.0, exp.(range(log(K / 100), log(facteur * K); length = npts)))
    O = map(grille) do C
        μB, μY = attendus(C, m, 1e6, :C, ctx.ws)
        observable(e, μB, μY, ctx)[1]
    end
    bas, haut, Kapp, napp = ajuster_hill(grille, O)
    return (K_app = Kapp, n_app = napp, bas = bas, haut = haut, C = grille, O = O)
end

"""Calibrations de toute une liste (clé, estimateur) sur le modèle `m` ; « 2H » reçoit la Hill."""
function calibrer_tous(ests, m::Modele, ctx::Contexte; C_max::Real = 150.0)
    cals = Dict{String,Calibration}()
    for (k, e) in ests
        cals[k] = calibrer(e, m, ctx; C_max, hill = k == "2H")
    end
    return cals
end
