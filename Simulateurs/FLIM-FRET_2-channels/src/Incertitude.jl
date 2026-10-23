# Incertitude.jl — propager l'incertitude des paramètres d'entrée (changement 2).
#
# Trois sources, par ordre d'importance attendue :
#   1. l'IRF et t₀ : la titration réajustée avec une IRF dilatée, ou un t₀ décalé,
#      donne d'autres paramètres ; chacun définit un « monde » cohérent (paramètres
#      ET IRF) dans lequel on refait les figures 1, 2 et 5 → bandes ;
#   2. la covariance des ajustements (global pour τ_D0, τ_DA ; isotherme pour
#      x₀, x_sat, K, n) : tirages de vecteurs de paramètres → variance paramétrique,
#      rapportée à côté de la variance photonique ;
#   3. les paramètres supposés (τ_A, Φ_A, η_YD) : balayages de l'étape 5, dont on
#      vérifie que le classement des estimateurs survit (`classement`).

"""
    monde_reajuste(m, irf_mat, r; s) -> (Modele, irf_mat)

Le monde défini par un réajustement de la titration (une ligne de
sorties/batterie/reel_perturbations.csv) : ses paramètres, l'IRF dilatée de `s`
avec laquelle il a été obtenu, et son t₀. x_sat est ramené à ≥ 0,001 (le
réajustement à IRF × 0,9 donne −0,003, sans sens physique).
"""
function monde_reajuste(m::Modele, irf_mat::Matrix{Float64}, r; s::Real = r.s_irf)
    p = remplacer(m.p; τ_D0 = r.tau_D0, τ_DA = r.tau_DA, x₀ = r.x0, x_sat = max(r.x_sat, 0.001), K_D = r.K, n = r.n)
    irf = s == 1 ? m.inst.irf : irf_dilatee(m.inst.irf, s)
    inst = remplacer(m.inst; irf = irf, t₀ = r.t0)
    mat = copy(irf_mat)
    s == 1 || (mat[:, 2] .= irf .* sum(irf_mat[:, 2]))
    return Modele(p, m.o, inst), mat
end

"""
    tirer_parametres(m, cov_τ, iso, J, rng) -> Vector{SensorParams}

J vecteurs de paramètres tirés de la loi normale des ajustements : (τ_D0, τ_DA)
de la covariance de l'ajustement global, (x₀, x_sat, K, n) de celle de
l'isotherme ; les deux blocs indépendants. Tirages hors domaine physique
(x_sat < 0, K ≤ 0, n ≤ 0, τ_DA ≥ τ_D0) rejetés et retirés.
"""
function tirer_parametres(m::Modele, cov_τ::AbstractMatrix, iso, J::Integer, rng::AbstractRNG)
    dτ = MvNormal([m.p.τ_D0, m.p.τ_DA], Symmetric(Matrix(cov_τ)))
    diso = MvNormal([iso.x₀, iso.x_sat, iso.K_D, iso.n], Symmetric(Matrix(iso.cov)))
    out = SensorParams[]
    while length(out) < J
        τ = rand(rng, dτ)
        θ = rand(rng, diso)
        (θ[2] >= 0 && θ[3] > 0 && θ[4] > 0 && 0 < τ[2] < τ[1] && θ[1] > θ[2]) || continue
        push!(out, remplacer(m.p; τ_D0 = τ[1], τ_DA = τ[2], x₀ = θ[1], x_sat = θ[2], K_D = θ[3], n = θ[4]))
    end
    return out
end

"""
    biais_sans_bruit(ests, cals, m_vrai, ctx, C) -> Dict(clé => Ĉ − C)

Lecture de chaque estimateur (calibré sur le monde nominal) sur les comptages
attendus du monde `m_vrai`, sans bruit : la limite N → ∞ du biais.
"""
function biais_sans_bruit(ests, cals, m_vrai::Modele, ctx::Contexte, C::Real)
    μB, μY = attendus(C, m_vrai, 1e6, :C, ctx.ws)
    out = Dict{String,Float64}()
    obs = Dict{Any,Tuple{Float64,Symbol}}()
    for (k, e) in ests
        o, st = get!(obs, e) do
            try
                observable(e, μB, μY, ctx)
            catch
                (NaN, :echec)
            end
        end
        out[k] = st === :echec ? NaN : inverser(cals[k], o)[1] - C
    end
    return out
end

"""Tau de Kendall (τ_b) entre deux classements donnés comme vecteurs de scores (plus petit = meilleur)."""
function kendall(a::AbstractVector, b::AbstractVector)
    n = length(a)
    c = d = ta = tb = 0
    for i in 1:n-1, j in i+1:n
        sa, sb = sign(a[i] - a[j]), sign(b[i] - b[j])
        if sa == 0 && sb == 0
            continue
        elseif sa == 0
            ta += 1
        elseif sb == 0
            tb += 1
        elseif sa == sb
            c += 1
        else
            d += 1
        end
    end
    return (c - d) / sqrt((c + d + ta) * (c + d + tb))
end
