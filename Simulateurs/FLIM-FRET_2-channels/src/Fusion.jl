# Fusion.jl — fusionner le ratio et la forme (changement 4).
#
# Le ratio est précis mais dérive avec le gain relatif des détecteurs (≈ 0,7 mM par
# % à 20 mM) ; la forme (5a, échelles bleue et jaune libres) est 6 à 8 fois plus
# bruitée mais insensible au gain. Les sensibilités sont orthogonales : on met
# ℓ = log g dans l'état d'un filtre de Kalman étendu, état ([Cl⁻], ℓ), avec deux
# mesures par cycle :
#   z₁ = log(N_J/N_B)     h₁ = ℓ + log R(C)    (R : courbe de calibration du ratio)
#   z₂ = Ĉ de la forme    h₂ = C
# Sur une fenêtre longue devant l'échantillonnage et courte devant la dérive,
# l'écart systématique entre les deux lectures est l'erreur de gain.

"""
    DeuxLectures(m, ctx; C_max, N_ref)

Estimateur 7 : renvoie les deux lectures (ratio et forme 5a) et leurs variances
issues de Fisher (bornes :totaux et :conjoint, tabulées en [Cl⁻] à N_ref photons
bleus puis remises à l'échelle du nombre de photons bleus observé). La forme est
aussi rendue en x, son espace naturel : 5a y est sans biais (test 6), alors que
l'inversion non linéaire de l'isotherme ajoute un biais de Jensen à N fini.
Son `observable` est la fusion statique (pondération inverse des variances),
pour figurer dans les grilles ; la fusion dynamique est `filtre_fusion`.
"""
struct DeuxLectures <: Estimateur
    forme::Conjoint
    p::SensorParams
    cal_ratio::Calibration
    cal_forme::Calibration
    C_tab::Vector{Float64}
    σ_ratio::Vector{Float64}      # mM, à N_ref photons bleus (normalisation C)
    σ_forme::Vector{Float64}      # mM
    σx_forme::Vector{Float64}     # sur x
    N_ref::Float64
    C_max::Float64
end
nom(::DeuxLectures) = "7 fusion ratio + forme"

function DeuxLectures(m::Modele, ctx::Contexte; C_max = 150.0, N_ref = 1e4)
    forme = Conjoint(m, ctx.ws)
    cr = calibrer(Ratio(), m, ctx; C_max)
    cf = calibrer(forme, m, ctx; C_max)
    C_tab = vcat(0.25, 0.5:0.5:4.5, 5.0:1.0:C_max)
    σr = [crlb_mM(:totaux, C, m, N_ref, :C, ctx.ws) for C in C_tab]
    σf = [crlb_mM(:conjoint, C, m, N_ref, :C, ctx.ws) for C in C_tab]
    σx = σf .* abs.(dfraction_dC.(C_tab, Ref(m.p)))
    return DeuxLectures(forme, m.p, cr, cf, C_tab, σr, σf, σx, N_ref, C_max)
end

"""σ tabulé à la concentration C (interpolation linéaire, bornée), remis à l'échelle de N_B photons bleus."""
function sigma_tab(C_tab, σ, C, N_ref, N_B)
    Cc = clamp(C, C_tab[1], C_tab[end])
    j = clamp(searchsortedlast(C_tab, Cc), 1, length(C_tab) - 1)
    w = (Cc - C_tab[j]) / (C_tab[j+1] - C_tab[j])
    return ((1 - w) * σ[j] + w * σ[j+1]) * sqrt(N_ref / max(N_B, 1.0))
end

"""
    lectures(e, nB, nY, ctx) -> (C_ratio, var_ratio, C_forme, var_forme, z_logR, var_logR, x_forme, N_B)

Les deux lectures et leurs variances de Fisher ; plus la mesure brute du ratio
(log N_J/N_B et sa variance de Poisson) et la lecture de forme en x, pour le
filtre à état augmenté.
"""
function lectures(e::DeuxLectures, nB, nY, ctx)
    NB, NY = sum(nB), sum(nY)
    Cr = (NB > 0 && NY > 0) ? inverser(e.cal_ratio, NY / NB)[1] : NaN
    xf, st = observable(e.forme, nB, nY, ctx)
    st === :echec && (xf = NaN)
    Cf = isnan(xf) ? NaN : inverser(e.cal_forme, xf)[1]
    vr = sigma_tab(e.C_tab, e.σ_ratio, isnan(Cr) ? 20.0 : Cr, e.N_ref, NB)^2
    vf = sigma_tab(e.C_tab, e.σ_forme, isnan(Cf) ? 20.0 : Cf, e.N_ref, NB)^2
    zl = (NB > 0 && NY > 0) ? log(NY / NB) : NaN
    return (Cr, vr, Cf, vf, zl, 1 / max(NB, 1) + 1 / max(NY, 1), xf, NB)
end

function observable(e::DeuxLectures, nB, nY, ctx)
    Cr, vr, Cf, vf = lectures(e, nB, nY, ctx)
    isnan(Cr) && isnan(Cf) && return (NaN, :echec)
    isnan(Cf) && return (Cr, :ok)
    isnan(Cr) && return (Cf, :ok)
    return ((Cr / vr + Cf / vf) / (1 / vr + 1 / vf), :ok)
end

"""Calibration identité (l'observable de l'estimateur 7 est déjà une concentration)."""
calibration_identite(C_max) = Calibration(collect(0.0:0.5:C_max), collect(0.0:0.5:C_max), true, Float64(C_max), nothing)

"""log R(C) et sa dérivée, à partir de la courbe de calibration du ratio."""
function logR(cal::Calibration, C)
    Cc = clamp(C, cal.C[1], cal.C[end])
    j = clamp(searchsortedlast(cal.C, Cc), 1, length(cal.C) - 1)
    w = (Cc - cal.C[j]) / (cal.C[j+1] - cal.C[j])
    R = (1 - w) * cal.O[j] + w * cal.O[j+1]
    dR = (cal.O[j+1] - cal.O[j]) / (cal.C[j+1] - cal.C[j])
    return log(R), dR / R
end

"""
    filtre_fusion(e, mesures, dt; q_C, q_ℓ, P0, mode) -> (Ĉ, ℓ̂, σ_C)

Filtre de Kalman étendu sur la suite de mesures (une `lectures` par pas).
`mode` : `:fusion` (état [Cl⁻], ℓ ; deux mesures), `:ratio` (état [Cl⁻] seul,
gain supposé constant et égal à 1, mesure z₁), `:forme` (état [Cl⁻], mesure z₂).
Mesures : z₁ = log(N_J/N_B), h₁ = ℓ + log R(C) ; z₂ = x̂ de la forme (5a),
h₂ = x(C), variance de Fisher de x remise à l'échelle des photons observés.
q_C : variance de la marche aléatoire de [Cl⁻] par pas (mM²) ; q_ℓ : celle de ℓ.
"""
function filtre_fusion(e::DeuxLectures, mes::Vector, dt::Real; q_C::Real, q_ℓ::Real, P0 = (100.0, 0.05^2), mode::Symbol = :fusion)
    n = length(mes)
    Ĉ, ℓ̂, σC = fill(NaN, n), fill(NaN, n), fill(NaN, n)
    # départ : première lecture valide de la mesure utilisée
    col = mode === :forme ? 3 : 1
    i0 = findfirst(t -> isfinite(t[col]), mes)
    x = [i0 === nothing ? 20.0 : mes[i0][col], 0.0]
    P = [P0[1] 0.0; 0.0 P0[2]]
    Q = [q_C 0.0; 0.0 (mode === :fusion ? q_ℓ : 0.0)]
    for k in 1:n
        P = P + Q
        Cr, vr, Cf, vf, zl, vl, xf, NB = mes[k]
        if mode !== :forme && isfinite(zl)
            lR, dlR = logR(e.cal_ratio, x[1])
            ℓ = mode === :fusion ? x[2] : 0.0
            Hm = mode === :fusion ? [dlR 1.0] : [dlR 0.0]
            S = (Hm*P*Hm')[1] + vl
            K = (P * Hm') ./ S
            x = x .+ vec(K) .* (zl - (ℓ + lR))
            P = (I - K * Hm) * P
        end
        if mode !== :ratio && isfinite(xf)
            Cc = clamp(x[1], 0.0, e.C_max)
            hx = fraction_active(Cc, e.p)
            dx = dfraction_dC(max(Cc, 1e-3), e.p)
            vx = sigma_tab(e.C_tab, e.σx_forme, Cc, e.N_ref, NB)^2
            Hm = [dx 0.0]
            S = (Hm*P*Hm')[1] + vx
            K = (P * Hm') ./ S
            x = x .+ vec(K) .* (xf - hx)
            P = (I - K * Hm) * P
        end
        x[1] = clamp(x[1], 0.0, e.C_max)
        Ĉ[k], ℓ̂[k], σC[k] = x[1], x[2], sqrt(max(P[1, 1], 0))
    end
    return Ĉ, ℓ̂, σC
end
