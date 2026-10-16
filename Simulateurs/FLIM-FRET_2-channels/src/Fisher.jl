# Fisher.jl — dérivées analytiques des deux canaux, information de Fisher,
# borne de Cramér-Rao en mM.
#
# Les comptages sont poissoniens et indépendants :
#   I_ij = Σ_canaux Σ_k (1/μ_k) ∂μ_k/∂θ_i ∂μ_k/∂θ_j.
# Les dérivées s'obtiennent en dérivant les amplitudes (Photophysics.jl) et les
# noyaux (Instrument.jl) avant convolution, puis en convoluant. Le test 7 les
# compare à la différentiation automatique.

const NOMS_PHYSIQUES = (:x, :τ_D0, :τ_DA, :τ_A, :t0)

"""
    formes_et_derivees(x, m, ws; t₀) -> (mB, mY, dB, dY)

Formes des deux canaux par excitation du donneur (fenêtre appliquée) et leurs
dérivées : `dB[:, j]`, `dY[:, j]` pour j ∈ (x, τ_D0, τ_DA, τ_A, t₀).
M2 lent seulement (les durées de vie sont constantes).
"""
function formes_et_derivees(x::Real, m::Modele, ws::EspaceTravail; t₀::Float64 = m.inst.t₀)
    m.p.mecanisme === :M2_lent || error("dérivées analytiques écrites pour :M2_lent seulement")
    (m.p.β_A == 0 && m.p.β_D == 0) || error("dérivées analytiques écrites sans blanchiment")
    α, dα_dx, dα_dτD0, dα_dτDA, dα_dτA = d_amplitudes_M2_lent(x, m.p, m.o)
    τs = durees_M2_lent(m.p)                  # (τ_DA, τ_D0, τ_A)
    n = ws.n
    B = [zeros(n) for _ in 1:3]
    dBτ = [zeros(n) for _ in 1:3]
    dBt = [zeros(n) for _ in 1:3]
    for i in 1:3
        bases_derivees!(B[i], dBτ[i], dBt[i], ws, τs[i], t₀)
    end
    # canal jaune décalé de Δt_jaune : ses propres bases
    Δ = m.inst.Δt_jaune
    BY, dBτY, dBtY = Δ == 0 ? (B, dBτ, dBt) : ([zeros(n) for _ in 1:3], [zeros(n) for _ in 1:3], [zeros(n) for _ in 1:3])
    if Δ != 0
        for i in 1:3
            bases_derivees!(BY[i], dBτY[i], dBtY[i], ws, τs[i], t₀ + Δ)
        end
    end
    sorties = map(1:2) do c
        Bc, dBτc, dBtc = c == 1 ? (B, dBτ, dBt) : (BY, dBτY, dBtY)
        mc = zeros(n)
        d = zeros(n, 5)
        for i in 1:3
            mc .+= α[c, i] .* Bc[i]
            d[:, 1] .+= dα_dx[c, i] .* Bc[i]
            d[:, 2] .+= dα_dτD0[c, i] .* Bc[i]
            d[:, 3] .+= dα_dτDA[c, i] .* Bc[i]
            d[:, 4] .+= dα_dτA[c, i] .* Bc[i]
            d[:, 5] .+= α[c, i] .* dBtc[i]
        end
        d[:, 2] .+= α[c, 2] .* dBτc[2]         # τ_D0 est le terme 2
        d[:, 3] .+= α[c, 1] .* dBτc[1]         # τ_DA est le terme 1
        d[:, 4] .+= α[c, 3] .* dBτc[3]         # τ_A est le terme 3
        masquer!(mc, m.inst.fenetre)
        for j in 1:5
            masquer!(view(d, :, j), m.inst.fenetre)
        end
        (mc, d)
    end
    return sorties[1][1], sorties[2][1], sorties[1][2], sorties[2][2]
end

"""
    information_fisher(C, m, N, norm, ws; params, canaux, echelles) -> (I, noms)

Information de Fisher pour θ = (params physiques..., échelles...).
- `params` ⊂ (:x, :τ_D0, :τ_DA, :τ_A, :t0) ; x doit être le premier ;
- `canaux` : (:B,) ou (:B, :Y) ;
- `echelles` : `:commune` (le rapport bleu/jaune est connu et informatif) ou
  `:separees` (une échelle libre par canal, comme l'estimateur 5).
Les échelles sont paramétrées en log (∂μ/∂ log S = μ).
"""
function information_fisher(C::Real, m::Modele, N::Real, norm::Symbol, ws::EspaceTravail;
                            params = (:x,), canaux = (:B, :Y), echelles::Symbol = :commune)
    x = fraction_active(C, m.p)
    mB, mY, dB, dY = formes_et_derivees(x, m, ws)
    mB_ref = norm === :A ? formes(C_REF_NORMALISATION_A, sans_blanchiment(m), ws)[1] : mB
    S = echelle(norm, N, mB, mY, mB_ref)
    W = m.inst.fenetre
    cols = [findfirst(==(p), NOMS_PHYSIQUES) for p in params]
    nech = echelles === :commune ? 1 : length(canaux)
    P = length(cols) + nech
    I = zeros(P, P)
    for (ic, c) in enumerate(canaux)
        mc, dc = c === :B ? (mB, dB) : (mY, dY)
        for k in W
            μ = S * mc[k]
            μ > 0 || continue
            g = zeros(P)
            for (j, col) in enumerate(cols)
                g[j] = S * dc[k, col]
            end
            g[length(cols)+(echelles === :commune ? 1 : ic)] = μ
            I .+= (g * g') ./ μ
        end
    end
    noms = vcat(collect(params), echelles === :commune ? [:logS] : [Symbol(:logS_, c) for c in canaux])
    return I, noms
end

"""
Information de Fisher des seuls totaux (N_B, N_Y) : ce que voit le ratio.
Paramètres (x, log S).
"""
function information_fisher_totaux(C::Real, m::Modele, N::Real, norm::Symbol, ws::EspaceTravail)
    x = fraction_active(C, m.p)
    mB, mY, dB, dY = formes_et_derivees(x, m, ws)
    mB_ref = norm === :A ? formes(C_REF_NORMALISATION_A, sans_blanchiment(m), ws)[1] : mB
    S = echelle(norm, N, mB, mY, mB_ref)
    I = zeros(2, 2)
    for (mc, dc) in ((mB, dB), (mY, dY))
        Nbar = S * sum(mc)
        g = [S * sum(view(dc, :, 1)), Nbar]
        I .+= (g * g') ./ Nbar
    end
    return I
end

"""Variance minimale du premier paramètre (x) : [I⁻¹]₁₁ ; Inf si I est singulière."""
function variance_crlb(I::AbstractMatrix)
    F = cholesky(Symmetric(I); check = false)
    issuccess(F) || return Inf
    return inv(F)[1, 1]
end

"""Matrice de corrélation des paramètres à la borne (diagnostic d'identifiabilité)."""
function correlations_crlb(I::AbstractMatrix)
    Σ = inv(Symmetric(I))
    d = sqrt.(diag(Σ))
    return Σ ./ (d * d')
end

"""σ[Cl⁻] (mM) à partir de la variance de x : σ_C = σ_x / |dx/dC|."""
sigma_C(var_x::Real, C::Real, p::SensorParams) = sqrt(var_x) / abs(dfraction_dC(C, p))

"""
    crlb_mM(definition, C, m, N, norm, ws) -> σ_C

Les bornes utilisées dans les sorties :
- `:absolue` : les deux histogrammes, τ et t₀ connus, rapport bleu/jaune
  connu (échelle commune) : toute l'information disponible sur x ;
- `:bleu` : canal bleu seul, τ connus ;
- `:bleu_tau_libres` : canal bleu seul, τ_D0 et τ_DA inconnus ;
- `:conjoint` : deux canaux, échelles séparées (estimateur 5a) ;
- `:conjoint_tauA_libre` : idem avec τ_A inconnu (estimateur 5c) ;
- `:totaux` : seulement N_B et N_Y (ratio).
"""
function crlb_mM(def::Symbol, C::Real, m::Modele, N::Real, norm::Symbol, ws::EspaceTravail)
    v = if def === :absolue
        variance_crlb(information_fisher(C, m, N, norm, ws; params = (:x,), canaux = (:B, :Y), echelles = :commune)[1])
    elseif def === :bleu
        variance_crlb(information_fisher(C, m, N, norm, ws; params = (:x,), canaux = (:B,), echelles = :separees)[1])
    elseif def === :bleu_tau_libres
        variance_crlb(information_fisher(C, m, N, norm, ws; params = (:x, :τ_D0, :τ_DA), canaux = (:B,), echelles = :separees)[1])
    elseif def === :conjoint
        variance_crlb(information_fisher(C, m, N, norm, ws; params = (:x,), canaux = (:B, :Y), echelles = :separees)[1])
    elseif def === :conjoint_tauA_libre
        variance_crlb(information_fisher(C, m, N, norm, ws; params = (:x, :τ_A), canaux = (:B, :Y), echelles = :separees)[1])
    elseif def === :totaux
        variance_crlb(information_fisher_totaux(C, m, N, norm, ws))
    else
        error("définition de borne inconnue : $def")
    end
    return sigma_C(v, C, m.p)
end
