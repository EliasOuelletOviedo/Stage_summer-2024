# Ajustement.jl — maximum de vraisemblance de Poisson.
#
# Deux solveurs servent à tous les estimateurs, à l'ajustement global de la
# titration et aux tests :
#   - `lm_poisson` : Levenberg-Marquardt (pas de Fisher scoring amorti) avec
#     bornes par ensemble actif, pour les modèles non linéaires ;
#   - `profil_lineaire` : modèles μ_c = A_c(u_c + x v_c), où les échelles A_c
#     se profilent analytiquement et il ne reste qu'une recherche 1-D en x.
# Poisson et pas moindres carrés : à 300 photons sur 229 canaux, la plupart
# des canaux contiennent 0, 1 ou 2 photons.

"""Demi-déviance de Poisson Σ[μ − n − n·log(μ/n)] (même minimum que −log L, mais ≈ m/2 à l'optimum)."""
function demi_deviance(n::AbstractVector, μ::AbstractVector)
    s = 0.0
    @inbounds for k in eachindex(n, μ)
        s += μ[k] - n[k]
        n[k] > 0 && (s -= n[k] * log(μ[k] / n[k]))
    end
    return s
end

"""
    ResultatLM

`θ` : estimation ; `converge` ; `borne` : paramètres collés à une borne ;
`iterations` ; `dev` : demi-déviance finale ; `fisher` : JᵀWJ à l'optimum
(information de Fisher attendue, W = 1/μ) ; `μ` : modèle à l'optimum.
"""
struct ResultatLM
    θ::Vector{Float64}
    converge::Bool
    borne::BitVector
    iterations::Int
    dev::Float64
    fisher::Matrix{Float64}
    μ::Vector{Float64}
end

"""
    lm_poisson(modele!, θ0, lo, hi, n; maxiter=200, tol_dec=1e-9, λ0=1e-3)

`modele!(μ, J, θ)` remplit les comptages attendus μ et le jacobien J (m × p).
Convergence quand le décrément de Newton gᵀH⁻¹g (en unités de log-vraisemblance)
passe sous `tol_dec` : c'est la baisse de déviance encore possible.
Les paramètres dont le gradient pousse contre une borne sont gelés (ensemble
actif) ; ils sont signalés dans `borne`.
"""
function lm_poisson(modele!, θ0::AbstractVector, lo::AbstractVector, hi::AbstractVector, n::AbstractVector;
                    maxiter::Int = 200, tol_dec::Float64 = 1e-9, λ0::Float64 = 1e-3)
    p, m = length(θ0), length(n)
    θ = clamp.(Float64.(θ0), lo, hi)
    μ, J = zeros(m), zeros(m, p)
    μn, Jn = zeros(m), zeros(m, p)
    modele!(μ, J, θ)
    echec() = ResultatLM(θ, false, falses(p), 0, Inf, zeros(p, p), μ)
    (all(isfinite, μ) && all(>(0), μ)) || return echec()
    dev = demi_deviance(n, μ)
    λ = λ0
    converge = false
    it = 0
    H = zeros(p, p)
    g = zeros(p)
    Jw = zeros(m, p)
    for outer it in 1:maxiter
        mul!(g, transpose(J), 1 .- n ./ μ)    # gradient de la demi-déviance
        Jw .= J ./ μ
        mul!(H, transpose(J), Jw)             # information de Fisher JᵀWJ
        # ensemble actif : bornes contre lesquelles le gradient pousse
        libre = [!((θ[j] <= lo[j] && g[j] > 0) || (θ[j] >= hi[j] && g[j] < 0)) for j in 1:p]
        F = findall(libre)
        if isempty(F)
            converge = true
            break
        end
        HF = Symmetric(H[F, F])
        gF = g[F]
        dec = try
            dot(gF, HF \ gF)
        catch
            Inf
        end
        if dec < tol_dec
            converge = true
            break
        end
        D = [max(H[j, j], 1e-12) for j in F]
        accepte = false
        while λ < 1e14
            δ = zeros(p)
            δF = try
                -((Matrix(HF) + λ * Diagonal(D)) \ gF)
            catch
                λ *= 10
                continue
            end
            δ[F] .= δF
            θn = clamp.(θ .+ δ, lo, hi)
            modele!(μn, Jn, θn)
            if all(isfinite, μn) && all(>(0), μn)
                devn = demi_deviance(n, μn)
                if devn < dev || (devn == dev && θn != θ)
                    θ, dev = θn, devn
                    μ, μn = μn, μ
                    J, Jn = Jn, J
                    λ = max(λ / 3, 1e-12)
                    accepte = true
                    break
                end
            end
            λ *= 4
        end
        if !accepte
            # plus aucun pas ne fait baisser la déviance : optimum numérique
            converge = dec < 1e3 * tol_dec
            break
        end
    end
    # information de Fisher à l'optimum
    Jw .= J ./ μ
    mul!(H, transpose(J), Jw)
    borne = BitVector([θ[j] <= lo[j] || θ[j] >= hi[j] for j in 1:p])
    return ResultatLM(θ, converge, borne, it, dev, copy(H), copy(μ))
end

# -----------------------------------------------------------------------------
# Modèles linéaires : μ_c = A_c (u_c + x v_c), échelles profilées
# -----------------------------------------------------------------------------

"""Un canal d'un modèle linéaire en x : comptages observés n, forme u + x·v (fenêtre seulement)."""
struct CanalLineaire
    n::Vector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
end

"""
Log-vraisemblance profilée ℓ(x) = Σ_c [Σ_k n_k log(u_k + x v_k) − N_c log(U_c + x V_c)]
et ses deux premières dérivées. Les échelles A_c = N_c/(U_c + x V_c) sont
réinjectées (vraisemblance de Poisson profilée = multinomiale de la forme).
"""
function profil_ℓ(canaux::Vector{CanalLineaire}, x::Float64)
    ℓ = d1 = d2 = 0.0
    for c in canaux
        N = U = V = 0.0
        @inbounds for k in eachindex(c.n)
            w = c.u[k] + x * c.v[k]
            U += c.u[k]
            V += c.v[k]
            nk = c.n[k]
            N += nk
            nk == 0 && continue
            r = c.v[k] / w
            ℓ += nk * log(w)
            d1 += nk * r
            d2 -= nk * r * r
        end
        S = U + x * V
        ℓ -= N * log(S)
        d1 -= N * V / S
        d2 += N * (V / S)^2
    end
    return ℓ, d1, d2
end

"""Intervalle de x où toutes les formes restent strictement positives."""
function domaine_positif(canaux::Vector{CanalLineaire})
    a, b = -Inf, Inf
    for c in canaux
        @inbounds for k in eachindex(c.u)
            u, v = c.u[k], c.v[k]
            if v > 0
                a = max(a, -u / v)
            elseif v < 0
                b = min(b, -u / v)
            elseif u <= 0
                return (NaN, NaN)
            end
        end
    end
    return a, b
end

"""
    profil_lineaire(canaux, lo, hi; tol=1e-12) -> (x̂, statut, ℓ'')

Maximise la vraisemblance profilée sur [lo, hi] ∩ domaine positif : balayage
grossier (la somme de deux canaux n'est pas forcément unimodale), puis Newton
protégé par un encadrement de ℓ'. `statut` ∈ (:ok, :borne, :echec).
"""
function profil_lineaire(canaux::Vector{CanalLineaire}, lo::Float64, hi::Float64; tol::Float64 = 1e-12, nbalayage::Int = 9)
    a, b = domaine_positif(canaux)
    isnan(a) && return (NaN, :echec, NaN)
    marge = 1e-9 * max(1.0, abs(hi - lo))
    lo2, hi2 = max(lo, a + marge), min(hi, b - marge)
    lo2 < hi2 || return (NaN, :echec, NaN)
    xs = range(lo2, hi2; length = nbalayage)
    ℓs = [profil_ℓ(canaux, x)[1] for x in xs]
    i = argmax(ℓs)
    # encadrement du maximum entre les voisins du meilleur point
    g, d = xs[max(i - 1, 1)], xs[min(i + 1, nbalayage)]
    x = xs[i]
    _, d1g, _ = profil_ℓ(canaux, g)
    _, d1d, _ = profil_ℓ(canaux, d)
    if i == 1 && d1g <= 0
        return (lo2, :borne, profil_ℓ(canaux, lo2)[3])
    elseif i == nbalayage && d1d >= 0
        return (hi2, :borne, profil_ℓ(canaux, hi2)[3])
    end
    # ℓ' change de signe dans [g, d] (ℓ'(g) ≥ 0 ≥ ℓ'(d) sauf aux bords)
    for _ in 1:100
        _, d1, d2 = profil_ℓ(canaux, x)
        if d1 > 0
            g = x
        else
            d = x
        end
        xn = d2 < 0 ? x - d1 / d2 : (g + d) / 2
        (xn <= g || xn >= d) && (xn = (g + d) / 2)
        if abs(xn - x) <= tol * max(1.0, abs(x)) || d - g <= tol
            x = xn
            break
        end
        x = xn
    end
    statut = (x <= lo2 + 10marge || x >= hi2 - 10marge) ? :borne : :ok
    return (x, statut, profil_ℓ(canaux, x)[3])
end

# -----------------------------------------------------------------------------
# Moindres carrés (ajustements de Hill et de Langmuir, sans bruit de Poisson)
# -----------------------------------------------------------------------------

"""
    lm_moindres_carres(f, θ0, lo, hi, x, y; w=ones) -> (θ, converge, cov)

Levenberg-Marquardt pondéré, jacobien par différences avant (petits problèmes
de 4 paramètres seulement). `cov` : (JᵀWJ)⁻¹·χ²_r.
"""
function lm_moindres_carres(f, θ0, lo, hi, x, y; w = ones(length(y)), maxiter = 500, tol = 1e-14)
    θ = clamp.(Float64.(θ0), lo, hi)
    p = length(θ)
    r(θ) = sqrt.(w) .* (f(x, θ) .- y)
    function jac(θ)
        r0 = r(θ)
        Jm = zeros(length(y), p)
        for j in 1:p
            h = 1e-7 * max(abs(θ[j]), 1e-3)
            θh = copy(θ)
            θh[j] = θ[j] + h > hi[j] ? θ[j] - h : θ[j] + h
            Jm[:, j] = (r(θh) .- r0) ./ (θh[j] - θ[j])
        end
        return Jm
    end
    ρ = r(θ)
    s = sum(abs2, ρ)
    λ = 1e-3
    converge = false
    for _ in 1:maxiter
        Jm = jac(θ)
        H = Jm' * Jm
        gr = Jm' * ρ
        D = Diagonal(max.(diag(H), 1e-12))
        accepte = false
        while λ < 1e16
            δ = -((H + λ * D) \ gr)
            θn = clamp.(θ .+ δ, lo, hi)
            ρn = r(θn)
            sn = sum(abs2, ρn)
            if sn <= s
                rel = (s - sn) / max(s, 1e-300)
                θ, ρ, s = θn, ρn, sn
                λ = max(λ / 3, 1e-15)
                accepte = true
                if rel < tol
                    converge = true
                end
                break
            end
            λ *= 4
        end
        (!accepte || converge) && (converge = true; break)
    end
    Jm = jac(θ)
    dof = max(length(y) - p, 1)
    cov = try
        inv(Jm' * Jm) .* (s / dof)
    catch
        fill(NaN, p, p)
    end
    return θ, converge, cov
end
