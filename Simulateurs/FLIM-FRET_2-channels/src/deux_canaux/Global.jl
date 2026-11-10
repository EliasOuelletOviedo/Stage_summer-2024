# Global.jl — étape 2 : l'ajustement global d'une longueur d'onde.
#
# Un ajustement par longueur d'onde, sur toutes les concentrations et les deux canaux à
# la fois (la structure qu'impose un mélange statique, avec moins de paramètres libres
# qu'un ajustement point par point) :
#   bleu  : diffusion + τ₂ + τ₃
#   jaune : diffusion + τ_J + τ_A, avec τ_J = τ₂ quand le lien est imposé
# Durées de vie communes aux concentrations et aux canaux ; amplitudes libres par
# concentration et par canal ; un t₀ par canal, lié entre concentrations ; fond figé.
#
# La « diffusion » est la composante qu'un ajustement libre pousse vers τ → 0 : une copie
# de l'IRF, dont l'intensité absolue ne dépend pas de la concentration (lumière diffusée
# par le milieu, ou résidu de déconvolution). Ce n'est pas une population : elle est
# rapportée en photons et exclue des fractions d'amplitude.

const TAU_DIFFUSION = 0.001      # ns
const SEUIL_SPEARMAN_8 = 0.738   # |ρ| de Spearman significatif à 5 % (bilatéral) pour 8 points

"""
Résultat de l'ajustement global d'une longueur d'onde. Amplitudes `aB[c, (diffusion,
τ₂, τ₃)]` et `aJ[c, (diffusion, τ_J, τ_A)]` en photons/ns à t = 0 (la diffusion en
photons par unité de sa base) ; `IB`, `IJ` : photons de chaque composante dans la
fenêtre ; `μB`, `μJ` : modèle lissé de chaque somme (256 canaux, fond compris, zéro hors
fenêtre).
"""
struct ModeleGlobal
    λ::Int
    lien::Bool
    cl::Vector{Float64}
    τ2::Float64
    τ3::Float64
    τA::Float64
    τJ::Float64
    σ_τ::Vector{Float64}          # (τ₂, τ₃, τ_A, τ_J), corrigés de la surdispersion
    t0B::Float64
    t0J::Float64
    aB::Matrix{Float64}
    aJ::Matrix{Float64}
    IB::Matrix{Float64}
    IJ::Matrix{Float64}
    fondB::Vector{Vector{Float64}}
    fondJ::Vector{Vector{Float64}}
    μB::Vector{Vector{Float64}}
    μJ::Vector{Vector{Float64}}
    χ2B::Vector{Float64}
    χ2J::Vector{Float64}
    χ2::Float64
    deviance::Float64
    nparam::Int
    converge::Bool
    duree_somme_s::Vector{Float64}
end

espaces(c::Cellule) = (S.EspaceTravail(S.instrument_banc(c.irf[BLEU, :], c.pas_ns)),
                       S.EspaceTravail(S.instrument_banc(c.irf[JAUNE, :], c.pas_ns)))

"""
    ajuster_global(cellules; lien, τ_init, t0_init, figer_tau) -> ModeleGlobal

Cellules d'une même longueur d'onde (sommes des acquisitions). `τ_init = (τ₂, τ₃, τ_A,
τ_J)` ; `figer_tau` fige les durées de vie (diagnostic de t₀ par concentration).
Maximum de vraisemblance de Poisson (Levenberg-Marquardt, dérivées analytiques).
"""
function ajuster_global(cellules::Vector{Cellule}; lien::Bool, τ_init = (0.7, 2.6, 3.2, 1.0), t0_init = (0.0, 0.0),
                        figer_tau::Bool = false, maxiter::Int = 500)
    C = length(cellules)
    wsB, wsJ = espaces(cellules[1])
    W = collect(FENETRE)
    m = length(W)
    K = lien ? 3 : 4                       # θ[1:K] = τ₂, τ₃, τ_A (, τ_J)
    kJ = lien ? 1 : 4
    P = K + 2 + 6C
    ia(c, j) = K + 2 + 6(c - 1) + j        # j = 1…3 bleu, 4…6 jaune
    hB = [somme(c, BLEU) for c in cellules]
    hJ = [somme(c, JAUNE) for c in cellules]
    durees = [n_acq(c) * c.duree_s for c in cellules]
    fB = [fond(c, BLEU, durees[i]) for (i, c) in enumerate(cellules)]
    fJ = [fond(c, JAUNE, durees[i]) for (i, c) in enumerate(cellules)]
    n = wsB.n
    nv() = zeros(n)
    DB, dDBτ, dDBt, DJ, dDJτ, dDJt = nv(), nv(), nv(), nv(), nv(), nv()
    BB, dBBτ, dBBt = [nv() for _ in 1:K], [nv() for _ in 1:K], [nv() for _ in 1:K]
    BJ, dBJτ, dBJt = [nv() for _ in 1:K], [nv() for _ in 1:K], [nv() for _ in 1:K]
    function bases!(θ)
        S.bases_derivees!(DB, dDBτ, dDBt, wsB, TAU_DIFFUSION, θ[K+1])
        S.bases_derivees!(DJ, dDJτ, dDJt, wsJ, TAU_DIFFUSION, θ[K+2])
        for k in 1:K
            S.bases_derivees!(BB[k], dBBτ[k], dBBt[k], wsB, θ[k], θ[K+1])
            S.bases_derivees!(BJ[k], dBJτ[k], dBJt[k], wsJ, θ[k], θ[K+2])
        end
    end
    function modele!(μ, J, θ)
        bases!(θ)
        fill!(J, 0.0)
        for c in 1:C
            r0 = (c - 1) * 2m
            aBd, aB2, aB3, aJd, aJJ, aJA = (θ[ia(c, j)] for j in 1:6)
            @inbounds for (i, b) in enumerate(W)
                r = r0 + i
                μ[r] = aBd * DB[b] + aB2 * BB[1][b] + aB3 * BB[2][b] + fB[c][b]
                J[r, 1] += aB2 * dBBτ[1][b]
                J[r, 2] += aB3 * dBBτ[2][b]
                J[r, K+1] = aBd * dDBt[b] + aB2 * dBBt[1][b] + aB3 * dBBt[2][b]
                J[r, ia(c, 1)] = DB[b]
                J[r, ia(c, 2)] = BB[1][b]
                J[r, ia(c, 3)] = BB[2][b]
                r = r0 + m + i
                μ[r] = aJd * DJ[b] + aJJ * BJ[kJ][b] + aJA * BJ[3][b] + fJ[c][b]
                J[r, kJ] += aJJ * dBJτ[kJ][b]
                J[r, 3] += aJA * dBJτ[3][b]
                J[r, K+2] = aJd * dDJt[b] + aJJ * dBJt[kJ][b] + aJA * dBJt[3][b]
                J[r, ia(c, 4)] = DJ[b]
                J[r, ia(c, 5)] = BJ[kJ][b]
                J[r, ia(c, 6)] = BJ[3][b]
            end
        end
    end
    y = reduce(vcat, [vcat(hB[c][W], hJ[c][W]) for c in 1:C])
    θ0 = zeros(P)
    θ0[1:3] .= τ_init[1:3]
    K == 4 && (θ0[4] = τ_init[4])
    θ0[K+1], θ0[K+2] = t0_init
    bases!(θ0)
    for c in 1:C
        for (canal, h, f, Bs) in ((1, hB[c], fB[c], (DB, BB[1], BB[2])), (2, hJ[c], fJ[c], (DJ, BJ[kJ], BJ[3])))
            X = reduce(hcat, [B[W] for B in Bs])
            w = 1 ./ max.(h[W], 1.0)
            β = (X' * (X .* w)) \ (X' * ((h[W] .- f[W]) .* w))
            θ0[ia(c, 3(canal - 1) + 1)] = β[1]
            θ0[ia(c, 3(canal - 1) + 2)] = canal == 1 ? max(β[2], 1.0) : β[2]
            θ0[ia(c, 3(canal - 1) + 3)] = max(β[3], 1.0)
        end
    end
    lo, hi = fill(-Inf, P), fill(Inf, P)
    lo[1:K] .= 0.02
    hi[1:K] .= 12.0
    lo[K+1:K+2] .= -1.5
    hi[K+1:K+2] .= 1.5
    for c in 1:C
        lo[ia(c, 2)] = lo[ia(c, 3)] = lo[ia(c, 6)] = 0.0     # τ₂, τ₃ du bleu et τ_A du jaune ≥ 0
    end
    if figer_tau
        lo[1:K] .= θ0[1:K]
        hi[1:K] .= θ0[1:K]
    end
    r = S.lm_poisson(modele!, θ0, lo, hi, y; maxiter, tol_dec = 1e-8)
    θ = r.θ
    bases!(θ)
    libres = findall(lo .< hi)
    resid = (y .- r.μ) .^ 2 ./ r.μ
    χ2 = sum(resid) / (length(y) - length(libres))
    σθ = fill(NaN, P)
    try
        σθ[libres] = sqrt.(max.(diag(inv(Symmetric(r.fisher[libres, libres]))), 0.0) .* max(χ2, 1.0))
    catch
    end
    aB = [θ[ia(c, j)] for c in 1:C, j in 1:3]
    aJ = [θ[ia(c, j)] for c in 1:C, j in 4:6]
    sW(v) = sum(@view v[W])
    IB = [aB[c, j] * sW((DB, BB[1], BB[2])[j]) for c in 1:C, j in 1:3]
    IJ = [aJ[c, j] * sW((DJ, BJ[kJ], BJ[3])[j]) for c in 1:C, j in 1:3]
    μB = [S.masquer!(aB[c, 1] .* DB .+ aB[c, 2] .* BB[1] .+ aB[c, 3] .* BB[2] .+ fB[c], FENETRE) for c in 1:C]
    μJ = [S.masquer!(aJ[c, 1] .* DJ .+ aJ[c, 2] .* BJ[kJ] .+ aJ[c, 3] .* BJ[3] .+ fJ[c], FENETRE) for c in 1:C]
    ddl = m - 3 - (K + 2) / (2C)
    χ2B = [sum(resid[(c-1)*2m+1:(c-1)*2m+m]) / ddl for c in 1:C]
    χ2J = [sum(resid[(c-1)*2m+m+1:c*2m]) / ddl for c in 1:C]
    τJ = lien ? θ[1] : θ[4]
    σ_τ = [σθ[1], σθ[2], σθ[3], lien ? σθ[1] : σθ[4]]
    return ModeleGlobal(cellules[1].λ, lien, [c.cl for c in cellules], θ[1], θ[2], θ[3], τJ, σ_τ, θ[K+1], θ[K+2],
        aB, aJ, IB, IJ, fB, fJ, μB, μJ, χ2B, χ2J, χ2, 2 * r.dev, length(libres), r.converge, durees)
end

"""
Le lien entre la composante rapide du jaune et τ₂ du bleu tient-il ? Écart de déviance
(lié − libre) divisé par le χ²ᵣ du modèle libre (surdispersion), comparé au χ² à 1 degré
de liberté à 5 % (3,84). Règle écrite avant de voir les données.
"""
function tester_lien(lie::ModeleGlobal, libre::ModeleGlobal)
    ΔD = lie.deviance - libre.deviance
    F = ΔD / libre.χ2
    return (lambda_nm = lie.λ, deviance_lie = lie.deviance, deviance_libre = libre.deviance, chi2r_lie = lie.χ2,
        chi2r_libre = libre.χ2, delta_deviance = ΔD, statistique = F, seuil = 3.84, lien_tient = F < 3.84,
        tau2_bleu = libre.τ2, tauJ_libre = libre.τJ, tau2_lie = lie.τ2)
end

"""Corrélation de rang de Spearman (rangs moyens en cas d'égalité)."""
function spearman(x::AbstractVector, y::AbstractVector)
    rang(v) = begin
        p = sortperm(v)
        r = zeros(length(v))
        i = 1
        while i <= length(v)
            j = i
            while j < length(v) && v[p[j+1]] == v[p[i]]
                j += 1
            end
            r[p[i:j]] .= (i + j) / 2
            i = j + 1
        end
        r
    end
    return cor(rang(x), rang(y))
end

"""
Forme d'un profil de χ²ᵣ contre la concentration (règle écrite d'avance) : « plat » si
l'écart max − min est sous 10 % de la médiane ; « monotone » si |ρ de Spearman| > 0,738 ;
« en cloche » si le maximum est intérieur et que les deux extrêmes sont sous la médiane ;
sinon « irrégulier ».
"""
function forme_profil(cl::AbstractVector, χ2::AbstractVector)
    med = median(χ2)
    (maximum(χ2) - minimum(χ2)) < 0.1med && return "plat"
    abs(spearman(cl, χ2)) > SEUIL_SPEARMAN_8 && return "monotone"
    k = argmax(χ2)
    (1 < k < length(χ2) && χ2[1] < med && χ2[end] < med) && return "en cloche"
    return "irrégulier"
end

"""
t₀ de chaque concentration : chaque somme réajustée seule, durées de vie figées aux
valeurs globales, t₀ des deux canaux et amplitudes libres. Une dérive monotone avec la
concentration signale un artefact corrélé au signal.
"""
function t0_par_cellule(cellules::Vector{Cellule}, mg::ModeleGlobal)
    τ = (mg.τ2, mg.τ3, mg.τA, mg.τJ)
    out = map(cellules) do c
        r = ajuster_global([c]; lien = mg.lien, τ_init = τ, t0_init = (mg.t0B, mg.t0J), figer_tau = true)
        (r.t0B, r.t0J)
    end
    return first.(out), last.(out)
end

"""
Gabarits du mélange statique à deux états : les modèles lissés des deux concentrations
extrêmes, fond retiré, les deux canaux d'un état normalisés par le même total bleu
(μ = A·[(1 − f)·état₀ + f·état₁] garde le rapport J/B de chaque état).
"""
function gabarits(mg::ModeleGlobal)
    T = map((1, length(mg.cl))) do s
        b = mg.μB[s] .- mg.fondB[s]
        j = mg.μJ[s] .- mg.fondJ[s]
        S.masquer!(b, FENETRE)
        S.masquer!(j, FENETRE)
        nb = sum(b)
        (b ./ nb, j ./ nb)
    end
    W = FENETRE
    return S.FormesBanc(T[1][1][W], (T[2][1].-T[1][1])[W], T[1][2][W], (T[2][2].-T[1][2])[W], W)
end

"""
Le mélange statique à deux états décrit-il les concentrations intermédiaires ? Chaque
somme est ajustée par μ = A_B·[(1 − f)·état₀ + f·état₁] (bleu) et A_J·[…] (jaune), fonds
figés : f, A_B, A_J libres. Rend f et le χ²ᵣ de chaque cellule (les deux canaux). Puis le
même mélange à une seule échelle (le J/B imposé par la forme) : `tension_rapport_forme`
= (χ²_commune − χ²_séparées) / χ²ᵣ_séparées, à comparer à 3,84 (χ² à 1 degré, 5 %).
Au-delà, le rapport J/B et la forme du déclin ne désignent pas le même point de la
titration.
"""
function tester_melange(cellules::Vector{Cellule}, mg::ModeleGlobal)
    G = gabarits(mg)
    W = FENETRE
    m = length(W)
    map(enumerate(cellules)) do (c, cel)
        nB, nJ = somme(cel, BLEU), somme(cel, JAUNE)
        bB, bJ = mg.fondB[c][W], mg.fondJ[c][W]
        y = vcat(nB[W], nJ[W])
        modele!(μ, J, θ) = begin
            f, AB, AJ = θ
            fill!(J, 0.0)
            @inbounds for i in 1:m
                gB = G.uB[i] + f * G.vB[i]
                gJ = G.uJ[i] + f * G.vJ[i]
                μ[i] = AB * gB + bB[i]
                μ[m+i] = AJ * gJ + bJ[i]
                J[i, 1] = AB * G.vB[i]
                J[m+i, 1] = AJ * G.vJ[i]
                J[i, 2] = gB
                J[m+i, 3] = gJ
            end
        end
        AB0 = max(sum(nB[W]) - sum(bB), 1.0)
        r = S.lm_poisson(modele!, [0.5, AB0, AB0 * sum(G.uJ .+ 0.5 .* G.vJ)], [-1.0, 0.0, 0.0], [2.0, Inf, Inf], y; maxiter = 300)
        χ2s = sum((y .- r.μ) .^ 2 ./ r.μ)
        # même mélange, une seule échelle : le rapport J/B est alors celui que la forme impose
        commun!(μ, J, θ) = begin
            f, A = θ
            @inbounds for i in 1:m
                gB = G.uB[i] + f * G.vB[i]
                gJ = G.uJ[i] + f * G.vJ[i]
                μ[i] = A * gB + bB[i]
                μ[m+i] = A * gJ + bJ[i]
                J[i, 1] = A * G.vB[i]
                J[m+i, 1] = A * G.vJ[i]
                J[i, 2] = gB
                J[m+i, 2] = gJ
            end
        end
        rc = S.lm_poisson(commun!, [r.θ[1], r.θ[2]], [-1.0, 0.0], [2.0, Inf], y; maxiter = 300)
        χ2c = sum((y .- rc.μ) .^ 2 ./ rc.μ)
        (f = r.θ[1], chi2r = χ2s / (2m - 3), chi2r_global = (mg.χ2B[c] + mg.χ2J[c]) / 2, f_commune = rc.θ[1],
            tension_rapport_forme = (χ2c - χ2s) / (χ2s / (2m - 3)))
    end
end

"""
Résidus de Pearson de l'ajustement global, canal par canal (concentrations × canaux
TCSPC), et leur premier vecteur singulier : la forme de résidu commune aux
concentrations, son poids par concentration et la part de la somme des carrés qu'elle porte.
"""
function residus_svd(cellules::Vector{Cellule}, mg::ModeleGlobal)
    W = FENETRE
    out = Dict{String,Any}()
    for (canal, nom, μs) in ((BLEU, "bleu", mg.μB), (JAUNE, "jaune", mg.μJ))
        R = reduce(vcat, [((somme(c, canal)[W] .- μs[i][W]) ./ sqrt.(μs[i][W]))' for (i, c) in enumerate(cellules)])
        F = svd(R)
        u, v = F.U[:, 1] .* F.S[1], F.V[:, 1]
        if sum(u) < 0
            u, v = -u, -v
        end
        out[nom] = (vecteur = v, poids = u, part = F.S[1]^2 / sum(abs2, F.S))
    end
    return out
end

"""
Une ligne par cellule : photons et fractions d'intensité des trois composantes de chaque
canal (diffusion comprise), fractions d'amplitude des composantes de durée de vie
(diffusion exclue), χ²ᵣ, coefficient de fuite λ = photons de la composante rapide du jaune
/ photons de la composante rapide du bleu (= f_J(rapide)/f_B(rapide)·J/B).
"""
function table_global(mg::ModeleGlobal)
    df = DataFrame()
    for c in eachindex(mg.cl)
        IB, IJ = mg.IB[c, :], mg.IJ[c, :]
        aB2, aB3 = mg.aB[c, 2], mg.aB[c, 3]
        aJJ, aJA = mg.aJ[c, 2], mg.aJ[c, 3]
        push!(df, (lambda_nm = mg.λ, cl_mM = mg.cl[c], modele = mg.lien ? "lie" : "libre",
            tau2_ns = mg.τ2, tau3_ns = mg.τ3, tauJ_ns = mg.τJ, tauA_ns = mg.τA,
            photons_diffusion_bleu = IB[1], photons_tau2_bleu = IB[2], photons_tau3_bleu = IB[3],
            photons_diffusion_jaune = IJ[1], photons_tauJ_jaune = IJ[2], photons_tauA_jaune = IJ[3],
            intensite_diffusion_bleu = IB[1] / sum(IB), intensite_tau2_bleu = IB[2] / sum(IB), intensite_tau3_bleu = IB[3] / sum(IB),
            intensite_diffusion_jaune = IJ[1] / sum(IJ), intensite_tauJ_jaune = IJ[2] / sum(IJ), intensite_tauA_jaune = IJ[3] / sum(IJ),
            amplitude_tau2_bleu = aB2 / (aB2 + aB3), amplitude_tau3_bleu = aB3 / (aB2 + aB3),
            amplitude_tauJ_jaune = aJJ / (aJJ + aJA), amplitude_tauA_jaune = aJA / (aJJ + aJA),
            diffusion_bleu_par_s = IB[1] / mg.duree_somme_s[c], diffusion_jaune_par_s = IJ[1] / mg.duree_somme_s[c],
            lambda_fuite = IJ[2] / IB[2], chi2r_bleu = mg.χ2B[c], chi2r_jaune = mg.χ2J[c]))
    end
    return df
end

"""Coefficient de fuite constant sur la titration ? Pente relative contre [Cl⁻] et ρ de Spearman."""
function tester_fuite(t::DataFrame)
    λf = t.lambda_fuite
    x = t.cl_mM
    b = cov(x, λf) / var(x)
    ρ = spearman(x, λf)
    return (lambda_fuite_moyen = mean(λf), lambda_fuite_min = minimum(λf), lambda_fuite_max = maximum(λf),
        variation_pct_sur_plage = 100 * b * (maximum(x) - minimum(x)) / mean(λf), spearman = ρ,
        constant = abs(ρ) <= SEUIL_SPEARMAN_8 && (maximum(λf) - minimum(λf)) < 0.1 * mean(λf))
end
