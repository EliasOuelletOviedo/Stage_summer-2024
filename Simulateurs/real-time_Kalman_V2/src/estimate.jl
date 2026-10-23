# estimate.jl — les lectures (J/B, 5a), leurs calibrations et leur covariance de Fisher.
#
# Une lecture reçoit le couple d'histogrammes, rien d'autre : ni C(t) ni g(t). Elle rend
# l'observable brute (pas la concentration), sa variance de Fisher à l'estimation courante,
# et sa linéarisation ∂z/∂n sur les canaux des deux fenêtres. Cette dernière donne la
# covariance entre deux lectures quelconques des mêmes photons : Σₖ aₖ·bₖ·μₖ.

"""Ce que chaque lecture renvoie au filtre (z, σ², ok), plus ∂z/∂n et μ̂ pour les covariances croisées."""
struct Lecture
    z::Float64
    σ²::Float64
    ok::Bool
    a::Vector{Float64}   # ∂z/∂nₖ, fenêtre bleue puis fenêtre jaune
    μ::Vector{Float64}   # comptages attendus estimés sur ces canaux
end
echec(m) = Lecture(NaN, NaN, false, fill(NaN, 2m), fill(NaN, 2m))

"""J/B : somme des comptages de chaque canal, fond soustrait, puis rapport. Aucun paramètre."""
struct CfgJB
    p::Params
end

"""
5a : un ajustement par canal, durées de vie figées, IRF du canal, fond figé, échelle libre
par canal (A_B et A_J indépendants absorbent g). μ = a_d·b_diffusion + a₁·b_rapide +
a₂·b_lente + fond : la diffusion (copie de l'IRF, τ_diff figé, base normalisée à un photon)
a son amplitude libre, comme dans le générateur. La fraction d'amplitude p = a₁/(a₁ + a₂)
(bleu), q de même (jaune), exclut la diffusion. Les deux canaux ne partagent aucun
paramètre : l'ajustement conjoint se sépare en deux. Observable : wₚ·p + w_q·q.
Sans `diffusion_5a`, la colonne de diffusion disparaît (le modèle à 4 paramètres d'avant).
"""
struct Cfg5a
    p::Params
    Xb::Matrix{Float64}
    Xj::Matrix{Float64}
    w::NTuple{2,Float64}
end

function Cfg5a(p::Params, w = (1.0, 0.0))
    W = p.fenetre
    base(τ, irf, t0) = convoluer([(1.0, τ)], irf, t0, p)[W]
    dif(irf, t0) = reglage(p, "lectures", "diffusion_5a") ? (b = base(p.τ.diff, irf, t0); b ./ sum(b)) : zeros(length(W), 0)
    return Cfg5a(p, hcat(dif(p.irf_b, p.t0[1]), base(p.τ.τ2, p.irf_b, p.t0[1]), base(p.τ.τ3, p.irf_b, p.t0[1])),
                 hcat(dif(p.irf_j, p.t0[2]), base(p.τ.τJ, p.irf_j, p.t0[2]), base(p.τ.τA, p.irf_j, p.t0[2])), w)
end

function estimer(hb, hj, cfg::CfgJB)
    p, W = cfg.p, cfg.p.fenetre
    m = length(W)
    B = sum(hb[W]) - sum(p.fond_b[W])
    J = sum(hj[W]) - sum(p.fond_j[W])
    B > 0 || return echec(m)
    z = J / B
    a = vcat(fill(-z / B, m), fill(1 / B, m))
    μ = vcat(hb[W], hj[W])
    return Lecture(z, sum(a .^ 2 .* μ), true, a, μ)
end

"""
Maximum de vraisemblance de Poisson pour μ = X·a + f (Fisher scoring, pas réduit de moitié
tant que la vraisemblance ne monte pas). Rend (a, information de Fisher, μ, convergé).
"""
function ajuster_canal(n, X, f, itmax)
    logL(μ) = sum(n .* log.(μ) .- μ)
    k = size(X, 2)
    a = fill(max(sum(n) - sum(f), 1.0) / k, k) ./ vec(sum(X; dims = 1))
    μ = X * a .+ f
    F = zeros(k, k)
    for _ in 1:itmax
        s = X' * (n ./ μ .- 1)
        F = X' * (X ./ μ)
        δ = F \ s
        dot(s, δ) < 1e-10 && return a, F, μ, true
        t, L0 = 1.0, logL(μ)
        while t > 1e-6
            μn = X * (a .+ t .* δ) .+ f
            if all(>(0), μn) && logL(μn) >= L0
                a, μ = a .+ t .* δ, μn
                break
            end
            t /= 2
        end
        t <= 1e-6 && return a, F, μ, false
    end
    return a, F, μ, false
end

"""
Fraction d'amplitude a₁/(a₁ + a₂) d'un canal (les deux dernières colonnes de X) et sa
linéarisation ∂/∂nₖ = ∇ᵀF⁻¹Xₖ/μₖ. F est l'information de Fisher complète, amplitude de la
diffusion comprise : la nuisance est inversée avec le reste, pas ignorée.
"""
function fraction(n, X, f, cfg::Cfg5a)
    a, F, μ, ok = ajuster_canal(n, X, f, reglage(cfg.p, "lectures", "iterations_max"))
    k = length(a)
    a1, a2 = a[k-1], a[k]
    s = a1 + a2
    x = a1 / s
    lo, hi = reglage(cfg.p, "lectures", "bornes_fraction")
    ok = ok && s > 0 && lo < x < hi
    ∇ = zeros(k)
    ∇[k-1], ∇[k] = a2 / s^2, -a1 / s^2
    return x, (X * (F \ ∇)) ./ μ, μ, ok
end

function estimer(hb, hj, cfg::Cfg5a)
    p, W = cfg.p, cfg.p.fenetre
    x, ab, μb, okb = fraction(hb[W], cfg.Xb, p.fond_b[W], cfg)
    y, aj, μj, okj = fraction(hj[W], cfg.Xj, p.fond_j[W], cfg)
    wp, wq = cfg.w
    a = vcat(wp .* ab, wq .* aj)
    μ = vcat(μb, μj)
    return Lecture(wp * x + wq * y, sum(a .^ 2 .* μ), okb && (wq == 0 || okj), a, μ)
end

"""Covariance de Fisher de deux lectures des mêmes histogrammes."""
covariance(l1::Lecture, l2::Lecture) = sum(l1.a .* l2.a .* (l1.μ .+ l2.μ) ./ 2)

"""R complet (2 × 2) des lectures fusionnées : variances propres, covariance croisée de Fisher."""
function R_jointe(l1::Lecture, l2::Lecture)
    c = covariance(l1, l2)
    return [l1.σ² c; c l2.σ²]
end

"""Une lecture et sa courbe de calibration (générateur sans bruit, g = 1, au budget de la course)."""
struct Lecteur{C}
    cfg::C
    courbe::Courbe
end

function grille(p::Params)
    a, b, h = reglage(p, "lectures", "grille_calibration_mM")
    return collect(a:h:b)
end

courbe(cfg, p::Params) = Courbe(grille(p), [estimer(histogrammes_c(c, 1.0, p)..., cfg).z for c in grille(p)])

"""
Les deux lecteurs d'un budget. L'observable de 5a est p seul, sauf si q réduit σ_C d'au
moins `gain_minimal_jaune` à C_reference ; on prend alors la combinaison de variance
minimale, w ∝ (p′/σₚ², q′/σ_q²), normalisée pour que z′ = 1 mM⁻¹ en C_reference.
"""
function lecteurs(p::Params)
    Cr = reglage(p, "lectures", "C_reference_mM")
    zp, zq = (courbe(Cfg5a(p, w), p) for w in ((1.0, 0.0), (0.0, 1.0)))
    lp, lq = (estimer(histogrammes_c(Cr, 1.0, p)..., Cfg5a(p, w)) for w in ((1.0, 0.0), (0.0, 1.0)))
    dp, dq = derivee(zp, Cr), derivee(zq, Cr)
    σC = (p = sqrt(lp.σ²) / abs(dp), q = sqrt(lq.σ²) / abs(dq), pq = 1 / sqrt(dp^2 / lp.σ² + dq^2 / lq.σ²))
    jaune = σC.pq < (1 - reglage(p, "lectures", "gain_minimal_jaune")) * σC.p
    w = jaune ? (dp / lp.σ², dq / lq.σ²) ./ (dp^2 / lp.σ² + dq^2 / lq.σ²) : (1.0, 0.0)
    cS = Courbe(zp.x, w[1] .* zp.y .+ w[2] .* zq.y)
    cfgR, cfgS = CfgJB(p), Cfg5a(p, w)
    # R de Fisher attendu sur la grille (sans bruit, g = 1) : le filtre l'évalue en son estimation
    Rg = [R_jointe(estimer(histogrammes_c(c, 1.0, p)..., cfgR), estimer(histogrammes_c(c, 1.0, p)..., cfgS)) for c in grille(p)]
    return (JB = Lecteur(cfgR, courbe(cfgR, p)), S = Lecteur(cfgS, cS), σC_reference = σC, w = w, C = grille(p), R = Rg)
end

"""
R attendu au point (C, g) : interpolation linéaire sur la grille, la variance du rapport
multipliée par g² et la covariance par g (z_R est proportionnel à g, z_S n'en dépend pas).
C'est ce R que reçoit le filtre. Celui de chaque tirage (`Lecture.σ²`) dépend du bruit du
tirage : des poids corrélés aux mesures biaisent la moyenne que fait le filtre.
"""
function R_attendue(lec, C, g)
    x = lec.C
    i = clamp(searchsortedlast(x, C), 1, length(x) - 1)
    u = clamp((C - x[i]) / (x[i+1] - x[i]), 0.0, 1.0)
    R = (1 - u) .* lec.R[i] .+ u .* lec.R[i+1]
    return [R[1, 1]*g^2 R[1, 2]*g; R[2, 1]*g R[2, 2]]
end
