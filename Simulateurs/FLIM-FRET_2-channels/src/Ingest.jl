# Ingest.jl — ta titration bleue → paramètres du modèle.
#
# Ajustement global : tous les histogrammes d'une protéine à la fois, τ_D0,
# τ_DA et le décalage t₀ liés entre toutes les concentrations, une fraction x
# et une intensité libres par histogramme, plus un fond constant par
# histogramme. Puis une isotherme à deux plateaux sur {x, C}.

const CHEMIN_SDTFILE = normpath(joinpath(@__DIR__, "..", "..", "..", "tests", "SdtFile.jl"))
include(CHEMIN_SDTFILE)

const DOSSIER_TITRATION = "/Users/eliasouellet-oviedo/Documents/Maîtrise/Données/2026-07-20_SCLM"
const FICHIER_IRF = "irf_20260710_16x_sclm.sdt"
const DOSSIERS_CONC = ["0 mM", "2.5 mM", "5 mM", "10 mM", "20 mM", "40 mM", "80 mM", "146 mM", "F"]
# « F » : concentration infinie, comme dans analyse_SCLM_2.jl (100000 mM là-bas).
const VALEURS_CONC = [0.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 146.0, Inf]

"""Un histogramme de la titration (canal 1 du .sdt, 256 canaux)."""
struct Histogramme
    proteine::String
    dossier::String
    C::Float64
    replicat::Int
    n::Vector{Float64}
    fichier::String
end
Histogramme(p, d, C, r, n) = Histogramme(p, d, C, r, n, "")

"""
Histogramme du canal 1 d'un .sdt, lu comme `read_sdt_frame` (tests/lifetime_analysis.jl) :
un bloc de plusieurs courbes est sommé jusqu'à `adc_re` canaux.
"""
function lire_sdt(chemin::AbstractString)
    brut = read(chemin)
    sdt = SdtFile.read_sdt(brut, basename(chemin))
    res = Int(sdt.measure_info[1].adc_re)
    plat = Float64.(vec(sdt.data[1]))
    h = length(plat) == res ? plat : vec(sum(reshape(plat, res, :), dims = 2))
    off = Int(sdt.header.meas_desc_block_offset)
    frame_time = reinterpret(Float32, brut[off+216:off+219])[1]
    return h, res, frame_time
end

"""
    lire_acquisition(chemin) -> (duree_s, arret_sur_temps)

Mode d'acquisition d'un .sdt : temps de collecte (col_t, s) et arrêt sur temps
(stopt). C'est ce qui tranche « temps fixé ou photons fixés » : à temps fixé,
l'intensité totale porte une information que le modèle doit reproduire.
"""
function lire_acquisition(chemin::AbstractString)
    sdt = SdtFile.read_sdt(read(chemin), basename(chemin))
    mi = sdt.measure_info[1]
    return (Float64(mi.col_t), mi.stopt != 0)
end

"""
IRF exactement comme `load_irf_from_sdt` : médiane soustraite, négatifs mis à
zéro, axe [t_ns comptages] (pas de 12,5/256 ns quand frame_time n'est pas un
pas de canal plausible, ce qui est le cas de ton fichier).
"""
function charger_irf_utilisateur(chemin::AbstractString = joinpath(DOSSIER_TITRATION, FICHIER_IRF))
    c, res, ft = lire_sdt(chemin)
    c .-= round(Int, round(median(c)))
    c[c.<=0] .= 0
    n = min(length(c), res)
    pas = Float64(ft) * 1e9
    (!isfinite(pas) || pas <= 0 || pas * n > 1000) && (pas = 12.5 / n)
    return hcat(collect(0:n-1) .* pas, c[1:n])
end

"""IRF du simulateur : la même, normalisée à somme 1."""
irf_normalisee(irf_mat::Matrix{Float64}) = irf_mat[:, 2] ./ sum(irf_mat[:, 2])

"""Tous les histogrammes d'une protéine (3 réplicats × 9 concentrations)."""
function lire_titration(proteine::AbstractString; dossier = DOSSIER_TITRATION)
    out = Histogramme[]
    for (d, C) in zip(DOSSIERS_CONC, VALEURS_CONC)
        rep = joinpath(dossier, proteine, d)
        for f in sort(filter(endswith(".sdt"), readdir(rep)))
            m = match(r"_c(\d+)\.sdt$", f)
            h, _, _ = lire_sdt(joinpath(rep, f))
            push!(out, Histogramme(proteine, d, C, m === nothing ? 0 : parse(Int, m[1]), h, joinpath(rep, f)))
        end
    end
    return out
end

# -----------------------------------------------------------------------------
# Ajustement global
# -----------------------------------------------------------------------------

"""
    ResultatGlobal

`τ_D0, τ_DA, t₀` liés ; par histogramme : `x`, `A` (intensité), `b` (fond par
canal) et leurs incertitudes ; `chi2` : χ²_r de Pearson par histogramme et
global ; `μ` : modèle par histogramme (fenêtre) ; `t0s` : le t₀ de chaque
histogramme (tous égaux, sauf avec `t0_par_histo`).
"""
struct ResultatGlobal
    τ_D0::Float64
    τ_DA::Float64
    t₀::Float64
    σ_glob::Vector{Float64}
    x::Vector{Float64}
    σ_x::Vector{Float64}
    A::Vector{Float64}
    b::Vector{Float64}
    chi2::Vector{Float64}
    chi2_global::Float64
    dev::Float64
    converge::Bool
    μ::Vector{Vector{Float64}}
    cov::Matrix{Float64}
    t0s::Vector{Float64}
    σ_t0s::Vector{Float64}
end

"""
    ajuster_global(histos, inst; init, t0_fixe, t0_par_histo, tau_D0_fixe, tau_DA_fixe) -> ResultatGlobal

Modèle pour l'histogramme h : μ_h = A_h·[x_h·B(τ_DA) + (1 − x_h)·B(τ_D0)] + b_h,
B = IRF ⊛ noyau (unité d'amplitude à t = 0), donc x_h est directement la
fraction d'amplitude du déclin court. Canal bleu pur donneur : la fuite de
l'accepteur dans le bleu (η_BA/η_BD ≈ 0,6 %) est négligée ici.

Options (changements 1 et 2) : `t0_fixe` impose le décalage de l'IRF (test d'un
canal de décalage) ; `t0_par_histo` libère un t₀ par histogramme (diagnostic :
un t₀ qui dérive avec [Cl⁻] serait un artefact corrélé au signal) ;
`tau_D0_fixe`, `tau_DA_fixe` imposent les durées de vie. L'IRF de l'ajustement
est celle de `inst`, qui peut différer de celle des données (hypothèse H5).
"""
function ajuster_global(histos::Vector{Histogramme}, inst::InstrumentParams;
                        τ_D0_init = 3.0, τ_DA_init = 0.8, t₀_init = 0.0, maxiter = 500,
                        t0_fixe::Union{Nothing,Real} = nothing, t0_par_histo::Bool = false,
                        tau_D0_fixe::Union{Nothing,Real} = nothing, tau_DA_fixe::Union{Nothing,Real} = nothing)
    ws = EspaceTravail(inst)
    W = collect(inst.fenetre)
    m = length(W)
    H = length(histos)
    nt = t0_par_histo ? H : 1
    P = 2 + nt + 3H
    it(h) = t0_par_histo ? 2 + h : 3
    iX(h) = 2 + nt + h
    iA(h) = 2 + nt + H + h
    ib(h) = 2 + nt + 2H + h
    n = vcat([h.n[W] for h in histos]...)

    nf = ws.n
    nouveaux() = [zeros(nf) for _ in 1:nt]
    BD0, dBD0τ, dBD0t, BDA, dBDAτ, dBDAt = nouveaux(), nouveaux(), nouveaux(), nouveaux(), nouveaux(), nouveaux()
    function modele!(μ, J, θ)
        τ_D0, τ_DA = θ[1], θ[2]
        for j in 1:nt
            bases_derivees!(BD0[j], dBD0τ[j], dBD0t[j], ws, τ_D0, θ[2+j])
            bases_derivees!(BDA[j], dBDAτ[j], dBDAt[j], ws, τ_DA, θ[2+j])
        end
        fill!(J, 0.0)
        for h in 1:H
            j = t0_par_histo ? h : 1
            x, A, b = θ[iX(h)], θ[iA(h)], θ[ib(h)]
            r0 = (h - 1) * m
            bD0, bDA = BD0[j], BDA[j]
            @inbounds for (i, k) in enumerate(W)
                r = r0 + i
                f = x * bDA[k] + (1 - x) * bD0[k]
                μ[r] = A * f + b
                J[r, 1] = A * (1 - x) * dBD0τ[j][k]
                J[r, 2] = A * x * dBDAτ[j][k]
                J[r, it(h)] = A * (x * dBDAt[j][k] + (1 - x) * dBD0t[j][k])
                J[r, iX(h)] = A * (bDA[k] - bD0[k])
                J[r, iA(h)] = f
                J[r, ib(h)] = 1.0
            end
        end
        return nothing
    end

    borne(v, a, b) = v === nothing ? (a, b) : (Float64(v), Float64(v))
    (lD0, hD0), (lDA, hDA), (lt, ht) = borne(tau_D0_fixe, 0.2, 12.0), borne(tau_DA_fixe, 0.01, 6.0), borne(t0_fixe, -1.5, 1.5)
    lo = vcat([lD0, lDA], fill(lt, nt), fill(-0.5, H), fill(0.0, H), fill(0.0, H))
    hi = vcat([hD0, hDA], fill(ht, nt), fill(1.5, H), fill(Inf, H), fill(Inf, H))
    t0_depart = t0_fixe === nothing ? t₀_init : Float64(t0_fixe)
    τD0_depart = tau_D0_fixe === nothing ? τ_D0_init : Float64(tau_D0_fixe)
    departs_τDA = tau_DA_fixe === nothing ? (τ_DA_init, 0.3, 0.6, 1.2, 1.8) : (Float64(tau_DA_fixe),)

    meilleur = nothing
    for τDA0 in departs_τDA
        τDA0 < τD0_depart || continue
        θ0 = zeros(P)
        θ0[1:2] = [τD0_depart, τDA0]
        θ0[3:2+nt] .= t0_depart
        base!(BD0[1], ws, τD0_depart, t0_depart)
        base!(BDA[1], ws, τDA0, t0_depart)
        for (h, hist) in enumerate(histos)
            θ0[iX(h)] = 0.5
            θ0[ib(h)] = 0.05 * minimum(hist.n[W])
            f = 0.5 .* (BDA[1][W] .+ BD0[1][W])
            θ0[iA(h)] = (sum(hist.n[W]) - m * θ0[ib(h)]) / sum(f)
        end
        r = lm_poisson(modele!, θ0, lo, hi, n; maxiter, tol_dec = 1e-8)
        (meilleur === nothing || r.dev < meilleur.dev) && (meilleur = r)
    end
    r = meilleur
    θ = r.θ
    libres = findall(lo .< hi)
    cov = fill(0.0, P, P)
    try
        cov[libres, libres] = inv(Symmetric(r.fisher[libres, libres]))
    catch
        cov .= NaN
    end
    μs = [r.μ[(h-1)*m+1:h*m] for h in 1:H]
    chi2 = [sum((histos[h].n[W] .- μs[h]) .^ 2 ./ μs[h]) / (m - 3 - (length(libres) - 3H) / H) for h in 1:H]
    chi2_g = sum(sum((histos[h].n[W] .- μs[h]) .^ 2 ./ μs[h]) for h in 1:H) / (H * m - length(libres))
    σ = sqrt.(max.(diag(cov), 0))
    t0s = t0_par_histo ? θ[3:2+H] : fill(θ[3], H)
    σt0s = t0_par_histo ? σ[3:2+H] : fill(σ[3], H)
    return ResultatGlobal(θ[1], θ[2], mean(t0s), [σ[1], σ[2], σ[3]], θ[iX(1):iX(H)], σ[iX(1):iX(H)],
        θ[iA(1):iA(H)], θ[ib(1):ib(H)], chi2, chi2_g, r.dev, r.converge, μs, cov, t0s, σt0s)
end

"""
    ajuster_libre(h, inst, t₀; τ_init) -> (θ, χ²_r)

Même modèle que l'ajustement global mais τ_D0 et τ_DA libres pour cet
histogramme seul (t₀ fixé à la valeur globale) : la différence de χ² avec
l'ajustement global est le coût du liage des durées de vie.
θ = (τ_D0, τ_DA, x, A, b).
"""
function ajuster_libre(h::Histogramme, inst::InstrumentParams, t₀::Float64; τ_init = (3.0, 0.5), x_init = 0.3)
    ws = EspaceTravail(inst)
    W = collect(inst.fenetre)
    n = h.n[W]
    nf = ws.n
    BD0, dBD0, BDA, dBDA, poub = zeros(nf), zeros(nf), zeros(nf), zeros(nf), zeros(nf)
    function modele!(μ, J, θ)
        τ_D0, τ_DA, x, A, b = θ
        bases_derivees!(BD0, dBD0, poub, ws, τ_D0, t₀)
        bases_derivees!(BDA, dBDA, poub, ws, τ_DA, t₀)
        @inbounds for (i, k) in enumerate(W)
            f = x * BDA[k] + (1 - x) * BD0[k]
            μ[i] = A * f + b
            J[i, 1] = A * (1 - x) * dBD0[k]
            J[i, 2] = A * x * dBDA[k]
            J[i, 3] = A * (BDA[k] - BD0[k])
            J[i, 4] = f
            J[i, 5] = 1.0
        end
    end
    base!(BD0, ws, τ_init[1], t₀)
    base!(BDA, ws, τ_init[2], t₀)
    A0 = sum(n) / sum(x_init .* BDA[W] .+ (1 - x_init) .* BD0[W])
    r = lm_poisson(modele!, [τ_init[1], τ_init[2], x_init, A0, 0.0], [0.2, 0.01, -0.5, 0.0, 0.0],
        [12.0, 6.0, 1.5, Inf, Inf], n; maxiter = 500, tol_dec = 1e-8)
    χ2 = sum((n .- r.μ) .^ 2 ./ r.μ) / (length(W) - 5)
    return r.θ, χ2
end

"""
χ²_r que ton ajustement mono-exponentiel devrait donner si M2 lent était
exact : ton ajustement appliqué au modèle global sans bruit (écart
systématique), plus 1 pour le bruit de Poisson.
"""
function chi2_mono_predit(inst_u::InstanceUtilisateur, μ_fenetre::Vector{Float64}, inst::InstrumentParams)
    h = zeros(inst.n_bins)
    h[inst.fenetre] .= μ_fenetre
    p, _ = ajuster_mono_utilisateur(inst_u, h)
    return 1 + chi2_utilisateur(inst_u, p, h), p[1]
end

# -----------------------------------------------------------------------------
# Isotherme à deux plateaux sur {x, C}
# -----------------------------------------------------------------------------

isotherme(C, θ) = [isinf(c) ? θ[2] : θ[2] + (θ[1] - θ[2]) * θ[3]^θ[4] / (θ[3]^θ[4] + c^θ[4]) for c in C]

"""
    ajuster_isotherme(C, x; n_libre=true) -> (x₀, x_sat, K_D, n, cov)

Moindres carrés non pondérés sur tous les réplicats (comme ton ajustement de
Hill). « F » entre comme C = ∞ et fixe donc x_sat.
"""
function ajuster_isotherme(C::Vector{Float64}, x::Vector{Float64}; n_libre::Bool = true)
    x0 = mean(x[C.==minimum(C)])
    xs = mean(x[C.==maximum(C)])
    θ0 = [x0, xs, 20.0, 1.0]
    lo = [-0.5, -0.5, 0.01, n_libre ? 0.2 : 1.0]
    hi = [1.5, 1.5, 1e4, n_libre ? 5.0 : 1.0]
    θ, ok, cov = lm_moindres_carres(isotherme, θ0, lo, hi, C, x)
    return (x₀ = θ[1], x_sat = θ[2], K_D = θ[3], n = θ[4], cov = cov, converge = ok)
end

# -----------------------------------------------------------------------------
# Phasor des données (diagnostic M2 : segment de droite)
# -----------------------------------------------------------------------------

"""
Phasor (première harmonique, ω = 2π/T) d'un histogramme restreint à la
fenêtre, corrigé de l'IRF (division par le phasor de l'IRF). Canaux centrés
en (k − ½)Δ.
"""
function phasor(h::AbstractVector, inst::InstrumentParams; fenetre = inst.fenetre)
    Δ = pas_temps(inst)
    ω = 2π / inst.T_rep
    z = sum(h[k] * cis(ω * (k - 0.5) * Δ) for k in fenetre) / sum(h[k] for k in fenetre)
    zi = sum(inst.irf[k] * cis(ω * (k - 0.5) * Δ) for k in eachindex(inst.irf)) / sum(inst.irf)
    return z / zi
end

"""Distance d'un point au segment (droite) [a, b] du plan complexe, relative à |b − a|."""
function distance_droite(z::Complex, a::Complex, b::Complex)
    d = b - a
    return abs(imag(conj(d) * (z - a))) / abs(d)^2
end
