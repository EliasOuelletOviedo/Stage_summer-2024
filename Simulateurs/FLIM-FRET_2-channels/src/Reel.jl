# Reel.jl — le « facteur 8 » sur les données réelles.
#
# Les histogrammes réels (≈ 10⁶ photons) sont amincis par tirage binomial à N
# photons : chaque photon est gardé avec la probabilité N/N_h. Le résultat est un
# vrai histogramme de N photons, avec l'IRF, le fond et les défauts des données,
# ce qu'aucune simulation ne garantit. Deux lectures y sont comparées :
#   2H : le pipeline actuel (mono-exponentiel du code de l'utilisateur + Hill) ;
#   4R : bi-exponentiel à τ_D0, τ_DA et t₀ figés (ajustement global), amplitude,
#        fraction courte et fond libres.
# Chacune est calibrée empiriquement sur les histogrammes complets, comme au banc.
# Le jumeau refait exactement la même chose sur des histogrammes tirés du modèle
# ajusté μ_h (fond compris) : si réel et jumeau concordent, le modèle décrit
# aussi le bruit des données, pas seulement leur forme moyenne.

"""Amincissement binomial : out[k] ~ Binomiale(n[k], p)."""
function amincir!(out::AbstractVector{Float64}, n::AbstractVector, p::Real, rng::AbstractRNG)
    @inbounds for k in eachindex(out, n)
        out[k] = n[k] > 0 ? Float64(rand(rng, Binomial(round(Int, n[k]), p))) : 0.0
    end
    return out
end

"""
    BiExpFigeFond(τ_D0, τ_DA, t₀, inst)

Estimateur 4R : μ = a·[B(τ_D0) + f·(B(τ_DA) − B(τ_D0))] + b sur la fenêtre, τ et
t₀ figés ; a, f et le fond b libres (Poisson, Levenberg-Marquardt). Renvoie f,
la fraction d'amplitude courte (le x de l'ajustement global).
"""
struct BiExpFigeFond <: Estimateur
    u::Vector{Float64}
    v::Vector{Float64}
    fenetre::UnitRange{Int}
end
nom(::BiExpFigeFond) = "4R biexp τ figés, fond libre"

function BiExpFigeFond(τ_D0::Real, τ_DA::Real, t₀::Real, inst::InstrumentParams)
    ws = EspaceTravail(inst)
    W = inst.fenetre
    BD0, BDA = zeros(ws.n), zeros(ws.n)
    base!(BD0, ws, Float64(τ_D0), Float64(t₀))
    base!(BDA, ws, Float64(τ_DA), Float64(t₀))
    return BiExpFigeFond(BD0[W], BDA[W] .- BD0[W], W)
end

function observable(e::BiExpFigeFond, nB, nY, ctx)
    n = nB[e.fenetre]
    N = sum(n)
    N > 0 || return (NaN, :echec)
    f0, st0, _ = profil_lineaire([CanalLineaire(n, e.u, e.v)], -0.5, 1.5)
    f0 = st0 === :echec ? 0.3 : clamp(f0, -0.4, 1.4)
    m = length(n)
    function modele!(μ, J, θ)
        a, f, b = θ
        @inbounds for i in 1:m
            g = e.u[i] + f * e.v[i]
            μ[i] = a * g + b
            J[i, 1] = g
            J[i, 2] = a * e.v[i]
            J[i, 3] = 1.0
        end
    end
    a0 = 0.98N / sum(e.u .+ f0 .* e.v)
    r = lm_poisson(modele!, [a0, f0, 0.02N / m], [0.0, -0.5, 0.0], [Inf, 1.5, Inf], n; maxiter = 200, tol_dec = 1e-8)
    r.converge || return (NaN, :echec)
    r.borne[2] && return (r.θ[2], :borne)
    return (r.θ[2], :ok)
end

"""Calibration par isotherme (4R) : x̂ → Ĉ, ramenée à [0, C_max]."""
function inverser_isotherme(iso, x::Real, C_max::Real)
    isnan(x) && return (NaN, false)
    hi, lo = iso.x₀, iso.x_sat                 # x décroît avec [Cl⁻]
    x >= hi && return (0.0, x > hi)
    x <= lo && return (Float64(C_max), true)
    C = iso.K_D * ((hi - x) / (x - lo))^(1 / iso.n)
    C > C_max && return (Float64(C_max), true)
    return (C, false)
end

"""
    lectures_amincies(source, hs, μs, Ns, R, inst, ctxs, e4, iso4, hill2; C_max, graine0) -> DataFrame

Pour chaque histogramme et chaque budget N : R tirages, chacun lu par 4R
(calibration `iso4`, isotherme) et 2H (calibration `hill2`, Hill).
`source = :reel` : amincissement binomial des histogrammes réels `hs` ;
`source = :jumeau` : tirages de Poisson du modèle ajusté `μs[h]` (fenêtre) remis
à N photons. Une ligne par (histogramme, N, tirage).
"""
function lectures_amincies(source::Symbol, hs::Vector{Histogramme}, μs, Ns, R::Int, inst::InstrumentParams,
                           ctxs::Vector{Contexte}, e4::BiExpFigeFond, iso4, hill2::NTuple{4,Float64}; C_max = 150.0, graine0 = 0)
    taches = [(h, N) for h in eachindex(hs) for N in Ns]
    lignes = Vector{Vector{NamedTuple}}(undef, length(taches))
    suivant = Threads.Atomic{Int}(1)
    W = inst.fenetre
    @sync for ctx in ctxs
        Threads.@spawn begin
            buf = zeros(inst.n_bins)
            while (i = Threads.atomic_add!(suivant, 1)) <= length(taches)
                h, N = taches[i]
                out = NamedTuple[]
                plein = hs[h].n
                Ntot = sum(plein[W])
                for r in 1:R
                    rng = StableRNG(graine(graine0 + 1000 * h + round(Int, log10(N) * 10), r))
                    if source === :reel
                        amincir!(buf, plein, N / Ntot, rng)
                    else
                        μ = zeros(inst.n_bins)
                        μ[W] .= μs[h] .* (N / sum(μs[h]))
                        tirer!(buf, μ, rng)
                    end
                    f, st4 = try
                        observable(e4, buf, buf, ctx)
                    catch
                        (NaN, :echec)
                    end
                    C4, hp4 = st4 === :echec ? (NaN, false) : inverser_isotherme(iso4, f, C_max)
                    p, _ = try
                        ajuster_mono_utilisateur(ctx.utilisateur, buf)
                    catch
                        ([NaN, NaN, NaN], NaN)
                    end
                    C2, hp2 = isnan(p[1]) ? (NaN, false) : inverser_hill(hill2, p[1], C_max)
                    push!(out, (source = String(source), histogramme = h, C_vrai = hs[h].C, dossier = hs[h].dossier,
                        replicat = hs[h].replicat, N = N, tirage = r, f_4R = f, C_4R = C4, hors_plage_4R = hp4,
                        tau_mono = p[1], C_2H = C2, hors_plage_2H = hp2))
                end
                lignes[i] = out
            end
        end
    end
    return DataFrame(reduce(vcat, lignes))
end

"""
Résumé par (source, concentration, N) : σ intra-histogramme (bruit de photons ;
moyenne quadratique sur les réplicats), σ robuste, biais, taux hors plage.
"""
function resumer_amincis(df::DataFrame)
    out = DataFrame()
    for g in groupby(df, [:source, :C_vrai, :N])
        isinf(g.C_vrai[1]) && continue
        ligne = Dict{Symbol,Any}(:source => g.source[1], :C_vrai => g.C_vrai[1], :N => g.N[1])
        for k in ("4R", "2H")
            Ĉ = g[!, "C_"*k]
            ok = isfinite.(Ĉ)
            parh = [Ĉ[ok.&(g.histogramme.==h)] for h in unique(g.histogramme)]
            σ_intra = sqrt(mean(var.(filter(v -> length(v) > 1, parh))))
            rob = sqrt(mean([(1.4826 * median(abs.(v .- median(v))))^2 for v in parh if length(v) > 1]))
            ligne[Symbol("sigma_", k)] = σ_intra
            ligne[Symbol("sigma_robuste_", k)] = rob
            ligne[Symbol("biais_", k)] = mean(Ĉ[ok]) - g.C_vrai[1]
            ligne[Symbol("hors_plage_", k)] = mean(g[!, "hors_plage_"*k][ok])
            ligne[Symbol("echec_", k)] = 1 - mean(ok)
        end
        ligne[:facteur] = ligne[:sigma_2H] / ligne[:sigma_4R]
        ligne[:facteur_robuste] = ligne[:sigma_robuste_2H] / ligne[:sigma_robuste_4R]
        push!(out, NamedTuple(ligne); cols = :union)
    end
    return sort!(out, [:source, :C_vrai, :N])
end
