# Lectures.jl — étape 3 : les neuf lectures d'un couple d'histogrammes (bleu, jaune).
#
# Chaque lecture rend une observable brute (un rapport, une durée de vie, une fraction,
# une coordonnée phasor) et un code d'échec. Le passage en mM se fait plus loin, par la
# courbe de calibration propre à chaque lecture (Calibrer.jl). Les nuisances connues (fond
# des détecteurs, diffusion) sont proportionnelles à la durée d'acquisition `t_s` de
# l'histogramme lu : 10 s pour une somme, 1 s pour une acquisition, p × 10 s pour un tirage.

const LECTURES = ["1", "1c", "2", "2J", "3", "4", "5a", "5d", "6"]
const NOMS_LECTURES = Dict(
    "1" => "1 rapport J/B", "1c" => "1c J/B corrigé de la fuite", "2" => "2 mono-exp., bleu (ton code)",
    "2J" => "2J mono-exp., jaune (ton code)", "3" => "3 bi-exp. τ libres, bleu", "4" => "4 bi-exp. τ figés, bleu",
    "5a" => "5a forme 2 canaux, échelles séparées", "5d" => "5d forme 2 canaux, échelle commune", "6" => "6 phasor, bleu")
"""Information exploitée par chaque lecture : la borne de Fisher qui lui sert de référence."""
const INFORMATION_LECTURE = Dict("1" => "rapport", "1c" => "toute", "2" => "bleu", "2J" => "jaune", "3" => "bleu",
    "4" => "bleu", "5a" => "forme", "5d" => "toute", "6" => "bleu")

# codes : 0 lecture valable ; 1 non-convergence ; 2 paramètre sur une borne ;
# 3 autre (pas de photons, composantes confondues) ; 4 hors de la plage calibrée et
# 5 lecture non monotone (attribués à l'inversion, Calibrer.jl)
const OK, NON_CONVERGENCE, BORNE, AUTRE, HORS_PLAGE, NON_MONOTONE = 0, 1, 2, 3, 4, 5

"""Tout ce qu'il faut pour lire un couple d'histogrammes à une longueur d'onde."""
struct Contexte
    λ::Int
    pas::Float64
    t0B::Float64
    t0J::Float64
    τ2::Float64
    τ3::Float64
    B2::Vector{Float64}             # bases (fenêtre, amplitude 1) aux durées de vie et t₀ globaux
    B3::Vector{Float64}
    JJ::Vector{Float64}
    JA::Vector{Float64}
    J3::Vector{Float64}             # jaune à τ₃ (fuite d'un donneur libéré, Blanchir.jl)
    DB::Vector{Float64}             # diffusion, normalisée à somme 1
    DJ::Vector{Float64}
    diffusion_par_s::NTuple{2,Float64}
    fond_par_s::Matrix{Float64}     # [canal, canal TCSPC], photons par seconde
    G::S.FormesBanc                 # gabarits du mélange à deux états (5a, 5d)
    P0::ComplexF64                  # axe du phasor : états extrêmes
    P1::ComplexF64
    ws::Vector{Any}                 # espaces de travail du bleu, un par fil (bi-exp. libre)
    monoB::Vector{Any}              # ton vec_to_lifetime, une instance par fil et par canal
    monoJ::Vector{Any}
end

"""
    contexte(mg, cellule; fils, mono) -> Contexte

`cellule` : n'importe quelle cellule de la longueur d'onde (IRF, fond, pas). `mono =
false` saute le chargement de ton code (lectures 2 et 2J rendues NaN).
"""
function contexte(mg::ModeleGlobal, cel::Cellule; fils::Int = Threads.nthreads(), mono::Bool = true)
    W = FENETRE
    instB = S.instrument_banc(cel.irf[BLEU, :], cel.pas_ns)
    instJ = S.instrument_banc(cel.irf[JAUNE, :], cel.pas_ns)
    wsB, wsJ = S.EspaceTravail(instB), S.EspaceTravail(instJ)
    base(ws, τ, t0) = (b = zeros(ws.n); S.base!(b, ws, τ, t0); b[W])
    DB = base(wsB, TAU_DIFFUSION, mg.t0B)
    DJ = base(wsJ, TAU_DIFFUSION, mg.t0J)
    diff = (mean(mg.IB[:, 1] ./ mg.duree_somme_s), mean(mg.IJ[:, 1] ./ mg.duree_somme_s))
    signal(μ, f) = S.masquer!(μ .- f, W)
    P0 = S.phasor_brut(signal(mg.μB[1], mg.fondB[1]), cel.pas_ns)
    P1 = S.phasor_brut(signal(mg.μB[end], mg.fondB[end]), cel.pas_ns)
    irf_u(canal) = hcat(collect(0:194) .* cel.pas_ns, cel.irf[canal, 1:195] .* 1e5)
    monoB = mono ? Any[S.charger_instance_utilisateur(irf_u(BLEU)) for _ in 1:fils] : Any[]
    monoJ = mono ? Any[S.charger_instance_utilisateur(irf_u(JAUNE)) for _ in 1:fils] : Any[]
    return Contexte(mg.λ, cel.pas_ns, mg.t0B, mg.t0J, mg.τ2, mg.τ3,
        base(wsB, mg.τ2, mg.t0B), base(wsB, mg.τ3, mg.t0B), base(wsJ, mg.τJ, mg.t0J), base(wsJ, mg.τA, mg.t0J),
        base(wsJ, mg.τ3, mg.t0J), DB ./ sum(DB), DJ ./ sum(DJ), diff, cel.fond ./ cel.duree_s, gabarits(mg),
        P0, P1, Any[S.EspaceTravail(instB) for _ in 1:fils], monoB, monoJ)
end

"""Nuisances connues d'un canal pour `t_s` secondes : fond seul, ou fond + diffusion."""
function connu(ctx::Contexte, canal::Int, t_s::Real; diffusion::Bool = true)
    k = ctx.fond_par_s[canal, FENETRE] .* t_s
    diffusion && (k .+= (canal == BLEU ? ctx.DB : ctx.DJ) .* (ctx.diffusion_par_s[canal] * t_s))
    return k
end

"""
Amplitudes d'un modèle linéaire μ = X·a + k (k connu), maximum de vraisemblance de
Poisson ; `lo` : bornes basses (0 ou −Inf). Rend (a, code).
"""
function amplitudes_poisson(n::AbstractVector, X::Matrix{Float64}, k::AbstractVector, lo::Vector{Float64})
    p = size(X, 2)
    modele!(μ, J, θ) = (mul!(μ, X, θ); μ .+= k; J .= X; nothing)
    w = 1 ./ max.(n, 1.0)
    θ0 = try
        (X' * (X .* w)) \ (X' * ((n .- k) .* w))
    catch
        fill(1.0, p)
    end
    tot = max(sum(n) - sum(k), 1.0)
    for j in 1:p
        (θ0[j] > 0 && isfinite(θ0[j])) || (θ0[j] = 0.1 * tot / max(sum(view(X, :, j)), 1e-12))
    end
    r = S.lm_poisson(modele!, θ0, lo, fill(Inf, p), n; maxiter = 100, tol_dec = 1e-9)
    r.converge || return (r.θ, NON_CONVERGENCE)
    any(r.borne) && return (r.θ, BORNE)
    return (r.θ, OK)
end

"""
Mélange à deux formes : μ = A·(u + x·v) + k par canal. `mode` : `:bleu` (bleu seul),
`:separees` (une échelle par canal, 5a) ou `:commune` (une seule échelle, 5d). x ∈
[−0,5 ; 1,5] ; sur une borne, la lecture échoue. Rend (x, code).
"""
function lire_melange(F::S.FormesBanc, nB, nJ, kB, kJ; mode::Symbol)
    W = F.fenetre
    b, y = nB[W], nJ[W]
    m = length(W)
    lo, hi = -0.5, 1.5
    if mode === :bleu
        mb!(μ, J, θ) = begin
            x, A = θ[1], exp(θ[2])
            @inbounds for i in 1:m
                g = F.uB[i] + x * F.vB[i]
                μ[i] = A * g + kB[i]
                J[i, 1] = A * F.vB[i]
                J[i, 2] = A * g
            end
        end
        A0 = max(sum(b) - sum(kB), 1.0) / sum(F.uB .+ 0.5 .* F.vB)
        r = S.lm_poisson(mb!, [0.5, log(A0)], [lo, -Inf], [hi, Inf], b; maxiter = 200, tol_dec = 1e-9)
    elseif mode === :commune
        mc!(μ, J, θ) = begin
            x, A = θ[1], exp(θ[2])
            @inbounds for i in 1:m
                gB = F.uB[i] + x * F.vB[i]
                gJ = F.uJ[i] + x * F.vJ[i]
                μ[i] = A * gB + kB[i]
                μ[m+i] = A * gJ + kJ[i]
                J[i, 1] = A * F.vB[i]
                J[m+i, 1] = A * F.vJ[i]
                J[i, 2] = A * gB
                J[m+i, 2] = A * gJ
            end
        end
        A0 = max(sum(b) + sum(y) - sum(kB) - sum(kJ), 1.0) / sum(F.uB .+ F.uJ .+ 0.5 .* (F.vB .+ F.vJ))
        r = S.lm_poisson(mc!, [0.5, log(A0)], [lo, -Inf], [hi, Inf], vcat(b, y); maxiter = 200, tol_dec = 1e-9)
    else
        ms!(μ, J, θ) = begin
            x, A, AJ = θ[1], exp(θ[2]), exp(θ[3])
            fill!(J, 0.0)
            @inbounds for i in 1:m
                gB = F.uB[i] + x * F.vB[i]
                gJ = F.uJ[i] + x * F.vJ[i]
                μ[i] = A * gB + kB[i]
                μ[m+i] = AJ * gJ + kJ[i]
                J[i, 1] = A * F.vB[i]
                J[m+i, 1] = AJ * F.vJ[i]
                J[i, 2] = A * gB
                J[m+i, 3] = AJ * gJ
            end
        end
        A0 = max(sum(b) - sum(kB), 1.0) / sum(F.uB .+ 0.5 .* F.vB)
        AJ0 = max(sum(y) - sum(kJ), 1.0) / sum(F.uJ .+ 0.5 .* F.vJ)
        r = S.lm_poisson(ms!, [0.5, log(A0), log(AJ0)], [lo, -Inf, -Inf], [hi, Inf, Inf], vcat(b, y); maxiter = 200, tol_dec = 1e-9)
    end
    r.converge || return (r.θ[1], NON_CONVERGENCE)
    r.borne[1] && return (r.θ[1], BORNE)
    return (r.θ[1], OK)
end

"""
Bi-exponentiel à τ libres sur le bleu (lecture 3) : μ = A₁B(τ₁) + A₂B(τ₂) + k, t₀ figé,
fond et diffusion connus. Rend la fraction d'amplitude de la composante courte. Échec :
non-convergence, τ ou amplitude sur une borne, composantes confondues (moins de 2 %).
"""
function lire_biexp_libre(nB, ws, ctx::Contexte, k::AbstractVector; bornes = (0.02, 20.0))
    W = FENETRE
    n = nB[W]
    sum(n) > 0 || return (NaN, AUTRE)
    N = ws.n
    B1, d1, B2, d2, poub = zeros(N), zeros(N), zeros(N), zeros(N), zeros(N)
    modele!(μ, J, θ) = begin
        S.bases_derivees!(B1, d1, poub, ws, θ[1], ctx.t0B)
        S.bases_derivees!(B2, d2, poub, ws, θ[2], ctx.t0B)
        @inbounds for (i, c) in enumerate(W)
            μ[i] = θ[3] * B1[c] + θ[4] * B2[c] + k[i]
            J[i, 1] = θ[3] * d1[c]
            J[i, 2] = θ[4] * d2[c]
            J[i, 3] = B1[c]
            J[i, 4] = B2[c]
        end
    end
    tot = max(sum(n) - sum(k), 1.0)
    θ0 = [ctx.τ2, ctx.τ3, 0.5tot / sum(ctx.B2), 0.5tot / sum(ctx.B3)]
    r = S.lm_poisson(modele!, θ0, [bornes[1], bornes[1], 0.0, 0.0], [bornes[2], bornes[2], Inf, Inf], n; maxiter = 200, tol_dec = 1e-8)
    τ1, τ2, A1, A2 = r.θ
    r.converge || return (NaN, NON_CONVERGENCE)
    any(r.borne) && return (NaN, BORNE)
    abs(τ1 - τ2) / max(τ1, τ2) < 0.02 && return (NaN, AUTRE)
    return ((τ1 < τ2 ? A1 : A2) / (A1 + A2), OK)
end

"""Ton ajustement mono-exponentiel (lectures 2 et 2J) : τ, ou un code d'échec."""
function lire_mono(inst, h::Vector{Float64}, pas::Real)
    sum(@view h[FENETRE]) > 0 || return (NaN, AUTRE)
    τ = S.ajuster_mono_banc(inst, copy(h); pas)
    isnan(τ) && return (NaN, NON_CONVERGENCE)
    (τ <= 2 * pas * 1.001 || τ >= 2 * 12.5 * 0.999) && return (τ, BORNE)
    return (τ, OK)
end

"""
    lire(ctx, hB, hJ, t_s, fil) -> (o, codes)

Les neuf lectures (ordre de `LECTURES`) d'un couple d'histogrammes de `t_s` secondes.
`fil` choisit l'espace de travail et les instances de ton code (un par tâche parallèle).
"""
function lire(ctx::Contexte, hB::AbstractVector, hJ::AbstractVector, t_s::Real, fil::Int)
    W = FENETRE
    o = fill(NaN, length(LECTURES))
    e = fill(Int8(OK), length(LECTURES))
    sB, sJ = sum(@view hB[W]), sum(@view hJ[W])
    # 1 : rapport J/B brut
    if sB > 0
        o[1] = sJ / sB
    else
        e[1] = AUTRE
    end
    # 1c : J/B moins la fuite estimée sur l'histogramme lui-même (composantes rapides)
    kB, kJ = connu(ctx, BLEU, t_s), connu(ctx, JAUNE, t_s)
    if sB > 0
        aB, cB = amplitudes_poisson(hB[W], hcat(ctx.B2, ctx.B3), kB, [0.0, 0.0])
        aJ, cJ = amplitudes_poisson(hJ[W], hcat(ctx.JJ, ctx.JA), kJ, [-Inf, 0.0])
        if cB == OK && cJ == OK
            o[2] = sJ / sB - (aJ[1] * sum(ctx.JJ)) / (aB[1] * sum(ctx.B2))
        else
            e[2] = max(cB, cJ)
        end
    else
        e[2] = AUTRE
    end
    # 2, 2J : ton code
    if isempty(ctx.monoB)
        e[3] = e[4] = AUTRE
    else
        o[3], e[3] = lire_mono(ctx.monoB[fil], Vector{Float64}(hB), ctx.pas)
        o[4], e[4] = lire_mono(ctx.monoJ[fil], Vector{Float64}(hJ), ctx.pas)
    end
    # 3 : bi-exponentiel à τ libres
    o[5], e[5] = lire_biexp_libre(hB, ctx.ws[fil], ctx, kB)
    # 4 : bi-exponentiel à τ figés (τ₂, τ₃) : fraction d'amplitude de τ₂ (population)
    F4 = S.FormesBanc(ctx.B3, ctx.B2 .- ctx.B3, Float64[], Float64[], W)
    o[6], e[6] = sB > 0 ? lire_melange(F4, hB, hJ, kB, kJ; mode = :bleu) : (NaN, AUTRE)
    # 5a, 5d : mélange des deux états extrêmes (gabarits mesurés ; la diffusion y est déjà)
    kB0, kJ0 = connu(ctx, BLEU, t_s; diffusion = false), connu(ctx, JAUNE, t_s; diffusion = false)
    o[7], e[7] = sB > 0 ? lire_melange(ctx.G, hB, hJ, kB0, kJ0; mode = :separees) : (NaN, AUTRE)
    o[8], e[8] = sB > 0 ? lire_melange(ctx.G, hB, hJ, kB0, kJ0; mode = :commune) : (NaN, AUTRE)
    # 6 : phasor brut du bleu projeté sur l'axe des deux états extrêmes
    if sB > 0
        o[9] = S.lecture_phasor(hB, ctx.P0, ctx.P1, ctx.pas)
    else
        e[9] = AUTRE
    end
    for i in eachindex(o)
        e[i] == OK || (o[i] = NaN)
    end
    return o, e
end

"""
Histogramme décalé de δ ns (plus tard si δ > 0), par interpolation linéaire entre canaux ;
zéro hors de la fenêtre. Sert à mesurer la sensibilité d'une lecture à un saut de t₀.
"""
function decaler(h::AbstractVector, δ::Real, pas::Real)
    s = δ / pas
    n = length(h)
    out = zeros(n)
    for k in 1:n
        x = k - s
        i = floor(Int, x)
        f = x - i
        a = 1 <= i <= n ? Float64(h[i]) : 0.0
        b = 1 <= i + 1 <= n ? Float64(h[i+1]) : 0.0
        out[k] = (1 - f) * a + f * b
    end
    out[setdiff(1:n, FENETRE)] .= 0.0
    return out
end
