# params.jl — lecture de params.toml : Langmuir, durées de vie, IRF, fond, réglages.
#
# Aucun nombre du modèle n'est écrit ailleurs que dans params.toml. Le bloc « step0 »
# (titration du 2026-10-16, 800 nm) est figé ; les réglages de la course (budget, dérive,
# filtre) se changent par `avec(p; ...)` sans relire le fichier.

using TOML, FFTW, LinearAlgebra

"""Langmuir y(C) = y₀ + (y∞ − y₀)·C/(K + C), ajustée par step0 sur les 8 points de la titration."""
struct Langmuir
    y0::Float64
    yinf::Float64
    K::Float64
end
(l::Langmuir)(c) = l.y0 + (l.yinf - l.y0) * c / (l.K + c)

struct Params
    toml::Dict{String,Any}
    pas::Float64                  # ns par canal TCSPC
    T::Float64                    # période laser (ns)
    n::Int
    fenetre::UnitRange{Int}
    t0::NTuple{2,Float64}         # décalage déclin/IRF, bleu et jaune (ns)
    irf_b::Vector{Float64}
    irf_j::Vector{Float64}
    τ::NamedTuple{(:diff, :τ2, :τ3, :τJ, :τA),NTuple{5,Float64}}
    lang::Dict{String,Langmuir}   # JB, d_bleu, p_bleu, d_jaune, q_jaune
    dt::Float64
    fond_b::Vector{Float64}       # photons de fond attendus par canal pour une acquisition de dt
    fond_j::Vector{Float64}
    N::Float64                    # photons bleus attendus par acquisition (fond exclu)
    derive::NamedTuple{(:forme, :alpha, :saut, :t_saut),Tuple{String,Float64,Float64,Float64}}
    plan::Any                     # rfft de longueur 2n : la période (12,5 ns) ne remplit pas les 16,4 ns du TAC
end

function charger_params(chemin = joinpath(@__DIR__, "..", "params.toml"))
    d = TOML.parsefile(chemin)
    haskey(d, "instrument") || error("params.toml sans bloc step0 : lancer step0_params.jl")
    ins, tv, tr = d["instrument"], d["durees_de_vie"], d["trajectoire"]
    n = ins["n_canaux"]
    fen = ins["fenetre"][1]:ins["fenetre"][2]
    dt = tr["dt_s"]
    fond(cps) = [k in fen ? cps * dt : 0.0 for k in 1:n]
    lang = Dict(k => Langmuir(v["y0"], v["yinf"], v["K_mM"]) for (k, v) in d["langmuir"])
    dv = d["derive"]
    return Params(d, ins["pas_ns"], ins["periode_ns"], n, fen, (ins["t0_bleu_ns"], ins["t0_jaune_ns"]),
        Float64.(ins["irf_bleu"]), Float64.(ins["irf_jaune"]),
        (diff = tv["tau_diffusion_ns"], τ2 = tv["tau2_ns"], τ3 = tv["tau3_ns"], τJ = tv["tauJ_ns"], τA = tv["tauA_ns"]),
        lang, dt, fond(ins["fond_bleu_cps_par_canal"]), fond(ins["fond_jaune_cps_par_canal"]),
        Float64(d["photons"]["N_bleu"]),
        (forme = dv["forme"], alpha = Float64(dv["alpha"]), saut = Float64(dv["saut"]), t_saut = Float64(dv["t_saut_s"])),
        plan_rfft(zeros(2n)))
end

"""Copie de `p` avec un autre budget ou une autre dérive (les scénarios d'une même course)."""
function avec(p::Params; N = p.N, forme = p.derive.forme, alpha = p.derive.alpha)
    d = (forme = forme, alpha = Float64(alpha), saut = p.derive.saut, t_saut = p.derive.t_saut)
    return Params((f === :N ? Float64(N) : f === :derive ? d : getfield(p, f) for f in fieldnames(Params))...)
end

reglage(p::Params, section, cle) = p.toml[section][cle]
jb(c, p::Params) = p.lang["JB"](c)

"""
Courbe monotone par morceaux (Hermite cubique de Fritsch–Carlson) : valeur, dérivée,
inverse. Hors des nœuds, prolongement linéaire pour le filtre (`h`, `dh`) ; l'inversion,
elle, refuse de sortir de la plage.
"""
struct Courbe
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
end

function Courbe(x::Vector{Float64}, y::Vector{Float64})
    s = diff(y) ./ diff(x)
    all(>(0), s) || all(<(0), s) || error("courbe de calibration non monotone")
    m = [s[1]; [s[i-1] * s[i] <= 0 ? 0.0 : 2 / (1 / s[i-1] + 1 / s[i]) for i in 2:length(s)]; s[end]]
    return Courbe(x, y, m)
end

function hermite(c::Courbe, x::Real)
    x <= c.x[1] && return c.y[1] + c.m[1] * (x - c.x[1]), c.m[1]
    x >= c.x[end] && return c.y[end] + c.m[end] * (x - c.x[end]), c.m[end]
    i = clamp(searchsortedlast(c.x, x), 1, length(c.x) - 1)
    h = c.x[i+1] - c.x[i]
    u = (x - c.x[i]) / h
    y0, y1, m0, m1 = c.y[i], c.y[i+1], c.m[i] * h, c.m[i+1] * h
    v = (2u^3 - 3u^2 + 1) * y0 + (u^3 - 2u^2 + u) * m0 + (-2u^3 + 3u^2) * y1 + (u^3 - u^2) * m1
    dv = ((6u^2 - 6u) * y0 + (3u^2 - 4u + 1) * m0 + (-6u^2 + 6u) * y1 + (3u^2 - 2u) * m1) / h
    return v, dv
end
valeur(c::Courbe, x) = hermite(c, x)[1]
derivee(c::Courbe, x) = hermite(c, x)[2]

"""C tel que la courbe vaille z ; NaN hors de la plage calibrée (un échec, pas une extrapolation)."""
function inverser(c::Courbe, z::Real)
    lo, hi = extrema(c.y)
    (isfinite(z) && lo <= z <= hi) || return NaN
    croissante = c.y[end] > c.y[1]
    a, b = c.x[1], c.x[end]
    for _ in 1:60
        m = (a + b) / 2
        (valeur(c, m) < z) == croissante ? (a = m) : (b = m)
    end
    return (a + b) / 2
end
