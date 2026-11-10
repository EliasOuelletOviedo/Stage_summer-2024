# Estimators.jl — les estimateurs, interface commune.
#
#   observable(e, nB, nY, ctx) -> (o, statut)
#
# `o` est l'observable propre à l'estimateur (ratio, τ apparent, fraction
# d'amplitude, x, coordonnée phasor) ; la concentration vient ensuite de SA
# courbe de calibration (Calibrate.jl). `statut` : :ok, :borne (valeur
# ramenée à une borne de recherche, gardée) ou :echec (non convergé, ou
# convergé sur une frontière d'un ajustement non linéaire : exclu des
# statistiques et compté dans le taux d'échec).
#
# Ce qui distingue ces estimateurs n'est pas la formule mais le nombre de
# paramètres libres et les données qu'ils consomment. Les estimateurs à τ
# figés sont construits à partir d'un modèle « cru » (les valeurs nominales de
# l'expérimentateur), qui peut différer du modèle vrai (balayages 3–5).

abstract type Estimateur end

"""Contexte d'un fil d'exécution : espace FFT et instance de ton code."""
struct Contexte
    ws::EspaceTravail
    utilisateur::Union{Nothing,InstanceUtilisateur}
end

# -----------------------------------------------------------------------------
# 1 — Ratio d'intensités (0 paramètre libre)
# -----------------------------------------------------------------------------

struct Ratio <: Estimateur end
nom(::Ratio) = "1 ratio J/B"
function observable(::Ratio, nB, nY, ctx)
    NB = sum(nB)
    NB > 0 || return (NaN, :echec)
    return (sum(nY) / NB, :ok)
end

# -----------------------------------------------------------------------------
# 2 — Mono-exponentiel MLE : ton code (2 paramètres utiles : τ, amplitude ;
#     ton modèle libère aussi le décalage d'IRF et un fond)
# -----------------------------------------------------------------------------

struct MonoUtilisateur <: Estimateur
    tentatives::Int
    bornes::Tuple{Float64,Float64}
    fixes::Union{Nothing,NTuple{2,Float64}}   # (décalage, fond) figés, ou libres (nothing)
end
MonoUtilisateur(t, b) = MonoUtilisateur(t, b, nothing)
nom(e::MonoUtilisateur) = e.fixes === nothing ? "2 mono (ton code)" : "2F mono, décalage et fond figés"

"""
2F : ton ajustement avec le décalage d'IRF et le fond figés à leurs valeurs sans
bruit à 20 mM (comme si on les avait mesurés une fois) ; seul τ reste libre.
"""
function MonoUtilisateurFige(m::Modele, ctx, irf_mat; C = C_REF_NORMALISATION_A)
    μB, _ = attendus(C, m, 1e6, :C, ctx.ws)
    p, _ = ajuster_mono_utilisateur(ctx.utilisateur, μB)
    return MonoUtilisateur(5, bornes_tau_utilisateur(irf_mat), (p[2], p[3]))
end

function observable(e::MonoUtilisateur, nB, nY, ctx)
    p, _ = ajuster_mono_utilisateur(ctx.utilisateur, nB; tentatives = e.tentatives, fixes = e.fixes)
    isnan(p[1]) && return (NaN, :echec)
    τ = p[1]
    (τ <= e.bornes[1] * (1 + 1e-6) || τ >= e.bornes[2] * (1 - 1e-6)) && return (τ, :echec)
    return (τ, :ok)
end

# -----------------------------------------------------------------------------
# 3 — Bi-exponentiel MLE, τ libres (4 paramètres : τ₁, τ₂, A₁, A₂)
# -----------------------------------------------------------------------------

struct BiExpLibre <: Estimateur
    τ_init::NTuple{2,Float64}     # (court, long) : valeurs nominales comme point de départ
    t₀::Float64
    fenetre::UnitRange{Int}
    bornes_τ::NTuple{2,Float64}
end
nom(::BiExpLibre) = "3 biexp τ libres"

function observable(e::BiExpLibre, nB, nY, ctx)
    ws = ctx.ws
    W = e.fenetre
    n = nB[W]
    sum(n) > 0 || return (NaN, :echec)
    nf = ws.n
    B1, d1, B2, d2, poub = zeros(nf), zeros(nf), zeros(nf), zeros(nf), zeros(nf)
    function modele!(μ, J, θ)
        bases_derivees!(B1, d1, poub, ws, θ[1], e.t₀)
        bases_derivees!(B2, d2, poub, ws, θ[2], e.t₀)
        @inbounds for (i, k) in enumerate(W)
            μ[i] = θ[3] * B1[k] + θ[4] * B2[k]
            J[i, 1] = θ[3] * d1[k]
            J[i, 2] = θ[4] * d2[k]
            J[i, 3] = B1[k]
            J[i, 4] = B2[k]
        end
    end
    # départ : τ nominaux, amplitudes par l'ajustement linéaire à τ figés
    base!(B1, ws, e.τ_init[1], e.t₀)
    base!(B2, ws, e.τ_init[2], e.t₀)
    f0, st0, _ = profil_lineaire([CanalLineaire(n, B2[W], B1[W] .- B2[W])], 0.0, 1.0)
    f0 = st0 === :echec ? 0.5 : clamp(f0, 0.02, 0.98)
    A0 = sum(n) / sum(f0 .* B1[W] .+ (1 - f0) .* B2[W])
    θ0 = [e.τ_init[1], e.τ_init[2], f0 * A0, (1 - f0) * A0]
    lo = [e.bornes_τ[1], e.bornes_τ[1], 0.0, 0.0]
    hi = [e.bornes_τ[2], e.bornes_τ[2], Inf, Inf]
    r = lm_poisson(modele!, θ0, lo, hi, n; maxiter = 200, tol_dec = 1e-8)
    r.converge || return (NaN, :echec)
    τ1, τ2, A1, A2 = r.θ
    (any(r.borne[1:2])) && return (NaN, :echec)
    abs(τ1 - τ2) / max(τ1, τ2) < 0.02 && return (NaN, :echec)    # composantes confondues
    Acourt = τ1 < τ2 ? A1 : A2
    return (Acourt / (A1 + A2), :ok)
end

# -----------------------------------------------------------------------------
# 4 — Bi-exponentiel MLE, τ figés (2 paramètres : A₁, A₂ ↔ fraction, échelle)
# -----------------------------------------------------------------------------

struct BiExpFige <: Estimateur
    u::Vector{Float64}       # B(τ_D0) sur la fenêtre
    v::Vector{Float64}       # B(τ_DA) − B(τ_D0)
    fenetre::UnitRange{Int}
end
nom(::BiExpFige) = "4 biexp τ figés"

function BiExpFige(m::Modele, ws::EspaceTravail)
    W = m.inst.fenetre
    BDA, BD0 = zeros(ws.n), zeros(ws.n)
    base!(BDA, ws, m.p.τ_DA, m.inst.t₀)
    base!(BD0, ws, m.p.τ_D0, m.inst.t₀)
    return BiExpFige(BD0[W], BDA[W] .- BD0[W], W)
end

function observable(e::BiExpFige, nB, nY, ctx)
    n = nB[e.fenetre]
    sum(n) > 0 || return (NaN, :echec)
    f, st, _ = profil_lineaire([CanalLineaire(n, e.u, e.v)], -0.5, 1.5)
    st === :echec && return (NaN, :echec)
    return (f, st)
end

# -----------------------------------------------------------------------------
# 5 — Conjoint deux canaux, τ figés : x directement
#   5a : échelles bleue et jaune libres et séparées (3 paramètres : x, S_B, S_J)
#   5b : idem, τ figés à des valeurs fausses (modèle cru ≠ modèle vrai)
#   5d : échelle commune (2 paramètres : x, S) — utilise aussi le rapport J/B
# -----------------------------------------------------------------------------

struct Conjoint <: Estimateur
    uB::Vector{Float64}
    vB::Vector{Float64}
    uY::Vector{Float64}
    vY::Vector{Float64}
    fenetre::UnitRange{Int}
    commune::Bool
    etiquette::String
end
nom(e::Conjoint) = e.etiquette

"""Formes des deux canaux affines en x : m(x) = u + x·v (exact sous M2 lent)."""
function Conjoint(m::Modele, ws::EspaceTravail; commune::Bool = false, etiquette = commune ? "5d conjoint échelle commune" : "5a conjoint τ figés")
    m.p.mecanisme === :M2_lent || error("l'estimateur conjoint suppose M2 lent")
    W = m.inst.fenetre
    mB0, mY0 = formes_x(0.0, m, ws)
    mB1, mY1 = formes_x(1.0, m, ws)
    return Conjoint(mB0[W], mB1[W] .- mB0[W], mY0[W], mY1[W] .- mY0[W], W, commune, etiquette)
end

function observable(e::Conjoint, nB, nY, ctx)
    W = e.fenetre
    b, y = nB[W], nY[W]
    sum(b) + sum(y) > 0 || return (NaN, :echec)
    canaux = e.commune ? [CanalLineaire(vcat(b, y), vcat(e.uB, e.uY), vcat(e.vB, e.vY))] :
             [CanalLineaire(b, e.uB, e.vB), CanalLineaire(y, e.uY, e.vY)]
    x, st, _ = profil_lineaire(canaux, -0.5, 1.5)
    st === :echec && return (NaN, :echec)
    return (x, st)
end

# -----------------------------------------------------------------------------
# 5c — Conjoint, τ_A libre (4 paramètres : x, τ_A, S_B, S_J)
# -----------------------------------------------------------------------------

struct ConjointTauALibre <: Estimateur
    m::Modele                 # modèle cru (τ_D0, τ_DA figés ; τ_A de départ)
    bornes_τA::NTuple{2,Float64}
    depart::Conjoint          # 5a pour le point de départ en x
end
nom(::ConjointTauALibre) = "5c conjoint τ_A libre"

function ConjointTauALibre(m::Modele, ws::EspaceTravail)
    return ConjointTauALibre(m, (1.05 * m.p.τ_DA, 10.0), Conjoint(m, ws))
end

function observable(e::ConjointTauALibre, nB, nY, ctx)
    W = e.m.inst.fenetre
    b, y = nB[W], nY[W]
    (sum(b) > 0 && sum(y) > 0) || return (NaN, :echec)
    x0, st0 = observable(e.depart, nB, nY, ctx)
    x0 = st0 === :echec ? 0.3 : clamp(x0, -0.2, 1.2)
    mW = length(W)
    n = vcat(b, y)
    function modele!(μ, J, θ)
        x, τA, lSB, lSY = θ
        mm = remplacer(e.m; p = remplacer(e.m.p; τ_A = τA))
        mB, mY, dB, dY = formes_et_derivees(x, mm, ctx.ws)
        SB, SY = exp(lSB), exp(lSY)
        @inbounds for (i, k) in enumerate(W)
            μ[i] = SB * mB[k]
            μ[mW+i] = SY * mY[k]
            J[i, 1] = SB * dB[k, 1]
            J[i, 2] = SB * dB[k, 4]
            J[i, 3] = μ[i]
            J[i, 4] = 0.0
            J[mW+i, 1] = SY * dY[k, 1]
            J[mW+i, 2] = SY * dY[k, 4]
            J[mW+i, 3] = 0.0
            J[mW+i, 4] = μ[mW+i]
        end
    end
    mB, mY = formes_x(x0, e.m, ctx.ws)
    θ0 = [x0, e.m.p.τ_A, log(sum(b) / sum(mB[W])), log(sum(y) / sum(mY[W]))]
    lo = [-0.5, e.bornes_τA[1], -Inf, -Inf]
    hi = [1.5, e.bornes_τA[2], Inf, Inf]
    r = lm_poisson(modele!, θ0, lo, hi, n; maxiter = 200, tol_dec = 1e-8)
    (r.converge && !r.borne[2]) || return (NaN, :echec)
    return (r.θ[1], r.borne[1] ? :borne : :ok)
end

# -----------------------------------------------------------------------------
# 6 — Phasor (canal bleu, 0 paramètre libre)
# -----------------------------------------------------------------------------

struct Phasor <: Estimateur
    P0::ComplexF64        # phasor du modèle cru à [Cl⁻] = 0
    P1::ComplexF64        # à [Cl⁻] = C_max : l'axe de projection
    inst::InstrumentParams
end
nom(::Phasor) = "6 phasor"

function Phasor(m::Modele, ws::EspaceTravail, C_max::Real)
    mB0, _ = formes(0.0, m, ws)
    mB1, _ = formes(C_max, m, ws)
    return Phasor(phasor(mB0, m.inst), phasor(mB1, m.inst), m.inst)
end

"""Coordonnée le long de l'axe P0 → P1 (0 à [Cl⁻] = 0, 1 à C_max) : sous M2 lent, la titration est un segment."""
function observable(e::Phasor, nB, nY, ctx)
    sum(nB[e.inst.fenetre]) > 0 || return (NaN, :echec)
    z = phasor(nB, e.inst)
    d = e.P1 - e.P0
    return (real((z - e.P0) * conj(d)) / abs2(d), :ok)
end

# -----------------------------------------------------------------------------
# Construction de la liste standard
# -----------------------------------------------------------------------------

"""
    estimateurs_standard(m_cru, ws, irf_mat; C_max, liste)

Les estimateurs construits sur le modèle cru `m_cru` (ce que l'expérimentateur
croit savoir : τ, η, ε). Clés : "1","2","2H","3","4","5a","5c","5d","6".
"""
function estimateurs_standard(m::Modele, ws::EspaceTravail, irf_mat::Matrix{Float64}; C_max = 150.0,
                              liste = ["1", "2", "2H", "2F", "3", "4", "5a", "5c", "5d", "6"], ctx = nothing)
    ctx === nothing && (liste = filter(!=("2F"), liste))
    mono = MonoUtilisateur(5, bornes_tau_utilisateur(irf_mat))   # 2 et 2H : même ajustement, deux calibrations
    d = Dict{String,Estimateur}(
        "1" => Ratio(),
        "2" => mono,
        "2H" => mono,
        "3" => BiExpLibre((m.p.τ_DA, m.p.τ_D0), m.inst.t₀, m.inst.fenetre, (0.02, 20.0)),
        "4" => BiExpFige(m, ws),
        "5a" => Conjoint(m, ws),
        "5b" => Conjoint(m, ws; etiquette = "5b conjoint τ faux"),
        "5c" => ConjointTauALibre(m, ws),
        "5d" => Conjoint(m, ws; commune = true),
        "6" => Phasor(m, ws, C_max),
    )
    "2F" in liste && (d["2F"] = MonoUtilisateurFige(m, ctx, irf_mat))
    return [(k, d[k]) for k in liste]
end

"""Noms affichés (sorties, figures)."""
const NOMS = Dict("1" => "1 ratio J/B", "2" => "2 mono (ton code)", "2H" => "2H mono + Hill (ton pipeline)",
    "2F" => "2F mono, décalage et fond figés",
    "3" => "3 biexp τ libres", "4" => "4 biexp τ figés", "5a" => "5a conjoint τ figés",
    "5b" => "5b conjoint τ faux", "5c" => "5c conjoint τ_A libre", "5d" => "5d conjoint échelle commune",
    "6" => "6 phasor")

"""Borne de Cramér-Rao « propre » de chaque estimateur (l'information des données qu'il consomme)."""
const BORNE_PROPRE = Dict("1" => :totaux, "2" => :bleu, "2H" => :bleu, "2F" => :bleu, "3" => :bleu_tau_libres, "4" => :bleu,
    "5a" => :conjoint, "5b" => :conjoint, "5c" => :conjoint_tauA_libre, "5d" => :absolue, "6" => :bleu)

"""Les observables qui ne dépendent que des données (pas du modèle cru) : mêmes valeurs quel que soit le modèle cru."""
const INDEPENDANTS_DU_MODELE = Set(["1", "2", "2H", "3"])   # 2F dépend du modèle (valeurs figées)
