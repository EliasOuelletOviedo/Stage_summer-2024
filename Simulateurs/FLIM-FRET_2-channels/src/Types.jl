# Types.jl — les structures de base du plan de conception.
#
# Tout canal se ramène à une somme d'exponentielles (ExpSum) ; le capteur,
# l'optique et l'instrument sont séparés pour qu'une valeur mesurée remplace
# une valeur inférée sans toucher au reste.

"""
    ExpSum(α, τ)

Un canal entièrement décrit : `Σᵢ αᵢ e^{-t/τᵢ}`, en photons/ns par excitation
du donneur. Certaines amplitudes sont négatives (la montée de l'émission
sensibilisée).
"""
struct ExpSum{T<:Real}
    α::Vector{T}
    τ::Vector{T}   # ns
end

"""
    SensorParams

Photophysique du capteur. Les trois durées de vie sont des constantes ; tout
l'effet du Cl⁻ passe par la fraction active x (isotherme à deux plateaux).
`K_D` est le K vrai, jamais un K apparent mesuré par un estimateur.
`mecanisme` : `:M2_lent` (mélange statique, défaut) ou `:M2_rapide` (échange
rapide, un seul τ du donneur qui glisse avec x).
"""
Base.@kwdef struct SensorParams
    τ_D0::Float64            # ns, donneur sans accepteur actif
    τ_DA::Float64            # ns, donneur avec accepteur actif, = τ_D0(1−E)
    τ_A::Float64 = 3.0       # ns, accepteur déprotoné
    x₀::Float64              # fraction active à [Cl⁻] = 0
    x_sat::Float64           # fraction active à [Cl⁻] saturant
    K_D::Float64             # mM
    n::Float64 = 1.0         # coefficient de Hill de l'isotherme
    Φ_D::Float64 = 0.87      # rendement quantique du donneur seul
    Φ_A::Float64 = 0.60      # rendement quantique de l'accepteur
    mecanisme::Symbol = :M2_lent
end

"""
    OpticsParams

`η[canal, émetteur]` : fraction des photons de chaque émetteur qui tombe dans
chaque canal (ligne 1 = bleu, ligne 2 = jaune ; colonne 1 = donneur,
colonne 2 = accepteur). `ε` : excitation directe de l'accepteur déprotoné,
relative à celle du donneur. `εΦ_dark` : produit ε_dark·Φ_dark de la forme
protonée de l'accepteur (risque à 800 nm, nul par défaut).
"""
Base.@kwdef struct OpticsParams
    η::SMatrix{2,2,Float64,4}
    ε::Float64 = 0.03
    εΦ_dark::Float64 = 0.0
end

"""
    InstrumentParams

`irf` : réponse instrumentale sur la grille des `n_bins` canaux, normalisée à
somme 1. `t₀` : décalage de l'IRF (ns), paramètre du modèle même s'il vaut 0.
`fenetre` : canaux valides du TAC (les autres ne reçoivent aucun photon).
"""
Base.@kwdef struct InstrumentParams
    n_bins::Int = 256
    T_rep::Float64 = 12.5
    irf::Vector{Float64}
    t₀::Float64 = 0.0
    fenetre::UnitRange{Int} = 14:242
end

pas_temps(inst::InstrumentParams) = inst.T_rep / inst.n_bins

"""Un monde simulé complet : capteur, optique, instrument."""
struct Modele
    p::SensorParams
    o::OpticsParams
    inst::InstrumentParams
end

"""
    remplacer(s; champs...)

Copie de la structure `s` avec quelques champs changés (balayages de nuisance,
valeurs volontairement fausses).
"""
function remplacer(s::T; kw...) where {T}
    T(; (f => get(kw, f, getfield(s, f)) for f in fieldnames(T))...)
end

remplacer(m::Modele; p = m.p, o = m.o, inst = m.inst) = Modele(p, o, inst)
