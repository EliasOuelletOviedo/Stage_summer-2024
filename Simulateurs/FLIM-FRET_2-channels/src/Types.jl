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
    Φ_A::Float64 = 0.60      # rendement quantique de l'accepteur, valeur de catalogue
    Φ_A_eff::Float64 = NaN   # rendement effectif dans le capteur (mesuré au banc) ; NaN : Φ_A
    mecanisme::Symbol = :M2_lent
    # Blanchiment pendant l'acquisition (changement 5) : fractions d'accepteurs (β_A) et de
    # donneurs (β_D) blanchies à la fin de l'acquisition, décroissance exponentielle de chacun ;
    # un donneur sans accepteur a la durée de vie τ_nu (mCerulean3 seul, FPbase).
    β_A::Float64 = 0.0
    β_D::Float64 = 0.0
    τ_nu::Float64 = 4.1
    # Produit du blanchiment de l'accepteur : false, il n'absorbe plus (le donneur devient nu,
    # τ_nu) ; true, il absorbe encore mais n'émet plus (extincteur sombre : le bleu est inchangé,
    # seule l'émission sensibilisée disparaît — ce que suggèrent les tranches : ratio instable,
    # durée de vie stable sous illumination répétée).
    blanchi_sombre::Bool = false
end

"""Rendement quantique de l'accepteur utilisé par le modèle : l'effectif s'il est connu."""
phi_A(p) = isnan(p.Φ_A_eff) ? p.Φ_A : p.Φ_A_eff

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
`Δt_jaune` : décalage temporel du détecteur jaune par rapport au bleu (0 = banc calibré).
"""
Base.@kwdef struct InstrumentParams
    n_bins::Int = 256
    T_rep::Float64 = 12.5
    irf::Vector{Float64}
    t₀::Float64 = 0.0
    fenetre::UnitRange{Int} = 14:242
    Δt_jaune::Float64 = 0.0   # ns, retard du canal jaune sur le bleu (changement 5)
    # Pas d'un canal (ns). NaN : la période remplit exactement les n_bins canaux (T_rep/n_bins,
    # titration du 2026-07-20). Sinon (banc du 2026-10-16 : TAC de 16,4 ns sur 256 canaux,
    # 64 ps, une période = 195,3 canaux), la convolution avec l'IRF est linéaire.
    pas::Float64 = NaN
end

pas_temps(inst::InstrumentParams) = isnan(inst.pas) ? inst.T_rep / inst.n_bins : inst.pas

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
