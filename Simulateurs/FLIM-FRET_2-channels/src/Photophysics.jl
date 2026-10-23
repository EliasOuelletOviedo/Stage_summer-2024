# Photophysics.jl — de [Cl⁻] aux sommes d'exponentielles des deux canaux.
# Aucune notion d'instrument ici.

"""
    fraction_active(C, p)

Isotherme à deux plateaux : x(C) = x_sat + (x₀ − x_sat)·K^n/(K^n + C^n).
`C = Inf` donne x_sat.
"""
function fraction_active(C::Real, p::SensorParams)
    isinf(C) && return p.x_sat
    Kn = p.K_D^p.n
    return p.x_sat + (p.x₀ - p.x_sat) * Kn / (Kn + C^p.n)
end

"""dx/dC de l'isotherme (mM⁻¹) ; c'est lui qui crée la zone de faible sensibilité."""
function dfraction_dC(C::Real, p::SensorParams)
    (isinf(C) || (C == 0 && p.n > 1)) && return 0.0
    C == 0 && p.n == 1 && return -(p.x₀ - p.x_sat) / p.K_D
    Kn = p.K_D^p.n
    Cn = C^p.n
    return -(p.x₀ - p.x_sat) * Kn * p.n * C^(p.n - 1) / (Kn + Cn)^2
end

"""Inverse analytique de l'isotherme (C ≥ 0) ; `nothing` hors de ]x_sat, x₀]."""
function concentration_de_x(x::Real, p::SensorParams)
    lo, hi = minmax(p.x_sat, p.x₀)
    (lo < x <= hi) || return nothing
    r = (p.x₀ - x) / (x - p.x_sat)
    r <= 0 && return 0.0
    return p.K_D * r^(1 / p.n)
end

# -----------------------------------------------------------------------------
# M2 lent : trois durées de vie constantes, amplitudes affines en x
# -----------------------------------------------------------------------------

"""
    amplitudes_M2_lent(x, τ_D0, τ_DA, τ_A, Φ_D, Φ_A, η, ε, εΦ_dark)

Amplitudes (photons/ns par excitation du donneur) des deux canaux sur les
trois durées de vie `(τ_DA, τ_D0, τ_A)` : matrice 2×3 (canal × terme).

- donneur : (Φ_D/τ_D0)[x e^{-t/τ_DA} + (1−x) e^{-t/τ_D0}]
- émission sensibilisée : s(e^{-t/τ_A} − e^{-t/τ_DA}),
  s = (Φ_A/τ_A)·x·k_t/(k_DA − k_A), k_t = 1/τ_DA − 1/τ_D0
- excitation directe : ε·x·(Φ_A/τ_A) e^{-t/τ_A}
- forme protonée (risque) : εΦ_dark·(1−x)/τ_A e^{-t/τ_A} ; on suppose qu'elle
  émet comme l'anion (transfert de proton à l'état excité), donc avec τ_A et
  le spectre de l'accepteur.

Générique en type pour la différentiation automatique (test 7).
"""
function amplitudes_M2_lent(x, τ_D0, τ_DA, τ_A, Φ_D, Φ_A, η, ε, εΦ_dark)
    k_D0, k_DA, k_A = 1 / τ_D0, 1 / τ_DA, 1 / τ_A
    k_t = k_DA - k_D0
    s = Φ_A * k_A * x * k_t / (k_DA - k_A)
    d = ε * x * Φ_A * k_A
    g = εΦ_dark * (1 - x) * k_A
    D_DA = Φ_D * k_D0 * x
    D_D0 = Φ_D * k_D0 * (1 - x)
    T = typeof(s + D_DA + η[1, 1])
    α = Matrix{T}(undef, 2, 3)
    for c in 1:2
        α[c, 1] = η[c, 1] * D_DA - η[c, 2] * s
        α[c, 2] = η[c, 1] * D_D0
        α[c, 3] = η[c, 2] * (s + d + g)
    end
    return α
end

amplitudes_M2_lent(x, p::SensorParams, o::OpticsParams) =
    amplitudes_M2_lent(x, p.τ_D0, p.τ_DA, p.τ_A, p.Φ_D, phi_A(p), o.η, o.ε, o.εΦ_dark)

durees_M2_lent(p::SensorParams) = (p.τ_DA, p.τ_D0, p.τ_A)

"""
    d_amplitudes_M2_lent(x, p, o) -> (α, dα_dx, dα_dτD0, dα_dτDA, dα_dτA)

Dérivées analytiques des amplitudes (matrices 2×3). Les dérivées en τ passent
par les taux : ∂/∂τ = −k² ∂/∂k.
"""
function d_amplitudes_M2_lent(x, p::SensorParams, o::OpticsParams)
    Φ_D, Φ_A, η, ε, εΦd = p.Φ_D, phi_A(p), o.η, o.ε, o.εΦ_dark
    k_D0, k_DA, k_A = 1 / p.τ_D0, 1 / p.τ_DA, 1 / p.τ_A
    k_t = k_DA - k_D0
    den = k_DA - k_A
    abs(den) > 1e-9 || error("τ_DA = τ_A : l'émission sensibilisée n'est plus une différence d'exponentielles")

    # dérivées de chaque composante (s, d, g, D_DA, D_D0)
    s_x, d_x, g_x, DDA_x, DD0_x = Φ_A * k_A * k_t / den, ε * Φ_A * k_A, -εΦd * k_A, Φ_D * k_D0, -Φ_D * k_D0
    s_kD0, DDA_kD0, DD0_kD0 = -Φ_A * k_A * x / den, Φ_D * x, Φ_D * (1 - x)
    s_kDA = Φ_A * k_A * x * (k_D0 - k_A) / den^2
    s_kA, d_kA, g_kA = Φ_A * x * k_t * k_DA / den^2, ε * x * Φ_A, εΦd * (1 - x)

    α = amplitudes_M2_lent(x, p, o)
    dα_dx = _assemble(η, s_x, d_x, g_x, DDA_x, DD0_x)
    # ∂/∂τ = −k² ∂/∂k
    dα_dτD0 = -k_D0^2 .* _assemble(η, s_kD0, 0.0, 0.0, DDA_kD0, DD0_kD0)
    dα_dτDA = -k_DA^2 .* _assemble(η, s_kDA, 0.0, 0.0, 0.0, 0.0)
    dα_dτA = -k_A^2 .* _assemble(η, s_kA, d_kA, g_kA, 0.0, 0.0)
    return α, dα_dx, dα_dτD0, dα_dτDA, dα_dτA
end

function _assemble(η, s, d, g, DDA, DD0)
    α = Matrix{Float64}(undef, 2, 3)
    for c in 1:2
        α[c, 1] = η[c, 1] * DDA - η[c, 2] * s
        α[c, 2] = η[c, 1] * DD0
        α[c, 3] = η[c, 2] * (s + d + g)
    end
    return α
end

# -----------------------------------------------------------------------------
# M2 rapide : échange rapide, le donneur voit un taux de transfert moyen x·k_t
# -----------------------------------------------------------------------------

"""
    amplitudes_M2_rapide(x, p, o) -> (τ, α)

Commutateur de mécanisme : un seul τ_eff = 1/(k_D0 + x k_t) qui glisse avec x.
La montée du canal jaune a alors la constante τ_eff, qui varie avec [Cl⁻] —
c'est le test de discrimination à faire sur le vrai canal jaune.
"""
function amplitudes_M2_rapide(x, p::SensorParams, o::OpticsParams)
    η = o.η
    k_D0, k_DA, k_A = 1 / p.τ_D0, 1 / p.τ_DA, 1 / p.τ_A
    k_t = k_DA - k_D0
    k_eff = k_D0 + x * k_t
    s = phi_A(p) * k_A * x * k_t / (k_eff - k_A)
    d = o.ε * x * phi_A(p) * k_A
    g = o.εΦ_dark * (1 - x) * k_A
    D = p.Φ_D * k_D0
    α = Matrix{Float64}(undef, 2, 2)
    for c in 1:2
        α[c, 1] = η[c, 1] * D - η[c, 2] * s
        α[c, 2] = η[c, 2] * (s + d + g)
    end
    return (1 / k_eff, p.τ_A), α
end

"""
    expsums(x, p, o) -> (ExpSum bleu, ExpSum jaune)

Les deux canaux pour une fraction active `x`, selon `p.mecanisme`.
"""
function expsums(x::Real, p::SensorParams, o::OpticsParams)
    if p.β_A > 0 || p.β_D > 0
        return canaux_etats(etats_blanchis(x, p)..., p, o)
    end
    if p.mecanisme === :M2_lent
        α = amplitudes_M2_lent(x, p, o)
        τ = collect(durees_M2_lent(p))
        return ExpSum(α[1, :], τ), ExpSum(α[2, :], τ)
    elseif p.mecanisme === :M2_rapide
        τ, α = amplitudes_M2_rapide(x, p, o)
        return ExpSum(α[1, :], collect(τ)), ExpSum(α[2, :], collect(τ))
    end
    error("mécanisme inconnu : $(p.mecanisme)")
end

"""Intégrale (photons par excitation) d'une somme d'exponentielles : Σ αᵢτᵢ."""
integrale(es::ExpSum) = sum(es.α .* es.τ)

"""Intégrales analytiques attendues (test 2) : donneur et émission sensibilisée."""
integrale_donneur(x, p::SensorParams) = p.Φ_D * (x * p.τ_DA / p.τ_D0 + 1 - x)
efficacite_FRET(p::SensorParams) = 1 - p.τ_DA / p.τ_D0
integrale_sensibilisee(x, p::SensorParams) = phi_A(p) * x * efficacite_FRET(p)

# -----------------------------------------------------------------------------
# Modèle à N états du donneur (changements 1 et 5)
# -----------------------------------------------------------------------------

"""
    EtatDonneur(τ, accepteur[, donneur])

Une sous-population de capteurs, dont le donneur a la durée de vie τ.
`accepteur` dit ce que reçoit l'énergie :
- `:actif`   accepteur déprotoné émetteur : émission sensibilisée et excitation directe (M2) ;
- `:sombre`  transfert vers un accepteur qui n'émet pas (hypothèse H2) ;
- `:protone` accepteur protoné, pas de transfert ; émission ε_dark·Φ_dark éventuelle ;
- `:nu`      pas d'accepteur (maturation incomplète, accepteur blanchi).
`donneur = false` : donneur blanchi ; seule l'excitation directe de l'accepteur reste.
"""
struct EtatDonneur
    τ::Float64
    accepteur::Symbol
    donneur::Bool
end
EtatDonneur(τ, accepteur) = EtatDonneur(τ, accepteur, true)

"""
    canaux_etats(etats, pops, p, o) -> (ExpSum bleu, ExpSum jaune)

Généralisation de M2 lent à un nombre quelconque d'états du donneur, de
populations `pops`. Le taux radiatif du donneur k_r = Φ_D/τ_D0 est le même dans
tous les états : l'amplitude du donneur à t = 0 vaut k_r × population, seule la
durée de vie change. Un état `:actif` de durée de vie τ transfère au taux
1/τ − 1/τ_D0 vers l'accepteur. Avec les deux états de M2 lent, redonne exactement
`expsums` (test 21).
"""
function canaux_etats(etats::AbstractVector{EtatDonneur}, pops::AbstractVector, p::SensorParams, o::OpticsParams)
    k_r, k_A, ΦA = p.Φ_D / p.τ_D0, 1 / p.τ_A, phi_A(p)
    τs, αB, αY = Float64[], Float64[], Float64[]
    ajoute!(τ, aD, aA) = (push!(τs, τ); push!(αB, o.η[1, 1] * aD + o.η[1, 2] * aA); push!(αY, o.η[2, 1] * aD + o.η[2, 2] * aA))
    for (e, q) in zip(etats, pops)
        q == 0 && continue
        if e.donneur
            ajoute!(e.τ, k_r * q, 0.0)
            if e.accepteur === :actif && e.τ < p.τ_D0
                k_t = 1 / e.τ - 1 / p.τ_D0
                s = ΦA * k_A * q * k_t / (1 / e.τ - k_A)
                ajoute!(e.τ, 0.0, -s)
                ajoute!(p.τ_A, 0.0, s)
            end
        end
        e.accepteur === :actif && o.ε > 0 && ajoute!(p.τ_A, 0.0, o.ε * q * ΦA * k_A)
        e.accepteur === :protone && o.εΦ_dark > 0 && ajoute!(p.τ_A, 0.0, o.εΦ_dark * q * k_A)
    end
    return ExpSum(αB, τs), ExpSum(αY, τs)
end

"""États et populations de M2 lent : (actif, τ_DA, x) et (protoné, τ_D0, 1 − x)."""
etats_M2(x, p::SensorParams) = ([EtatDonneur(p.τ_DA, :actif), EtatDonneur(p.τ_D0, :protone)], [x, 1 - x])

"""Moyenne de e^{−κu} sur u ∈ [0, 1] : fraction moyenne de molécules intactes pendant l'acquisition."""
moyenne_survie(κ) = κ < 1e-12 ? 1.0 : -expm1(-κ) / κ

"""
    etats_blanchis(x, p) -> (etats, pops)

Blanchiment pendant l'acquisition (changement 5) : le nombre de donneurs et
d'accepteurs intacts décroît exponentiellement, chacun à son taux, et
indépendamment l'un de l'autre. `β_D`, `β_A` sont les fractions blanchies à la
fin de l'acquisition (taux × durée = −ln(1 − β)). L'histogramme accumulé est
celui des populations moyennées sur l'acquisition :
- capteur intact (actif x, protoné 1 − x) : moyenne de e^{−(κ_A+κ_D)u} ;
- accepteur blanchi, donneur intact : donneur nu, durée de vie τ_nu (ou, si
  `blanchi_sombre`, extincteur sombre : transfert inchangé, aucune émission) ;
- donneur blanchi, accepteur intact : excitation directe de l'accepteur seule.
Simplification : les taux ne dépendent pas de l'état (un donneur éteint par FRET
passe en réalité moins de temps excité, et l'accepteur, excité surtout par
transfert, blanchit plus vite quand x est grand).
"""
function etats_blanchis(x, p::SensorParams)
    κA, κD = -log1p(-p.β_A), -log1p(-p.β_D)
    mAD, mD, mA = moyenne_survie(κA + κD), moyenne_survie(κD), moyenne_survie(κA)
    if p.blanchi_sombre
        # accepteur blanchi en extincteur sombre : il transfère encore (fraction x), sans émettre
        etats = [EtatDonneur(p.τ_DA, :actif), EtatDonneur(p.τ_D0, :protone), EtatDonneur(p.τ_DA, :sombre),
                 EtatDonneur(p.τ_D0, :nu), EtatDonneur(p.τ_DA, :actif, false), EtatDonneur(p.τ_D0, :protone, false)]
        pops = [x * mAD, (1 - x) * mAD, x * (mD - mAD), (1 - x) * (mD - mAD), x * (mA - mAD), (1 - x) * (mA - mAD)]
        return etats, pops
    end
    etats = [EtatDonneur(p.τ_DA, :actif), EtatDonneur(p.τ_D0, :protone), EtatDonneur(p.τ_nu, :nu),
             EtatDonneur(p.τ_DA, :actif, false), EtatDonneur(p.τ_D0, :protone, false)]
    pops = [x * mAD, (1 - x) * mAD, mD - mAD, x * (mA - mAD), (1 - x) * (mA - mAD)]
    return etats, pops
end
