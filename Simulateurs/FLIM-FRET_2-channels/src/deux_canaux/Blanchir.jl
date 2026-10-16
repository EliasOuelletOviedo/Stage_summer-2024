# Blanchir.jl — étape 6 : le photoblanchiment.
#
# Ce qui compte n'est pas le taux, c'est ce que devient la molécule blanchie :
#   A  l'accepteur blanchi cesse d'absorber : le FRET disparaît (le donneur se libère) ;
#   B  l'accepteur blanchi reste un extincteur sombre : le FRET persiste ;
#   C  le donneur blanchit.
# D'abord mesurer : les 10 acquisitions de 1 s de chaque cellule portent la signature
# (tendances de cinq observables ; la paire décisive est durée de vie du bleu × intensité
# bleue). Ensuite injecter : chaque scénario, à taux contrôlé, dans les tirages.
#
# Le blanchiment est une décroissance exponentielle de la population active pendant
# l'acquisition, pas un facteur appliqué à la fin. L'histogramme d'une acquisition est la
# somme dans le temps des histogrammes instantanés, donc un mélange intact/blanchi dont
# le poids intact est la survie MOYENNE s̄ = β / (−ln(1 − β)) quand une fraction β est
# blanchie à la fin (s̄ = 0,949 pour β = 10 %, et non 0,90).
#
# Décomposition utilisée (modèle global de la cellule) :
#   bleu  = diffusion + F (τ₂, état FRET) + D₀ (τ₃, donneur sans FRET) + fond
#   jaune = diffusion + L (τ_J, composante rapide) + Acc (τ_A, émission de l'accepteur) + fond
# A : F → s̄F ; la population convertie (1 − s̄)·a₂ rejoint D₀ (même amplitude par
#     molécule : même taux radiatif) et fuit dans le jaune avec la forme de τ₃, au
#     coefficient κ = a_J/a₂ de la composante rapide ; Acc → s̄·Acc.
# B : bleu inchangé ; L inchangé ; Acc → s̄·Acc.
# C : F, D₀, L et Acc → × s̄ (toute l'émission de l'accepteur est supposée sensibilisée :
#     l'excitation directe n'est pas séparable dans ces données).
# La diffusion et le fond ne blanchissent pas.

const SCENARIOS = (:A, :B, :C)
const TAUX_BLANCHIMENT = (0.01, 0.02, 0.05, 0.10)

"""Survie moyenne de la population active quand une fraction β est blanchie en fin d'acquisition."""
survie_moyenne(β::Real) = β <= 0 ? 1.0 : β / -log1p(-β)

pleine(v::AbstractVector) = (x = zeros(256); x[FENETRE] .= v; x)

"""
    modele_blanchi(mg, ctx, i, scenario, β) -> (μB, μJ)

Somme attendue de la cellule `i` (256 canaux, fond compris) quand une fraction β des
molécules est blanchie pendant l'acquisition, selon le scénario.
"""
function modele_blanchi(mg::ModeleGlobal, ctx::Contexte, i::Int, scenario::Symbol, β::Real)
    W = FENETRE
    s = survie_moyenne(β)
    dif = mg.IB[i, 1] .* ctx.DB
    F = mg.aB[i, 2] .* ctx.B2
    D0 = mg.aB[i, 3] .* ctx.B3
    difJ = mg.IJ[i, 1] .* ctx.DJ
    L = mg.aJ[i, 2] .* ctx.JJ
    Acc = mg.aJ[i, 3] .* ctx.JA
    fB, fJ = mg.fondB[i][W], mg.fondJ[i][W]
    if scenario === :A
        κ = mg.aJ[i, 2] / mg.aB[i, 2]
        μB = dif .+ s .* F .+ D0 .+ (1 - s) * mg.aB[i, 2] .* ctx.B3 .+ fB
        μJ = difJ .+ s .* L .+ (1 - s) * κ * mg.aB[i, 2] .* ctx.J3 .+ s .* Acc .+ fJ
    elseif scenario === :B
        μB = dif .+ F .+ D0 .+ fB
        μJ = difJ .+ L .+ s .* Acc .+ fJ
    elseif scenario === :C
        μB = dif .+ s .* (F .+ D0) .+ fB
        μJ = difJ .+ s .* (L .+ Acc) .+ fJ
    else
        error("scénario inconnu : $scenario")
    end
    return pleine(μB), pleine(μJ)
end

"""
Facteurs d'injection d'un scénario dans une cellule : par canal TCSPC, la fraction des
photons qui survit (r = min(μ_β/μ, 1)) et l'excès attendu (e = max(μ_β − μ, 0)) de la
somme. Un tirage blanchi = Binomial(tirage intact, r) + Poisson(p·e) : son espérance est
p·μ_β, et il reste couplé au tirage intact.
"""
function facteurs_injection(mg::ModeleGlobal, ctx::Contexte, i::Int, scenario::Symbol, β::Real)
    μβB, μβJ = modele_blanchi(mg, ctx, i, scenario, β)
    μB, μJ = mg.μB[i], mg.μJ[i]
    r(a, b) = [b[k] > 0 ? min(a[k] / b[k], 1.0) : 1.0 for k in eachindex(a)]
    e(a, b) = max.(a .- b, 0.0)
    return (rB = r(μβB, μB), rJ = r(μβJ, μJ), eB = e(μβB, μB), eJ = e(μβJ, μJ))
end

function injecter!(bB, bJ, inj, p::Real, rng::AbstractRNG)
    @inbounds for k in eachindex(bB)
        if bB[k] > 0 && inj.rB[k] < 1
            bB[k] = Float64(rand(rng, Binomial(round(Int, bB[k]), inj.rB[k])))
        end
        inj.eB[k] > 0 && (bB[k] += Float64(rand(rng, Poisson(p * inj.eB[k]))))
        if bJ[k] > 0 && inj.rJ[k] < 1
            bJ[k] = Float64(rand(rng, Binomial(round(Int, bJ[k]), inj.rJ[k])))
        end
        inj.eJ[k] > 0 && (bJ[k] += Float64(rand(rng, Poisson(p * inj.eJ[k]))))
    end
    return nothing
end

const SIGNATURES = Dict(
    :A => (bleu = :monte, jaune = :baisse, JB = :baisse, tau_bleu = :monte, amplitude_FRET = :baisse),
    :B => (bleu = :stable, jaune = :baisse, JB = :baisse, tau_bleu = :stable, amplitude_FRET = :stable),
    :C => (bleu = :baisse, jaune = :baisse, JB = :stable, tau_bleu = :stable, amplitude_FRET = :stable))

sens_tendance(t::Real) = t > 2 ? :monte : t < -2 ? :baisse : :stable

"""
    signature_blanchiment(cel, ctx) -> NamedTuple

Tendances sur les acquisitions successives de la cellule : intensités bleue et jaune, J/B,
durée de vie moyenne du bleu (pondérée en intensité sur τ₂ et τ₃) et fraction d'amplitude
de la composante FRET (τ₂). Pente relative en % par acquisition, et t = pente/erreur.
Règle écrite d'avance (une tendance compte si |t| > 2) :
  toutes stables                                                  → « aucun »
  durée de vie du bleu et intensité bleue montent                 → A
  durée de vie et intensité bleue stables, J/B baisse             → B
  durée de vie stable, intensités bleue et jaune baissent, J/B stable → C
  durée de vie, J/B et amplitude stables, les deux intensités montent → « excitation »
    (dérive de puissance ou de mise au point, pas du blanchiment)
  autre                                                           → « indéterminé »
Le taux mesuré (fraction blanchie par seconde) vient de la pente de l'intensité jaune,
qui baisse dans les trois scénarios ; il n'est retenu que pour A, B ou C.
"""
function signature_blanchiment(cel::Cellule, ctx::Contexte)
    n = n_acq(cel)
    W = FENETRE
    obs = Dict(k => Float64[] for k in (:bleu, :jaune, :JB, :tau_bleu, :amplitude_FRET))
    for i in 1:n
        hB, hJ = acquisition(cel, BLEU, i), acquisition(cel, JAUNE, i)
        nb, nj = photons(hB), photons(hJ)
        a, _ = amplitudes_poisson(hB[W], hcat(ctx.B2, ctx.B3), connu(ctx, BLEU, cel.duree_s), [0.0, 0.0])
        I2, I3 = a[1] * sum(ctx.B2), a[2] * sum(ctx.B3)
        push!(obs[:bleu], nb)
        push!(obs[:jaune], nj)
        push!(obs[:JB], nj / nb)
        push!(obs[:tau_bleu], (I2 * ctx.τ2 + I3 * ctx.τ3) / (I2 + I3))
        push!(obs[:amplitude_FRET], a[1] / (a[1] + a[2]))
    end
    t = Dict(k => tendance(v) for (k, v) in obs)
    sens = Dict(k => sens_tendance(t[k].t) for k in keys(t))
    τb, b, j, r, a = sens[:tau_bleu], sens[:bleu], sens[:jaune], sens[:JB], sens[:amplitude_FRET]
    scenario = if all(==(:stable), values(sens))
        "aucun"
    elseif τb === :monte && b === :monte
        "A"
    elseif τb === :stable && b === :stable && r === :baisse
        "B"
    elseif τb === :stable && b === :baisse && j === :baisse && r === :stable
        "C"
    elseif τb === :stable && r === :stable && a === :stable && b === :monte && j === :monte
        "excitation"
    else
        "indéterminé"
    end
    accord = scenario in ("A", "B", "C") ? count(k -> getfield(SIGNATURES[Symbol(scenario)], k) === sens[k], keys(sens)) : 0
    rel(k) = 100 * t[k].pente / mean(obs[k])
    tj = t[:jaune]
    taux = (scenario in ("A", "B", "C") && tj.t < -2) ? -tj.pente / mean(obs[:jaune]) / cel.duree_s : 0.0
    borne = max(0.0, -(tj.pente - 2tj.se) / mean(obs[:jaune]) / cel.duree_s)
    return (lambda_nm = cel.λ, cl_mM = cel.cl,
        pente_bleu_pct = rel(:bleu), t_bleu = t[:bleu].t, pente_jaune_pct = rel(:jaune), t_jaune = t[:jaune].t,
        pente_JB_pct = rel(:JB), t_JB = t[:JB].t, pente_tau_bleu_pct = rel(:tau_bleu), t_tau_bleu = t[:tau_bleu].t,
        pente_amplitude_FRET_pct = rel(:amplitude_FRET), t_amplitude_FRET = t[:amplitude_FRET].t,
        tau_bleu_ns = mean(obs[:tau_bleu]), amplitude_FRET = mean(obs[:amplitude_FRET]),
        sens = join(["$(k)=$(sens[k])" for k in (:bleu, :jaune, :JB, :tau_bleu, :amplitude_FRET)], " "),
        scenario = scenario, accord_sur_5 = accord, taux_mesure_par_s = taux, borne_2sigma_par_s = borne)
end
