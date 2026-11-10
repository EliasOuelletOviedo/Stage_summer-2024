# Batterie.jl — tester le module d'ingestion (changement 1).
#
# Des titrations synthétiques dont on connaît la vérité, générées sous des
# hypothèses différentes, passent dans le même ajustement global que les données
# réelles (Ingest.jl). Même IRF, mêmes photons et même fond par histogramme, même
# grille de concentrations que la titration du 2026-07-20 : la comparaison porte
# sur le mécanisme seul. Pour chaque jeu : le χ²ᵣ, les paramètres retrouvés contre
# les paramètres générés, et le biais qu'aurait une lecture par le ratio calibrée
# sur les paramètres retrouvés.

"""La titration réelle et son ajustement : la grille commune à tous les jeux."""
struct Reference
    hs::Vector{Histogramme}
    g::ResultatGlobal
    N::Vector{Float64}          # photons dans la fenêtre, par histogramme
    C::Vector{Float64}
end

function Reference(hs::Vector{Histogramme}, inst::InstrumentParams; kw...)
    g = ajuster_global(hs, inst; kw...)
    return Reference(hs, g, [sum(h.n[inst.fenetre]) for h in hs], [h.C for h in hs])
end

"""
Une hypothèse : `etats(x, p)` donne les sous-populations du donneur pour une
fraction active x (voir `canaux_etats`) ; `x_de(x)` transforme la fraction ajustée
sur les données réelles en fraction générée (H4 la dilate pour garder x·E) ;
`irf` : IRF de génération (`nothing` : la mesurée) ; `t0s` : un t₀ par histogramme
(gigue), sinon le t₀ global de l'ajustement réel.
"""
Base.@kwdef struct Hypothese
    cle::String
    description::String
    p::SensorParams
    etats::Function = etats_M2
    x_de::Function = identity
    irf::Union{Nothing,Vector{Float64}} = nothing
    t0s::Union{Nothing,Vector{Float64}} = nothing
end

"""
    titration_synthetique(hyp, ref, m, rng; bruit) -> Vector{Histogramme}

Canal bleu complet (donneur et fuite de l'accepteur dans le bleu, η_BA), mis à
l'échelle pour que chaque histogramme ait les photons de son homologue réel,
fond compris. `bruit = false` : comptages attendus.
"""
function titration_synthetique(hyp::Hypothese, ref::Reference, m::Modele, rng; bruit::Bool = true)
    inst = hyp.irf === nothing ? m.inst : remplacer(m.inst; irf = hyp.irf)
    ws = EspaceTravail(inst)
    W = inst.fenetre
    out = Histogramme[]
    for (h, hr) in enumerate(ref.hs)
        x = hyp.x_de(ref.g.x[h])
        eB, _ = canaux_etats(hyp.etats(x, hyp.p)..., hyp.p, m.o)
        t0 = hyp.t0s === nothing ? ref.g.t₀ : hyp.t0s[h]
        μ = forward(eB, remplacer(inst; t₀ = t0), ws)
        masquer!(μ, W)
        b = ref.g.b[h]
        S = (ref.N[h] - b * length(W)) / sum(μ[W])
        μ[W] .= S .* μ[W] .+ b
        n = bruit ? tirer!(zeros(length(μ)), μ, rng) : copy(μ)
        push!(out, Histogramme(hr.proteine, hr.dossier, hr.C, hr.replicat, n))
    end
    return out
end

# -----------------------------------------------------------------------------
# Vérité d'une hypothèse : ce que le canal bleu « devrait » rendre
# -----------------------------------------------------------------------------

"""
Fraction d'amplitude des états qui transfèrent (ce que mesure x dans l'ajustement)
et déficit d'intensité bleue relatif au donneur seul (« x·E » : 1 − Σ q·τ/τ_D0).
"""
function verite_bleue(etats, pops, p::SensorParams)
    don = [e.donneur for e in etats]
    transf = [e.donneur && e.accepteur in (:actif, :sombre) && e.τ < p.τ_D0 for e in etats]
    x_amp = sum(pops[transf]) / sum(pops[don])
    deficit = 1 - sum(q * e.τ / p.τ_D0 for (e, q) in zip(etats, pops) if e.donneur) / sum(pops)
    return x_amp, deficit
end

"""Rapport J/B (totaux de la fenêtre, sans bruit ni fond) d'un monde décrit par des états."""
function ratio_etats(etats, pops, p::SensorParams, m::Modele, ws::EspaceTravail)
    eB, eY = canaux_etats(etats, pops, p, m.o)
    μB = forward(eB, m.inst, ws)
    μY = forward(eY, m.inst.Δt_jaune == 0 ? m.inst : remplacer(m.inst; t₀ = m.inst.t₀ + m.inst.Δt_jaune), ws)
    W = m.inst.fenetre
    return sum(μY[W]) / sum(μB[W])
end

"""
    biais_ratio(hyp, iso_vrai, m_ajuste, m, ws; C) -> Vector{(Ĉ, hors_plage)}

Lecture par le ratio calibrée sur le modèle reconstruit à partir des paramètres
retrouvés (`m_ajuste`), appliquée au monde vrai de l'hypothèse (isotherme
`iso_vrai` sur la fraction générée, états de l'hypothèse) : Ĉ − C est le biais.
"""
function biais_ratio(hyp::Hypothese, iso_vrai, m_ajuste::Modele, m::Modele, ws::EspaceTravail; C = (5.0, 20.0, 60.0), C_max = 150.0)
    grille = grille_calibration(C_max)
    O = [begin
        μB, μY = formes(c, m_ajuste, ws)
        sum(μY) / sum(μB)
    end for c in grille]
    croiss = O[end] > O[1]
    all(>(0), croiss ? diff(O) : -diff(O)) || return [(NaN, true) for _ in C]
    cal = Calibration(collect(grille), O, croiss, C_max, nothing)
    pv = remplacer(hyp.p; x₀ = iso_vrai.x₀, x_sat = iso_vrai.x_sat, K_D = iso_vrai.K_D, n = iso_vrai.n)
    return map(C) do c
        R = ratio_etats(hyp.etats(fraction_active(c, pv), hyp.p)..., hyp.p, m, ws)
        inverser(cal, R)
    end
end

# -----------------------------------------------------------------------------
# Signatures comparables aux données réelles
# -----------------------------------------------------------------------------

"""
    signatures(hs, g, inst, inst_u) -> NamedTuple

Ce que l'on peut comparer entre une titration (réelle ou synthétique) et son
ajustement global, en plus du χ²ᵣ :
- distance des phasors (fond retiré, IRF de l'ajustement) au segment τ_D0–τ_DA
  ajusté, en longueur de segment (médiane, max) ;
- écart du τ mono (code de l'utilisateur, décalage et fond libres) entre les
  données et le modèle ajusté du même histogramme : à 0 mM, à saturation (« F »),
  et rms sur toute la titration. Il mesure ce que le mono voit et que le modèle à
  deux états ne reproduit pas.
"""
function signatures(hs::Vector{Histogramme}, g::ResultatGlobal, inst::InstrumentParams, inst_u)
    W = inst.fenetre
    ws = EspaceTravail(inst)
    BD0, BDA = zeros(inst.n_bins), zeros(inst.n_bins)
    base!(BD0, ws, g.τ_D0, g.t₀)
    base!(BDA, ws, g.τ_DA, g.t₀)
    P0, P1 = phasor(BD0, inst), phasor(BDA, inst)
    dist = map(eachindex(hs)) do h
        n = copy(hs[h].n)
        n[W] .-= g.b[h]
        distance_droite(phasor(n, inst), P0, P1)
    end
    δτ = fill(NaN, length(hs))
    if inst_u !== nothing
        for h in eachindex(hs)
            μ = zeros(inst.n_bins)
            μ[W] .= g.μ[h]
            τ_d = ajuster_mono_utilisateur(inst_u, hs[h].n)[1][1]
            τ_m = ajuster_mono_utilisateur(inst_u, μ)[1][1]
            δτ[h] = τ_d - τ_m
        end
    end
    C = [h.C for h in hs]
    return (phasor_mediane = median(dist), phasor_max = maximum(dist),
        dtau_mono_0 = mean(δτ[C.==minimum(C)]), dtau_mono_sat = mean(δτ[C.==maximum(C)]), dtau_mono_rms = sqrt(mean(δτ .^ 2)))
end

# -----------------------------------------------------------------------------
# Un passage de la batterie
# -----------------------------------------------------------------------------

"""
    passer(hyp, ref, m; graine, s_irf, ajustement...) -> NamedTuple

Génère la titration de `hyp`, l'ajuste avec une IRF dilatée de `s_irf` (1 : la
même) et les options de `ajuster_global`, puis compare à la vérité.
"""
function passer(hyp::Hypothese, ref::Reference, m::Modele; graine::Integer = 1, s_irf::Real = 1.0, bruit::Bool = true,
                etiquette = hyp.cle, inst_u = nothing, ajustement...)
    rng = StableRNG(graine)
    hs = titration_synthetique(hyp, ref, m, rng; bruit)
    inst_fit = s_irf == 1 ? m.inst : remplacer(m.inst; irf = irf_dilatee(m.inst.irf, s_irf))
    g = ajuster_global(hs, inst_fit; τ_D0_init = ref.g.τ_D0, τ_DA_init = ref.g.τ_DA, t₀_init = ref.g.t₀, ajustement...)
    C = [h.C for h in hs]
    iso = ajuster_isotherme(C, g.x)
    xg = hyp.x_de.(ref.g.x)
    iso_v = ajuster_isotherme(C, xg)
    # vérité aux deux plateaux
    v0 = verite_bleue(hyp.etats(iso_v.x₀, hyp.p)..., hyp.p)
    vs = verite_bleue(hyp.etats(iso_v.x_sat, hyp.p)..., hyp.p)
    E = 1 - g.τ_DA / g.τ_D0
    # biais du ratio : modèle reconstruit (paramètres retrouvés, optique nominale, IRF de l'ajustement)
    ws = EspaceTravail(inst_fit)
    m_aj = Modele(remplacer(m.p; τ_D0 = g.τ_D0, τ_DA = min(g.τ_DA, 0.999 * g.τ_D0), x₀ = iso.x₀, x_sat = iso.x_sat,
        K_D = iso.K_D, n = iso.n), m.o, remplacer(inst_fit; t₀ = g.t₀))
    b = try
        biais_ratio(hyp, iso_v, m_aj, m, ws)
    catch
        [(NaN, true) for _ in 1:3]
    end
    dof = length(hs) * length(m.inst.fenetre) - (3 + 3length(hs))
    sg = signatures(hs, g, inst_fit, inst_u)
    return (jeu = etiquette, hypothese = hyp.cle, description = hyp.description, s_irf_ajustement = s_irf, graine = graine,
        chi2_r = g.chi2_global, sigma_chi2_poisson = sqrt(2 / dof), converge = g.converge,
        tau_D0_vrai = hyp.p.τ_D0, tau_D0 = g.τ_D0, tau_DA_vrai = hyp.p.τ_DA, tau_DA = g.τ_DA,
        E_vrai = 1 - hyp.p.τ_DA / hyp.p.τ_D0, E = E, t0 = g.t₀,
        x0_vrai = v0[1], x0 = iso.x₀, x_sat_vrai = vs[1], x_sat = iso.x_sat, K_vrai = iso_v.K_D, K = iso.K_D,
        n_vrai = iso_v.n, n = iso.n, x0E_vrai = v0[2], x0E = iso.x₀ * E, xsatE_vrai = vs[2], xsatE = iso.x_sat * E,
        biais_ratio_5mM = b[1][1] - 5.0, biais_ratio_20mM = b[2][1] - 20.0, biais_ratio_60mM = b[3][1] - 60.0,
        hors_plage_20mM = b[2][2], sg...)
end

"""
Verdict sur une ligne de la batterie, relatif au contrôle H1 (même IRF) et au
χ²ᵣ observé sur les données réelles :
- « indiscernable » : χ²ᵣ dans le bruit de Poisson du contrôle (< contrôle + 3σ) ;
- « compatible » : discernable, mais un χ²ᵣ inférieur à celui des données réelles ;
- « exclu » : χ²ᵣ supérieur à celui des données réelles, ou paramètres retrouvés
  incompatibles avec ceux des données (renseigné à part).
"""
function verdict(chi2, chi2_controle, σ, chi2_reel)
    chi2 < chi2_controle + 3σ && return "indiscernable"
    chi2 <= chi2_reel && return "compatible"
    return "exclu"
end
