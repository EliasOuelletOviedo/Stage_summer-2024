# Étape 12 — accueillir le canal jaune mesuré (changement 3).
#
# Sans banc :
#   1. η_YD tiré des images déjà acquises : le K apparent du ratio TIFF (44,2 mM,
#      prédit 40,6) et sa plage dynamique (×2,99) dépendent de η_YD, le seul paramètre
#      purement supposé de ce calcul ; on inverse les deux, dans le monde nominal et
#      dans les mondes IRF/t₀ de l'étape 11 ;
#   2. le module d'ingestion du jaune (src/Jaune.jl) testé en aller-retour sur un
#      banc simulé : trois concentrations (0, 30 mM, saturant), temps fixé, valeurs
#      vraies volontairement différentes des valeurs supposées ;
#   3. les prédictions que la mesure confrontera : J/B à chaque concentration et sa
#      plage dynamique (indépendante du gain), sous le modèle et sous H2/H3a ; la
#      montée de l'émission sensibilisée est-elle observable ?
# Avec le banc : lire les .sdt des deux canaux, puis `ajuster_jaune` et
# `parametres_mesures` (voir le README, § Changement 3).
#
#   julia --project=. -t 8 step12_canal_jaune.jl
#
# Sorties : sorties/canal_jaune/*.csv, sorties/calibration/eta_YD_inverse.toml,
#           sorties/canal_jaune/parametres_mesures_exemple.toml, sorties/figures/fig12_canal_jaune.png

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, CSV, DataFrames, StableRNGs, TOML, Dates, StaticArrays

m, irf_mat, cfg = S.charger_configuration()
println(S.decrire(m))
dossier = S.dossier_sorties("canal_jaune")
mkpath(dossier)
ctx = S.Contexte(S.EspaceTravail(m.inst), S.charger_instance_utilisateur(irf_mat))
ws = ctx.ws
prot = cfg["capteur"]["proteine"]

# -----------------------------------------------------------------------------
# 1. η_YD depuis les images TIFF (K apparent et plage dynamique de B/J)
# -----------------------------------------------------------------------------
mes = CSV.read(normpath(joinpath(@__DIR__, "..", "..", "tests", "csv_exports", "hill_fit_summary.csv")), DataFrame)
r = mes[(mes.protein.==prot).&(mes.metric.=="fret"), :]
K_mes, plage_mes = r.kd_mM[1], r.asymptote[1] / r.baseline[1]
mondes = [("nominal", m)]
fp = S.dossier_sorties("batterie", "reel_perturbations.csv")
if isfile(fp)
    reel = CSV.read(fp, DataFrame)
    for (prefixe, nomm) in (("IRF ×0.90", "IRF × 0,9"), ("IRF ×1.10", "IRF × 1,1"), ("t₀ libre par histogramme", "t₀ par histogramme"))
        i = findfirst(x -> startswith(x, prefixe), reel.perturbation)
        i === nothing || push!(mondes, (nomm, S.monde_reajuste(m, irf_mat, reel[i, :])[1]))
    end
end
inversions = DataFrame()
courbes = DataFrame()
for (nomm, mm) in mondes
    c = S.courbe_fuite(mm, S.Contexte(S.EspaceTravail(mm.inst), nothing))
    c.monde .= nomm
    global courbes = vcat(courbes, c)
    fK = S.inverser_fuite(c, :K_app, K_mes)
    fP = S.inverser_fuite(c, :plage, plage_mes)
    i1 = argmin(abs.(c.facteur .- 1))
    push!(inversions, (monde = nomm, K_app_predit = c.K_app[i1], plage_predite = c.plage[i1], facteur_par_K = fK, eta_YD_par_K = fK * m.o.η[2, 1],
        facteur_par_plage = fP, eta_YD_par_plage = fP * m.o.η[2, 1]))
    @printf("%-20s K app B/J prédit %.1f mM (mesuré %.1f) → η_YD × %.2f ; plage prédite ×%.2f (mesurée ×%.2f) → η_YD × %.2f\n",
        nomm, c.K_app[i1], K_mes, fK, c.plage[i1], plage_mes, fP)
end
CSV.write(joinpath(dossier, "eta_YD_courbes.csv"), courbes)
CSV.write(joinpath(dossier, "eta_YD_inversion.csv"), inversions)
open(S.dossier_sorties("calibration", "eta_YD_inverse.toml"), "w") do io
    n = inversions[inversions.monde.=="nominal", :]
    TOML.print(io, Dict("date" => string(Dates.today()), "eta_YD_spectres" => m.o.η[2, 1],
        "par_K_apparent" => Dict("K_mesure_TIFF" => K_mes, "facteur" => n.facteur_par_K[1], "eta_YD" => n.eta_YD_par_K[1],
            "plage_mondes_IRF" => [minimum(inversions.facteur_par_K), maximum(inversions.facteur_par_K)]),
        "par_plage_dynamique" => Dict("plage_mesuree_TIFF" => plage_mes, "facteur" => n.facteur_par_plage[1], "eta_YD" => n.eta_YD_par_plage[1],
            "plage_mondes_IRF" => [minimum(inversions.facteur_par_plage), maximum(inversions.facteur_par_plage)]),
        "reserve" => "Images TIFF non photoniques (Bliq) : un décalage d'intensité change K apparent et plage ; à confirmer au banc (step12)");
        sorted = true)
end

# -----------------------------------------------------------------------------
# 2. Aller-retour : banc simulé → ingestion du bleu → ingestion du jaune
# -----------------------------------------------------------------------------
banc = S.remplacer(S.avec_fuite(m, 0.8); p = S.remplacer(m.p; τ_A = 2.8, Φ_A_eff = 0.50),
    o = S.remplacer(S.avec_fuite(m, 0.8).o; ε = 0.04), inst = S.remplacer(m.inst; Δt_jaune = 0.030))
vrai = S.verite_jaune(banc)
Cs = [0.0, 30.0, Inf]
W = m.inst.fenetre
fond_B, fond_J = 50.0, 30.0
ar = DataFrame()
exemple = nothing
for graine in 1:20
    rng = StableRNG(12_000 + graine)
    hB, hY = S.Histogramme[], Vector{Float64}[]
    for (i, C) in enumerate(Cs)
        μB, μY = S.attendus(C, banc, 1e6, :A, ws)       # temps fixé : 10⁶ photons bleus à 20 mM
        μB[W] .+= fond_B
        μY[W] .+= fond_J
        nB, nY = zeros(length(μB)), zeros(length(μY))
        S.tirer_deux!(nB, nY, μB, μY, rng)
        push!(hB, S.Histogramme(prot, string(C), C, 1, nB))
        push!(hY, nY)
    end
    # Deux procédures : en deux temps (bleu puis jaune, A et x du bleu tenus pour exacts)
    # et conjointe (recommandée). τ_D0 et τ_DA figés à leurs valeurs de la titration.
    gB = S.ajuster_global(hB, m.inst; τ_D0_init = m.p.τ_D0, τ_DA_init = m.p.τ_DA, tau_D0_fixe = m.p.τ_D0, tau_DA_fixe = m.p.τ_DA)
    r2 = S.ajuster_jaune(hY, gB, m.inst; Δt_libre = true)
    rc = S.ajuster_bleu_jaune([h.n for h in hB], hY, m.inst, m; τ_fixes = true, Δt_libre = true)
    graine == 1 && (global exemple = rc)
    for (protocole, rj, χJ, χB) in (("deux temps", r2, r2.chi2, gB.chi2_global), ("conjoint", rc, rc.chi2_jaune, rc.chi2_bleu))
        for (k, v) in pairs(vrai)
            est = getfield(rj, k)
            push!(ar, (protocole = protocole, graine = graine, parametre = String(k), vrai = v, estime = est, sigma = rj.σ[k],
                ecart_norme = (est - v) / rj.σ[k], chi2_r_jaune = χJ, chi2_r_bleu = χB))
        end
    end
end
CSV.write(joinpath(dossier, "aller_retour.csv"), ar)
println("\nAller-retour du canal jaune (20 bancs simulés, valeurs vraies ≠ valeurs supposées) :")
for g in groupby(ar, [:protocole, :parametre]; sort = false)
    @printf("  %-20s %-4s vrai %8.4f  estimé %8.4f ± %.4f (σ de Fisher %.4f)  écart normé %+.2f ± %.2f\n", g.protocole[1], g.parametre[1], g.vrai[1],
        mean(g.estime), std(g.estime), mean(g.sigma), mean(g.ecart_norme), std(g.ecart_norme))
end
open(joinpath(dossier, "parametres_mesures_exemple.toml"), "w") do io
    println(io, "# Exemple du bloc écrit par `parametres_mesures` après la mesure du canal jaune,")
    println(io, "# ici sur le banc SIMULÉ de l'aller-retour (graine 1) : ne pas recopier tel quel.\n")
    println(io, S.parametres_mesures(exemple, m; date = "AAAA-MM-JJ"))
end

# -----------------------------------------------------------------------------
# 3. Prédictions que la mesure confrontera
# -----------------------------------------------------------------------------
avec_p(mm; kw...) = S.remplacer(mm; p = S.remplacer(mm.p; kw...))
avec_o(mm; kw...) = S.remplacer(mm; o = S.remplacer(mm.o; kw...))
"""J/B (totaux, sans bruit ni fond) d'un monde décrit par ses états."""
function jb_etats(C, mm, etats)
    eB, eY = S.canaux_etats(etats(S.fraction_active(C, mm.p), mm.p)..., mm.p, mm.o)
    μB, μY = S.forward(eB, mm.inst, ws), S.forward(eY, mm.inst, ws)
    return sum(μY[W]) / sum(μB[W])
end
Cg = vcat(0.0, exp.(range(log(0.5), log(1000), length = 60)), Inf)
pred = DataFrame()
scenarios = Tuple{String,Any,Any}[("nominal", m, S.etats_M2)]
for fη in (0.5, 1.5), fΦ in (0.7, 1.3), fτ in (0.8, 1.2), ε in (0.01, 0.05)
    mm = avec_o(avec_p(S.avec_fuite(m, fη); Φ_A_eff = fΦ * S.phi_A(m.p), τ_A = fτ * m.p.τ_A); ε = ε)
    push!(scenarios, (@sprintf("supposés η×%.1f Φ×%.1f τ_A×%.1f ε=%.2f", fη, fΦ, fτ, ε), mm, S.etats_M2))
end
push!(scenarios, ("H2 accepteur non émissif", m, (x, q) -> ([S.EtatDonneur(q.τ_DA, :sombre), S.EtatDonneur(q.τ_D0, :protone)], [x, 1 - x])))
push!(scenarios, ("H3a τ=1,5 ns 20 %", m, (x, q) -> ([S.EtatDonneur(q.τ_DA, :actif), S.EtatDonneur(1.5, :actif), S.EtatDonneur(q.τ_D0, :protone)],
    [0.8x, 0.2x, 1 - x])))
for (nomm, mm) in mondes[2:end]
    push!(scenarios, ("monde " * nomm, mm, S.etats_M2))
end
for (nom_s, mm, et) in scenarios, C in Cg
    push!(pred, (scenario = nom_s, C_mM = C, JB = jb_etats(C, mm, et)))
end
CSV.write(joinpath(dossier, "predictions_JB.csv"), pred)
println("\nJ/B prédit (totaux de la fenêtre, gains égaux) :")
for nom_s in unique(pred.scenario)
    d = pred[pred.scenario.==nom_s, :]
    j0, js = d.JB[d.C_mM.==0][1], d.JB[isinf.(d.C_mM)][1]
    i30 = argmin(abs.(d.C_mM .- 30))
    (startswith(nom_s, "supposés") && !(nom_s in ("supposés η×0.5 Φ×1.3 τ_A×1.2 ε=0.05", "supposés η×1.5 Φ×0.7 τ_A×0.8 ε=0.01"))) && continue
    @printf("  %-40s J/B(0) = %.3f  J/B(30) = %.3f  J/B(sat) = %.3f  plage ×%.2f\n", nom_s, j0, d.JB[i30], js, j0 / js)
end
sup = combine(groupby(pred[startswith.(pred.scenario, "supposés"), :], :C_mM), :JB => minimum => :min, :JB => maximum => :max)
plages = [begin
    d = pred[pred.scenario.==s, :]
    d.JB[d.C_mM.==0][1] / d.JB[isinf.(d.C_mM)][1]
end for s in unique(pred.scenario[startswith.(pred.scenario, "supposés")])]
@printf("  enveloppe des paramètres supposés : J/B(0) %.2f–%.2f, J/B(sat) %.2f–%.2f, plage ×%.2f–%.2f\n",
    sup.min[sup.C_mM.==0][1], sup.max[sup.C_mM.==0][1], sup.min[isinf.(sup.C_mM)][1], sup.max[isinf.(sup.C_mM)][1], extrema(plages)...)

# montée observable ?
mt = DataFrame()
for τ in (0.05, 0.10, 0.15, 0.30, 0.50, 1.0, 2.1), N in (1e5, 1e6), Δl in (false, true)
    σ = S.sigma_montee(τ, m, N, ws; ε_libre = true, Δt_libre = Δl)
    push!(mt, (tau_DA = τ, N_jaune = N, epsilon_libre = true, dt_libre = Δl, sigma_tau_DA = σ, z = τ / σ))
end
CSV.write(joinpath(dossier, "montee.csv"), mt)
println("\nMontée de l'émission sensibilisée, canal jaune seul à 0 mM (borne de Cramér-Rao, ε libre) :")
for r in eachrow(mt)
    @printf("  τ_DA = %.2f ns, %g photons jaunes, Δt %-24s : σ(τ_DA) = %.3f ns (τ_DA/σ = %.1f)\n", r.tau_DA, r.N_jaune,
        r.dt_libre ? "inconnu" : "mesuré au diffuseur", r.sigma_tau_DA, r.z)
end
# constante de temps du jaune : moyenne pondérée en intensité, 0 mM contre saturant
τmoy(es) = sum(es.α .* es.τ .^ 2) / sum(es.α .* es.τ)
for C in (0.0, Inf)
    _, eY = S.expsums(S.fraction_active(C, m.p), m.p, m.o)
    @printf("τ moyen (pondéré en intensité) du jaune à %s : %.2f ns (τ_D0 = %.2f, τ_A = %.2f)\n", isinf(C) ? "saturation" : "0 mM", τmoy(eY), m.p.τ_D0, m.p.τ_A)
end

include(joinpath(@__DIR__, "src", "Figures.jl"))
using .Figures, CairoMakie
set_theme!(Figures.theme_sim())
if isdefined(Figures, :figure_canal_jaune)
    Figures.figure_canal_jaune(dossier, S.dossier_sorties("figures", "fig12_canal_jaune.png"); K_mes, plage_mes)
end
