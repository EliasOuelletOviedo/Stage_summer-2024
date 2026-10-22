# Étape 16 — déclins purs (déconvolués de l'IRF) de la titration SCLM du 2026-10-16, 800 nm.
#
# Pour chaque concentration (0 à 131 mM) : somme des 10 acquisitions de 1 s, puis le même
# ajustement global que l'étape 15 (Global.jl) : IRF mesurée de chaque détecteur
# (IRF_800_16X.sdt, fond retiré), t₀ par canal, fond des détecteurs figé, composante de
# diffusion (copie de l'IRF), durées de vie communes, amplitudes libres par concentration.
# Le lien τ_J = τ₂ est testé comme à l'étape 15 et le modèle retenu de la même façon.
#
# Le déclin pur est la partie fluorescente du modèle, sans IRF : la réponse à UNE impulsion
# laser, t = 0 à l'excitation (t₀ de chaque canal retiré), sans fond, sans diffusion, sans
# la queue de l'impulsion précédente, sans coupure du TAC. Comptes par canal TCSPC :
#   n_k = Σ_j a_j·τ_j·(e^{−t_k/τ_j} − e^{−t_{k+1}/τ_j}),  t_k = k·Δ, Δ = pas du TAC.
# Sur une période (12,5 ns), leur somme vaut Σ a_j·τ_j·(1 − e^{−T/τ_j}) : les photons
# fluorescents émis par la somme des 10 acquisitions, ramenés à un TAC qui verrait tout.
#
# Pourquoi un modèle plutôt qu'une déconvolution directe (Fourier, Richardson-Lucy) : la
# déconvolution directe amplifie le bruit de Poisson canal par canal et dépend d'un
# réglage arbitraire (nombre d'itérations, filtre). Ici, la déconvolution est exacte pour
# le modèle, et le χ²ᵣ dit ce que le modèle ne capte pas.
#
#   julia --project=. step16_declins_purs.jl
#   julia --project=. step16_declins_purs.jl --donnees Data/2026-10-16/SCLM --lambda 800

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
include(joinpath(@__DIR__, "src", "deux_canaux", "DeuxCanaux.jl"))
using .SimFLIMFRET, .DeuxCanaux
const S = SimFLIMFRET
const D = DeuxCanaux
using Printf, Statistics, CSV, DataFrames

args = copy(ARGS)
valeur(nom, defaut) = (i = findfirst(==(nom), args); i === nothing ? defaut : args[i+1])
donnees = normpath(joinpath(@__DIR__, valeur("--donnees", "Data/2026-10-16/SCLM")))
λ = parse(Int, valeur("--lambda", "800"))
jeu = basename(donnees)
date = basename(dirname(donnees))
dossier = joinpath(@__DIR__, "sorties", "declins_purs_$(date)_$(jeu)")
mkpath(dossier)

# -----------------------------------------------------------------------------
# 1. Cellules (somme des 10 acquisitions) et ajustement global
# -----------------------------------------------------------------------------
cellules, sources = D.charger(donnees; λs = (λ,))
println("$(jeu), $(λ) nm : $(length(cellules)) concentrations ($(join([@sprintf("%g", c.cl) for c in cellules], ", ")) mM), ",
    "$(D.n_acq(cellules[1])) acquisitions de $(cellules[1].duree_s) s sommées ; fond : $(sources[λ])")
lie = D.ajuster_global(cellules; lien = true)
libre = D.ajuster_global(cellules; lien = false, τ_init = (lie.τ2, lie.τ3, lie.τA, 1.0), t0_init = (lie.t0B, lie.t0J))
test = D.tester_lien(lie, libre)
mg = test.lien_tient ? lie : libre
@printf("Modèle retenu : %s (Δdéviance/χ²ᵣ = %.0f, seuil 3,84) ; χ²ᵣ %.2f, convergé : %s\n",
    mg.lien ? "lié" : "libre", test.statistique, mg.χ2, mg.converge)
@printf("  bleu  : τ₂ %.3f, τ₃ %.3f ns ; t₀ %+.1f ps\n  jaune : τ_J %.3f, τ_A %.3f ns ; t₀ %+.1f ps\n",
    mg.τ2, mg.τ3, 1000mg.t0B, mg.τJ, mg.τA, 1000mg.t0J)

# -----------------------------------------------------------------------------
# 2. Déclins purs
# -----------------------------------------------------------------------------
Δ = cellules[1].pas_ns
T = 12.5                                   # période du laser (ns)
nb = floor(Int, T / Δ + 1e-9)              # canaux d'une période complète
t = (0:nb) .* Δ                            # bords des canaux
bin_exp(a, τ) = [a * τ * (exp(-t[k] / τ) - exp(-t[k+1] / τ)) for k in 1:nb]
etiquette(c) = replace(@sprintf("%g", c), "." => "p")
csv = DataFrame(canal = 0:nb-1, t_debut_ns = t[1:nb], t_fin_ns = t[2:nb+1])
resume = DataFrame()
for (c, cl) in enumerate(mg.cl)
    pB = bin_exp(mg.aB[c, 2], mg.τ2) .+ bin_exp(mg.aB[c, 3], mg.τ3)
    pJ = bin_exp(mg.aJ[c, 2], mg.τJ) .+ bin_exp(mg.aJ[c, 3], mg.τA)
    csv[!, "bleu_$(etiquette(cl))mM"] = pB
    csv[!, "jaune_$(etiquette(cl))mM"] = pJ
    hB, hJ = D.somme(cellules[c], D.BLEU), D.somme(cellules[c], D.JAUNE)
    push!(resume, (cl_mM = cl,
        photons_mesures_bleu = D.photons(hB), photons_mesures_jaune = D.photons(hJ),
        photons_purs_bleu = sum(pB), photons_purs_jaune = sum(pJ),
        a_tau2_bleu = mg.aB[c, 2], a_tau3_bleu = mg.aB[c, 3], a_tauJ_jaune = mg.aJ[c, 2], a_tauA_jaune = mg.aJ[c, 3],
        photons_diffusion_bleu = mg.IB[c, 1], photons_diffusion_jaune = mg.IJ[c, 1],
        chi2r_bleu = mg.χ2B[c], chi2r_jaune = mg.χ2J[c], canaux_negatifs_jaune = count(<(0), pJ)))
end
# colonnes : canal, temps, puis tout le bleu, puis tout le jaune (concentrations croissantes)
select!(csv, :canal, :t_debut_ns, :t_fin_ns, r"^bleu_", r"^jaune_")
CSV.write(joinpath(dossier, "declins_purs_$(λ)nm.csv"), csv)
param = DataFrame(lambda_nm = λ, modele = mg.lien ? "lie" : "libre", tau2_bleu_ns = mg.τ2, tau3_bleu_ns = mg.τ3,
    tauJ_jaune_ns = mg.τJ, tauA_jaune_ns = mg.τA, t0_bleu_ps = 1000mg.t0B, t0_jaune_ps = 1000mg.t0J,
    pas_ns = Δ, periode_ns = T, chi2r = mg.χ2, delta_deviance_sur_chi2r_lien = test.statistique)
CSV.write(joinpath(dossier, "parametres_$(λ)nm.csv"), param)
CSV.write(joinpath(dossier, "resume_$(λ)nm.csv"), resume)
for r in eachrow(resume)
    @printf("  %6.1f mM : bleu %.3g mesurés → %.3g purs ; jaune %.3g → %.3g ; χ²ᵣ %.2f / %.2f%s\n", r.cl_mM,
        r.photons_mesures_bleu, r.photons_purs_bleu, r.photons_mesures_jaune, r.photons_purs_jaune, r.chi2r_bleu, r.chi2r_jaune,
        r.canaux_negatifs_jaune > 0 ? " ; $(r.canaux_negatifs_jaune) canaux jaunes négatifs" : "")
end

# -----------------------------------------------------------------------------
# 3. Tout sur l'axe du TAC : IRF, histogrammes bruts (non déconvolués), déclins purs
# -----------------------------------------------------------------------------
# Une ligne par canal TCSPC (256, temps du TAC). IRF : comptes bruts du fichier, et l'IRF
# utilisée (fond retiré, somme 1). Histogrammes : somme brute des 10 acquisitions, fond
# compris ; hors de la fenêtre (canaux 13–194 en base 0), le TAC coupe. Déclin pur : le
# même qu'au point 2, recalé sur le TAC : son t = 0 est placé au pic de l'IRF du canal
# plus le t₀ ajusté, puis intégré sur chaque canal (le premier canal touché l'est en partie).
# C'est un repère pour comparer à l'œil : le modèle exact est données ≈ IRF ⊛ pur + diffusion + fond.
chemin_irf = joinpath(dirname(donnees), "IRF_$(λ)_16X.sdt")
irf_brute, _, duree_irf, _ = S.lire_sdt_voies(chemin_irf)
n = size(irf_brute, 2)
tt = (0:n) .* Δ
pic = Dict(D.BLEU => argmax(cellules[1].irf[D.BLEU, :]), D.JAUNE => argmax(cellules[1].irf[D.JAUNE, :]))
debut = Dict(D.BLEU => (pic[D.BLEU] - 1) * Δ + mg.t0B, D.JAUNE => (pic[D.JAUNE] - 1) * Δ + mg.t0J)
bin_tac(a, τ, ts) = [a * τ * (exp(-max(tt[k] - ts, 0) / τ) - exp(-max(tt[k+1] - ts, 0) / τ)) for k in 1:n]
tout = DataFrame(canal = 0:n-1, t_debut_ns = tt[1:n], t_fin_ns = tt[2:n+1],
    dans_fenetre = [k in D.FENETRE for k in 1:n],
    t_apres_excitation_bleu_ns = tt[1:n] .- debut[D.BLEU], t_apres_excitation_jaune_ns = tt[1:n] .- debut[D.JAUNE],
    irf_bleu_brute = irf_brute[D.BLEU, :], irf_jaune_brute = irf_brute[D.JAUNE, :],
    irf_bleu_utilisee = cellules[1].irf[D.BLEU, :], irf_jaune_utilisee = cellules[1].irf[D.JAUNE, :])
for (canal, nom) in ((D.BLEU, "bleu"), (D.JAUNE, "jaune")), (c, cl) in enumerate(mg.cl)
    tout[!, "$(nom)_brut_$(etiquette(cl))mM"] = D.somme(cellules[c], canal)
end
for (c, cl) in enumerate(mg.cl)
    tout[!, "bleu_pur_$(etiquette(cl))mM"] = bin_tac(mg.aB[c, 2], mg.τ2, debut[D.BLEU]) .+ bin_tac(mg.aB[c, 3], mg.τ3, debut[D.BLEU])
end
for (c, cl) in enumerate(mg.cl)
    tout[!, "jaune_pur_$(etiquette(cl))mM"] = bin_tac(mg.aJ[c, 2], mg.τJ, debut[D.JAUNE]) .+ bin_tac(mg.aJ[c, 3], mg.τA, debut[D.JAUNE])
end
CSV.write(joinpath(dossier, "irf_bruts_purs_$(λ)nm.csv"), tout)
@printf("IRF : %s (%.0f s), pic au canal %d (bleu) et %d (jaune) en base 0 ; t = 0 des déclins purs à %.3f ns (bleu) et %.3f ns (jaune) du TAC\n",
    basename(chemin_irf), duree_irf, pic[D.BLEU] - 1, pic[D.JAUNE] - 1, debut[D.BLEU], debut[D.JAUNE])

# Échantillon inconnu (dossier `unknown`) : histogrammes bruts sommés et IRF seulement
inc = D.charger_inconnue(donnees; λs = (λ,))
if !isempty(inc)
    u = inc[1]
    CSV.write(joinpath(dossier, "inconnu_irf_bruts_$(λ)nm.csv"), DataFrame(canal = 0:n-1, t_debut_ns = tt[1:n], t_fin_ns = tt[2:n+1],
        dans_fenetre = [k in D.FENETRE for k in 1:n],
        irf_bleu_brute = irf_brute[D.BLEU, :], irf_jaune_brute = irf_brute[D.JAUNE, :],
        irf_bleu_utilisee = u.irf[D.BLEU, :], irf_jaune_utilisee = u.irf[D.JAUNE, :],
        bleu_brut_inconnu = D.somme(u, D.BLEU), jaune_brut_inconnu = D.somme(u, D.JAUNE)))
    @printf("Inconnu : %d acquisitions sommées ; %.3g photons bleus, %.3g jaunes dans la fenêtre\n",
        D.n_acq(u), D.photons(D.somme(u, D.BLEU)), D.photons(D.somme(u, D.JAUNE)))
end

# -----------------------------------------------------------------------------
# 4. Figures : les deux canaux côte à côte, toutes les concentrations
# -----------------------------------------------------------------------------
include(joinpath(@__DIR__, "src", "Figures.jl"))
using CairoMakie
const F = Figures
set_theme!(F.theme_sim())
couleurs = F.RAMPE[round.(Int, range(1, length(F.RAMPE), length(mg.cl)))]
tc = (t[1:nb] .+ t[2:nb+1]) ./ 2
fig = Figure(size = (1150, 520))
Label(fig[0, 1:2], "$(jeu) · $(λ) nm · déclins purs (IRF déconvoluée, fond et diffusion retirés), somme des 10 acquisitions";
    halign = :left, fontsize = 15, color = F.ENCRE)
haut = maximum(maximum(csv[!, n]) for n in names(csv, r"mM$"))
for (j, (canal, nom, τs)) in enumerate((("bleu", "Canal bleu (1)", @sprintf("τ₂ %.2f ns, τ₃ %.2f ns", mg.τ2, mg.τ3)),
                                       ("jaune", "Canal jaune (2)", @sprintf("τ_J %.2f ns, τ_A %.2f ns", mg.τJ, mg.τA))))
    ax = Axis(fig[1, j]; yscale = log10, xlabel = "temps après l'excitation (ns)", ylabel = j == 1 ? "photons par canal TCSPC" : "",
        title = nom, subtitle = F.virg(τs))
    for (c, cl) in enumerate(mg.cl)
        y = csv[!, "$(canal)_$(etiquette(cl))mM"]
        ok = y .> 0
        lines!(ax, tc[ok], y[ok]; color = couleurs[c], label = F.virg(@sprintf("%g mM", cl)))
    end
    xlims!(ax, 0, T)
    ylims!(ax, haut / 1e3, haut * 1.5)
    j == 2 && Legend(fig[1, 3], ax, "[Cl⁻]"; framevisible = false)
end
save(joinpath(dossier, "declins_purs_$(λ)nm.png"), fig; px_per_unit = 2)

# Bruts seuls, chaque courbe divisée par son maximum (les déclins purs restent dans le CSV)
W = collect(D.FENETRE)
tcc = (tt[1:n] .+ tt[2:n+1]) ./ 2
fig2 = Figure(size = (1150, 540))
Label(fig2[0, 1:2], "$(jeu) · $(λ) nm · histogrammes bruts (somme des 10 acquisitions), normalisés à 1 au maximum";
    halign = :left, fontsize = 15, color = F.ENCRE)
# échelle des y partagée : du plus petit compte normalisé des deux canaux jusqu'à 1
ymin = minimum(minimum(b[ok]) / maximum(b[ok]) for nom in ("bleu", "jaune") for cl in mg.cl
    for b in (tout[!, "$(nom)_brut_$(etiquette(cl))mM"],) for ok in ((b .> 0) .& tout.dans_fenetre,))
for (j, (canal, nom, titre)) in enumerate(((D.BLEU, "bleu", "Canal bleu (1)"), (D.JAUNE, "jaune", "Canal jaune (2)")))
    ax = Axis(fig2[1, j]; yscale = log10, xlabel = "temps du TAC (ns)", ylabel = j == 1 ? "comptes normalisés (max = 1)" : "",
        title = titre, subtitle = F.virg(@sprintf("fenêtre utile %.2f–%.2f ns", tt[first(W)], tt[last(W)+1])))
    for (c, cl) in enumerate(mg.cl)
        b = tout[!, "$(nom)_brut_$(etiquette(cl))mM"]
        okb = (b .> 0) .& tout.dans_fenetre
        lines!(ax, tcc[okb], b[okb] ./ maximum(b[okb]); color = couleurs[c], linewidth = 1.5, label = F.virg(@sprintf("%g mM", cl)))
    end
    xlims!(ax, 0, 12.5)
    ylims!(ax, ymin, 1.0)
    g = filter(v -> ymin <= v <= 1, [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0])
    ax.yticks = (g, [F.virg(@sprintf("%g", v)) for v in g])
    j == 2 && Legend(fig2[1, 3], ax, "[Cl⁻]"; framevisible = false)
end
save(joinpath(dossier, "irf_bruts_purs_$(λ)nm.png"), fig2; px_per_unit = 2)
println("Sorties : $(dossier)")
