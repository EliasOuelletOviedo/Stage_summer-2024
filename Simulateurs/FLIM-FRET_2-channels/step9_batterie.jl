# Étape 9 — la batterie H1–H5 : le module d'ingestion retrouve-t-il la vérité ?
#
# Changement 1 de la note « Simulateur FLIM-FRET : changements à apporter ».
# Des titrations synthétiques (même IRF, mêmes photons et même fond par
# histogramme, même grille de concentrations que le 2026-07-20) passent dans
# l'ajustement global du canal bleu :
#   H1  le modèle ajusté actuel (contrôle), avec la fuite de l'accepteur dans le bleu ;
#   H2  transfert vers un accepteur qui n'émet pas (zéro émission sensibilisée) ;
#   H2b extinction statique stricte : complexe non émissif ;
#   H3a un troisième état du donneur, transfert intermédiaire (τ = 1,5 ns) ;
#   H3b un troisième état, donneur sans accepteur (τ_nu = 4,1 ns) ;
#   H4  τ_DA = 1,0 ns, x dilaté pour garder x·E ;
#   H5  H4 ajusté avec une IRF fausse de ±10 % ;
#   plus H1 ajusté avec une IRF fausse (±5, ±10 %) et H1 avec la gigue de t₀ observée.
# Sur les données réelles : réajustement avec l'IRF dilatée (±5, ±10 %) ou t₀
# décalé d'un canal, profil de χ²ᵣ en τ_DA, et t₀ libre par histogramme.
#
#   julia --project=. step9_batterie.jl
#
# Sorties : sorties/batterie/{batterie,reel_perturbations,profil_tau_DA,t0_par_histogramme}.csv,
#           sorties/figures/fig9_batterie_ingestion.png

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, Random, CSV, DataFrames, StableRNGs, TOML, Dates

m, irf_mat, cfg = S.charger_configuration()
println(S.decrire(m))
dossier = S.dossier_sorties("batterie")
mkpath(dossier)
prot = cfg["capteur"]["proteine"]
Δ = S.pas_temps(m.inst)

hs = S.lire_titration(prot)
inst_u = S.charger_instance_utilisateur(irf_mat)
ref = S.Reference(hs, m.inst; τ_D0_init = m.p.τ_D0, τ_DA_init = m.p.τ_DA)
sig_reel = S.signatures(hs, ref.g, m.inst, inst_u)
H = length(hs)
C = ref.C
@printf("Données réelles : χ²ᵣ = %.3f, τ_D0 = %.4f, τ_DA = %.4f, t₀ = %.4f ns ; phasors à %.4f (max %.4f) du segment ; τ mono données − modèle : %+.3f ns à 0 mM, %+.3f ns à saturation (rms %.3f)\n",
    ref.g.chi2_global, ref.g.τ_D0, ref.g.τ_DA, ref.g.t₀, sig_reel.phasor_mediane, sig_reel.phasor_max, sig_reel.dtau_mono_0,
    sig_reel.dtau_mono_sat, sig_reel.dtau_mono_rms)

"""Largeur à mi-hauteur d'une IRF (interpolation linéaire), en ns."""
function fwhm(irf, Δ)
    k = argmax(irf)
    h = irf[k] / 2
    i = k
    while i > 1 && irf[i] > h
        i -= 1
    end
    j = k
    while j < length(irf) && irf[j] > h
        j += 1
    end
    a = i + (h - irf[i]) / (irf[i+1] - irf[i])
    b = j - 1 + (irf[j-1] - h) / (irf[j-1] - irf[j])
    return (b - a) * Δ
end
ecart_type(irf) = (t = (1:length(irf)) .- 0.5; μ = sum(t .* irf) / sum(irf); sqrt(sum((t .- μ) .^ 2 .* irf) / sum(irf) + 1 / 12) * Δ)
@printf("IRF mesurée : largeur à mi-hauteur %.0f ps, écart-type %.0f ps (×0,9 : %.0f ps ; ×1,1 : %.0f ps)\n",
    1000fwhm(m.inst.irf, Δ), 1000ecart_type(m.inst.irf), 1000ecart_type(S.irf_dilatee(m.inst.irf, 0.9)), 1000ecart_type(S.irf_dilatee(m.inst.irf, 1.1)))

# -----------------------------------------------------------------------------
# 1. Données réelles : perturbations de l'IRF et de t₀
# -----------------------------------------------------------------------------
function resume_reel(etiquette, g; s = 1.0, t0_mode = "libre")
    iso = S.ajuster_isotherme(C, g.x)
    E = 1 - g.τ_DA / g.τ_D0
    return (perturbation = etiquette, s_irf = s, t0_mode = t0_mode, chi2_r = g.chi2_global, tau_D0 = g.τ_D0,
        tau_DA = g.τ_DA, E = E, t0 = g.t₀, x0 = iso.x₀, x_sat = iso.x_sat, K = iso.K_D, n = iso.n,
        x0E = iso.x₀ * E, xsatE = iso.x_sat * E, sigma_x0 = sqrt(max(iso.cov[1, 1], 0)))
end
inst_s(s) = s == 1 ? m.inst : S.remplacer(m.inst; irf = S.irf_dilatee(m.inst.irf, s))
ajuster(inst; kw...) = S.ajuster_global(hs, inst; τ_D0_init = ref.g.τ_D0, τ_DA_init = ref.g.τ_DA, t₀_init = ref.g.t₀, kw...)

reel = [resume_reel("aucune", ref.g)]
for s in (0.8, 0.85, 0.9, 0.95, 1.02, 1.05, 1.1)
    push!(reel, resume_reel(@sprintf("IRF ×%.2f", s), ajuster(inst_s(s)); s))
end
for d in (-1, 1)
    push!(reel, resume_reel(@sprintf("t₀ %+d canal (%+.1f ps)", d, 1000d * Δ), ajuster(m.inst; t0_fixe = ref.g.t₀ + d * Δ);
        t0_mode = @sprintf("fixé %+d canal", d)))
end
g_t0 = ajuster(m.inst; t0_par_histo = true)
push!(reel, resume_reel("t₀ libre par histogramme", g_t0; t0_mode = "par histogramme"))
σ_gigue = std(g_t0.t0s)
for d in (-1, 1)
    push!(reel, resume_reel(@sprintf("t₀ %+.1f ps (gigue observée)", 1000d * σ_gigue), ajuster(m.inst; t0_fixe = ref.g.t₀ + d * σ_gigue);
        t0_mode = @sprintf("fixé %+.4f ns", d * σ_gigue)))
end
reel = DataFrame(reel)
reel.delta_x0 = reel.x0 .- reel.x0[1]
reel.delta_x0E = reel.x0E .- reel.x0E[1]
CSV.write(joinpath(dossier, "reel_perturbations.csv"), reel)
println("\nDonnées réelles réajustées :")
for r in eachrow(reel)
    @printf("  %-30s χ²ᵣ=%7.3f  τ_D0=%.3f τ_DA=%.4f E=%.3f  x₀=%.4f (Δ %+.4f) x_sat=%.4f K=%.1f n=%.2f  x₀E=%.4f (Δ %+.4f)\n",
        r.perturbation, r.chi2_r, r.tau_D0, r.tau_DA, r.E, r.x0, r.delta_x0, r.x_sat, r.K, r.n, r.x0E, r.delta_x0E)
end

# -----------------------------------------------------------------------------
# 2. Profil de χ²ᵣ en τ_DA (données réelles, trois IRF)
# -----------------------------------------------------------------------------
profil = DataFrame()
for s in (0.9, 1.0, 1.1), τ in (0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.30, 0.50, 0.75, 1.0, 1.5, 2.1)
    g = ajuster(inst_s(s); tau_DA_fixe = τ)
    push!(profil, (s_irf = s, tau_DA = τ, chi2_r = g.chi2_global, tau_D0 = g.τ_D0, x0_brut = mean(g.x[C.==0])))
end
CSV.write(joinpath(dossier, "profil_tau_DA.csv"), profil)
println("\nProfil χ²ᵣ(τ_DA imposé), IRF ×1 :")
for r in eachrow(profil[profil.s_irf.==1.0, :])
    @printf("  τ_DA = %.2f ns : χ²ᵣ = %8.3f\n", r.tau_DA, r.chi2_r)
end

# -----------------------------------------------------------------------------
# 3. Hypothèses
# -----------------------------------------------------------------------------
p = m.p
E1 = 1 - p.τ_DA / p.τ_D0
τ_H4 = 1.0
E4 = 1 - τ_H4 / p.τ_D0
τ_m, τ_nu = 1.5, 4.1
H1 = S.Hypothese(cle = "H1", description = "modèle ajusté actuel (contrôle)", p = p)
H1g = S.Hypothese(cle = "H1-gigue", description = "H1 avec la gigue de t₀ observée entre échantillons", p = p, t0s = g_t0.t0s)
H2 = S.Hypothese(cle = "H2", description = "transfert vers un accepteur non émissif", p = p,
    etats = (x, q) -> ([S.EtatDonneur(q.τ_DA, :sombre), S.EtatDonneur(q.τ_D0, :protone)], [x, 1 - x]))
H2b = S.Hypothese(cle = "H2b", description = "extinction statique stricte (complexe non émissif)", p = p,
    etats = (x, q) -> ([S.EtatDonneur(q.τ_DA, :sombre, false), S.EtatDonneur(q.τ_D0, :protone)], [x, 1 - x]))
H3a(f; τm = τ_m, t0s = nothing) = S.Hypothese(cle = "H3a", description = @sprintf("3 états : %.0f %% de x à τ = %.1f ns", 100f, τm), p = p,
    etats = (x, q) -> ([S.EtatDonneur(q.τ_DA, :actif), S.EtatDonneur(τm, :actif), S.EtatDonneur(q.τ_D0, :protone)],
        [(1 - f) * x, f * x, 1 - x]), t0s = t0s)
H3b(f) = S.Hypothese(cle = "H3b", description = @sprintf("3 états : %.0f %% de donneur nu (τ = %.1f ns)", 100f, τ_nu), p = p,
    etats = (x, q) -> ([S.EtatDonneur(q.τ_DA, :actif), S.EtatDonneur(q.τ_D0, :protone), S.EtatDonneur(τ_nu, :nu)],
        [(1 - f) * x, (1 - f) * (1 - x), f]))
H4 = S.Hypothese(cle = "H4", description = @sprintf("τ_DA = %.1f ns, x × %.2f (même x·E)", τ_H4, E1 / E4),
    p = S.remplacer(p; τ_DA = τ_H4), x_de = x -> x * E1 / E4)

lignes = []
function ajoute!(hyp; graines = 1:3, etiquette = hyp.cle, kw...)
    for gr in graines
        t = @elapsed r = S.passer(hyp, ref, m; graine = gr, etiquette, inst_u, kw...)
        push!(lignes, r)
        @printf("  %-12s graine %d  s=%.2f  χ²ᵣ=%7.3f  τ_D0=%.3f (%.3f)  τ_DA=%.4f (%.4f)  x₀=%.3f (%.3f)  x₀E=%.3f (%.3f)  K=%.1f (%.1f)  biais ratio 20 mM %+.2f mM  [%.1f s]\n",
            etiquette, gr, r.s_irf_ajustement, r.chi2_r, r.tau_D0, r.tau_D0_vrai, r.tau_DA, r.tau_DA_vrai, r.x0, r.x0_vrai,
            r.x0E, r.x0E_vrai, r.K, r.K_vrai, r.biais_ratio_20mM, t)
    end
end
println("\nBatterie (valeur retrouvée, vraie entre parenthèses) :")
ajoute!(H1)
ajoute!(H1g)
ajoute!(H2)
ajoute!(H2b)
ajoute!(H3a(0.2))
ajoute!(H3b(0.1))
ajoute!(H4)
ajoute!(H4; etiquette = "H5 (IRF ×0,9)", s_irf = 0.9)
ajoute!(H4; etiquette = "H5 (IRF ×1,1)", s_irf = 1.1)
for s in (0.9, 0.95, 1.05, 1.1)
    ajoute!(H1; etiquette = @sprintf("H1 (IRF ×%.2f)", s), s_irf = s)
end
# rejet croisé : H4 ajusté avec τ_DA imposé à la valeur réelle, H1 avec τ_DA = 1,0 ns
ajoute!(H4; graines = 1:1, etiquette = "H4, τ_DA imposé 0,15", tau_DA_fixe = p.τ_DA)
ajoute!(H1; graines = 1:1, etiquette = "H1, τ_DA imposé 1,0", tau_DA_fixe = τ_H4)
# calibration du χ²ᵣ : durée de vie et fraction du troisième état
for τm in (0.8, 1.5, 2.5), f in (0.1, 0.2, 0.4)
    (τm == τ_m && f == 0.2) && continue
    ajoute!(H3a(f; τm); graines = 1:1, etiquette = @sprintf("H3a τ=%.1f %.0f %%", τm, 100f))
end
for f in (0.02, 0.05, 0.2)
    ajoute!(H3b(f); graines = 1:1, etiquette = @sprintf("H3b %.0f %%", 100f))
end
# combinaison : troisième état et gigue de t₀
ajoute!(H3a(0.3; t0s = g_t0.t0s); graines = 1:1, etiquette = "H3a 30 % + gigue")
bat = DataFrame(lignes)
χ2_ctrl = mean(bat.chi2_r[bat.jeu.=="H1"])
σχ = bat.sigma_chi2_poisson[1]
bat.verdict = [S.verdict(c, χ2_ctrl, σχ, ref.g.chi2_global) for c in bat.chi2_r]
# la ligne des données réelles, pour les figures (pas de vérité)
iso_r = S.ajuster_isotherme(C, ref.g.x)
E_r = 1 - ref.g.τ_DA / ref.g.τ_D0
reel_ligne = (jeu = "réel", hypothese = "réel", description = "titration du 2026-07-20", s_irf_ajustement = 1.0, graine = 0,
    chi2_r = ref.g.chi2_global, sigma_chi2_poisson = σχ, converge = ref.g.converge, tau_D0_vrai = NaN, tau_D0 = ref.g.τ_D0,
    tau_DA_vrai = NaN, tau_DA = ref.g.τ_DA, E_vrai = NaN, E = E_r, t0 = ref.g.t₀, x0_vrai = NaN, x0 = iso_r.x₀,
    x_sat_vrai = NaN, x_sat = iso_r.x_sat, K_vrai = NaN, K = iso_r.K_D, n_vrai = NaN, n = iso_r.n, x0E_vrai = NaN,
    x0E = iso_r.x₀ * E_r, xsatE_vrai = NaN, xsatE = iso_r.x_sat * E_r, biais_ratio_5mM = NaN, biais_ratio_20mM = NaN,
    biais_ratio_60mM = NaN, hors_plage_20mM = false, sig_reel..., verdict = "")
bat = vcat(DataFrame([reel_ligne]), bat; cols = :union)
bat.chi2_r_reel .= ref.g.chi2_global
CSV.write(joinpath(dossier, "batterie.csv"), bat)

println("\nRésumé (moyenne des graines) :")
res = combine(groupby(bat, :jeu, sort = false), :chi2_r => mean => :chi2_r, :tau_DA => mean => :tau_DA, :tau_DA_vrai => first => :tau_DA_vrai,
    :x0 => mean => :x0, :x0_vrai => first => :x0_vrai, :x0E => mean => :x0E, :x0E_vrai => first => :x0E_vrai,
    :biais_ratio_20mM => mean => :biais_20, :phasor_mediane => mean => :ph, :dtau_mono_0 => mean => :dt0,
    :dtau_mono_sat => mean => :dts, :verdict => first => :verdict)
for r in eachrow(res)
    @printf("  %-22s χ²ᵣ=%7.3f  τ_DA %.3f/%.3f  x₀ %.3f/%.3f  x₀E %.3f/%.3f  biais ratio 20 mM %+7.2f mM  phasor %.4f  δτ mono %+.3f/%+.3f  %s\n",
        r.jeu, r.chi2_r, r.tau_DA, r.tau_DA_vrai, r.x0, r.x0_vrai, r.x0E, r.x0E_vrai, r.biais_20, r.ph, r.dt0, r.dts, r.verdict)
end

# -----------------------------------------------------------------------------
# 4. t₀ libre par histogramme : données réelles et contrôle H1
# -----------------------------------------------------------------------------
hs_h1 = S.titration_synthetique(H1, ref, m, StableRNG(1))
g_h1 = S.ajuster_global(hs_h1, m.inst; τ_D0_init = ref.g.τ_D0, τ_DA_init = ref.g.τ_DA, t0_par_histo = true)
t0df = vcat(DataFrame(jeu = "réel", C_mM = C, dossier = [h.dossier for h in hs], replicat = [h.replicat for h in hs],
        t0_ns = g_t0.t0s, sigma_t0_ns = g_t0.σ_t0s, x = g_t0.x),
    DataFrame(jeu = "H1 (contrôle)", C_mM = C, dossier = [h.dossier for h in hs], replicat = [h.replicat for h in hs],
        t0_ns = g_h1.t0s, sigma_t0_ns = g_h1.σ_t0s, x = g_h1.x))
CSV.write(joinpath(dossier, "t0_par_histogramme.csv"), t0df)

"""Corrélation de Spearman (rangs moyens pour les ex æquo)."""
function spearman(a, b)
    rang(v) = (p = sortperm(v); r = similar(v, Float64); r[p] = 1:length(v);
        for u in unique(v); i = findall(==(u), v); r[i] .= mean(r[i]); end; r)
    return cor(rang(a), rang(b))
end
Cr = [isinf(c) ? 1e6 : c for c in C]
resume_t0 = Dict{String,Any}()
for (nomj, g) in (("réel", g_t0), ("H1 (contrôle)", g_h1))
    ρ = spearman(Cr, g.t0s)
    # gigue entre concentrations contre bruit entre réplicats
    groupes = [g.t0s[C.==c] for c in unique(C)]
    intra = sqrt(mean(var.(groupes)))
    inter = std(mean.(groupes))
    # permutations : probabilité d'un |ρ| au moins aussi grand sans tendance (concentrations permutées en bloc)
    rng = StableRNG(42)
    uC = unique(Cr)
    nperm = 20000
    plus = count(1:nperm) do _
        perm = Dict(zip(uC, uC[randperm(rng, length(uC))]))
        abs(spearman([perm[c] for c in Cr], g.t0s)) >= abs(ρ)
    end
    @printf("t₀ par histogramme, %-14s : χ²ᵣ %.3f → %.3f ; ρ de Spearman (t₀, [Cl⁻]) = %+.2f (p ≈ %.2f, permutations des 9 concentrations) ; dispersion entre concentrations %.1f ps, entre réplicats %.1f ps, σ statistique %.1f ps\n",
        nomj, nomj == "réel" ? ref.g.chi2_global : NaN, g.chi2_global, ρ, plus / nperm, 1000inter, 1000intra, 1000median(g.σ_t0s))
    resume_t0[nomj == "réel" ? "reel" : "H1"] = Dict("rho_spearman" => ρ, "p_permutation" => plus / nperm, "chi2_r" => g.chi2_global,
        "dispersion_inter_ps" => 1000inter, "dispersion_intra_ps" => 1000intra, "sigma_stat_ps" => 1000median(g.σ_t0s))
end
iso_t0 = S.ajuster_isotherme(C, g_t0.x)
resume = Dict{String,Any}(
    "date" => string(Dates.today()), "chi2_r_reel" => ref.g.chi2_global, "sigma_chi2_poisson" => σχ,
    "t0_par_histogramme" => resume_t0,
    "isotherme_t0_par_histogramme" => Dict("x0" => iso_t0.x₀, "x_sat" => iso_t0.x_sat, "K_D" => iso_t0.K_D, "n" => iso_t0.n,
        "tau_D0" => g_t0.τ_D0, "tau_DA" => g_t0.τ_DA, "chi2_r" => g_t0.chi2_global),
    "signatures_reelles" => Dict(String(k) => v for (k, v) in pairs(sig_reel)),
    "irf_fwhm_ps" => 1000fwhm(m.inst.irf, Δ),
)
open(io -> TOML.print(io, resume; sorted = true), joinpath(dossier, "batterie_resume.toml"), "w")

include(joinpath(@__DIR__, "src", "Figures.jl"))
using .Figures, CairoMakie
set_theme!(Figures.theme_sim())
if isdefined(Figures, :figure_batterie)
    Figures.figure_batterie(dossier, S.dossier_sorties("figures", "fig9_batterie_ingestion.png"))
    println("figure : ", S.dossier_sorties("figures", "fig9_batterie_ingestion.png"))
end
