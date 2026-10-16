# Étape 1 — ta titration bleue → paramètres du modèle, et le test de M2 lent.
#
# Pour chaque protéine de 2026-07-20_SCLM (3 réplicats × 9 concentrations,
# émission 470/40, IRF irf_20260710_16x_sclm.sdt comme dans analyse_SCLM_2.jl) :
#   1. ton ajustement mono-exponentiel (vec_to_lifetime) et son χ²_r, point par point ;
#   2. l'ajustement global M2 : τ_D0, τ_DA, t₀ liés sur les 27 histogrammes,
#      x, intensité et fond libres par histogramme ;
#   3. le même modèle avec τ libres par histogramme (coût du liage) ;
#   4. le χ²_r mono que M2 lent prédit, à comparer au tien ;
#   5. l'isotherme à deux plateaux sur les x (n libre, n = 1, sans « F ») ;
#   6. les phasors des données (sous M2 lent : un segment de droite).
# Sorties : sorties/titration/titration_<protéine>.csv (un histogramme par
# ligne) et parametres_<protéine>.toml (lu par parametres.toml).
#
#   julia --project=. step1_titration.jl [protéines...]

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, DataFrames, CSV, TOML, Dates

proteines = isempty(ARGS) ? ["mCer3_Topaz", "SCLM", "M2-SCLM"] : ARGS
dossier = S.dossier_sorties("titration")
mkpath(dossier)

cfg = TOML.parsefile(S.FICHIER_PARAMETRES)
f = cfg["instrument"]["fenetre"]
irf_mat = S.charger_irf_utilisateur(S.valeur_cfg(cfg["instrument"]["irf"]))
inst = S.InstrumentParams(irf = S.irf_normalisee(irf_mat), fenetre = f[1]:f[2])
inst_u = S.charger_instance_utilisateur(irf_mat)

for prot in proteines
    hs = S.lire_titration(prot)
    H = length(hs)
    println("\n===== $prot ($H histogrammes)")

    # 1. ton ajustement mono-exponentiel
    mono = [S.ajuster_mono_utilisateur(inst_u, h.n)[1] for h in hs]
    χ2_mono = [S.chi2_utilisateur(inst_u, mono[i], hs[i].n) for i in 1:H]

    # 2. ajustement global M2
    g = S.ajuster_global(hs, inst; τ_D0_init = 1.2 * maximum(p[1] for p in mono), τ_DA_init = 0.5)
    @printf("global : τ_D0 = %.4f ± %.4f ns, τ_DA = %.4f ± %.4f ns, t₀ = %.4f ns, χ²_r = %.3f (convergé : %s)\n",
        g.τ_D0, g.σ_glob[1], g.τ_DA, g.σ_glob[2], g.t₀, g.chi2_global, g.converge)

    # 3. τ libres par histogramme
    libres = [S.ajuster_libre(hs[i], inst, g.t₀; τ_init = (g.τ_D0, g.τ_DA), x_init = g.x[i]) for i in 1:H]

    # 4. χ² mono prédit par M2 lent
    pred = [S.chi2_mono_predit(inst_u, g.μ[i], inst) for i in 1:H]

    # 5. isothermes
    C = [h.C for h in hs]
    iso = S.ajuster_isotherme(C, g.x)
    iso1 = S.ajuster_isotherme(C, g.x; n_libre = false)
    sansF = isfinite.(C)
    isoF = S.ajuster_isotherme(C[sansF], g.x[sansF])
    for (nomiso, r) in (("n libre", iso), ("n = 1", iso1), ("sans F", isoF))
        @printf("isotherme %-8s : x₀ = %.4f, x_sat = %.4f, K_D = %.2f mM, n = %.3f\n", nomiso, r.x₀, r.x_sat, r.K_D, r.n)
    end

    # 6. phasors (fond ajusté retiré) et distance au segment du modèle
    Δ = S.pas_temps(inst)
    ws = S.EspaceTravail(inst)
    BD0, BDA = zeros(inst.n_bins), zeros(inst.n_bins)
    S.base!(BD0, ws, g.τ_D0, g.t₀)
    S.base!(BDA, ws, g.τ_DA, g.t₀)
    P_D0, P_DA = S.phasor(BD0, inst), S.phasor(BDA, inst)
    phas = map(1:H) do i
        h = copy(hs[i].n)
        h[inst.fenetre] .-= g.b[i]
        S.phasor(h, inst)
    end
    dist = [S.distance_droite(z, P_D0, P_DA) for z in phas]

    # 7. le test d'intensité (la « porte » de la note de changements) : acquisitions à
    #    temps fixé ? Si oui, les comptages totaux, fond retiré, doivent monter avec [Cl⁻]
    #    comme (1 − E·x(C))/(1 − E·x₀), l'amplitude du donneur à t = 0 restant constante.
    acq = [S.lire_acquisition(h.fichier) for h in hs]
    m_fen = length(inst.fenetre)
    N_net = [sum(h.n[inst.fenetre]) for h in hs] .- g.b .* m_fen
    E = 1 - g.τ_DA / g.τ_D0
    x_iso = S.isotherme(C, [iso.x₀, iso.x_sat, iso.K_D, iso.n])
    I_pred = (1 .- E .* x_iso) ./ (1 - E * iso.x₀)
    N0 = mean(N_net[C.==0])
    rapport_mes = mean(N_net[isinf.(C)]) / N0
    rapport_pred = (1 - E * iso.x_sat) / (1 - E * iso.x₀)
    @printf("acquisitions : durée %s s, arrêt sur temps : %s\n", join(unique(round.(first.(acq), digits = 3)), "/"), join(unique(last.(acq)), "/"))
    @printf("test d'intensité : comptages nets saturant / 0 mM mesurés ×%.3f, prédits (1 − E·x_sat)/(1 − E·x₀) = ×%.3f\n", rapport_mes, rapport_pred)

    df = DataFrame(proteine = prot, concentration = [h.dossier for h in hs], C_mM = C, replicat = [h.replicat for h in hs],
        duree_acquisition_s = first.(acq), arret_sur_temps = last.(acq),
        N_fenetre = [sum(h.n[inst.fenetre]) for h in hs], N_net = N_net, intensite_nette_sur_0mM = N_net ./ N0,
        intensite_predite_sur_0mM = I_pred,
        tau_mono_ns = [p[1] for p in mono], chi2_mono = χ2_mono, chi2_mono_predit_M2 = first.(pred),
        x = g.x, sigma_x = g.σ_x, intensite_relative = g.A ./ mean(g.A), fond_par_canal = g.b, chi2_global = g.chi2,
        tau_D0_libre = [l[1][1] for l in libres], tau_DA_libre = [l[1][2] for l in libres], x_libre = [l[1][3] for l in libres],
        chi2_libre = last.(libres), phasor_g = real.(phas), phasor_s = imag.(phas), distance_segment = dist)
    CSV.write(joinpath(dossier, "titration_$(prot).csv"), df)
    # données et modèle global, canal par canal, pour la figure de titration
    courbes = DataFrame(canal = collect(inst.fenetre), t_ns = (collect(inst.fenetre) .- 0.5) .* Δ)
    for i in 1:H
        etiq = replace(hs[i].dossier, " " => "") * "_c$(hs[i].replicat)"
        courbes[!, "donnees_" * etiq] = hs[i].n[inst.fenetre]
        courbes[!, "modele_" * etiq] = g.μ[i]
    end
    CSV.write(joinpath(dossier, "ajustement_global_$(prot).csv"), courbes)

    for i in 1:H
        @printf("  %-7s c%d  χ²mono=%6.2f (M2 prédit %6.2f)  x=%.3f  A/Ā=%.3f  χ²glob=%.2f  χ²libre=%.2f  τ libres=(%.3f, %.3f)\n",
            hs[i].dossier, hs[i].replicat, χ2_mono[i], pred[i][1], g.x[i], g.A[i] / mean(g.A), g.chi2[i], libres[i][2],
            libres[i][1][1], libres[i][1][2])
    end
    @printf("χ²_r moyen : global lié %.3f ; τ libres %.3f ; mono (ton code) %.3f\n", mean(g.chi2), mean(last.(libres)), mean(χ2_mono))
    @printf("distance des phasors au segment τ_D0–τ_DA : médiane %.4f, max %.4f (en longueur de segment)\n", median(dist), maximum(dist))

    out = Dict{String,Any}(
        "proteine" => prot, "date" => string(Dates.today()),
        "source" => "$(S.DOSSIER_TITRATION)/$prot, IRF $(S.valeur_cfg(cfg["instrument"]["irf"])), fenêtre $(f[1])–$(f[2])",
        "tau_D0" => g.τ_D0, "tau_D0_sigma" => g.σ_glob[1], "tau_DA" => g.τ_DA, "tau_DA_sigma" => g.σ_glob[2],
        "t0" => g.t₀, "chi2_r_global" => g.chi2_global, "chi2_r_libres_moyen" => mean(last.(libres)),
        "x0" => iso.x₀, "x_sat" => iso.x_sat, "K_D" => iso.K_D, "n" => iso.n,
        "isotherme_n1" => Dict("x0" => iso1.x₀, "x_sat" => iso1.x_sat, "K_D" => iso1.K_D),
        "isotherme_sans_F" => Dict("x0" => isoF.x₀, "x_sat" => isoF.x_sat, "K_D" => isoF.K_D, "n" => isoF.n),
        "intensite_relative_cv" => std(g.A) / mean(g.A),
        "acquisition_duree_s" => unique(first.(acq)), "acquisition_arret_sur_temps" => all(last.(acq)),
        "test_intensite" => Dict("mesure_sat_sur_0" => rapport_mes, "predit_sat_sur_0" => rapport_pred, "E" => E),
        "phasor_distance_mediane" => median(dist),
        "phasor_tau_D0" => [real(P_D0), imag(P_D0)], "phasor_tau_DA" => [real(P_DA), imag(P_DA)],
    )
    open(io -> TOML.print(io, out; sorted = true), joinpath(dossier, "parametres_$(prot).toml"), "w")
end
