# Étape 8 — figures de synthèse pour présenter les résultats (format de tes slides :
# 16:9, 2560 × 1440 px, français et anglais, fond clair et fond noir).
#
#   julia --project=. step8_synthese.jl
#
# Deux diapositives :
#   synthese_{FR,EN}{,_dark}.png           ce que dit le simulateur (étapes 1, 3, 4, 5, 7), avec les
#                                          bandes d'incertitude de l'étape 11 ;
#   synthese_revision_{FR,EN}{,_dark}.png  ce que la révision du 2026-10-15 a vérifié et changé
#                                          (étapes 1, 9, 10, 11, 13 et les nouveaux balayages de l'étape 5).
#
# Diapositive 1, six panneaux titrés par leur conclusion :
#   A  la titration : M2 lent, K vrai ;
#   B  le modèle retrouve tes K apparents (τ mono et ratio TIFF), plage des mondes IRF/t₀ ;
#   C  la borne de Cramér-Rao selon la donnée utilisée : l'information est dans le rapport J/B ;
#   D  σ[Cl⁻] selon le budget de photons, à 20 mM, bandes IRF/t₀ ;
#   E  biais ajouté par une dérive de gain du jaune et par un τ_D0 décalé ;
#   F  erreur statique en boucle fermée.
# Diapositive 2 :
#   A  τ_DA = 0,15 ns tient (profil de χ²ᵣ, H4 et H5) ;
#   B  à temps fixé, l'intensité suit le modèle ;
#   C  le facteur 8 sur les données réelles amincies ;
#   D  l'incertitude des paramètres ne change pas le classement ;
#   E  la fusion corrige la dérive de gain ;
#   F  chaque lecture craint une perturbation différente.

using CairoMakie, CSV, DataFrames, TOML, Printf, Statistics

const RACINE = joinpath(@__DIR__, "sorties")
lire(f) = CSV.read(joinpath(RACINE, f), DataFrame; types = Dict(:estimateur => String, :normalisation => String,
    :mode_calibration => String))
sel(df; kw...) = df[reduce(.&, [df[!, k] .== v for (k, v) in kw]; init = trues(nrow(df))), :]

prot = TOML.parsefile(joinpath(@__DIR__, "parametres.toml"))["capteur"]["proteine"]
tit = TOML.parsefile(joinpath(RACINE, "titration", "parametres_$(prot).toml"))
titr = CSV.read(joinpath(RACINE, "titration", "titration_$(prot).csv"), DataFrame)
kapp = CSV.read(joinpath(RACINE, "calibration", "K_apparents.csv"), DataFrame; types = Dict(:estimateur => String))
bornes = CSV.read(joinpath(RACINE, "calibration", "bornes.csv"), DataFrame; types = Dict(:normalisation => String))
grille = lire(joinpath("resultats", "resultats_principaux.csv"))
nuis = lire(joinpath("resultats", "nuisances.csv"))
mobile = CSV.read(joinpath(RACINE, "resultats", "cible_mobile.csv"), DataFrame)
# révision du 2026-10-15
mondes = lire(joinpath("incertitude", "mondes_grille.csv"))
mondes_K = CSV.read(joinpath(RACINE, "incertitude", "mondes_K_apparents.csv"), DataFrame; types = Dict(:estimateur => String, :monde => String))
profil = CSV.read(joinpath(RACINE, "batterie", "profil_tau_DA.csv"), DataFrame)
batterie = CSV.read(joinpath(RACINE, "batterie", "batterie.csv"), DataFrame; types = Dict(:jeu => String))
f8 = CSV.read(joinpath(RACINE, "reel", "facteur8.csv"), DataFrame; types = Dict(:source => String))
fusion = CSV.read(joinpath(RACINE, "resultats", "fusion_resume.csv"), DataFrame; types = Dict(:capteur => String, :filtre => String))

# -----------------------------------------------------------------------------
# Thèmes (mêmes valeurs que tes slides : slide_resultats.py) et textes
# -----------------------------------------------------------------------------
const THEMES = Dict(
    "" => (SURFACE = "#fcfcfb", INK = "#0b0b0b", INK_2 = "#52514e", MUET = "#898781", GRID = "#e6e6e2", AXE = "#c3c2b7",
        ZONE = "#f0efec", GRIS = "#c3c2b7",
        c = Dict("1" => "#2a78d6", "5d" => "#008300", "5a" => "#e87ba4", "4" => "#eda100", "2H" => "#eb6834"),
        N3 = "#86b6ef", N4 = "#184f95"),
    "_dark" => (SURFACE = "#000000", INK = "#ffffff", INK_2 = "#c3c2b7", MUET = "#898781", GRID = "#2f2f2c", AXE = "#52514e",
        ZONE = "#1a1a19", GRIS = "#52514e",
        c = Dict("1" => "#3987e5", "5d" => "#008300", "5a" => "#d55181", "4" => "#c98500", "2H" => "#d95926"),
        N3 = "#86b6ef", N4 = "#3987e5"),
)

const TEXTES = Dict(
    "FR" => Dict(
        :titre => "Simulateur FLIM-FRET à deux canaux : lire [Cl⁻] sur mCer3–Topaz",
        :sous_titre => "Modèle calé sur la titration du 2026-07-20 · Monte-Carlo de Poisson avec l'IRF mesurée, l'ajustement mono-exponentiel actuel et ses fenêtres · 28/28 tests de validation",
        :A => "A   La titration confirme le mélange statique", :A_sub => "fraction FRET-active x, 27 histogrammes, τ liés (modèle M2 lent)",
        :A_x => "[Cl⁻] (mM)", :A_y => "fraction active x",
        :A_txt => "τ_D0 = %.2f ns, τ_DA = %.2f ns (E = %.2f)\nK vrai = %.1f mM, n = %.2f ; χ²ᵣ = %.2f",
        :B => "B   Le modèle retrouve les K apparents mesurés", :B_sub => "prédits depuis le K vrai, contre les ajustements de Hill expérimentaux",
        :B_x => "K apparent (mM)", :B_rows => ["x (lecture directe)", "τ mono-exponentiel", "ratio TIFF C2/C1 (B/J)"],
        :B_pred => "prédit", :B_mes => "mesuré", :B_vrai => "K vrai",
        :C => "C   L'information est dans le rapport J/B", :C_sub => "σ[Cl⁻] minimal selon la donnée utilisée (20 mM, 10⁴ photons)",
        :C_x => "borne de Cramér-Rao (mM)",
        :C_rows => ["tout (J/B + forme)", "totaux J/B seuls", "forme seule, 2 canaux", "canal bleu, τ connus", "canal bleu, τ libres"],
        :C_qui => ["1, 5d", "1", "5a", "4", "3"], :C_atteint => "atteinte par ",
        :D => "D   Le ratio atteint la borne", :D_sub => "σ[Cl⁻] à 20 mM selon le budget (temps d'acquisition fixé)",
        :D_x => "photons bleus", :D_y => "σ[Cl⁻] (mM)", :D_pente => "pente −½",
        :E => "E   Chaque lecture craint une perturbation différente", :E_sub => "biais ajouté si la tranche diffère de la calibration (20 mM, 10⁴ photons)",
        :E1_x => "gain du canal jaune (%)", :E2_x => "τ_D0 tranche / calibration (%)", :E_y => "biais (mM)",
        :F => "F   En boucle fermée, il reste le biais", :F_sub => "[Cl⁻] vraie − consigne (20 mM), régime établi, 10 Hz",
        :F_y => "erreur statique (mM)", :F_N => "photons / échantillon",
        :leg => Dict("1" => "1  ratio J/B", "5d" => "5d  conjoint, échelle commune", "5a" => "5a  conjoint, forme seule",
            "4" => "4  bi-exp., τ figés (bleu)", "2H" => "pipeline actuel (mono + Hill)", "borne" => "borne de Cramér-Rao"),
        :F_cats => ["1", "5d", "5a", "4", "pipeline\nactuel"],
        :pied => "mCer3_Topaz · normalisation A (N = photons bleus attendus à 20 mM ; le jaune s'y ajoute) · 1 000 à 2 000 réalisations par point · estimations ramenées à [0, 150] mM · barres et bandes grises : titration réajustée avec l'IRF ×0,9 et ×1,1, t₀ ± 3,5 ps, t₀ par histogramme",
    ),
    "EN" => Dict(
        :titre => "Two-channel FLIM-FRET simulator: reading [Cl⁻] from mCer3–Topaz",
        :sous_titre => "Model fitted to the 2026-07-20 titration · Poisson Monte Carlo with the measured IRF, the lab's mono-exponential fit and gates · 28/28 validation tests",
        :A => "A   The titration supports a static mixture", :A_sub => "FRET-active fraction x, 27 histograms, linked lifetimes (slow-exchange M2)",
        :A_x => "[Cl⁻] (mM)", :A_y => "active fraction x",
        :A_txt => "τ_D0 = %.2f ns, τ_DA = %.2f ns (E = %.2f)\ntrue K = %.1f mM, n = %.2f; χ²ᵣ = %.2f",
        :B => "B   The model predicts the measured apparent K", :B_sub => "predicted from the true K, vs. the lab's Hill fits",
        :B_x => "apparent K (mM)", :B_rows => ["x (direct readout)", "mono-exp. τ", "TIFF ratio C2/C1 (B/Y)"],
        :B_pred => "predicted", :B_mes => "measured", :B_vrai => "true K",
        :C => "C   The information is in the Y/B ratio", :C_sub => "minimal σ[Cl⁻] by data used (20 mM, 10⁴ photons)",
        :C_x => "Cramér-Rao bound (mM)",
        :C_rows => ["all (Y/B + shape)", "Y/B totals only", "shape only, 2 channels", "blue channel, known τ", "blue channel, free τ"],
        :C_qui => ["1, 5d", "1", "5a", "4", "3"], :C_atteint => "reached by ",
        :D => "D   The ratio reaches the bound", :D_sub => "σ[Cl⁻] at 20 mM vs. photon budget (fixed acquisition time)",
        :D_x => "blue photons", :D_y => "σ[Cl⁻] (mM)", :D_pente => "slope −½",
        :E => "E   Each readout fears a different perturbation", :E_sub => "added bias when the slice differs from the calibration (20 mM, 10⁴ photons)",
        :E1_x => "yellow-channel gain (%)", :E2_x => "τ_D0 slice / calibration (%)", :E_y => "bias (mM)",
        :F => "F   In closed loop, only the bias remains", :F_sub => "true [Cl⁻] − set point (20 mM), steady state, 10 Hz",
        :F_y => "steady-state error (mM)", :F_N => "photons / sample",
        :leg => Dict("1" => "1  Y/B ratio", "5d" => "5d  joint, common scale", "5a" => "5a  joint, shape only",
            "4" => "4  bi-exp., fixed τ (blue)", "2H" => "current pipeline (mono + Hill)", "borne" => "Cramér-Rao bound"),
        :F_cats => ["1", "5d", "5a", "4", "current\npipeline"],
        :pied => "mCer3_Topaz · normalization A (N = expected blue photons at 20 mM; the yellow channel comes on top) · 1,000–2,000 realizations per point · estimates clipped to [0, 150] mM · grey bars and bands: titration refitted with the IRF ×0.9 and ×1.1, t₀ ± 3.5 ps, t₀ per histogram",
    ),
)

hexa(s) = Makie.to_color(s)
virgule(x; d = 1) = replace(@sprintf("%.*f", d, x), "." => ",")

function graduations_log(lo, hi, langue)
    v = Float64[]
    for e in floor(Int, log10(lo)):ceil(Int, log10(hi)), m in (1, 2, 5)
        x = m * 10.0^e
        lo <= x <= hi && push!(v, x)
    end
    fmt(x) = x >= 1 ? @sprintf("%.0f", x) : (langue == "FR" ? replace(@sprintf("%g", x), "." => ",") : @sprintf("%g", x))
    return (v, fmt.(v))
end
puissances = ([1e2, 1e3, 1e4, 1e5], ["10²", "10³", "10⁴", "10⁵"])

# -----------------------------------------------------------------------------
# Données des panneaux
# -----------------------------------------------------------------------------
isotherme(C) = tit["x_sat"] + (tit["x0"] - tit["x_sat"]) * tit["K_D"]^tit["n"] / (tit["K_D"]^tit["n"] + C^tit["n"])
kval(k) = kapp[findfirst(==(k), kapp.estimateur), :K_app]
b20 = sel(bornes; normalisation = "A", N_photons = 1e4, C_mM = 20.0)
bornes_C = [b20.absolue[1], b20.totaux[1], b20.conjoint[1], b20.bleu[1], b20.bleu_tau_libres[1]]
d20 = sel(grille; C_vrai = 20.0, normalisation = "A", mode_calibration = "exacte")
nb(p, v, e) = (r = sel(nuis; parametre_nuisance = p, C_vrai = 20.0, N_photons = 1e4, mode_calibration = "nominale", estimateur = e);
               r = r[abs.(r.valeur_nuisance .- v).<1e-9, :]; nrow(r) == 0 ? NaN : r.biais_mM[1])
base(e) = nb("eta_YD_relatif", 1.0, e)          # monde nominal, mêmes germes
function serie_nuisance(p, vals, e; e_tau = e)
    y = [v == 1.0 ? 0.0 : nb(p, v, e_tau) - base(e) for v in vals]
    return 100 .* (vals .- 1), y
end
function erreur_statique(N, k)
    d = mobile[(mobile.boucle.=="fermée").&(mobile.N.==N).&isfinite.(mobile[!, "est_"*k]).&(mobile.t.>=120), :]
    return mean(d.C_vrai) - 20.0
end

# -----------------------------------------------------------------------------
function synthese(langue, suffixe)
    T, X = THEMES[suffixe], TEXTES[langue]
    C = Dict(k => hexa(v) for (k, v) in T.c)
    INK, INK2, MUET, GRID, AXE = hexa.((T.INK, T.INK_2, T.MUET, T.GRID, T.AXE))
    set_theme!(Theme(fontsize = 15, backgroundcolor = hexa(T.SURFACE), textcolor = INK,
        Axis = (backgroundcolor = hexa(T.SURFACE), xgridcolor = GRID, ygridcolor = GRID, leftspinecolor = AXE,
            bottomspinecolor = AXE, topspinevisible = false, rightspinevisible = false, xtickcolor = AXE, ytickcolor = AXE,
            xticklabelcolor = INK2, yticklabelcolor = INK2, xlabelcolor = INK2, ylabelcolor = INK2, titlecolor = INK,
            titlealign = :left, titlesize = 18, titlefont = :bold, subtitlecolor = INK2, subtitlesize = 13,
            xminorticksvisible = false, yminorticksvisible = false),
        Legend = (framevisible = false, labelcolor = INK2, backgroundcolor = hexa(T.SURFACE)),
        Lines = (linewidth = 2.5,), Scatter = (markersize = 11, strokewidth = 1.5, strokecolor = hexa(T.SURFACE))))
    fig = Figure(size = (1600, 900))
    vf(x; d = 1) = langue == "FR" ? virgule(x; d) : @sprintf("%.*f", d, x)

    # en-tête
    Label(fig[0, 1:3], X[:titre]; fontsize = 26, font = :bold, color = INK, halign = :left, tellwidth = false)
    Label(fig[1, 1:3], X[:sous_titre]; fontsize = 14, color = INK2, halign = :left, tellwidth = false)

    # --- A : titration ---------------------------------------------------------
    pos(Cv) = Cv == 0 ? log10(0.8) : isinf(Cv) ? log10(400.0) : log10(Cv)
    axA = Axis(fig[2, 1]; title = X[:A], subtitle = X[:A_sub], xlabel = X[:A_x], ylabel = X[:A_y],
        xticks = (pos.([0, 2.5, 10, 40, 146, Inf]), ["0", langue == "FR" ? "2,5" : "2.5", "10", "40", "146", "F"]))
    Cf = exp10.(range(log10(0.8), log10(400); length = 300))
    lines!(axA, log10.(Cf), isotherme.(Cf); color = C["1"])
    scatter!(axA, pos.(titr.C_mM), titr.x; color = INK, markersize = 7, strokewidth = 0)
    vlines!(axA, [log10(tit["K_D"])]; color = AXE, linewidth = 1)
    texte_A = Printf.format(Printf.Format(X[:A_txt]), tit["tau_D0"], tit["tau_DA"], 1 - tit["tau_DA"] / tit["tau_D0"],
        tit["K_D"], tit["n"], tit["chi2_r_global"])
    langue == "FR" && (texte_A = replace(texte_A, r"(\d)\.(\d)" => s"\1,\2"))
    text!(axA, log10(0.85), 0.07; text = texte_A,
        color = INK2, fontsize = 13, align = (:left, :bottom))

    # --- B : K apparents ---------------------------------------------------------
    axB = Axis(fig[2, 2]; title = X[:B], subtitle = X[:B_sub], xlabel = X[:B_x], yticks = (3:-1:1, X[:B_rows]),
        ygridvisible = false)
    yB = [3, 2, 1]
    Kv = tit["K_D"]
    pred = [Kv, kval("2"), kval("1inv")]
    mes = [NaN, kval("mesure_tau"), kval("mesure_fret")]
    vlines!(axB, [Kv]; color = AXE, linewidth = 1)
    for i in 1:3
        isnan(mes[i]) || lines!(axB, [pred[i], mes[i]], [yB[i], yB[i]]; color = GRID, linewidth = 6)
    end
    # plage des mondes IRF/t₀ (chaque monde a son propre K vrai)
    for (i, k) in ((1, ""), (2, "2H"), (3, "1inv"))
        isempty(k) && continue
        d = mondes_K[mondes_K.estimateur.==k, :]
        lo, hi = min(minimum(d.K_app), pred[i]), max(maximum(d.K_app), pred[i])
        lines!(axB, [lo, hi], [yB[i] + 0.22, yB[i] + 0.22]; color = hexa(T.GRIS), linewidth = 4)
    end
    dK = vcat(Kv, unique(mondes_K.K_vrai))
    lines!(axB, [minimum(dK), maximum(dK)], [yB[1] + 0.22, yB[1] + 0.22]; color = hexa(T.GRIS), linewidth = 4)
    scatter!(axB, pred, yB; color = hexa(T.SURFACE), strokecolor = C["1"], strokewidth = 2.5, markersize = 16, label = X[:B_pred])
    scatter!(axB, mes, yB; color = INK, marker = :diamond, markersize = 15, strokewidth = 0, label = X[:B_mes])
    for i in 1:3
        txt = isnan(mes[i]) ? vf(pred[i]) : vf(pred[i]) * " / " * vf(mes[i])
        text!(axB, max(pred[i], isnan(mes[i]) ? 0 : mes[i]), yB[i]; text = "  " * txt * " mM", color = INK2, fontsize = 13,
            align = (:left, :center))
    end
    text!(axB, Kv, 3.55; text = " " * X[:B_vrai], color = MUET, fontsize = 12, align = (:left, :center))
    xlims!(axB, 0, 70)
    ylims!(axB, 0.4, 3.7)
    axislegend(axB; position = :rt, labelsize = 13)

    # --- C : bornes selon la donnée ----------------------------------------------
    axC = Axis(fig[2, 3]; title = X[:C], subtitle = X[:C_sub], xlabel = X[:C_x], yticks = (1:5, X[:C_rows]),
        ygridvisible = false, yreversed = true)
    couleurs_C = [C["1"], C["1"], hexa(T.GRIS), hexa(T.GRIS), hexa(T.GRIS)]
    barplot!(axC, 1:5, bornes_C; direction = :x, color = couleurs_C, gap = 0.35, strokewidth = 0)
    for (i, v) in enumerate(bornes_C)
        text!(axC, v, i; text = "  " * vf(v; d = 2) * " mM · " * X[:C_atteint] * X[:C_qui][i], color = INK2, fontsize = 12,
            align = (:left, :center))
    end
    xlims!(axC, 0, 1.75 * maximum(bornes_C))

    # --- D : σ selon N ----------------------------------------------------------------
    axD = Axis(fig[3, 1]; title = X[:D], subtitle = X[:D_sub], xlabel = X[:D_x], ylabel = X[:D_y], xscale = log10, yscale = log10,
        xticks = puissances, yticks = graduations_log(0.2, 100, langue))
    for k in ("2H", "4", "5a", "1")
        dm = mondes[(mondes.C_vrai.==20.0).&(mondes.normalisation.=="A").&(mondes.estimateur.==k), :]
        s0 = sel(d20; estimateur = k)
        g = combine(groupby(vcat(select(dm, :N_photons, :sigma_mM), select(s0, :N_photons, :sigma_mM)), :N_photons),
            :sigma_mM => minimum => :lo, :sigma_mM => maximum => :hi)
        sort!(g, :N_photons)
        band!(axD, g.N_photons, g.lo, g.hi; color = (C[k], 0.18))
    end
    for k in ("2H", "4", "5a", "5d", "1")
        s = sort(sel(d20; estimateur = k), :N_photons)
        lines!(axD, s.N_photons, s.sigma_mM; color = C[k], label = X[:leg][k])
        scatter!(axD, s.N_photons, s.sigma_mM; color = C[k], label = X[:leg][k], markersize = 9)
    end
    b = sort(sel(d20; estimateur = "1"), :N_photons)
    lines!(axD, b.N_photons, b.CRLB_mM; color = INK, linewidth = 1.5, linestyle = :dash, label = X[:leg]["borne"])
    ylims!(axD, 0.12, 120)

    # --- E : robustesse -------------------------------------------------------------
    gE = fig[3, 2] = GridLayout()
    axE1 = Axis(gE[1, 1]; xlabel = X[:E1_x], ylabel = X[:E_y], title = X[:E], subtitle = X[:E_sub])
    axE2 = Axis(gE[1, 2]; xlabel = X[:E2_x], title = " ", subtitle = " ")
    gain = [0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]
    tauD = [0.90, 0.95, 1.0, 1.05, 1.10]
    for ax in (axE1, axE2)
        hlines!(ax, [0.0]; color = AXE, linewidth = 1)
        vlines!(ax, [0.0]; color = AXE, linewidth = 1)
    end
    for (k, kτ) in (("4", "4"), ("5a", "5b"), ("1", "1"))
        x, y = serie_nuisance("gain_jaune_relatif", gain, k)
        k == "4" || (lines!(axE1, x, y; color = C[k]); scatter!(axE1, x, y; color = C[k], markersize = 8))
        x, y = serie_nuisance("tau_D0_relatif", tauD, k; e_tau = kτ)
        lines!(axE2, x, y; color = C[k])
        scatter!(axE2, x, y; color = C[k], markersize = 8)
    end
    linkyaxes!(axE1, axE2)
    hideydecorations!(axE2; grid = false, ticks = false)
    ylims!(axE1, -22, 47)
    colgap!(gE, 12)

    # --- F : boucle fermée --------------------------------------------------------------
    axF = Axis(fig[3, 3]; title = X[:F], subtitle = X[:F_sub], ylabel = X[:F_y],
        xticks = (1:5, X[:F_cats]), xgridvisible = false)
    hlines!(axF, [0.0]; color = AXE, linewidth = 1)
    ks = ["1", "5d", "5a", "4", "2H"]
    for (j, (N, col, lab)) in enumerate(((1e3, T.N3, "10³"), (1e4, T.N4, "10⁴")))
        v = [erreur_statique(N, k) for k in ks]
        barplot!(axF, collect(1:5) .+ (j - 1.5) * 0.36, v; width = 0.32, color = hexa(col), strokewidth = 0,
            label = lab * " " * X[:F_N])
    end
    axislegend(axF; position = :lb, labelsize = 13)

    # légende commune (D, E) et pied
    elems = [[LineElement(color = C[k], linewidth = 2.5), MarkerElement(color = C[k], marker = :circle, markersize = 9,
        strokecolor = hexa(T.SURFACE), strokewidth = 1.5)] for k in ("1", "5d", "5a", "4", "2H")]
    push!(elems, [LineElement(color = INK, linewidth = 1.5, linestyle = :dash)])
    Legend(fig[4, 1:3], elems, [X[:leg][k] for k in ("1", "5d", "5a", "4", "2H", "borne")];
        orientation = :horizontal, tellwidth = false, labelsize = 14, colgap = 26)
    Label(fig[5, 1:3], X[:pied]; fontsize = 12, color = MUET, halign = :left, tellwidth = false)

    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 34)
    chemin = joinpath(RACINE, "figures", "synthese_$(langue)$(suffixe).png")
    save(chemin, fig; px_per_unit = 1.6)
    return chemin
end

# -----------------------------------------------------------------------------
# Diapositive 2 — la révision du 2026-10-15
# -----------------------------------------------------------------------------
const TEXTES2 = Dict(
    "FR" => Dict(
        :titre => "Révision du simulateur : ce qui tient, ce qui change",
        :sous_titre => "Hypothèses concurrentes passées dans l'ingestion, test d'intensité, facteur 8 sur les données réelles, incertitude des paramètres, fusion, décalage et blanchiment · 28/28 tests",
        :A => "A   τ_DA = 0,15 ns tient", :A_sub => "χ²ᵣ des 27 histogrammes réels, τ_DA imposé",
        :A_x => "τ_DA imposé (ns)", :A_y => "χ²ᵣ", :A_irf => ["IRF mesurée", "IRF × 0,9", "IRF × 1,1"],
        :A_txt => "τ_DA = 1 ns : χ²ᵣ = %.1f\nsimulé à 1 ns, retrouvé à %.2f ns\n(%.2f–%.2f avec une IRF fausse de 10 %%)",
        :B => "B   À temps fixé, l'intensité suit le modèle", :B_sub => "relatif à 0 mM, acquisitions de 10 s",
        :B_x => "[Cl⁻] (mM)", :B_y => "relatif à 0 mM",
        :B_leg => ["comptages nets mesurés", "prédits : (1 − E·x)/(1 − E·x₀)", "amplitude du donneur à t = 0"],
        :B_txt => "×%.2f mesuré, ×%.2f prédit",
        :C => "C   Le facteur 8 se vérifie sur les données", :C_sub => "σ[Cl⁻] à 20 mM · ● réel aminci   ○ jumeau simulé   - - prédiction",
        :C_x => "photons par histogramme", :C_y => "σ[Cl⁻] (mM)",
        :C_leg => ["4R  bi-exp., τ figés", "pipeline actuel (mono + Hill)", "données réelles", "jumeau simulé", "prédiction du simulateur"],
        :C_txt => "×%.1f à 10⁴ photons\n(jumeau ×%.1f, prédit ×%.1f)",
        :D => "D   L'incertitude ne change pas le classement", :D_sub => "σ[Cl⁻] à 20 mM, 10⁴ photons",
        :D_x => "σ[Cl⁻] (mM)", :D_leg => ["bruit de photons (nominal)", "mondes IRF ±10 %, t₀", "photons + paramètres ajustés"],
        :E => "E   La fusion corrige la dérive de gain", :E_sub => "RMSE sur [Cl⁻], 10 Hz, 10⁴ photons par échantillon",
        :E_x => "échelle de temps de la dérive de gain", :E_y => "RMSE (mM)",
        :E_leg => ["ratio seul", "forme seule (5a)", "fusion ([Cl⁻], gain)"],
        :F => "F   Aucune lecture ne résiste à tout", :F_sub => "biais ajouté à 20 mM, 10⁴ photons, calibration non corrigée",
        :F_x => "biais ajouté (mM)",
        :F_rows => ["gain du jaune +5 %", "retard du jaune +20 ps", "accepteur blanchi 10 %,\nnon absorbant", "accepteur blanchi 10 %,\nextincteur sombre",
            "η_YD réel = 0,75 × supposé", "τ_D0 en tranche +5 %"],
        :F_leg => ["1  ratio J/B", "5a  forme, 2 canaux", "4  bi-exp., τ figés (bleu)"],
        :pied => "mCer3_Topaz · titration du 2026-07-20 · normalisation A · H2 (accepteur non émissif) est indiscernable par le canal bleu : seul le canal jaune le tranchera · la gigue de t₀ entre échantillons (3,6 ps) et un état FRET intermédiaire expliquent l'essentiel du χ²ᵣ",
        :tg => ["1 min", "5 min", "20 min", "1 h"],
    ),
    "EN" => Dict(
        :titre => "Simulator revision: what holds, what changes",
        :sous_titre => "Competing hypotheses through the ingestion, intensity test, factor 8 on real data, parameter uncertainty, fusion, channel delay and photobleaching · 28/28 tests",
        :A => "A   τ_DA = 0.15 ns holds", :A_sub => "χ²ᵣ of the 27 real histograms, τ_DA imposed",
        :A_x => "imposed τ_DA (ns)", :A_y => "χ²ᵣ", :A_irf => ["measured IRF", "IRF × 0.9", "IRF × 1.1"],
        :A_txt => "τ_DA = 1 ns: χ²ᵣ = %.1f\nsimulated at 1 ns, recovered at %.2f ns\n(%.2f–%.2f with a 10%% wrong IRF)",
        :B => "B   At fixed time, intensity follows the model", :B_sub => "relative to 0 mM, 10 s acquisitions",
        :B_x => "[Cl⁻] (mM)", :B_y => "relative to 0 mM",
        :B_leg => ["measured net counts", "predicted: (1 − E·x)/(1 − E·x₀)", "donor amplitude at t = 0"],
        :B_txt => "×%.2f measured, ×%.2f predicted",
        :C => "C   The factor 8 holds on real data", :C_sub => "σ[Cl⁻] at 20 mM · ● thinned real   ○ simulated twin   - - prediction",
        :C_x => "photons per histogram", :C_y => "σ[Cl⁻] (mM)",
        :C_leg => ["4R  bi-exp., fixed τ", "current pipeline (mono + Hill)", "real data", "simulated twin", "simulator prediction"],
        :C_txt => "×%.1f at 10⁴ photons\n(twin ×%.1f, predicted ×%.1f)",
        :D => "D   Uncertainty does not change the ranking", :D_sub => "σ[Cl⁻] at 20 mM, 10⁴ photons",
        :D_x => "σ[Cl⁻] (mM)", :D_leg => ["photon noise (nominal)", "IRF ±10%, t₀ worlds", "photons + fitted parameters"],
        :E => "E   Fusion corrects the gain drift", :E_sub => "RMSE on [Cl⁻], 10 Hz, 10⁴ photons per sample",
        :E_x => "gain drift time scale", :E_y => "RMSE (mM)",
        :E_leg => ["ratio only", "shape only (5a)", "fusion ([Cl⁻], gain)"],
        :F => "F   No readout resists everything", :F_sub => "added bias at 20 mM, 10⁴ photons, uncorrected calibration",
        :F_x => "added bias (mM)",
        :F_rows => ["yellow gain +5%", "yellow delay +20 ps", "10% acceptors bleached,\nnon-absorbing", "10% acceptors bleached,\ndark quencher",
            "true η_YD = 0.75 × assumed", "slice τ_D0 +5%"],
        :F_leg => ["1  Y/B ratio", "5a  shape, 2 channels", "4  bi-exp., fixed τ (blue)"],
        :pied => "mCer3_Topaz · 2026-07-20 titration · normalization A · H2 (non-emissive acceptor) is indistinguishable in the blue channel: only the yellow channel will settle it · sample-to-sample t₀ jitter (3.6 ps) and an intermediate FRET state explain most of χ²ᵣ",
        :tg => ["1 min", "5 min", "20 min", "1 h"],
    ),
)

nref(e) = (r = sel(nuis; parametre_nuisance = "reference", C_vrai = 20.0, N_photons = 1e4, mode_calibration = "nominale", estimateur = e);
           nrow(r) == 0 ? NaN : r.biais_mM[1])
ajoute(p, v, e) = nb(p, v, e) - nref(e == "5b" ? "5a" : e)

function synthese_revision(langue, suffixe)
    T, X = THEMES[suffixe], TEXTES2[langue]
    C = Dict(k => hexa(v) for (k, v) in T.c)
    INK, INK2, MUET, GRID, AXE = hexa.((T.INK, T.INK_2, T.MUET, T.GRID, T.AXE))
    SURF = hexa(T.SURFACE)
    set_theme!(Theme(fontsize = 15, backgroundcolor = SURF, textcolor = INK,
        Axis = (backgroundcolor = SURF, xgridcolor = GRID, ygridcolor = GRID, leftspinecolor = AXE,
            bottomspinecolor = AXE, topspinevisible = false, rightspinevisible = false, xtickcolor = AXE, ytickcolor = AXE,
            xticklabelcolor = INK2, yticklabelcolor = INK2, xlabelcolor = INK2, ylabelcolor = INK2, titlecolor = INK,
            titlealign = :left, titlesize = 18, titlefont = :bold, subtitlecolor = INK2, subtitlesize = 13,
            xminorticksvisible = false, yminorticksvisible = false),
        Legend = (framevisible = false, labelcolor = INK2, backgroundcolor = SURF),
        Lines = (linewidth = 2.5,), Scatter = (markersize = 11, strokewidth = 1.5, strokecolor = SURF)))
    fig = Figure(size = (1600, 900))
    vf(x; d = 1) = langue == "FR" ? virgule(x; d) : @sprintf("%.*f", d, x)
    fr(s) = langue == "FR" ? replace(s, r"(\d)\.(\d)" => s"\1,\2") : s
    Label(fig[0, 1:3], X[:titre]; fontsize = 26, font = :bold, color = INK, halign = :left, tellwidth = false)
    Label(fig[1, 1:3], X[:sous_titre]; fontsize = 14, color = INK2, halign = :left, tellwidth = false)

    # --- A : profil de χ²ᵣ en τ_DA ------------------------------------------------
    tA = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    axA = Axis(fig[2, 1]; title = X[:A], subtitle = X[:A_sub], xlabel = X[:A_x], ylabel = X[:A_y], xscale = log10,
        xticks = (tA, [fr(@sprintf("%g", t)) for t in tA]))
    for (sv, st, lw, lab) in ((0.9, :dash, 1.5, X[:A_irf][2]), (1.1, :dot, 1.5, X[:A_irf][3]), (1.0, :solid, 2.5, X[:A_irf][1]))
        d = sort(profil[profil.s_irf.==sv, :], :tau_DA)
        lines!(axA, d.tau_DA, d.chi2_r; color = sv == 1 ? C["1"] : MUET, linestyle = st, linewidth = lw, label = lab)
    end
    h4 = batterie[batterie.jeu.=="H4", :]
    h5 = batterie[startswith.(batterie.jeu, "H5"), :]
    p1 = profil[(profil.s_irf.==1.0).&(profil.tau_DA.==1.0), :chi2_r][1]
    text!(axA, 0.05, 11.5; text = fr(Printf.format(Printf.Format(X[:A_txt]), p1, mean(h4.tau_DA), minimum(h5.tau_DA), maximum(h5.tau_DA))),
        color = INK2, fontsize = 13, align = (:left, :top))
    ylims!(axA, 1.5, 12)
    axislegend(axA; position = :rb, labelsize = 12)

    # --- B : test d'intensité --------------------------------------------------------
    pos(Cv) = Cv == 0 ? log10(0.8) : isinf(Cv) ? log10(400.0) : log10(Cv)
    axB = Axis(fig[2, 2]; title = X[:B], subtitle = X[:B_sub], xlabel = X[:B_x], ylabel = X[:B_y],
        xticks = (pos.([0, 2.5, 10, 40, 146, Inf]), ["0", fr("2.5"), "10", "40", "146", "F"]))
    o = sortperm(pos.(titr.C_mM))
    a0 = mean(titr.intensite_relative[titr.C_mM.==0])
    hlines!(axB, [1.0]; color = AXE, linewidth = 1)
    scatter!(axB, pos.(titr.C_mM), titr.intensite_relative ./ a0; color = SURF, strokecolor = MUET, strokewidth = 1.5, markersize = 9,
        label = X[:B_leg][3])
    lines!(axB, pos.(titr.C_mM)[o], titr.intensite_predite_sur_0mM[o]; color = C["1"], label = X[:B_leg][2])
    scatter!(axB, pos.(titr.C_mM), titr.intensite_nette_sur_0mM; color = INK, markersize = 8, strokewidth = 0, label = X[:B_leg][1])
    ti = tit["test_intensite"]
    text!(axB, pos(0.0), 1.85; text = fr(Printf.format(Printf.Format(X[:B_txt]), ti["mesure_sat_sur_0"], ti["predit_sat_sur_0"])),
        color = INK2, fontsize = 14, align = (:left, :top))
    ylims!(axB, 0.75, 1.95)
    axislegend(axB; position = (0.02, 0.82), labelsize = 12)

    # --- C : facteur 8 -------------------------------------------------------------------
    Nc = [1e3, 3e3, 1e4, 3e4]
    axC = Axis(fig[2, 3]; title = X[:C], subtitle = X[:C_sub], xlabel = X[:C_x], ylabel = X[:C_y], xscale = log10, yscale = log10,
        xticks = (Nc, ["10³", "3·10³", "10⁴", "3·10⁴"]), yticks = graduations_log(2, 100, langue))
    for (k, col) in (("4R", C["4"]), ("2H", C["2H"]))
        for (src, st, plein) in (("simulateur", :dash, nothing), ("jumeau", :solid, false), ("reel", :solid, true))
            d = sort(f8[(f8.source.==src).&(f8.C_vrai.==20.0), :], :N)
            y = d[!, "sigma_"*k]
            lines!(axC, d.N, y; color = col, linestyle = st, linewidth = src == "simulateur" ? 1.5 : 2.5)
            plein === nothing || scatter!(axC, d.N, y; color = plein ? col : SURF, strokecolor = col, strokewidth = 2, markersize = 11)
        end
    end
    fv(src) = f8.facteur[(f8.source.==src).&(f8.C_vrai.==20.0).&(f8.N.==1e4)][1]
    text!(axC, 1.05e3, 3.0; text = fr(Printf.format(Printf.Format(X[:C_txt]), fv("reel"), fv("jumeau"), fv("simulateur"))),
        color = INK2, fontsize = 13, align = (:left, :bottom))
    text!(axC, 3e3, 30.0; text = X[:C_leg][1], color = C["4"], fontsize = 13, align = (:left, :bottom), offset = (6, 2))
    text!(axC, 1e3, 75.0; text = X[:C_leg][2], color = C["2H"], fontsize = 13, align = (:left, :bottom))
    ylims!(axC, 2.5, 110)

    # --- D : classement et incertitude ----------------------------------------------------
    ks = ["1", "5d", "5a", "4", "2H"]
    axD = Axis(fig[3, 1]; title = X[:D], subtitle = X[:D_sub], xlabel = X[:D_x], xscale = log10, yreversed = true,
        yticks = (1:5, [TEXTES[langue][:leg][k] for k in ks]), xticks = graduations_log(0.5, 100, langue), ygridvisible = false)
    for (i, k) in enumerate(ks)
        s0 = sel(d20; estimateur = k, N_photons = 1e4)
        dm = mondes[(mondes.C_vrai.==20.0).&(mondes.N_photons.==1e4).&(mondes.normalisation.=="A").&(mondes.estimateur.==k), :]
        lo, hi = min(minimum(dm.sigma_mM), s0.sigma_mM[1]), max(maximum(dm.sigma_mM), s0.sigma_mM[1])
        lines!(axD, [lo, hi], [i, i]; color = hexa(T.GRIS), linewidth = 7)
        scatter!(axD, [s0.sigma_total_mM[1]], [i]; color = SURF, strokecolor = C[k], strokewidth = 2.2, markersize = 15)
        scatter!(axD, [s0.sigma_mM[1]], [i]; color = C[k], markersize = 13)
    end
    xlims!(axD, 0.5, 100)
    elD = [[MarkerElement(color = INK2, marker = :circle, markersize = 12)], [LineElement(color = hexa(T.GRIS), linewidth = 7)],
        [MarkerElement(color = SURF, strokecolor = INK2, strokewidth = 2, marker = :circle, markersize = 14)]]
    axislegend(axD, elD, X[:D_leg]; position = :rt, labelsize = 12)

    # --- E : fusion -----------------------------------------------------------------------------
    τs = [60.0, 300.0, 1200.0, 3600.0]
    axE = Axis(fig[3, 2]; title = X[:E], subtitle = X[:E_sub], xlabel = X[:E_x], ylabel = X[:E_y], xscale = log10, yscale = log10,
        xticks = (τs, X[:tg]), yticks = graduations_log(0.2, 50, langue))
    dE = fusion[startswith.(fusion.capteur, "mCer3").&(fusion.N.==1e4), :]
    for (f, col, lab) in (("ratio", C["1"], X[:E_leg][1]), ("forme", C["5a"], X[:E_leg][2]), ("fusion", INK, X[:E_leg][3]))
        g = sort(combine(groupby(dE[dE.filtre.==f, :], :tau_g_s), :RMSE_mM => mean => :r), :tau_g_s)
        lines!(axE, g.tau_g_s, g.r; color = col, label = lab)
        scatter!(axE, g.tau_g_s, g.r; color = col, markersize = 10, label = lab)
    end
    axislegend(axE; position = :rt, labelsize = 12, merge = true)
    ylims!(axE, 0.3, 40)

    # --- F : matrice de robustesse -------------------------------------------------------------------
    scen = [("gain_jaune_relatif", 1.05), ("decalage_jaune_ps", 20.0), ("blanchiment_accepteur", 0.1), ("blanchiment_accepteur_sombre", 0.1),
        ("eta_YD_relatif", 0.75), ("tau_D0_relatif", 1.05)]
    axF = Axis(fig[3, 3]; title = X[:F], subtitle = X[:F_sub], xlabel = X[:F_x], yreversed = true,
        yticks = (1:length(scen), X[:F_rows]), xscale = Makie.Symlog10(1.0), xticks = [-10, -1, 0, 1, 10], ygridvisible = false,
        yticklabelsize = 12)
    vlines!(axF, [0.0]; color = AXE, linewidth = 1)
    for (j, (k, dy)) in enumerate((("1", -0.25), ("5a", 0.0), ("4", 0.25)))
        v = [ajoute(p, val, (startswith(p, "tau") && k == "5a") ? "5b" : k) for (p, val) in scen]
        barplot!(axF, collect(1:length(scen)) .+ dy, v; direction = :x, width = 0.24, color = C[k], strokewidth = 0)
    end
    xlims!(axF, -20, 20)
    axislegend(axF, [PolyElement(color = C[k]) for k in ("1", "5a", "4")], X[:F_leg]; position = :rt, labelsize = 12)

    Label(fig[4, 1:3], X[:pied]; fontsize = 12, color = MUET, halign = :left, tellwidth = false)
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 34)
    chemin = joinpath(RACINE, "figures", "synthese_revision_$(langue)$(suffixe).png")
    save(chemin, fig; px_per_unit = 1.6)
    return chemin
end

for langue in ("FR", "EN"), suffixe in ("", "_dark")
    println(synthese(langue, suffixe))
    println(synthese_revision(langue, suffixe))
end
