# Dessins.jl — figures de l'analyse à deux canaux (étape 15), dans le module Figures.
#
# Trois figures demandées par le document (RMSE contre concentration ; RMSE contre budget à
# 20 mM ; biais par point de pourcentage blanchi), la carte finale, les diagnostics de
# l'ajustement global et les courbes de calibration. Chaque en-tête rappelle que la
# référence est la somme des 10 acquisitions (reproductibilité, pas justesse).

const LECTURES_DC = ["1", "1c", "5d", "5a", "6", "4", "2", "2J", "3"]
const NOM_DC = Dict("1" => "1 rapport J/B", "1c" => "1c J/B corrigé de la fuite", "2" => "2 mono-exp., bleu (ton code)",
    "2J" => "2J mono-exp., jaune (ton code)", "3" => "3 bi-exp. τ libres, bleu", "4" => "4 bi-exp. τ figés, bleu",
    "5a" => "5a forme 2 canaux, échelles séparées", "5d" => "5d forme 2 canaux, échelle commune", "6" => "6 phasor, bleu")
couleur_dc(k) = k == "1c" ? COULEUR["1"] : k == "2J" ? COULEUR["2"] : COULEUR[k]
marqueur_dc(k) = k == "1c" ? :xcross : k == "2J" ? :star5 : MARQUEUR[k]
style_dc(k) = k in ("1c", "2J") ? :dash : :solid
const COULEUR_LAMBDA = Dict(750 => RAMPE[2], 800 => RAMPE[5], 850 => RAMPE[9])
# grille des concentrations et nom du jeu (SCLM, KIRIN…), fixés par `figures_deux_canaux` d'après les données
const CONC_DC = Ref([0.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0, 131.0])
const JEU_DC = Ref("")
pos_dc(c) = something(findfirst(==(c), CONC_DC[]), NaN)
ticks_dc() = (1:length(CONC_DC[]), [virg(@sprintf("%g", c)) for c in CONC_DC[]])
cmax_dc() = virg(@sprintf("%g", last(CONC_DC[])))
jeu_dc() = isempty(JEU_DC[]) ? "" : JEU_DC[] * " · "
const REFERENCE_DC = "Référence : la somme des 10 acquisitions de chaque cellule. Elle mesure la reproductibilité, pas la justesse : " *
                     "un biais commun à la référence et aux tirages est invisible ici."

"""Lit un CSV des sorties (lignes « # » ignorées) ; les colonnes d'étiquettes restent du texte."""
function lire_dc(f)
    df = CSV.read(f, DataFrame; comment = "#", stripwhitespace = true)
    for c in (:lecture, :gagnante, :scenario_blanchiment, :canal, :modele)
        c in propertynames(df) && (df[!, c] = string.(df[!, c]))
    end
    return df
end

"""Bande de la fenêtre de clamping (5–35 mM : nœuds 5, 10, 20) sur l'axe catégoriel des concentrations."""
clamping_dc!(ax) = (p = [pos_dc(c) for c in CONC_DC[] if 5 <= c <= 35]; vspan!(ax, minimum(p) - 0.4, maximum(p) + 0.4; color = (ZONE, 1.0)))

function legende_dc!(pos, lectures; titre = "lecture")
    el = [[LineElement(color = couleur_dc(k), linestyle = style_dc(k), linewidth = 2),
           MarkerElement(color = couleur_dc(k), marker = marqueur_dc(k), markersize = 10, strokecolor = SURFACE, strokewidth = 1)] for k in lectures]
    Legend(pos, el, [NOM_DC[k] for k in lectures], titre; framevisible = false, labelsize = 11, rowgap = 2)
end

"""Fig. 17 — RMSE (vraie) contre concentration, par lecture : budgets 10³ et visé, une colonne par λ."""
function figure_dc_rmse_concentration(flat, chemin; budgets = (1e3, 1e4), essai = false)
    f = flat[flat.scenario_blanchiment.=="aucun", :]
    λs = sort(unique(f.lambda_nm))
    fig = Figure(size = (430 * length(λs) + 300, 820))
    for (i, N) in enumerate(budgets), (j, λ) in enumerate(λs)
        ax = Axis(fig[i, j]; yscale = log10, yticks = graduations_log(0.1, 1000), xticks = ticks_dc(),
            xlabel = i == length(budgets) ? "[Cl⁻] (mM)" : "", ylabel = j == 1 ? "RMSE vraie (mM)" : "",
            title = virg("$(λ) nm, $(fmtN(N)) photons bleus"))
        clamping_dc!(ax)
        for k in LECTURES_DC
            d = sort(f[(f.lambda_nm.==λ).&(f.budget_photons.==N).&(f.lecture.==k), :], :cl_mM)
            nrow(d) == 0 && continue
            x = pos_dc.(d.cl_mM)
            y = d.rmse_total_mM
            ok = isfinite.(y) .& (y .> 0)
            lines!(ax, x[ok], y[ok]; color = couleur_dc(k), linestyle = style_dc(k), linewidth = 2)
            borne = (d.cl_mM .== first(CONC_DC[])) .| (d.cl_mM .== last(CONC_DC[])) .| (d.taux_echec .>= 0.01)
            p = ok .& .!borne
            scatter!(ax, x[p], y[p]; color = couleur_dc(k), marker = marqueur_dc(k), markersize = 10)
            q = ok .& borne
            scatter!(ax, x[q], y[q]; color = SURFACE, strokecolor = couleur_dc(k), strokewidth = 1.6, marker = marqueur_dc(k), markersize = 10)
        end
        ylims!(ax, 0.1, 1000)
    end
    legende_dc!(fig[1:length(budgets), length(λs)+1], LECTURES_DC)
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 17 — Titration du 2026-10-16 : RMSE de chaque lecture selon la concentration",
        "RMSE vraie = √(biais² + σ² + σ_calibration²) en mM, σ corrigé de l'amincissement (√(1 − p)). Bande grise : fenêtre de clamping " *
        "(5–35 mM). Marqueur vide : 1 % d'échecs ou plus, ou borne de la calibration (0 et $(cmax_dc()) mM, rapportées mais non jugées). " * REFERENCE_DC;
        largeur = length(λs) + 1)
    save(chemin, fig; px_per_unit = 2)
end

"""Fig. 18 — RMSE contre budget à 20 mM, bornes de Fisher, et accord réel/jumeau."""
function figure_dc_rmse_budget(flat, fisher, jum, chemin; C = 20.0, essai = false)
    f = flat[(flat.scenario_blanchiment.=="aucun").&(flat.cl_mM.==C), :]
    λs = sort(unique(f.lambda_nm))
    Ns = sort(unique(f.budget_photons))
    tN = (Ns, fmtN.(Ns))
    fig = Figure(size = (430 * (length(λs) + 1) + 300, 560))
    for (j, λ) in enumerate(λs)
        ax = Axis(fig[1, j]; xscale = log10, yscale = log10, xticks = tN, yticks = graduations_log(0.05, 2000),
            xlabel = "photons bleus", ylabel = j == 1 ? "RMSE vraie (mM)" : "", title = virg(@sprintf("%d nm, %g mM", λ, C)))
        for k in LECTURES_DC
            d = sort(f[(f.lambda_nm.==λ).&(f.lecture.==k), :], :budget_photons)
            ok = isfinite.(d.rmse_total_mM) .& (d.rmse_total_mM .> 0)
            lines!(ax, d.budget_photons[ok], d.rmse_total_mM[ok]; color = couleur_dc(k), linestyle = style_dc(k), linewidth = 2)
            c = ok .& (d.taux_echec .>= 0.01)
            p = ok .& .!c
            scatter!(ax, d.budget_photons[p], d.rmse_total_mM[p]; color = couleur_dc(k), marker = marqueur_dc(k), markersize = 9)
            scatter!(ax, d.budget_photons[c], d.rmse_total_mM[c]; color = SURFACE, strokecolor = couleur_dc(k), strokewidth = 1.6,
                marker = marqueur_dc(k), markersize = 9)
        end
        b = sort(fisher[(fisher.lambda_nm.==λ).&(fisher.cl_mM.==C), :], :budget_photons)
        lines!(ax, b.budget_photons, b.toute; color = ENCRE, linestyle = :dash, linewidth = 1.2)
        lines!(ax, b.budget_photons, b.bleu; color = MUET, linestyle = :dot, linewidth = 1.4)
        ylims!(ax, 0.05, 2000)
    end
    ax = Axis(fig[1, length(λs)+1]; xscale = log10, xticks = tN, xlabel = "photons bleus", ylabel = "σ réel / σ jumeau",
        title = virg(@sprintf("Réel contre jumeau, %g mM (toutes λ)", C)))
    hspan!(ax, 0.9, 1.1; color = (ZONE, 1.0))
    hlines!(ax, [1.0]; color = AXE, linewidth = 1)
    jj = jum[jum.cl_mM.==C, :]
    for k in LECTURES_DC
        d = sort(jj[jj.lecture.==k, :], :budget_photons)
        nrow(d) == 0 && continue
        g = combine(groupby(d, :budget_photons), :rapport_sigma => (v -> (w = filter(isfinite, v); isempty(w) ? NaN : median(w))) => :r)
        g = g[isfinite.(g.r), :]
        lines!(ax, g.budget_photons, g.r; color = couleur_dc(k), linestyle = style_dc(k), linewidth = 1.6)
        scatter!(ax, g.budget_photons, g.r; color = couleur_dc(k), marker = marqueur_dc(k), markersize = 8)
    end
    ylims!(ax, 0.5, 2.0)
    el = [LineElement(color = ENCRE, linestyle = :dash), LineElement(color = MUET, linestyle = :dot)]
    g = GridLayout(fig[1, length(λs)+2])
    legende_dc!(g[1, 1], LECTURES_DC)
    Legend(g[2, 1], el, ["toute l'information (2 canaux)", "bleu seul"], "borne de Fisher"; framevisible = false, labelsize = 11)
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 18 — RMSE selon le budget de photons à $(Int(C)) mM, bornes de Fisher et jumeau",
        "Sommes réelles amincies (400 tirages sous 10³ photons, 200 au-dessus). Borne de Fisher au modèle global ajusté, fonds figés, " *
        "échelle libre. À droite : σ des tirages réels sur σ du jumeau (Poisson du modèle lissé), médiane des trois longueurs d'onde. " * REFERENCE_DC;
        largeur = length(λs) + 2)
    save(chemin, fig; px_per_unit = 2)
end

"""Fig. 19 — biais en mM par point de pourcentage blanchi, et taux tolérable, par lecture et scénario."""
function figure_dc_blanchiment(pentes, tol, chemin; budget = 1e4, essai = false)
    p = pentes[pentes.budget_photons.==budget, :]
    fig = Figure(size = (1500, 860))
    noms_sc = Dict("A" => "A — l'accepteur blanchi n'absorbe plus (FRET perdu)", "B" => "B — l'accepteur blanchi éteint encore (FRET gardé)",
        "C" => "C — le donneur blanchit")
    xt = (1:length(LECTURES_DC), LECTURES_DC)
    for (j, sc) in enumerate(("A", "B", "C"))
        ax = Axis(fig[1, j]; xticks = xt, yscale = Makie.Symlog10(0.1), xlabel = "", ylabel = j == 1 ? "biais (mM) par % blanchi" : "",
            yticks = ([-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100], ["−100", "−10", "−1", "−0,1", "0", "0,1", "1", "10", "100"]),
            title = noms_sc[sc])
        hlines!(ax, [0.0]; color = AXE, linewidth = 1)
        for (i, k) in enumerate(LECTURES_DC), (λ, dx) in ((750, -0.22), (800, 0.0), (850, 0.22))
            d = p[(p.lecture.==k).&(p.scenario_blanchiment.==sc).&(p.lambda_nm.==λ).&(5 .<= p.cl_mM .<= 35), :]
            v = filter(isfinite, d.pente_mM_par_pct)
            isempty(v) && continue
            rangebars!(ax, [i + dx], [minimum(v)], [maximum(v)]; color = (COULEUR_LAMBDA[λ], 0.6), whiskerwidth = 4)
            scatter!(ax, [i + dx], [median(v)]; color = COULEUR_LAMBDA[λ], marker = marqueur_dc(k), markersize = 10)
        end
        ax2 = Axis(fig[2, j]; xticks = xt, yscale = log10, yticks = graduations_log(0.01, 1000), xlabel = "lecture",
            ylabel = j == 1 ? "taux tolérable (% blanchi) pour |biais| < 1 mM" : "", title = "Taux tolérable sur 5–35 mM")
        for (i, k) in enumerate(LECTURES_DC), (λ, dx) in ((750, -0.22), (800, 0.0), (850, 0.22))
            d = tol[(tol.lecture.==k).&(tol.scenario_blanchiment.==sc).&(tol.lambda_nm.==λ), :]
            nrow(d) == 1 || continue
            t = d.taux_tolerable_5_35_pct[1]
            if isinf(t) || t > 1000         # pente nulle : la lecture ne voit pas ce blanchiment
                scatter!(ax2, [i + dx], [800.0]; color = SURFACE, strokecolor = COULEUR_LAMBDA[λ], strokewidth = 1.6, marker = :utriangle, markersize = 11)
            elseif isfinite(t)
                scatter!(ax2, [i + dx], [max(t, 0.01)]; color = COULEUR_LAMBDA[λ], marker = marqueur_dc(k), markersize = 10)
            end
        end
        hlines!(ax2, [1.0, 10.0]; color = AXE, linewidth = 1, linestyle = :dash)
        ylims!(ax2, 0.01, 1000)
    end
    Legend(fig[1:2, 4], [MarkerElement(color = COULEUR_LAMBDA[λ], marker = :circle, markersize = 10) for λ in (750, 800, 850)],
        ["750 nm", "800 nm", "850 nm"], "longueur d'onde"; framevisible = false)
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 19 — Photoblanchiment injecté : biais par point de pourcentage blanchi",
        "Pente du biais contre le taux blanchi (0, 1, 2, 5, 10 % des molécules pendant l'acquisition, décroissance exponentielle), à " *
        "$(fmtN(budget)) photons bleus ; point = médiane sur 5–35 mM, barre = étendue. En bas : 1 mM / la plus forte pente sur 5–35 mM ; " *
        "triangle vide en haut : pente nulle, la lecture ne voit pas ce blanchiment. " * REFERENCE_DC;
        largeur = 4)
    save(chemin, fig; px_per_unit = 2)
end

"""Fig. 20 — la carte finale : lecture gagnante par (concentration, budget), une carte par λ."""
function figure_dc_carte(carte, chemin; essai = false)
    λs = sort(unique(carte.lambda_nm))
    Ns = sort(unique(carte.budget_photons))
    fig = Figure(size = (420 * length(λs) + 260, 620))
    for (j, λ) in enumerate(λs)
        ax = Axis(fig[1, j]; xticks = (1:length(Ns), fmtN.(Ns)), yticks = ticks_dc(), xlabel = "photons bleus",
            ylabel = j == 1 ? "[Cl⁻] (mM)" : "", title = "$(λ) nm", xgridvisible = false, ygridvisible = false)
        for r in eachrow(carte[carte.lambda_nm.==λ, :])
            x, y = findfirst(==(r.budget_photons), Ns), pos_dc(r.cl_mM)
            isfinite(y) || continue
            col = r.gagnante == "aucune" ? GRILLE : couleur_dc(r.gagnante)
            poly!(ax, Rect2f(x - 0.48, y - 0.48, 0.96, 0.96); color = r.gagnante == "aucune" ? SURFACE : (col, 0.85), strokecolor = GRILLE, strokewidth = 1)
            text!(ax, x, y; text = r.gagnante == "aucune" ? "—" : r.gagnante, align = (:center, :center), fontsize = 11,
                color = r.gagnante == "aucune" ? MUET : SURFACE, font = :bold)
        end
        xlims!(ax, 0.5, length(Ns) + 0.5)
        ylims!(ax, 0.5, 8.5)
    end
    gag = sort(unique(filter(!=("aucune"), carte.gagnante)); by = k -> findfirst(==(k), LECTURES_DC))
    Legend(fig[1, length(λs)+1], [PolyElement(color = (couleur_dc(k), 0.85)) for k in gag], [NOM_DC[k] for k in gag], "gagnante";
        framevisible = false, labelsize = 11)
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 20 — Carte finale : la meilleure lecture selon la concentration et le budget",
        "Dans chaque case, la lecture de plus petite RMSE vraie parmi celles qui échouent moins d'une fois sur cent et dont le biais reste " *
        "sous 1 mM à ce point (aux bornes de la calibration, seuls les échecs d'ajustement comptent). « — » : aucune lecture ne passe. " *
        "1 et 5d diffèrent de moins de 3 % partout : là où l'un gagne, l'autre est à égalité statistique (σ connu à 5 %). " * REFERENCE_DC;
        largeur = length(λs) + 1)
    save(chemin, fig; px_per_unit = 2)
end

"""Fig. 21 — l'ajustement global : χ²ᵣ par point, mélange à deux états, t₀, fuite, populations, résidus."""
function figure_dc_global(gc, lien, svd, chemin; pas = 0.06400390634153652, essai = false)
    λs = sort(unique(gc.lambda_nm))
    fig = Figure(size = (1550, 900))
    a1, a2, a3, a4, a5 = (Axis(fig[i, j]; xticks = ticks_dc(), xlabel = "[Cl⁻] (mM)") for (i, j) in ((1, 1), (1, 2), (1, 3), (2, 1), (2, 2)))
    a6 = Axis(fig[2, 3])
    a1.title = "χ²ᵣ par concentration (plein : bleu ; tireté : jaune)"
    a2.title = "Mélange à deux états : χ²ᵣ en excès du global"
    a3.title = "t₀ par concentration, relatif à la moyenne (ps)"
    a4.title = "Coefficient de fuite λ = I_J(rapide) / I_B(τ₂)"
    a5.title = "Composante τ₂ du bleu : amplitude (plein), intensité (tireté)"
    a6.title = "Résidus : premier vecteur singulier (plein : bleu ; tireté : jaune)"
    a6.xlabel = "temps (ns)"
    for λ in λs
        d = sort(gc[gc.lambda_nm.==λ, :], :cl_mM)
        x = pos_dc.(d.cl_mM)
        c = COULEUR_LAMBDA[λ]
        lines!(a1, x, d.chi2r_bleu; color = c)
        lines!(a1, x, d.chi2r_jaune; color = c, linestyle = :dash)
        scatterlines!(a2, x, d.exces_melange; color = c, markersize = 7, label = "$(λ) nm")
        lines!(a3, x, d.t0_bleu_par_cellule_ps .- mean(d.t0_bleu_par_cellule_ps); color = c)
        lines!(a3, x, d.t0_jaune_par_cellule_ps .- mean(d.t0_jaune_par_cellule_ps); color = c, linestyle = :dash)
        scatterlines!(a4, x, d.lambda_fuite; color = c, markersize = 7)
        lines!(a5, x, d.amplitude_tau2_bleu; color = c)
        lines!(a5, x, d.intensite_tau2_bleu; color = c, linestyle = :dash)
        for (canal, st) in (("bleu", :solid), ("jaune", :dash))
            s = sort(svd[(svd.lambda_nm.==λ).&(svd.canal.==canal), :], :canal_tcspc)
            nrow(s) == 0 && continue
            lines!(a6, (s.canal_tcspc .- s.canal_tcspc[1]) .* pas, s.vecteur; color = c, linestyle = st, linewidth = 1.4)
        end
    end
    hlines!(a2, [0.0]; color = AXE, linewidth = 1)
    hlines!(a3, [0.0]; color = AXE, linewidth = 1)
    axislegend(a2; position = :lt, framevisible = false)
    lignes = [virg(@sprintf("%d nm : Δdéviance/χ²ᵣ = %.0f (seuil 3,84) → le lien %s ; τ₂ bleu %.2f ns, τ_J libre %.2f ns", r.lambda_nm,
        r.statistique, r.lien_tient ? "tient" : "force", r.tau2_bleu, r.tauJ_libre)) for r in eachrow(lien)]
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 21 — Ajustement global par longueur d'onde (8 concentrations × 2 canaux)",
        "Bleu : diffusion + τ₂ + τ₃ ; jaune : diffusion + τ_J + τ_A ; durées de vie communes, amplitudes libres, t₀ par canal, fond figé. " *
        "Lien de la composante rapide du jaune à τ₂ : " * join(lignes, " ; ") * "."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

"""Fig. 22 — les courbes de calibration de chaque lecture (normalisées), et leur K apparent."""
function figure_dc_calibration(refs, cal, chemin; essai = false)
    λs = sort(unique(refs.lambda_nm))
    fig = Figure(size = (1400, 980))
    for (q, k) in enumerate(LECTURES_DC)
        i, j = fldmod1(q, 3)
        ax = Axis(fig[i, j]; xticks = ticks_dc(), xlabel = i == 3 ? "[Cl⁻] (mM)" : "", ylabel = j == 1 ? "(o − o₀) / (o_max − o₀)" : "",
            title = NOM_DC[k])
        hlines!(ax, [0.0, 1.0]; color = AXE, linewidth = 1)
        txt = String[]
        for λ in λs
            d = sort(refs[(refs.lambda_nm.==λ).&(refs.lecture.==k), :], :cl_mM)
            nrow(d) == 0 && continue
            y = (d.o .- d.o[1]) ./ (d.o[end] - d.o[1])
            σ = abs.(d.sigma_ref ./ (d.o[end] - d.o[1]))
            x = pos_dc.(d.cl_mM)
            errorbars!(ax, x, y, σ; color = COULEUR_LAMBDA[λ], whiskerwidth = 4)
            scatterlines!(ax, x, y; color = COULEUR_LAMBDA[λ], marker = marqueur_dc(k), markersize = 8)
            c = cal[(cal.lambda_nm.==λ).&(cal.lecture.==k), :]
            nrow(c) == 1 && push!(txt, virg(@sprintf("%d nm : K %s%s", λ, isfinite(c.K_apparent_mM[1]) ? @sprintf("%.0f mM", c.K_apparent_mM[1]) : "—",
                c.monotone[1] ? "" : " (non monotone)")))
        end
        text!(ax, 0.98, 0.04; text = join(txt, "\n"), space = :relative, align = (:right, :bottom), fontsize = 10, color = ENCRE2)
    end
    Legend(fig[1:3, 4], [MarkerElement(color = COULEUR_LAMBDA[λ], marker = :circle, markersize = 9) for λ in λs], ["$(λ) nm" for λ in λs],
        "longueur d'onde"; framevisible = false)
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 22 — Courbes de calibration des neuf lectures (8 points à pleins photons)",
        "Chaque lecture ramenée sur l'écart entre 0 et $(cmax_dc()) mM (o_max) ; barres : incertitude de chaque point (dispersion des 10 acquisitions / √10). " *
        "K apparent : Langmuir ajustée sur la courbe. Une lecture non monotone n'est pas inversée. " * REFERENCE_DC; largeur = 4)
    save(chemin, fig; px_per_unit = 2)
end

"""Toutes les figures de l'étape 15, depuis les sorties persistées de `dossier`."""
function figures_deux_canaux(dossier, chemin_fig; budget = 1e4, essai = false, jeu = "")
    set_theme!(theme_sim())
    flat = lire_dc(joinpath(dossier, "resultats_deux_canaux.csv"))
    CONC_DC[] = sort(unique(flat.cl_mM))
    JEU_DC[] = jeu
    fisher = lire_dc(joinpath(dossier, "fisher.csv"))
    jum = lire_dc(joinpath(dossier, "jumeau_contre_reel.csv"))
    figure_dc_rmse_concentration(flat, chemin_fig("fig17_deux_canaux_rmse_concentration.png"); budgets = (1e3, budget), essai)
    figure_dc_rmse_budget(flat, fisher, jum, chemin_fig("fig18_deux_canaux_rmse_budget_20mM.png"); essai)
    fp = joinpath(dossier, "pentes_blanchiment.csv")
    if isfile(fp)
        figure_dc_blanchiment(lire_dc(fp), lire_dc(joinpath(dossier, "blanchiment_tolerable.csv")),
            chemin_fig("fig19_deux_canaux_blanchiment.png"); budget, essai)
    end
    figure_dc_carte(lire_dc(joinpath(dossier, "carte_finale.csv")), chemin_fig("fig20_deux_canaux_carte_finale.png"); essai)
    figure_dc_global(lire_dc(joinpath(dossier, "global_cellules.csv")), lire_dc(joinpath(dossier, "global_lien.csv")),
        lire_dc(joinpath(dossier, "global_residus_svd.csv")), chemin_fig("fig21_deux_canaux_ajustement_global.png"); essai)
    figure_dc_calibration(lire_dc(joinpath(dossier, "lectures_pleins_photons.csv")), lire_dc(joinpath(dossier, "calibration.csv")),
        chemin_fig("fig22_deux_canaux_calibration.png"); essai)
    return nothing
end

"""
Fig. 23 — l'échantillon inconnu : la concentration lue par chaque lecture, à chaque longueur
d'onde (barre fine : ±σ total ; barre épaisse : ±σ statistique), et la moyenne des
longueurs d'onde (losange noir, ±σ gonflé par √χ²ᵣ). Bande : J/B, la lecture recommandée.
"""
function figure_dc_inconnue(tab, comb, chemin; essai = false, jeu = "")
    JEU_DC[] = jeu
    fig = Figure(size = (1350, 820))
    ny = length(LECTURES_DC)
    ax = Axis(fig[1, 1]; yticks = (1:ny, [NOM_DC[k] for k in LECTURES_DC]), yreversed = true, xlabel = "[Cl⁻] prédite (mM)",
        title = "Concentration de l'échantillon inconnu, par lecture et par longueur d'onde", ygridvisible = false)
    r1 = comb[comb.lecture.=="1", :]
    if nrow(r1) == 1 && isfinite(r1.c_mM[1])
        vspan!(ax, r1.c_mM[1] - r1.sigma_gonfle_mM[1], r1.c_mM[1] + r1.sigma_gonfle_mM[1]; color = (ZONE, 1.0))
        vlines!(ax, [r1.c_mM[1]]; color = AXE, linewidth = 1)
    end
    for (i, k) in enumerate(LECTURES_DC)
        for (λ, dy) in ((750, -0.27), (800, -0.09), (850, 0.09))
            d = tab[(tab.lecture.==k).&(tab.lambda_nm.==λ), :]
            (nrow(d) == 1 && isfinite(d.c_mM[1])) || continue
            c, σt, σs = d.c_mM[1], d.sigma_total_mM[1], d.sigma_stat_mM[1]
            col = COULEUR_LAMBDA[λ]
            isfinite(σt) && rangebars!(ax, [i + dy], [c - σt], [c + σt]; direction = :x, color = (col, 0.55), linewidth = 1.2, whiskerwidth = 6)
            isfinite(σs) && rangebars!(ax, [i + dy], [c - σs], [c + σs]; direction = :x, color = col, linewidth = 3.5, whiskerwidth = 0)
            if d.monotone[1]
                scatter!(ax, [c], [i + dy]; color = col, marker = marqueur_dc(k), markersize = 11)
            else
                scatter!(ax, [c], [i + dy]; color = SURFACE, strokecolor = col, strokewidth = 1.6, marker = marqueur_dc(k), markersize = 11)
            end
        end
        r = comb[comb.lecture.==k, :]
        if nrow(r) == 1 && isfinite(r.c_mM[1])
            rangebars!(ax, [i + 0.27], [r.c_mM[1] - r.sigma_gonfle_mM[1]], [r.c_mM[1] + r.sigma_gonfle_mM[1]]; direction = :x, color = ENCRE,
                linewidth = 1.6, whiskerwidth = 6)
            scatter!(ax, [r.c_mM[1]], [i + 0.27]; color = ENCRE, marker = :diamond, markersize = 11)
            text!(ax, r.c_mM[1] + r.sigma_gonfle_mM[1], i + 0.27; text = virg(@sprintf("  %.1f ± %.1f", r.c_mM[1], r.sigma_gonfle_mM[1])),
                align = (:left, :center), fontsize = 10, color = ENCRE2)
        end
    end
    el = vcat([MarkerElement(color = COULEUR_LAMBDA[λ], marker = :circle, markersize = 10) for λ in (750, 800, 850)],
        [MarkerElement(color = ENCRE, marker = :diamond, markersize = 10),
         MarkerElement(color = SURFACE, strokecolor = ENCRE2, strokewidth = 1.5, marker = :circle, markersize = 10)])
    Legend(fig[1, 2], el, ["750 nm", "800 nm", "850 nm", "moyenne des λ (±σ gonflé)", "courbe non monotone (indicatif)"];
        framevisible = false, labelsize = 11)
    entete!(fig, (essai ? "ESSAI — " : "") * jeu_dc() * "Fig. 23 — Échantillon inconnu : la concentration prédite par chaque lecture",
        "Somme des 10 acquisitions lue par chaque lecture, inversée par sa courbe de calibration (Hermite monotone sur les 8 points). " *
        "Barre épaisse : σ statistique (Monte-Carlo : mesure de l'inconnu + bruit des 8 points) ; barre fine : σ total, avec la sensibilité au " *
        "saut de t₀ observé entre cellules. Bande grise : J/B moyenné sur les trois longueurs d'onde. " * REFERENCE_DC; largeur = 2)
    save(chemin, fig; px_per_unit = 2)
end
