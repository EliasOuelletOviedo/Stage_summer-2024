# figures.jl — le rendu, après la simulation : les vidéos, le résumé, le réglage du filtre,
# Fisher contre Monte-Carlo, le test d'amplitudes.

using CairoMakie

const SURFACE = colorant"#fcfcfb"
const COULEURS = Dict("J/B" => colorant"#d9822b", "5a" => colorant"#7b5ea7", "Kalman" => colorant"#1b7f79",
                      "ancien" => colorant"#9aa5a4", "bleu" => colorant"#3b6fd6", "jaune" => colorant"#c99a06")
# une couleur par réglage du filtre (marche aléatoire : dégradé de vert ; modèle à vitesse : rouges ; OU : brun)
const PALETTE_DERIVE = [colorant"#222222", colorant"#4c78a8", colorant"#1b7f79", colorant"#d9822b", colorant"#b8336a"]
couleur_filtre(nom) = startswith(nom, "MA") ? COULEURS["Kalman"] : startswith(nom, "OU") ? colorant"#8c564b" : colorant"#c03d3e"
axe(pos; kw...) = Axis(pos; backgroundcolor = SURFACE, kw...)

"""
Trois panneaux empilés : [Cl⁻] (vérité, J/B, 5a, Kalman ± 2√P₁₁), ĝ (vérité, estimation),
les deux histogrammes de l'instant courant en log. Une image par pas de 0,1 s ; à 30 ips,
la vidéo va trois fois plus vite que la course. Garde la dernière image en PNG.
"""
function video(T, H, p::Params, chemin)
    t, n, W = T[:t], length(T[:t]), p.fenetre
    ips = p.toml["video"]["ips"]
    fig = Figure(size = (900, 980), backgroundcolor = SURFACE, fontsize = 15)
    ax1 = Axis(fig[1, 1], ylabel = "[Cl⁻] (mM)", backgroundcolor = SURFACE)
    ax2 = Axis(fig[2, 1], ylabel = "g, gain du canal jaune", xlabel = "temps (s)", backgroundcolor = SURFACE)
    ax3 = Axis(fig[3, 1], yscale = log10, xlabel = "temps TCSPC (ns)", ylabel = "comptages (0,1 s)", backgroundcolor = SURFACE)
    linkxaxes!(ax1, ax2)
    k = Observable(1)
    pts(v) = @lift Point2f.(t[1:$k], v[1:$k])
    σ = 2 .* sqrt.(T[:P11])
    band!(ax1, @lift(t[1:$k]), @lift(T[:C_kalman][1:$k] .- σ[1:$k]), @lift(T[:C_kalman][1:$k] .+ σ[1:$k]),
          color = (COULEURS["Kalman"], 0.25))
    lines!(ax1, t, T[:C_vrai], color = :black, linewidth = 2, label = "vérité")
    scatter!(ax1, pts(T[:C_ratio]), markersize = 4, color = COULEURS["J/B"], label = "J/B")
    scatter!(ax1, pts(T[:C_5a]), markersize = 4, color = COULEURS["5a"], label = "5a")
    lines!(ax1, pts(T[:C_kalman]), color = COULEURS["Kalman"], linewidth = 2.5, label = "Kalman ± 2√P₁₁")
    axislegend(ax1, position = :rt, framevisible = false)
    lines!(ax2, t, T[:g_vrai], color = :black, linewidth = 2, label = "vérité")
    lines!(ax2, pts(T[:g_chapeau]), color = COULEURS["Kalman"], linewidth = 2.5, label = "ĝ")
    axislegend(ax2, position = :lt, framevisible = false)
    x = (W .- 1) .* p.pas
    for (i, nom) in ((1, "bleu"), (2, "jaune"))
        stairs!(ax3, x, @lift(max.(H[$k][i][W], 0.7)), color = COULEURS[nom], step = :center, label = nom)
    end
    axislegend(ax3, position = :rt, framevisible = false)
    xlims!(ax1, t[1], t[end] + p.dt)
    ylims!(ax1, -5, 50)
    ylims!(ax2, 0.9, 1 + max(0.25, 1.2 * (maximum(T[:g_vrai]) - 1)))
    ylims!(ax3, 0.7, 1.5 * maximum(maximum(h[2][W]) for h in H))
    record(fig, chemin * ".mp4", 1:n; framerate = ips) do i
        k[] = i
    end
    save(chemin * ".png", fig)
end

"""
Résumé au réglage de référence (N = 300 retiré : ce n'est pas un régime de fonctionnement) :
(a) biais des 30 dernières secondes contre la dérive ; (b) le même biais pour chaque modèle de g ;
(c) RMSE contre le budget ; (d) ĝ autour du saut, moyenne des répétitions, pour chaque modèle de g.
"""
function figure_resume(resume, saut, p::Params, chemin)
    ref, anc = reference(p), nom_ma(minimum(p.toml["kalman"]["balayage_sigma_wg"]))
    Nv = p.toml["video"]["N_bleu"]
    Ns = p.toml["criteres"]["budgets"]
    αs = p.toml["derive"]["balayage_alpha"]
    nomα(α) = nom_derive("rampe", α, p)
    lect = [("J/B", "J/B", COULEURS["J/B"], :solid), ("5a", "5a", COULEURS["5a"], :solid),
            (ref, "Kalman, " * ref * " (référence)", COULEURS["Kalman"], :solid), (anc, "Kalman, " * anc * " (ancien défaut)", COULEURS["ancien"], :dash)]
    fig = Figure(size = (1300, 950), backgroundcolor = SURFACE, fontsize = 14)
    a = axe(fig[1, 1], title = @sprintf("(a) biais sur les 30 dernières secondes, N = %g", Nv), xlabel = "dérive α du gain jaune (%)", ylabel = "biais moyen (mM)")
    for (f, lab, col, st) in lect
        r = [ligne(resume, Nv, nomα(α), f) for α in αs]
        errorbars!(a, 100 .* αs, getfield.(r, :biais_30s), getfield.(r, :biais_30s_sd), color = col)
        scatterlines!(a, 100 .* αs, getfield.(r, :biais_30s), color = col, linestyle = st, label = lab)
    end
    hlines!(a, [0.0], color = :black, linewidth = 0.8)
    axislegend(a, position = :lb, framevisible = false)
    b = axe(fig[1, 2], title = @sprintf("(b) modèles de g : biais sur les 30 dernières secondes, N = %g", Nv), xlabel = "dérive α du gain jaune (%)", ylabel = "biais moyen (mM)")
    for nom in unique(r.filtre for r in resume if r.lecture == "Kalman" && (r.filtre == ref || !startswith(r.filtre, "MA")))
        r = [ligne(resume, Nv, nomα(α), nom) for α in αs]
        scatterlines!(b, 100 .* αs, getfield.(r, :biais_30s), color = couleur_filtre(nom), label = nom,
                      linestyle = nom == ref ? :solid : occursin("1e-05", nom) ? :dot : occursin("1e-03", nom) ? :dashdot : :dash)
    end
    hlines!(b, [0.0], color = :black, linewidth = 0.8)
    axislegend(b, position = :lb, framevisible = false)
    c = axe(fig[2, 1], title = "(c) RMSE sur toute la course, rampe 10 %", xlabel = "photons bleus par acquisition", ylabel = "RMSE (mM)",
            xscale = log10, yscale = log10, xticks = (Ns, [@sprintf("%g", N) for N in Ns]))
    for (f, lab, col, st) in lect
        scatterlines!(c, Ns, [ligne(resume, N, "rampe 10 %", f).rmse for N in Ns], color = col, linestyle = st, label = lab)
    end
    axislegend(c, position = :rt, framevisible = false)
    d = axe(fig[2, 2], title = @sprintf("(d) ĝ autour du saut de %g %% à %g s, N = %g (moyenne de %d courses)", 100p.derive.saut,
            p.derive.t_saut, Nv, p.toml["course"]["repetitions"]), xlabel = "temps (s)", ylabel = "g")
    for (nom, s) in sort(collect(saut); by = first)
        (nom == ref || nom == anc || !startswith(nom, "MA")) || continue
        lines!(d, s.t, s.gm, color = nom == anc ? COULEURS["ancien"] : couleur_filtre(nom), label = nom,
               linestyle = nom in (ref, anc) ? :solid : occursin("1e-05", nom) ? :dot : occursin("1e-03", nom) ? :dashdot : :dash)
    end
    s = saut[ref]
    lines!(d, s.t, s.g, color = :black, linewidth = 2, label = "vérité")
    xlims!(d, p.derive.t_saut - 10, s.t[end])
    axislegend(d, position = :rb, framevisible = false)
    save(chemin, fig)
end

"""
Le critère de réglage (correctif 4) : P₁₁/err² contre σ_w,g, une courbe par dérive, avec la
bande des critères de sortie ; en dessous, ce que coûte chaque réglage en biais et en RMSE.
"""
function figure_reglage(resume, p::Params, chemin)
    σs = p.toml["kalman"]["balayage_sigma_wg"]
    ders = [[nom_derive("rampe", α, p) for α in p.toml["derive"]["balayage_alpha"]]; nom_derive("saut", 0, p)]
    Ns = p.toml["criteres"]["budgets"]
    bande = p.toml["criteres"]["P11_sur_err2"]
    fig = Figure(size = (1300, 950), backgroundcolor = SURFACE, fontsize = 14)
    for (j, N) in enumerate(Ns)
        a = axe(fig[1, j], title = @sprintf("P₁₁ annoncé ÷ erreur quadratique réelle, N = %g", N), xlabel = "σ_w,g (par pas)",
                ylabel = "P₁₁ / err²", xscale = log10, yscale = log10)
        hspan!(a, bande..., color = (:black, 0.06))
        hlines!(a, [1.0], color = :black, linewidth = 0.8)
        for (d, col) in zip(ders, PALETTE_DERIVE)
            scatterlines!(a, σs, [ligne(resume, N, d, nom_ma(σ)).P11_sur_err2 for σ in σs], color = col, label = d)
        end
        axislegend(a, position = :rb, framevisible = false)
    end
    for (j, (champ, titre)) in enumerate(((:biais_30s, "biais sur les 30 dernières secondes (mM)"), (:rmse, "RMSE sur toute la course (mM)")))
        a = axe(fig[2, j], title = titre * " ; trait plein N = $(Int(Ns[1])), tirets N = $(Int(Ns[2]))", xlabel = "σ_w,g (par pas)", xscale = log10)
        for (N, st) in zip(Ns, (:solid, :dash)), (d, col) in zip(ders, PALETTE_DERIVE)
            lines!(a, σs, [getfield(ligne(resume, N, d, nom_ma(σ)), champ) for σ in σs], color = col, linestyle = st, label = N == Ns[1] ? d : nothing)
            scatter!(a, σs, [getfield(ligne(resume, N, d, nom_ma(σ)), champ) for σ in σs], color = col, markersize = 6)
        end
        if champ === :biais_30s
            hlines!(a, [0.0], color = :black, linewidth = 0.8)
            hlines!(a, [-0.3, 0.3], color = :black, linewidth = 0.6, linestyle = :dot)
        end
        axislegend(a, position = champ === :biais_30s ? :rb : :rt, framevisible = false)
    end
    save(chemin, fig)
end

"""Fisher contre Monte-Carlo : var(J/B), var(5a), covariance et corrélation, aux trois budgets."""
function figure_fisher(verif, chemin)
    fig = Figure(size = (1500, 420), backgroundcolor = SURFACE, fontsize = 14)
    N = [v.N_photons_bleus for v in verif]
    M = [v.tirages for v in verif]
    for (j, (reel, fisher, titre, log)) in enumerate(((:var_ratio_reelle, :var_ratio_fisher, "var(J/B)", true),
            (:var_5a_reelle, :var_5a_fisher, "var(z_5a) (mM²-équivalent)", true), (:cov_reelle, :cov_fisher, "cov(J/B, z_5a)", false),
            (:rho_reel, :rho_fisher, "corrélation ρ", false)))
        a = axe(fig[1, j], title = titre, xlabel = "photons bleus par acquisition", xscale = log10, yscale = log ? log10 : identity,
                xticks = (N, [@sprintf("%g", n) for n in N]))
        r, f = [getfield(v, reel) for v in verif], [getfield(v, fisher) for v in verif]
        e = j <= 2 ? r .* sqrt.(2 ./ (M .- 1)) :
            j == 3 ? [sqrt((v.var_ratio_reelle * v.var_5a_reelle + v.cov_reelle^2) / (v.tirages - 1)) for v in verif] :
            (1 .- r .^ 2) ./ sqrt.(M .- 1)
        scatterlines!(a, N, f, color = COULEURS["Kalman"], label = "Fisher")
        errorbars!(a, N, r, e, color = :black, whiskerwidth = 8)
        scatter!(a, N, r, color = :black, marker = :circle, markersize = 7, label = "Monte-Carlo (±1σ)")
        j >= 3 && hlines!(a, [0.0], color = :black, linewidth = 0.5)
        j == 1 && axislegend(a, position = :rt, framevisible = false)
    end
    save(chemin, fig)
end

"""Test d'amplitudes sans bruit : p et q retrouvés par 5a contre injectés, avec et sans la diffusion du générateur."""
function figure_amplitudes(amp, chemin)
    C = [a.C_mM for a in amp]
    fig = Figure(size = (1200, 800), backgroundcolor = SURFACE, fontsize = 14)
    for (j, (x, nom, col)) in enumerate(((:p, "p (fraction d'amplitude de τ₂, bleu)", COULEURS["bleu"]),
                                          (:q, "q (fraction d'amplitude de τ_J, jaune)", COULEURS["jaune"])))
        inj = [getfield(a, Symbol(x, :_injecte)) for a in amp]
        av = [getfield(a, Symbol(x, :_avec_diffusion)) for a in amp]
        sa = [getfield(a, Symbol(x, :_sans_diffusion)) for a in amp]
        a = axe(fig[1, j], title = nom, ylabel = string(x))
        lines!(a, C, inj, color = :black, linewidth = 2, label = "injecté (Langmuir de la titration)")
        scatter!(a, C, av, color = col, label = "5a, générateur avec diffusion")
        scatter!(a, C, sa, color = :black, marker = :xcross, markersize = 9, label = "5a, générateur sans diffusion")
        axislegend(a, position = x === :p ? :rt : :lt, framevisible = false)
        b = axe(fig[2, j], title = "écart relatif 5a − injecté", xlabel = "[Cl⁻] (mM)", ylabel = "écart (%)")
        scatterlines!(b, C, 100 .* (av .- inj) ./ inj, color = col, label = "avec diffusion")
        scatterlines!(b, C, 100 .* (sa .- inj) ./ inj, color = :black, marker = :xcross, label = "sans diffusion")
        d = [getfield(a, x === :p ? :d_bleu : :d_jaune) for a in amp]
        lines!(b, C, d ./ maximum(d) .* maximum(100 .* (av .- inj) ./ inj), color = (col, 0.5), linestyle = :dash,
               label = "fraction de diffusion d(C), échelle arbitraire")
        axislegend(b, position = :rt, framevisible = false)
    end
    save(chemin, fig)
end
