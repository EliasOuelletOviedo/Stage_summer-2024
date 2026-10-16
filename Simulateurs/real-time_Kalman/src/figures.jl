# figures.jl — le rendu, après la simulation : la vidéo d'une course et la figure de résumé.

using CairoMakie

const SURFACE = colorant"#fcfcfb"
const COULEURS = Dict("J/B" => colorant"#d9822b", "5a" => colorant"#7b5ea7", "Kalman R complet" => colorant"#1b7f79",
                      "Kalman R diagonal" => colorant"#9aa5a4", "bleu" => colorant"#3b6fd6", "jaune" => colorant"#c99a06")

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
          color = (COULEURS["Kalman R complet"], 0.25))
    lines!(ax1, t, T[:C_vrai], color = :black, linewidth = 2, label = "vérité")
    scatter!(ax1, pts(T[:C_ratio]), markersize = 4, color = COULEURS["J/B"], label = "J/B")
    scatter!(ax1, pts(T[:C_5a]), markersize = 4, color = COULEURS["5a"], label = "5a")
    lines!(ax1, pts(T[:C_kalman]), color = COULEURS["Kalman R complet"], linewidth = 2.5, label = "Kalman ± 2√P₁₁")
    axislegend(ax1, position = :rt, framevisible = false)
    lines!(ax2, t, T[:g_vrai], color = :black, linewidth = 2, label = "vérité")
    lines!(ax2, pts(T[:g_chapeau]), color = COULEURS["Kalman R complet"], linewidth = 2.5, label = "ĝ")
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
Résumé : (a) biais des 30 dernières secondes contre la dérive, (b) P₁₁ annoncé ÷ erreur
réelle du Kalman avec R complet et diagonal, (c) RMSE contre le budget, (d) ĝ après le saut.
"""
function figure_resume(resume, saut, p::Params, chemin)
    Nv = p.toml["video"]["N_bleu"]
    αs = p.toml["derive"]["balayage_alpha"]
    lect = ["J/B", "5a", "Kalman R complet", "Kalman R diagonal"]
    σ0 = p.toml["kalman"]["sigma_wg"]
    ligne(N, nom, l) = resume[findfirst(r -> r[1] == N && r[2] == nom && r[3] == σ0 && r[4] == l, resume)]
    nomα(α) = α == 0 ? "aucune" : @sprintf("rampe %g %%", 100α)
    fig = Figure(size = (1200, 900), backgroundcolor = SURFACE, fontsize = 14)
    ax = [Axis(fig[i, j], backgroundcolor = SURFACE, xscale = (i, j) == (2, 1) ? log10 : identity,
               yscale = (i, j) == (2, 1) ? log10 : identity) for i in 1:2, j in 1:2]
    ax[1, 1].title = @sprintf("(a) biais sur les 30 dernières secondes, N = %g", Nv)
    ax[1, 1].xlabel, ax[1, 1].ylabel = "dérive α du gain jaune (%)", "biais moyen (mM)"
    for l in lect
        r = [ligne(Nv, nomα(α), l) for α in αs]
        errorbars!(ax[1, 1], 100 .* αs, getindex.(r, 7), getindex.(r, 8), color = COULEURS[l])
        scatterlines!(ax[1, 1], 100 .* αs, getindex.(r, 7), color = COULEURS[l], label = l)
    end
    hlines!(ax[1, 1], [0.0], color = :black, linewidth = 0.8)
    axislegend(ax[1, 1], position = :lb, framevisible = false)
    ax[1, 2].title = "(b) calibration du filtre : P₁₁ annoncé ÷ erreur quadratique réelle"
    ax[1, 2].xlabel, ax[1, 2].ylabel = "dérive α du gain jaune (%)", "P₁₁ / err²"
    for (N, st) in zip(p.toml["photons"]["balayage"], (:dot, :dash, :solid)), l in lect[3:4]
        lines!(ax[1, 2], 100 .* αs, [ligne(N, nomα(α), l)[12] for α in αs], color = COULEURS[l], linestyle = st,
               label = @sprintf("%s, N = %g", replace(l, "Kalman " => ""), N))
    end
    hlines!(ax[1, 2], [1.0], color = :black, linewidth = 0.8)
    axislegend(ax[1, 2], position = :rb, framevisible = false, nbanks = 2)
    ax[2, 1].title = @sprintf("(c) RMSE sur toute la course, rampe %g %%", 100p.toml["video"]["alpha"])
    ax[2, 1].xlabel, ax[2, 1].ylabel = "photons bleus par acquisition", "RMSE (mM)"
    Ns = p.toml["photons"]["balayage"]
    for l in lect
        scatterlines!(ax[2, 1], Ns, [ligne(N, nomα(p.toml["video"]["alpha"]), l)[5] for N in Ns], color = COULEURS[l], label = l)
    end
    axislegend(ax[2, 1], position = :rt, framevisible = false)
    ax[2, 2].title = @sprintf("(d) ĝ autour d'un saut de %g %% à %g s, N = %g (R complet)", 100p.derive.saut, p.derive.t_saut, Nv)
    ax[2, 2].xlabel, ax[2, 2].ylabel = "temps (s)", "g"
    σs = p.toml["kalman"]["balayage_sigma_wg"]
    for (i, σ) in enumerate(σs)
        T = saut[(Nv, σ)]
        lines!(ax[2, 2], T[:t], T[:g_chapeau], color = (COULEURS["Kalman R complet"], 0.35 + 0.65(i - 1) / max(length(σs) - 1, 1)),
               label = @sprintf("ĝ, σ_w,g = %.0e", σ))
    end
    lines!(ax[2, 2], saut[(Nv, σ0)][:t], saut[(Nv, σ0)][:g_vrai], color = :black, linewidth = 2, label = "vérité")
    axislegend(ax[2, 2], position = :lt, framevisible = false)
    save(chemin, fig)
end

"""
Distribution de l'écart Ĉ − C (mesuré − attendu) sur les 900 pas d'une course, un panneau par
technique : histogramme normalisé, gaussienne de même moyenne et même écart-type, et l'erreur
RMS. Pour J/B et 5a, les pas hors de la plage calibrée (échecs) sont exclus et comptés.
"""
function figure_ecarts(T, p::Params, chemin, titre)
    fig = Figure(size = (1500, 520), backgroundcolor = SURFACE, fontsize = 15)
    Label(fig[0, 1:3], titre, fontsize = 17, font = :bold)
    for (j, (nom, col, cle)) in enumerate((("J/B", "J/B", :C_ratio), ("5a", "5a", :C_5a), ("Kalman", "Kalman R complet", :C_kalman)))
        e = T[cle] .- T[:C_vrai]
        v = filter(!isnan, e)
        rms, m, s = sqrt(mean(v .^ 2)), mean(v), std(v)
        ax = Axis(fig[1, j], title = nom, xlabel = "écart Ĉ − C (mM)", ylabel = j == 1 ? "densité" : "", backgroundcolor = SURFACE)
        hist!(ax, v, bins = 40, normalization = :pdf, color = (COULEURS[col], 0.55), strokecolor = COULEURS[col], strokewidth = 0.5)
        x = range(m - 4s, m + 4s; length = 300)
        lines!(ax, x, exp.(-(x .- m) .^ 2 ./ (2s^2)) ./ (s * sqrt(2π)), color = :black, linewidth = 1.5)
        vlines!(ax, [0.0], color = :black, linestyle = :dash, linewidth = 0.8)
        txt = @sprintf("RMS = %.2f mM\nmoyenne = %+.2f mM\nσ = %.2f mM", rms, m, s)
        length(v) < length(e) && (txt *= @sprintf("\néchecs : %d / %d pas", length(e) - length(v), length(e)))
        text!(ax, 0.03, 0.97, text = txt, space = :relative, align = (:left, :top), fontsize = 14)
    end
    save(chemin, fig)
end
