# figures.jl — le rendu, après la simulation : les quatre figures du plan, la vérification de
# Fisher, une vidéo de la boucle fermée.

using CairoMakie

const SURFACE = colorant"#fcfcfb"
const COUL = Dict("fusion" => colorant"#1b7f79", "ratio" => colorant"#d9822b", "ratio_libre" => colorant"#b8336a",
                  "5a" => colorant"#7b5ea7", "consigne" => colorant"#888888")
axe(pos; kw...) = Axis(pos; backgroundcolor = SURFACE, kw...)

"""Moyenne glissante centrée sur `w` points (NaN ignorés)."""
lisse(x, w) = [moy(x[max(1, i - w ÷ 2):min(end, i + w ÷ 2)]) for i in eachindex(x)]

"""Bandes grises pour les phases de régime établi."""
bandes!(ax, p) = foreach(w -> vspan!(ax, w[1], w[2], color = (:black, 0.05)), p.toml["scenario"]["etabli"])

"""
Figure 1 — la figure de l'article. Rampe de gain de 20 %, N = 1000, moyenne des courses, lissée
sur 1 s. En haut, l'erreur de régulation que voit l'expérimentateur (consigne − Ĉ) ; en bas,
l'erreur réelle (consigne − C vrai).
"""
function figure_article(traces, p, chemin)
    fig = Figure(size = (1100, 760), backgroundcolor = SURFACE, fontsize = 15)
    Label(fig[0, 1], "Rampe de gain de 20 %, 1 000 photons bleus : le panneau du haut est ce qu'on observe, celui du bas est ce qui se passe",
          font = :bold, fontsize = 15)
    a1 = axe(fig[1, 1], ylabel = "consigne − Ĉ (mM)", title = "erreur de régulation affichée")
    a2 = axe(fig[2, 1], ylabel = "consigne − C vrai (mM)", xlabel = "temps (s)", title = "erreur réelle")
    linkxaxes!(a1, a2)
    w = round(Int, 1.0 / p.dt)
    for (nom, lab) in (("fusion", "fusion J/B + 5a"), ("ratio", "J/B seul, g figé"), ("ratio_libre", "J/B seul, g libre"))
        haskey(traces, nom) || continue
        res = traces[nom]
        t = res[1][1][:t]
        reg = mean([c[:consigne] .- T[nom][:C_hat] for (c, T) in res])
        vra = mean([c[:consigne] .- c[:C_vrai] for (c, T) in res])
        lines!(a1, t, lisse(reg, w), color = COUL[nom], linewidth = 2, label = lab, linestyle = nom == "ratio_libre" ? :dash : :solid)
        lines!(a2, t, lisse(vra, w), color = COUL[nom], linewidth = 2, label = lab, linestyle = nom == "ratio_libre" ? :dash : :solid)
    end
    for a in (a1, a2)
        bandes!(a, p)
        hlines!(a, [0.0], color = :black, linewidth = 0.8)
        ylims!(a, -14, 9)
    end
    axislegend(a2, position = :lb, framevisible = false)
    save(chemin, fig)
end

"""Figure 2 — la blancheur par phase : innovation moyenne du rapport contre β, transitoires et plateaux, deux budgets."""
function figure_blancheur(BO, p, β_ret, chemin)
    K = p.toml["kalman"]
    βs, σs, Ns = K["balayage_beta"], K["balayage_sigma_C0_diagnostic"], p.toml["photons"]["balayage"]
    fig = Figure(size = (1200, 480), backgroundcolor = SURFACE, fontsize = 14)
    x = [max(β, 3e-3) for β in βs]
    for (j, N) in enumerate(Ns)
        a = axe(fig[1, j], title = @sprintf("N = %g, boucle ouverte sans dérive", N), xlabel = "β (β = 0 placé à 3·10⁻³)",
                ylabel = "innovation moyenne du rapport (σ)", xscale = log10)
        hspan!(a, -0.05, 0.05, color = (:black, 0.07))
        hlines!(a, [-0.274], color = :black, linestyle = :dot, linewidth = 1)
        for (σ, st) in zip(σs, (:solid, :dash, :dot)), (cle, col, lab) in ((:innov_transitoires, COUL["ratio"], "transitoires"), (:innov_plateaux, COUL["fusion"], "plateaux"))
            v = [getfield(BO[findfirst(l -> l.N_photons_bleus == N && l.derive == "aucune" && l.filtre == @sprintf("fusion β=%g σC0=%g", β, σ), BO)], cle) for β in βs]
            scatterlines!(a, x, v, color = col, linestyle = st, label = @sprintf("%s, σ_C,0 = %g", lab, σ))
        end
        vlines!(a, [max(β_ret, 3e-3)], color = :black, linewidth = 0.8)
        j == 1 && axislegend(a, position = :lb, framevisible = false, labelsize = 11)
    end
    save(chemin, fig)
end

"""Figure 3 — le coude de couplage : les quatre signatures contre Kp, β retenu (plein) et β = 0 (tirets)."""
function figure_couplage(SIG, Kn, Kc, β_ret, p, chemin)
    cr = p.toml["criteres"]
    fig = Figure(size = (1300, 820), backgroundcolor = SURFACE, fontsize = 14)
    champs = ((:agitation, "agitation de la commande, σ(uₖ − uₖ₋₁)"), (:correlation_P11_agitation, "corrélation P₁₁ – |Δu| (max sur ±1 s)"),
              (:respiration_P11, "respiration de P₁₁ (crête à crête / moyenne)"), (:rho12_etabli, "ρ₁₂ en régime établi"))
    for (i, (c, titre)) in enumerate(champs)
        a = axe(fig[(i-1)÷2+1, (i-1)%2+1], title = titre, xlabel = "Kp / Kp nominal", xscale = log10)
        for β in sort(unique(x.beta for x in SIG if x.essai == "Kp"))
            L = filter(x -> x.essai == "Kp" && x.beta == β, SIG)
            lab = @sprintf("β = %g", β) * (β == β_ret ? " (retenu)" : "") * (β == 0 ? " (Q constant)" : "")
            scatterlines!(a, [x.Kp_sur_nominal for x in L], [getfield(x, c) for x in L], color = β == 0 ? COUL["fusion"] : COUL["ratio"],
                          linestyle = β == 0 ? :solid : :dash, label = lab)
        end
        vlines!(a, [1.0], color = :black, linewidth = 0.8)
        isnan(Kc) || vlines!(a, [Kc / Kn], color = COUL["ratio"], linewidth = 1.5)
        c === :correlation_P11_agitation && hlines!(a, [-cr["correlation_P11_agitation_max"], cr["correlation_P11_agitation_max"]], color = :black, linestyle = :dot)
        c === :respiration_P11 && hlines!(a, [cr["respiration_P11_max"]], color = :black, linestyle = :dot)
        i == 1 && axislegend(a, position = :lt, framevisible = false)
    end
    save(chemin, fig)
end

"""Figure 4 — l'observabilité conditionnelle de k_pompe : √P₄₄ et |C − C_repos| contre le temps ; k̂ contre la vérité."""
function figure_k_pompe(traces7, ph0, p, chemin)
    fig = Figure(size = (1200, 760), backgroundcolor = SURFACE, fontsize = 14)
    c, T = traces7[ph0.k]
    t = c[:t]
    a1 = axe(fig[1, 1], ylabel = "√P₄₄ (s⁻¹)", title = @sprintf("k_pompe vrai = %g s⁻¹ : l'incertitude ne baisse que loin du repos", ph0.k))
    a1b = Axis(fig[1, 1], yaxisposition = :right, ylabel = "|C − C_repos| (mM)", ygridvisible = false, xgridvisible = false)
    hidespines!(a1b); hidexdecorations!(a1b)
    lines!(a1b, t, abs.(c[:C_vrai] .- ph0.C_repos), color = (COUL["consigne"], 0.6), linewidth = 2)
    lines!(a1, t, sqrt.(T["fusion4"][:P44]), color = COUL["fusion"], linewidth = 2.5)
    a2 = axe(fig[2, 1], ylabel = "k̂_pompe (s⁻¹)", xlabel = "temps (s)", yscale = log10)
    linkxaxes!(a1, a2)
    for (k, (ck, Tk)) in sort(collect(traces7); by = first)
        kh, s = Tk["fusion4"][:k_hat], sqrt.(Tk["fusion4"][:P44])
        band!(a2, t, max.(kh .- s, 1e-3), kh .+ s, color = (COUL["fusion"], 0.15))
        lines!(a2, t, max.(kh, 1e-3), color = COUL["fusion"], linewidth = 2)
        hlines!(a2, [k], color = :black, linestyle = :dash)
    end
    save(chemin, fig)
end

"""Fisher contre Monte-Carlo, aux budgets de l'étape 2."""
function figure_fisher(verif, chemin)
    fig = Figure(size = (1200, 400), backgroundcolor = SURFACE, fontsize = 14)
    N = [v.N_photons_bleus for v in verif]
    for (j, (r, f, titre)) in enumerate(((:var_ratio_reelle, :var_ratio_fisher, "var(J/B)"), (:var_5a_reelle, :var_5a_fisher, "var(z_5a)"),
                                         (:cov_reelle, :cov_fisher, "cov(J/B, z_5a)")))
        a = axe(fig[1, j], title = titre, xlabel = "photons bleus par acquisition", xscale = log10, xticks = (N, [@sprintf("%g", n) for n in N]))
        scatterlines!(a, N, [getfield(v, f) for v in verif], color = COUL["fusion"], label = "Fisher")
        scatter!(a, N, [getfield(v, r) for v in verif], color = :black, label = "Monte-Carlo")
        j == 1 && axislegend(a, position = :rt, framevisible = false)
    end
    save(chemin, fig)
end

"""Vidéo de la boucle fermée sous rampe de 20 % : fusion et J/B seul côte à côte (répétition 1)."""
function video_boucle(traces, p, chemin)
    (cf, Tf), (cr_, Tr) = traces["fusion"][1], traces["ratio"][1]
    t, n = cf[:t], length(cf[:t])
    fig = Figure(size = (1100, 900), backgroundcolor = SURFACE, fontsize = 15)
    k = Observable(1)
    pts(x, v) = @lift Point2f.(x[1:$k], v[1:$k])
    for (i, (c, T, nom, titre)) in enumerate(((cf, Tf, "fusion", "fusion J/B + 5a"), (cr_, Tr, "ratio", "J/B seul, g figé")))
        a = axe(fig[i, 1], ylabel = "[Cl⁻] (mM)", title = titre)
        lines!(a, t, c[:consigne], color = COUL["consigne"], linestyle = :dash, label = "consigne")
        lines!(a, pts(t, c[:C_vrai]), color = :black, linewidth = 2, label = "C vrai")
        lines!(a, pts(t, T[nom][:C_hat]), color = COUL[nom], linewidth = 2, label = "Ĉ (ce que voit le régulateur)")
        xlims!(a, 0, t[end] + p.dt); ylims!(a, 0, 40)
        axislegend(a, position = :rt, framevisible = false)
    end
    a3 = axe(fig[3, 1], ylabel = "commande u", xlabel = "temps (s)")
    lines!(a3, pts(t, cf[:u]), color = COUL["fusion"], label = "fusion")
    lines!(a3, pts(t, cr_[:u]), color = COUL["ratio"], label = "J/B seul")
    xlims!(a3, 0, t[end] + p.dt); ylims!(a3, 0, 1)
    axislegend(a3, position = :rt, framevisible = false)
    record(fig, chemin * ".mp4", 1:n; framerate = p.toml["video"]["ips"]) do i
        k[] = i
    end
    save(chemin * ".png", fig)
end
