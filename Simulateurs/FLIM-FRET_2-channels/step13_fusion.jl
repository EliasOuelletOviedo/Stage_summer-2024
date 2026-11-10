# Étape 13 — fusionner le ratio et la forme (changement 4) ; extension de la figure 7.
#
# Le gain du canal jaune relatif au bleu dérive (marche aléatoire sur log g, 5 %
# d'écart-type après τ_g, τ_g de 1 min à 1 h). [Cl⁻] varie lentement autour de
# 20 mM (20 + 8·sin(2πt/300 s)), échantillonnée à 10 Hz pendant 30 min, flux de
# photons fixe (normalisation A). Trois filtres de Kalman sur les mêmes mesures :
#   ratio seul   état [Cl⁻], gain supposé égal à 1 (ce qu'on fait aujourd'hui) ;
#   forme seule  état [Cl⁻], lecture 5a (insensible au gain) ;
#   fusion       état ([Cl⁻], log g), les deux lectures (estimateur 7).
# Plus une réponse à un saut de gain de 5 % ([Cl⁻] constante) : combien de temps
# la forme met-elle à discipliner le ratio ? Deux capteurs : mCer3_Topaz (E = 0,96)
# et le capteur hypothétique à τ_DA = 2,1 ns (E = 0,43) de la figure 8.
#
#   julia --project=. -t 8 step13_fusion.jl
#
# Sorties : sorties/resultats/fusion_{resume,trajectoire,saut}.csv, sorties/figures/fig7b_fusion.png

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, CSV, DataFrames, StableRNGs

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
ctxs = S.creer_contextes(m, irf_mat; avec_utilisateur = false)
dt = 0.1
T = 1800.0
n = round(Int, T / dt)
temps = collect(0:n-1) .* dt
C_traj(t) = 20.0 + 8.0 * sin(2π * t / 300)
# Bruit de processus de [Cl⁻] : 1 mM/s (l'échelon de 5 à 25 mM en 20 s de la figure 7). Il doit
# dépasser, en unités de log R, celui du gain (q_ℓ) : sinon le filtre attribue au gain toute
# variation du ratio et jette sa précision (avec 0,2 mM/s, la fusion ne faisait pas mieux que la forme).
q_C = (1.0 * dt)^2
capteurs = [("mCer3_Topaz (E = 0,96)", m), ("hypothétique, τ_DA = 2,1 ns (E = 0,43)", S.remplacer(m; p = S.remplacer(m.p; τ_DA = 2.1)))]
e7s = [S.DeuxLectures(mc, ctxs[1]; C_max) for (_, mc) in capteurs]
for ((nomc, mc), e7) in zip(capteurs, e7s)
    i = findfirst(==(20.0), e7.C_tab)
    @printf("%s : à 20 mM et 10⁴ photons bleus, σ ratio %.2f mM, σ forme (5a) %.2f mM (rapport %.1f)\n", nomc, e7.σ_ratio[i], e7.σ_forme[i], e7.σ_forme[i] / e7.σ_ratio[i])
end

"""Suite des lectures (estimateur 7) le long d'une trajectoire, gain g(t) appliqué au jaune."""
function mesures(e7, mc, N, Cs, g, rng, ctx)
    nb = mc.inst.n_bins
    nB, nY = zeros(nb), zeros(nb)
    out = Vector{NTuple{8,Float64}}(undef, length(Cs))
    for k in eachindex(Cs)
        μB, μY = S.attendus(Cs[k], mc, N, :A, ctx.ws)
        μY .*= g[k]
        S.tirer_deux!(nB, nY, μB, μY, rng)
        out[k] = S.lectures(e7, nB, nY, ctx)
    end
    return out
end

taches = [(ic, τg, N, r) for ic in eachindex(capteurs) for τg in (60.0, 300.0, 1200.0, 3600.0) for N in (1e3, 1e4) for r in 1:6]
resultats = Vector{Any}(undef, length(taches))
trajectoire = Ref{Any}(nothing)
suivant = Threads.Atomic{Int}(1)
t_calc = @elapsed @sync for ctx in ctxs
    Threads.@spawn while (i = Threads.atomic_add!(suivant, 1)) <= length(taches)
        ic, τg, N, r = taches[i]
        mc, e7 = capteurs[ic][2], e7s[ic]
        rng = StableRNG(S.graine(13_000_000 + 1000ic + round(Int, τg) + round(Int, log10(N)), r))
        g = S.derive_gain(rng, n, dt; τ_g = τg)
        Cs = C_traj.(temps)
        mes = mesures(e7, mc, N, Cs, g, rng, ctx)
        q_ℓ = 0.05^2 / τg * dt
        lignes = NamedTuple[]
        deux = n ÷ 2 + 1:n                                   # seconde moitié : filtres stabilisés
        sorties = Dict{Symbol,Any}()
        for mode in (:ratio, :forme, :fusion)
            Ĉ, ℓ̂, _ = S.filtre_fusion(e7, mes, dt; q_C, q_ℓ, mode)
            sorties[mode] = (Ĉ, ℓ̂)
            err = Ĉ[deux] .- Cs[deux]
            push!(lignes, (capteur = capteurs[ic][1], tau_g_s = τg, N = N, realisation = r, filtre = String(mode),
                RMSE_mM = sqrt(mean(err .^ 2)), biais_mM = mean(err), RMSE_gain_pct = mode === :fusion ?
                100 * sqrt(mean((ℓ̂[deux] .- log.(g[deux])) .^ 2)) : NaN, derive_gain_pct = 100 * std(log.(g))))
        end
        # lectures brutes (sans filtre) pour référence
        Cr = [t[1] for t in mes]
        push!(lignes, (capteur = capteurs[ic][1], tau_g_s = τg, N = N, realisation = r, filtre = "ratio brut",
            RMSE_mM = sqrt(mean((Cr[deux] .- Cs[deux]) .^ 2)), biais_mM = mean(Cr[deux] .- Cs[deux]), RMSE_gain_pct = NaN,
            derive_gain_pct = 100 * std(log.(g))))
        resultats[i] = lignes
        if ic == 1 && τg == 300.0 && N == 1e4 && r == 1
            k = 1:10:n
            trajectoire[] = DataFrame(t = temps[k], C_vrai = Cs[k], g_vrai = g[k], ratio = sorties[:ratio][1][k],
                forme = sorties[:forme][1][k], fusion = sorties[:fusion][1][k], g_fusion = exp.(sorties[:fusion][2][k]),
                lecture_ratio = Cr[k], lecture_forme = [t[3] for t in mes[k]])
        end
    end
end
res = DataFrame(reduce(vcat, resultats))
CSV.write(S.dossier_sorties("resultats", "fusion_resume.csv"), res)
CSV.write(S.dossier_sorties("resultats", "fusion_trajectoire.csv"), trajectoire[])
@printf("%d séries de 30 min en %.0f s\n", length(taches), t_calc)
println("\nRMSE sur [Cl⁻] (mM, seconde moitié ; moyenne des réalisations) :")
for g in groupby(res, [:capteur, :N, :tau_g_s]; sort = true)
    @printf("  %-40s N = %-6g τ_g = %5.0f s :", g.capteur[1], g.N[1], g.tau_g_s[1])
    for f in ("ratio brut", "ratio", "forme", "fusion")
        d = g[g.filtre.==f, :]
        @printf("  %s %6.2f (biais %+5.2f)", f, mean(d.RMSE_mM), mean(d.biais_mM))
    end
    @printf("  | gain : dérive %.1f %%, erreur fusion %.2f %%\n", mean(g.derive_gain_pct), mean(filter(isfinite, g.RMSE_gain_pct)))
end

# --- réponse à un saut de gain de 5 % ([Cl⁻] = 20 mM constante) ----------------------
saut = DataFrame()
Tsaut, t_saut = 240.0, 60.0
ns = round(Int, Tsaut / dt)
ts = collect(0:ns-1) .* dt
for (ic, (nomc, mc)) in enumerate(capteurs), N in (1e3, 1e4)
    e7 = e7s[ic]
    acc = zeros(ns)
    R = 40
    for r in 1:R
        rng = StableRNG(S.graine(14_000_000 + 10ic + round(Int, log10(N)), r))
        g = [t < t_saut ? 1.0 : 1.05 for t in ts]
        mes = mesures(e7, mc, N, fill(20.0, ns), g, rng, ctxs[1])
        Ĉ, _, _ = S.filtre_fusion(e7, mes, dt; q_C = q_C, q_ℓ = 0.05^2 / 300 * dt, mode = :fusion)
        acc .+= (Ĉ .- 20.0) ./ R
    end
    # biais du ratio seul juste après le saut, et temps pour que la fusion en revienne à 1/e
    br = mean(begin
        rng = StableRNG(S.graine(15_000_000 + 10ic + round(Int, log10(N)), 1))
        mes = mesures(e7, mc, N, fill(20.0, ns), [t < t_saut ? 1.0 : 1.05 for t in ts], rng, ctxs[1])
        Ĉr, _, _ = S.filtre_fusion(e7, mes, dt; q_C = q_C, q_ℓ = 0.0, mode = :ratio)
        Ĉr[ts.>=t_saut+30] .- 20.0
    end)
    # moyenne glissante d'une seconde, pic cherché après le saut, retour à 1/e du pic
    lisse = [mean(acc[max(1, i - 5):min(ns, i + 5)]) for i in 1:ns]
    apres = findall(ts .>= t_saut)
    pic = apres[argmax(abs.(lisse[apres]))]
    seuil = abs(lisse[pic]) / ℯ
    k = findfirst(i -> i > pic && abs(lisse[i]) <= seuil, eachindex(lisse))
    τr = k === nothing ? NaN : ts[k] - t_saut
    acc = lisse
    for (t, b) in zip(ts[1:5:end], acc[1:5:end])
        push!(saut, (capteur = nomc, N = N, t = t, biais_fusion_mM = b, biais_ratio_seul_mM = t >= t_saut ? br : 0.0, temps_retour_s = τr))
    end
    @printf("saut de gain de 5 %% — %-40s N = %-6g : biais du ratio seul %+.2f mM ; fusion : pic %+.2f mM, retour à 1/e en %.1f s (%s mesures)\n",
        nomc, N, br, acc[pic], τr, isfinite(τr) ? string(round(Int, τr / dt)) : "—")
end
CSV.write(S.dossier_sorties("resultats", "fusion_saut.csv"), saut)

include(joinpath(@__DIR__, "src", "Figures.jl"))
using .Figures, CairoMakie
set_theme!(Figures.theme_sim())
if isdefined(Figures, :figure_fusion)
    Figures.figure_fusion(res, trajectoire[], saut, S.dossier_sorties("figures", "fig7b_fusion.png"))
end
