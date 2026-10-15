# Étape 7 (optionnelle) — la cible mobile.
#
# A. Boucle ouverte : [Cl⁻] suit un échelon de 5 à 25 mM avec une constante de
#    temps de 20 s, échantillonné à 10 Hz (normalisation A : le flux de photons
#    est fixe, 10⁴ photons bleus par échantillon à 20 mM, ton point de
#    fonctionnement). Chaque estimateur suit la trajectoire ; un filtre de
#    Kalman scalaire (marche aléatoire) lisse le ratio avec R tiré de la borne
#    de Cramér-Rao évaluée à l'estimation courante (usage 2 du plan).
# B. Boucle fermée : un intégrateur proportionnel-intégral commande un
#    actionneur qui fait entrer du Cl⁻ (eNpHR, puissance ≥ 0) pour tenir
#    20 mM ; la cellule relaxe vers 5 mM avec 20 s de constante de temps. Un
#    intégrateur amène l'ESTIMATION à la consigne : l'erreur statique sur la
#    vraie [Cl⁻] est l'opposé du biais de l'estimateur à la consigne.
#
#   julia --project=. -t 8 step7_cible_mobile.jl
#
# Sorties : sorties/resultats/cible_mobile.csv, sorties/figures/fig7_cible_mobile.png

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf, Statistics, CSV, DataFrames, StableRNGs

m, irf_mat, cfg = S.charger_configuration()
C_max = Float64(cfg["experience"]["C_max"])
ctxs = S.creer_contextes(m, irf_mat)
ws = ctxs[1].ws
ests = S.estimateurs_standard(m, ws, irf_mat; C_max, liste = ["1", "2", "2H", "4", "5a", "5d"])
cals = S.calibrer_tous(ests, m, ctxs[1]; C_max)

dt, T_tot = 0.1, 120.0
temps = collect(0:dt:T_tot-dt)
C_traj(t) = t < 10 ? 5.0 : 25.0 - 20.0 * exp(-(t - 10) / 20)

"""Estimation de chaque estimateur sur un histogramme tiré à [Cl⁻] = C (germe fixé)."""
function estimer_tous(C, N, graine_, ctx)
    μB, μY = S.attendus(C, m, N, :A, ctx.ws)
    nB, nY = zeros(length(μB)), zeros(length(μY))
    S.tirer_deux!(nB, nY, μB, μY, StableRNG(graine_))
    out = Dict{String,Float64}()
    obs = Dict{Any,Tuple{Float64,Symbol}}()
    for (k, e) in ests
        o, st = get!(() -> S.observable(e, nB, nY, ctx), obs, e)
        out[k] = st === :echec ? NaN : S.inverser(cals[k], o)[1]
    end
    return out
end

# --- A. boucle ouverte --------------------------------------------------------
lignes = []
for N in (1e3, 1e4)
    res = Vector{Dict{String,Float64}}(undef, length(temps))
    suivant = Threads.Atomic{Int}(1)
    @sync for ctx in ctxs
        Threads.@spawn while (i = Threads.atomic_add!(suivant, 1)) <= length(temps)
            res[i] = estimer_tous(C_traj(temps[i]), N, S.graine(8_000_000 + round(Int, log10(N)), i), ctx)
        end
    end
    # Kalman scalaire sur le ratio : marche aléatoire, q adapté à la pente maximale de la trajectoire
    q = (20.0 / 20.0 * dt)^2
    x̂, P = res[1]["1"], 100.0
    for (i, t) in enumerate(temps)
        z = res[i]["1"]
        R = S.crlb_mM(:totaux, clamp(x̂, 0.5, C_max), m, N, :A, ws)^2
        P += q
        Kg = P / (P + R)
        x̂ += Kg * (z - x̂)
        P *= (1 - Kg)
        push!(lignes, (boucle = "ouverte", N = N, t = t, C_vrai = C_traj(t), consigne = NaN, kalman_ratio = x̂,
            (Symbol("est_" * k) => res[i][k] for (k, _) in ests)...))
    end
end

# --- B. boucle fermée -------------------------------------------------------------
C_repos, τ_p, gain, u_max = 5.0, 20.0, 1.75, 1.0     # C stationnaire = 5 + 20·1,75·u : 40 mM à pleine puissance
Kp, Ki = 0.01, 0.004                                  # par mM et par mM·s
consigne = 20.0
T_fb = 240.0
temps_fb = collect(0:dt:T_fb-dt)
for N in (1e3, 1e4), (k, _) in ests
    C, integ = C_repos, 0.0
    for (i, t) in enumerate(temps_fb)
        Ĉ = estimer_tous(C, N, S.graine(9_000_000 + round(Int, log10(N)), i), ctxs[1])[k]
        isnan(Ĉ) && (Ĉ = consigne)        # un échec ne commande rien
        e = consigne - Ĉ
        u = Kp * e + Ki * integ
        if 0 <= u <= u_max                 # anti-emballement : on n'intègre pas en saturation
            integ += e * dt
        end
        u = clamp(u, 0, u_max)
        push!(lignes, (boucle = "fermée", N = N, t = t, C_vrai = C, consigne = consigne, kalman_ratio = NaN,
            (Symbol("est_" * kk) => (kk == k ? Ĉ : NaN) for (kk, _) in ests)...))
        C += dt * (-(C - C_repos) / τ_p + gain * u)
    end
end
df = DataFrame(lignes)
CSV.write(S.dossier_sorties("resultats", "cible_mobile.csv"), df)

# --- résumé -------------------------------------------------------------------------
println("Boucle ouverte, RMSE sur toute la trajectoire (mM) :")
for N in (1e3, 1e4)
    d = df[(df.boucle.=="ouverte").&(df.N.==N), :]
    @printf("  N = %g :", N)
    for (k, _) in ests
        v = d[!, "est_"*k]
        ok = isfinite.(v)
        @printf("  %s %.2f", k, sqrt(mean((v[ok] .- d.C_vrai[ok]) .^ 2)))
    end
    @printf("  | Kalman(ratio) %.2f\n", sqrt(mean((d.kalman_ratio .- d.C_vrai) .^ 2)))
end
println("Boucle fermée, consigne 20 mM : [Cl⁻] vraie sur les 120 dernières secondes (moyenne ± écart-type) :")
for N in (1e3, 1e4), (k, _) in ests
    d = df[(df.boucle.=="fermée").&(df.N.==N).&isfinite.(df[!, "est_"*k]).&(df.t.>=120), :]
    @printf("  N = %-6g %-3s  %6.2f ± %5.2f mM (erreur statique %+6.2f mM)\n", N, k, mean(d.C_vrai), std(d.C_vrai), mean(d.C_vrai) - consigne)
end

include(joinpath(@__DIR__, "src", "Figures.jl"))
using .Figures, CairoMakie
set_theme!(Figures.theme_sim())
Figures.figure_cible_mobile(df, S.dossier_sorties("figures", "fig7_cible_mobile.png"))
