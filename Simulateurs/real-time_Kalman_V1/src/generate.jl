# generate.jl — trajectoire du Cl⁻, dérive du gain jaune, histogrammes des deux canaux.
#
# Six étapes, de l'instant t aux deux histogrammes entiers : concentration, amplitudes
# (Langmuir), périodisation analytique, convolution avec l'IRF du canal, mise à l'échelle
# (la dérive g n'entre qu'ici, et seulement sur le jaune), bruit de Poisson.

using Distributions, StableRNGs

"""C(t) : repos, montée exponentielle de t_a à t_b, puis extrusion par KCC2 vers le repos."""
function concentration(t, p::Params)
    r = p.toml["trajectoire"]
    Cr, ΔC, ta, tb = r["C_repos_mM"], r["delta_C_mM"], r["t_debut_montee_s"], r["t_fin_montee_s"]
    t < ta && return Cr
    t < tb && return Cr + ΔC * (1 - exp(-(t - ta) / r["tau_montee_s"]))
    return Cr + ΔC * (1 - exp(-(tb - ta) / r["tau_montee_s"])) * exp(-(t - tb) / r["tau_KCC2_s"])
end

temps(p::Params) = (0:round(Int, p.toml["trajectoire"]["duree_s"] / p.dt)-1) .* p.dt

"""g(t), facteur multiplicatif du canal jaune, inconnu de toutes les lectures."""
function derive(t, p::Params)
    d = p.derive
    d.forme == "rampe" && return 1 + d.alpha * t / p.toml["trajectoire"]["duree_s"]
    d.forme == "saut" && return t >= d.t_saut ? 1 + d.saut : 1.0
    error("dérive inconnue : $(d.forme)")
end

"""
Composantes de chaque canal en C : la fraction d'intensité `d` de la diffusion (copie de
l'IRF) et les couples (α, τ) des deux populations, α en fractions d'amplitude.
"""
function expsums(c, p::Params)
    L, τ = p.lang, p.τ
    pb, qj = L["p_bleu"](c), L["q_jaune"](c)
    return (d = L["d_bleu"](c), es = [(pb, τ.τ2), (1 - pb, τ.τ3)]),
           (d = L["d_jaune"](c), es = [(qj, τ.τJ), (1 - qj, τ.τA)])
end

"""
Photons d'un déclin périodisé Σα·e^{-t/τ}/(1 − e^{-T/τ}) (série géométrique des impulsions
précédentes sommée analytiquement) reçus entre a et a + Δ après une impulsion, Δ < T.
"""
function integrale_periodique(a, α, τ, Δ, T)
    a = mod(a, T)
    q = α * τ / -expm1(-T / τ)
    a + Δ <= T && return q * (exp(-a / τ) - exp(-(a + Δ) / τ))
    return q * (exp(-a / τ) - exp(-T / τ) + -expm1(-(a + Δ - T) / τ))
end

"""
IRF ⊛ déclin périodisé, par FFT sur 2n points (plan calculé une fois). Le noyau couvre les
retards −(n−1)…(n−1) et il est périodique en T : un retard négatif donne la fin de
l'impulsion précédente. La période (12,5 ns) ne remplit pas les n canaux du TAC (16,4 ns),
d'où 2n et non n : une convolution circulaire sur n replierait à 16,4 ns. Rend la fenêtre.
"""
function convoluer(es, irf, t0, p::Params)
    n = p.n
    f = zeros(2n)
    for i in 0:2n-1
        m = i < n ? i : i - 2n
        f[i+1] = sum(integrale_periodique(m * p.pas - t0, α, τ, p.pas, p.T) for (α, τ) in es)
    end
    μ = irfft((p.plan * vcat(irf, zeros(n))) .* (p.plan * f), 2n)[1:n]
    μ[setdiff(1:n, p.fenetre)] .= 0.0
    return μ
end

"""Un canal : N photons dans la fenêtre, dont la fraction d en diffusion."""
function echelle(canal, irf, t0, N, p::Params)
    pop = convoluer(canal.es, irf, t0, p)
    dif = convoluer([(1.0, p.τ.diff)], irf, t0, p)
    W = p.fenetre
    return N .* ((1 - canal.d) .* pop ./ sum(pop[W]) .+ canal.d .* dif ./ sum(dif[W]))
end

"""
    histogrammes_c(c, g, p; rng) -> (hb, hj)

Les deux histogrammes à la concentration c et au gain jaune g. Sans `rng` : les comptages
attendus (fond compris), sans tirage — le mode qui construit les calibrations et sert aux tests.
"""
function histogrammes_c(c, g, p::Params; rng = nothing)
    eb, ej = expsums(c, p)
    μb = echelle(eb, p.irf_b, p.t0[1], p.N, p) .+ p.fond_b
    μj = echelle(ej, p.irf_j, p.t0[2], p.N * jb(c, p) * g, p) .+ p.fond_j
    (any(<(0), μb) || any(<(0), μj)) && error("comptage attendu négatif à C = $c mM : paramètres incohérents")
    rng === nothing && return μb, μj
    return Float64.(rand.(rng, Poisson.(μb))), Float64.(rand.(rng, Poisson.(μj)))
end

histogrammes(t, p::Params; rng = nothing) = histogrammes_c(concentration(t, p), derive(t, p), p; rng)

"""Le générateur du pas k de la répétition r : la graine ne dépend que de (r, k), un pas se rejoue seul."""
rng_pas(k, r, p::Params) = StableRNG(p.toml["course"]["graine"] + 10_000r + k)
