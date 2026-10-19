# boucle.jl — une course de 90 s : la cellule, les deux histogrammes, les lectures, les filtres,
# le régulateur. Simule tout d'abord, rend ensuite.
#
# Un filtre « pilote » fournit Ĉ au régulateur (boucle fermée) ; `pilote = :oracle` fait réguler
# sur le vrai C, ce qui donne la boucle ouverte des courses 1 et 2 : la commande ne dépend alors
# d'aucun filtre, et tous les filtres « passifs » voient exactement les mêmes histogrammes.

"""Spécification d'un filtre : nom, modèle (:fusion, :fusion4, :ratio, :ratio_libre), β du Q(u), σ_C,0 (défaut : params)."""
spec(nom, modele, β) = (nom = nom, modele = modele, β = Float64(β))
spec(nom, modele, β, σC0) = (nom = nom, modele = modele, β = Float64(β), σC0 = Float64(σC0))

function filtre(p::Params, lec, ph::Physio, s)
    k = p.toml["kalman"]
    cR, cS = lec.JB.courbe, lec.S.courbe
    x0 = [k["C0_mM"], k["g0"]]
    σ0 = [k["sigma_C0_init_mM"], k["sigma_g0"]]
    σC0 = get(s, :σC0, k["sigma_C0_mM"])
    rat = s.modele in (:ratio, :ratio_libre)
    s.modele === :ratio && (σ0[2] = 0.0)                     # g figé à 1 : la calibration unique
    σw = [s.modele === :ratio ? 0.0 : s.modele === :ratio_libre ? k["sigma_wg_ratio_libre"] : k["sigma_wg"]]
    if !rat
        push!(x0, 0.0); push!(σ0, k["sigma_gp0"]); push!(σw, k["sigma_wgp"])
    end
    if s.modele === :fusion4
        push!(x0, k["k_pompe_init_par_s"]); push!(σ0, k["sigma_k_init_par_s"]); push!(σw, k["sigma_wk_par_s"])
    end
    return EKF(s.modele, x0, Matrix(Diagonal(σ0 .^ 2)), ph, σC0, s.β, σw, 1 - p.dt / k["tau_ou_s"],
               Function[c -> valeur(cR, c), c -> valeur(cS, c)], Function[c -> derivee(cR, c), c -> derivee(cS, c)],
               rat ? [true, false] : [true, true])
end

const TRACES = [:C_hat, :g_hat, :P11, :P22, :rho12, :k_hat, :P44, :innov_ratio, :innov_5a]

"""
    course(p, lec, ph, r; pilote, passifs, reg, ph_filtre) -> (commun, par filtre)

`pilote` : une spécification de filtre, ou `:oracle`. `passifs` : d'autres filtres qui
regardent sans agir. `reg` : le régulateur (PI). `ph_filtre` : la physiologie que les filtres
croient (par défaut, la vraie).
"""
function course(p::Params, lec, ph::Physio, r; pilote, passifs = [], reg::PI, ph_filtre::Physio = ph)
    ts = collect(temps(p))
    n = length(ts)
    sc = p.toml["scenario"]
    specs = pilote === :oracle ? passifs : [pilote; passifs]
    F = Dict(s.nom => filtre(p, lec, ph_filtre, s) for s in specs)
    T = Dict(s.nom => Dict(c => fill(NaN, n) for c in TRACES) for s in specs)
    com = Dict(c => fill(NaN, n) for c in (:t, :C_vrai, :g_vrai, :consigne, :u, :C_ratio, :C_5a, :z_ratio, :z_5a))
    C, u = ph.C_repos, 0.0
    for (k, t) in enumerate(ts)
        g = derive(t, p)
        hb, hj = histogrammes_c(max(C, 0.0), g, p; rng = rng_pas(k, r, p))
        lR, lS = estimer(hb, hj, lec.JB.cfg), estimer(hb, hj, lec.S.cfg)
        for s in specs
            f = F[s.nom]
            predire!(f, u)
            ν = mettre_a_jour!(f, [lR.z, lS.z], R_attendue(lec, f.x[1], f.x[2]), [lR.ok, lS.ok])
            d = diagnostics(f)
            for (c, v) in zip(TRACES, (d.C, d.g, d.P11, d.P22, d.rho12, d.k, d.P44, ν[1], ν[2]))
                T[s.nom][c][k] = v
            end
        end
        cons = consigne(t, p)
        Ĉ, P11 = pilote === :oracle ? (C, 0.0) : (T[pilote.nom][:C_hat][k], T[pilote.nom][:P11][k])
        u = t < sc["t_activation_s"] - 1e-9 ? 0.0 : commande!(reg, cons - Ĉ, P11, p.dt)
        for (c, v) in zip((:t, :C_vrai, :g_vrai, :consigne, :u, :z_ratio, :z_5a), (t, C, g, cons, u, lR.z, lS.z))
            com[c][k] = v
        end
        com[:C_ratio][k] = lR.ok ? inverser(lec.JB.courbe, lR.z) : NaN
        com[:C_5a][k] = lS.ok ? inverser(lec.S.courbe, lS.z) : NaN
        C = pas_cellule(C, u, ph, rng_cellule(k, r, p))
    end
    return com, T
end
