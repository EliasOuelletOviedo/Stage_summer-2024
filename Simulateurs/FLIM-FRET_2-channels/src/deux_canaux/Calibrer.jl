# Calibrer.jl — étape 4 : une courbe de calibration par lecture, et la borne de Fisher.
#
# Chaque lecture reçoit sa propre courbe, construite sur les 8 points à pleins photons
# (les sommes), puis inversée pour chaque tirage : c'est ce que l'on ferait au banc.
#   - monotonie vérifiée avant d'inverser : une lecture non monotone est marquée, pas
#     interpolée ;
#   - un tirage dont l'observable sort de l'intervalle calibré compte comme un échec, pas
#     comme une valeur extrapolée ;
#   - l'incertitude des 8 points eux-mêmes (photons finis) est mesurée sur la dispersion
#     des 10 acquisitions de chaque cellule et propagée en mM par la pente de la courbe.
# Interpolation monotone par morceaux (Hermite cubique de Fritsch–Carlson, « PCHIP »),
# exacte aux nœuds. Une Langmuir ajustée sur chaque courbe donne le K apparent.

"""Courbe de calibration d'une lecture à une longueur d'onde."""
struct Courbe
    lecture::String
    λ::Int
    c::Vector{Float64}       # nœuds (mM)
    o::Vector{Float64}       # observable à pleins photons
    d::Vector{Float64}       # pentes de Hermite aux nœuds
    monotone::Bool
    σ_ref::Vector{Float64}   # incertitude de chaque nœud (dispersion des 10 acquisitions / √10)
end

"""Pentes de Fritsch–Carlson (moyenne harmonique pondérée ; 0 aux extrema locaux)."""
function pentes_pchip(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    δ = diff(y) ./ h
    d = zeros(n)
    n == 2 && return fill(δ[1], 2)
    for i in 2:n-1
        if δ[i-1] * δ[i] <= 0
            d[i] = 0.0
        else
            w1, w2 = 2h[i] + h[i-1], h[i] + 2h[i-1]
            d[i] = (w1 + w2) / (w1 / δ[i-1] + w2 / δ[i])
        end
    end
    # bords : formule à trois points, bornée pour garder la monotonie
    bord(h1, h2, δ1, δ2) = begin
        v = ((2h1 + h2) * δ1 - h1 * δ2) / (h1 + h2)
        sign(v) != sign(δ1) ? 0.0 : (sign(δ1) != sign(δ2) && abs(v) > 3abs(δ1)) ? 3δ1 : v
    end
    d[1] = bord(h[1], h[2], δ[1], δ[2])
    d[n] = bord(h[n-1], h[n-2], δ[n-1], δ[n-2])
    return d
end

"""
    courbe(lecture, λ, c, o, σ_ref) -> Courbe

Monotone si les 7 écarts successifs ont tous le même signe (aucun nul ni manquant).
"""
function courbe(lecture::String, λ::Int, c::Vector{Float64}, o::Vector{Float64}, σ_ref::Vector{Float64})
    p = sortperm(c)
    c, o, σ_ref = c[p], o[p], σ_ref[p]
    Δ = diff(o)
    mono = all(isfinite, o) && (all(>(0), Δ) || all(<(0), Δ))
    d = all(isfinite, o) ? pentes_pchip(c, o) : fill(NaN, length(c))   # aussi pour une courbe non monotone (prédiction indicative)
    return Courbe(lecture, λ, c, o, d, mono, σ_ref)
end

"""
Plus forte violation de la monotonie, en écarts-types : pour chaque pas de sens opposé au
sens global (o₁₃₁ − o₀), |Δo| / √(σᵢ² + σᵢ₊₁²). Sous 2, la violation reste dans le bruit
de la calibration. Rend (z, « cᵢ–cᵢ₊₁ »), ou (0, "") pour une courbe monotone.
"""
function violation_monotonie(C::Courbe)
    s = sign(C.o[end] - C.o[1])
    z, ou = 0.0, ""
    for i in 1:length(C.c)-1
        Δ = C.o[i+1] - C.o[i]
        if sign(Δ) != s
            v = abs(Δ) / sqrt(C.σ_ref[i]^2 + C.σ_ref[i+1]^2)
            if v > z
                z, ou = v, "$(C.c[i])–$(C.c[i+1])"
            end
        end
    end
    return z, ou
end

"""Valeur et pente de la courbe en x (dans la plage des nœuds)."""
function evaluer(C::Courbe, x::Real)
    i = clamp(searchsortedlast(C.c, x), 1, length(C.c) - 1)
    h = C.c[i+1] - C.c[i]
    t = (x - C.c[i]) / h
    h00, h10, h01, h11 = 2t^3 - 3t^2 + 1, t^3 - 2t^2 + t, -2t^3 + 3t^2, t^3 - t^2
    v = h00 * C.o[i] + h10 * h * C.d[i] + h01 * C.o[i+1] + h11 * h * C.d[i+1]
    dv = ((6t^2 - 6t) * C.o[i] + (3t^2 - 4t + 1) * h * C.d[i] + (-6t^2 + 6t) * C.o[i+1] + (3t^2 - 2t) * h * C.d[i+1]) / h
    return v, dv
end

"""
    inverser(C, o) -> (ĉ en mM, code)

Code `OK`, `HORS_PLAGE` (observable hors de l'intervalle des 8 points) ou `NON_MONOTONE`.
Dichotomie dans l'intervalle de nœuds qui encadre o (la courbe y est monotone).
"""
function inverser(C::Courbe, o::Real)
    C.monotone || return (NaN, NON_MONOTONE)
    isfinite(o) || return (NaN, AUTRE)
    lo, hi = extrema(C.o)
    (o < lo || o > hi) && return (NaN, HORS_PLAGE)
    sens = C.o[end] > C.o[1] ? 1 : -1
    i = 1
    while i < length(C.c) - 1 && sens * (o - C.o[i+1]) > 0
        i += 1
    end
    a, b = C.c[i], C.c[i+1]
    for _ in 1:60
        m = (a + b) / 2
        v, _ = evaluer(C, m)
        if sens * (v - o) < 0
            a = m
        else
            b = m
        end
    end
    return ((a + b) / 2, OK)
end

"""
Incertitude de calibration en mM à chaque nœud : σ_ref / |pente|. Pente de Hermite aux
nœuds intérieurs ; pente du segment voisin aux deux bornes (la formule de bord de
Fritsch–Carlson peut y valoir zéro).
"""
function sigma_calibration(C::Courbe)
    C.monotone || return fill(NaN, length(C.c))
    n = length(C.c)
    pente(i) = i == 1 ? (C.o[2] - C.o[1]) / (C.c[2] - C.c[1]) : i == n ? (C.o[n] - C.o[n-1]) / (C.c[n] - C.c[n-1]) : C.d[i]
    return [C.σ_ref[i] / abs(pente(i)) for i in 1:n]
end

"""
    langmuir(c, o, σ) -> NamedTuple

o(c) = o₀ + Δ·c/(K + c), moindres carrés pondérés par 1/σ² (non pondérés si σ manque).
K est profilé (o₀ et Δ linéaires) ; intervalle à 68 % où χ² ≤ χ²_min + χ²ᵣ (le χ²ᵣ
absorbe l'écart au modèle).
"""
function langmuir(c::Vector{Float64}, o::Vector{Float64}, σ::Vector{Float64})
    ok = isfinite.(o)
    c, o, σ = c[ok], o[ok], σ[ok]
    length(c) >= 4 || return (K = NaN, K_bas = NaN, K_haut = NaN, o0 = NaN, delta = NaN, chi2r = NaN)
    w = all(x -> isfinite(x) && x > 0, σ) ? 1 ./ σ .^ 2 : ones(length(c))
    function χ2(K)
        X = hcat(ones(length(c)), c ./ (K .+ c))
        β = (X' * (X .* w)) \ (X' * (o .* w))
        return sum(w .* (o .- X * β) .^ 2), β
    end
    Ks = exp10.(range(-1, 4; length = 501))
    vals = [χ2(K)[1] for K in Ks]
    i = argmin(vals)
    a, b = Ks[max(i - 1, 1)], Ks[min(i + 1, length(Ks))]
    for _ in 1:80                                  # section dorée en log K
        x1, x2 = exp(log(a) + 0.382 * (log(b) - log(a))), exp(log(a) + 0.618 * (log(b) - log(a)))
        χ2(x1)[1] < χ2(x2)[1] ? (b = x2) : (a = x1)
    end
    K = sqrt(a * b)
    v, β = χ2(K)
    χ2r = v / max(length(c) - 3, 1)
    seuil = v + max(χ2r, 1.0)
    dedans = Ks[vals.<=seuil]
    return (K = K, K_bas = isempty(dedans) ? NaN : minimum(dedans), K_haut = isempty(dedans) ? NaN : maximum(dedans),
        o0 = β[1], delta = β[2], chi2r = χ2r)
end

# -----------------------------------------------------------------------------
# Borne de Fisher sur [Cl⁻], au modèle ajusté
# -----------------------------------------------------------------------------

"""Interpolation de Hermite monotone (PCHIP) d'une table et sa dérivée en x."""
function pchip_et_derivee(x::Vector{Float64}, y::Vector{Float64}, d::Vector{Float64}, t::Real)
    C = Courbe("", 0, x, y, d, true, zeros(length(x)))
    return evaluer(C, t)
end

"""
    bornes_fisher(mg, ctx, budgets) -> DataFrame

Borne de Cramér-Rao sur [Cl⁻] (mM) à chaque concentration et chaque budget de photons
bleus, sur le modèle global ajusté. Les composantes de chaque canal (diffusion, deux
exponentielles) et le fond sont ramenés au total bleu de leur cellule, puis interpolés
entre les concentrations (PCHIP, dérivées analytiques) : à budget de photons bleus fixé,
seule la forme relative porte de l'information, et l'intensité absolue (qui change d'un
échantillon à l'autre avec la quantité de capteur et l'excitation) n'en porte pas. Bases
et fonds figés ; une échelle d'intensité libre joue la nuisance. Cinq informations :
`rapport` (J/B seul), `toute` (les deux canaux, une échelle : forme et rapport), `forme`
(une échelle par canal : la forme seule), `bleu`, `jaune`.
"""
function bornes_fisher(mg::ModeleGlobal, ctx::Contexte, budgets)
    cl = mg.cl
    W = FENETRE
    NB = [sum(mg.μB[i][W]) for i in eachindex(cl)]
    # formes normalisées (somme 1) : bleu diffusion, τ₂, τ₃, fond ; jaune diffusion, τ_J, τ_A, fond
    formeB = (ctx.DB, ctx.B2 ./ sum(ctx.B2), ctx.B3 ./ sum(ctx.B3), mg.fondB[1][W] ./ sum(mg.fondB[1][W]))
    formeJ = (ctx.DJ, ctx.JJ ./ sum(ctx.JJ), ctx.JA ./ sum(ctx.JA), mg.fondJ[1][W] ./ sum(mg.fondJ[1][W]))
    # poids de chaque forme par photon bleu de la cellule, et leurs pentes de Hermite
    PB = [vcat([mg.IB[:, j] for j in 1:3], [[sum(mg.fondB[i][W]) for i in eachindex(cl)]])[q] ./ NB for q in 1:4]
    PJ = [vcat([mg.IJ[:, j] for j in 1:3], [[sum(mg.fondJ[i][W]) for i in eachindex(cl)]])[q] ./ NB for q in 1:4]
    # pentes de Hermite ; aux deux bornes, la pente du segment voisin (comme sigma_calibration)
    bords!(d, v) = (d[1] = (v[2] - v[1]) / (cl[2] - cl[1]); d[end] = (v[end] - v[end-1]) / (cl[end] - cl[end-1]); d)
    dPB = [bords!(pentes_pchip(cl, v), v) for v in PB]
    dPJ = [bords!(pentes_pchip(cl, v), v) for v in PJ]
    m = length(W)
    df = DataFrame()
    for c in cl
        μB, μJ, dB, dJ = zeros(m), zeros(m), zeros(m), zeros(m)
        for q in 1:4
            a, da = pchip_et_derivee(cl, PB[q], dPB[q], c)
            μB .+= a .* formeB[q]
            dB .+= da .* formeB[q]
            a, da = pchip_et_derivee(cl, PJ[q], dPJ[q], c)
            μJ .+= a .* formeJ[q]
            dJ .+= da .* formeJ[q]
        end
        for N in budgets
            mB, mJ, gB, gJ = N .* μB, N .* μJ, N .* dB, N .* dJ
            info(cols, μ) = (D = reduce(hcat, cols); D' * (D ./ μ))
            # σ² = 1 / (I_cc − I_cn I_nn⁻¹ I_nc) : complément de Schur des nuisances ; Inf si la
            # dérivée en concentration se confond avec les échelles (aucune information)
            σc(I) = (Ic = I[1, 1] - dot(I[1, 2:end], Symmetric(I[2:end, 2:end]) \ I[2:end, 1]); Ic > 1e-12 * I[1, 1] ? 1 / sqrt(Ic) : Inf)
            z = zeros(m)
            I_toute = info([vcat(gB, gJ), vcat(mB, mJ)], vcat(mB, mJ))
            I_forme = info([vcat(gB, gJ), vcat(mB, z), vcat(z, mJ)], vcat(mB, mJ))
            I_bleu = info([gB, mB], mB)
            I_jaune = info([gJ, mJ], mJ)
            NBt, NJt = sum(mB), sum(mJ)
            R = NJt / NBt
            dR = (sum(gJ) * NBt - NJt * sum(gB)) / NBt^2
            σR = R * sqrt(1 / NBt + 1 / NJt)
            push!(df, (lambda_nm = mg.λ, cl_mM = c, budget_photons = N, rapport = σR / abs(dR), toute = σc(I_toute),
                forme = σc(I_forme), bleu = σc(I_bleu), jaune = σc(I_jaune)))
        end
    end
    return df
end

# -----------------------------------------------------------------------------
# Prédiction d'une concentration inconnue
# -----------------------------------------------------------------------------

"""
Toutes les concentrations où la courbe vaut `o` : une par intervalle de nœuds qui
l'encadre (la courbe de Hermite est monotone dans chaque intervalle, même quand elle ne
l'est pas sur toute la plage). Vide si `o` sort de la plage calibrée.
"""
function inverser_toutes(C::Courbe, o::Real)
    sols = Float64[]
    isfinite(o) || return sols
    for i in 1:length(C.c)-1
        a, b = C.o[i], C.o[i+1]
        (min(a, b) <= o <= max(a, b) && a != b) || continue
        sens = b > a ? 1 : -1
        lo, hi = C.c[i], C.c[i+1]
        for _ in 1:60
            m = (lo + hi) / 2
            sens * (evaluer(C, m)[1] - o) < 0 ? (lo = m) : (hi = m)
        end
        x = (lo + hi) / 2
        (isempty(sols) || abs(x - sols[end]) > 1e-6) && push!(sols, x)
    end
    return sols
end

"""Inversion par interpolation linéaire entre les nœuds (première solution)."""
function inverser_lineaire(C::Courbe, o::Real)
    for i in 1:length(C.c)-1
        a, b = C.o[i], C.o[i+1]
        (min(a, b) <= o <= max(a, b) && a != b) && return C.c[i] + (o - a) / (b - a) * (C.c[i+1] - C.c[i])
    end
    return NaN
end

"""
    predire(C, o, σ_o; n, rng) -> NamedTuple

Concentration de l'échantillon inconnu lue sur la courbe `C`, et ses incertitudes
statistiques par Monte-Carlo (`n` tirages) : l'observable de l'inconnu varie de σ_o
(dispersion de ses 10 acquisitions / √10), les 8 nœuds de σ_ref ; chaque tirage
reconstruit la courbe et l'inverse. `sigma_mesure` et `sigma_calibration` isolent chaque
source, `sigma` les combine ; les intervalles (quantiles) gardent l'asymétrie de la
courbe. `solutions` : toutes les racines de la courbe centrale (plusieurs si elle n'est
pas monotone).
"""
function predire(C::Courbe, o::Real, σ_o::Real; n::Int = 20_000, rng::AbstractRNG = StableRNG(2026_10_16))
    sols = inverser_toutes(C, o)
    centrale = length(sols) == 1 ? sols[1] : NaN
    σn = replace(C.σ_ref, NaN => 0.0)
    σu = isfinite(σ_o) ? σ_o : 0.0
    function tirer(mesure::Bool, noeuds::Bool)
        v = Float64[]
        hors = 0
        for _ in 1:n
            oo = noeuds ? C.o .+ σn .* randn(rng, length(C.o)) : C.o
            ou = mesure ? o + σu * randn(rng) : o
            s = inverser_toutes(Courbe(C.lecture, C.λ, C.c, oo, pentes_pchip(C.c, oo), true, C.σ_ref), ou)
            if isempty(s)
                hors += 1
                continue
            end
            ref = isfinite(centrale) ? centrale : isempty(sols) ? s[1] : sols[1]
            push!(v, s[argmin(abs.(s .- ref))])
        end
        return v, hors / n
    end
    tot, hp = tirer(true, true)
    mes, _ = tirer(true, false)
    cal, _ = tirer(false, true)
    q(v, p) = length(v) < 10 ? NaN : quantile(v, p)
    sd(v) = length(v) < 10 ? NaN : std(v)
    return (c = centrale, solutions = sols, mediane = q(tot, 0.5), bas68 = q(tot, 0.16), haut68 = q(tot, 0.84),
        bas95 = q(tot, 0.025), haut95 = q(tot, 0.975), sigma = sd(tot), sigma_mesure = sd(mes), sigma_calibration = sd(cal),
        hors_plage_mc = hp, c_lineaire = inverser_lineaire(C, o))
end

"""
    hill(c, o, σ) -> NamedTuple

o(c) = o₀ + Δ·cⁿ/(Kⁿ + cⁿ), moindres carrés pondérés par 1/σ² : une forme paramétrique de
la courbe, pour vérifier l'interpolation entre les nœuds. K et n sur une grille, puis une
grille fine autour du meilleur point ; o₀ et Δ linéaires.
"""
function hill(c::Vector{Float64}, o::Vector{Float64}, σ::Vector{Float64})
    w = all(x -> isfinite(x) && x > 0, σ) ? 1 ./ σ .^ 2 : ones(length(c))
    function χ2(K, m)
        X = hcat(ones(length(c)), c .^ m ./ (K^m .+ c .^ m))
        β = (X' * (X .* w)) \ (X' * (o .* w))
        return sum(w .* (o .- X * β) .^ 2), β
    end
    meilleur = (Inf, 1.0, 1.0)
    for K in exp10.(range(0, 3; length = 121)), m in range(0.4, 3.5; length = 63)
        v = χ2(K, m)[1]
        v < meilleur[1] && (meilleur = (v, K, m))
    end
    for K in meilleur[2] .* exp10.(range(-0.03, 0.03; length = 41)), m in meilleur[3] .+ range(-0.06, 0.06; length = 41)
        m > 0 || continue
        v = χ2(K, m)[1]
        v < meilleur[1] && (meilleur = (v, K, m))
    end
    v, K, m = meilleur
    _, β = χ2(K, m)
    return (o0 = β[1], delta = β[2], K = K, n = m, chi2r = v / max(length(c) - 4, 1))
end

"""Inversion de la Hill : c = K·(f/(1 − f))^(1/n), f = (o − o₀)/Δ ; NaN hors de (0, 1)."""
inverser_hill(h, o::Real) = (f = (o - h.o0) / h.delta; 0 < f < 1 ? h.K * (f / (1 - f))^(1 / h.n) : NaN)
