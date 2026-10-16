# Banc.jl — les acquisitions de banc à deux canaux (2026-10-16 et suivantes).
#
# Format : .sdt de Becker & Hickl, un bloc de trois courbes de 256 canaux par fichier
# (voie 1 = bleu, donneur ; voie 2 = jaune, accepteur ; voie 3 vide), comme ch1 = CFP
# et ch2 = YFP dans FLIMApp. TAC de 16,4 ns sur 256 canaux : pas de 64,004 ps, une
# période laser (12,5 ns) = 195,3 canaux ; le TAC coupe les canaux 2 à 13 et le dernier
# canal est partiel, d'où la fenêtre 14–195. Une IRF par longueur d'onde, vue par les
# deux détecteurs à la fois (diffuseur) : chaque voie a donc sa propre IRF, et le
# retard entre détecteurs est dans ces IRF.

const FENETRE_BANC = 14:195

"""Un histogramme à deux canaux du banc (une acquisition)."""
struct HistoBanc
    echantillon::String
    C::Float64
    λ::Int
    cycle::Int
    bleu::Vector{Float64}
    jaune::Vector{Float64}
    duree_s::Float64
end

"""
    lire_sdt_voies(chemin) -> (d, pas_ns, duree_s, arret_sur_temps)

Toutes les voies d'un .sdt : `d[voie, canal]`. Pas = plage du TAC / (gain × canaux).
"""
function lire_sdt_voies(chemin::AbstractString)
    s = SdtFile.read_sdt(read(chemin), basename(chemin))
    mi = s.measure_info[1]
    d = Float64.(s.data[1])
    size(d, 2) == mi.adc_re || (d = permutedims(d))
    pas = mi.tac_r / (mi.tac_g * mi.adc_re) * 1e9
    return d, pas, Float64(mi.col_t), mi.stopt != 0
end

"""
    irf_banc(chemin, voie; fenetre) -> (irf normalisée, pas_ns, fond par canal)

IRF d'une voie : fond retiré (moyenne des canaux valides avant le pic, plutôt que la
médiane, qui inclut ici les canaux vides coupés par le TAC), négatifs mis à zéro, hors
fenêtre à zéro, normalisée à somme 1.
"""
function irf_banc(chemin::AbstractString, voie::Int; fenetre = FENETRE_BANC)
    d, pas, _, _ = lire_sdt_voies(chemin)
    c = d[voie, :]
    k = argmax(c)
    fond = mean(c[first(fenetre):k-4])
    irf = max.(c .- fond, 0.0)
    irf[setdiff(1:length(c), fenetre)] .= 0.0
    return irf ./ sum(irf), pas, fond
end

"""Toutes les acquisitions d'un dossier de banc : échantillon => [Cl⁻] (mM)."""
function lire_banc(dossier::AbstractString; echantillons = ["test_0" => 0.0, "test_131" => 131.0], λs = (750, 800, 850))
    out = HistoBanc[]
    for (ech, C) in echantillons, λ in λs
        rep = joinpath(dossier, ech)
        for f in sort(filter(x -> startswith(x, "$(λ)_") && endswith(x, ".sdt"), readdir(rep)))
            m = match(r"_c(\d+)\.sdt$", f)
            d, _, duree, _ = lire_sdt_voies(joinpath(rep, f))
            push!(out, HistoBanc(ech, C, λ, m === nothing ? 0 : parse(Int, m[1]), d[1, :], d[2, :], duree))
        end
    end
    return out
end

"""Instrument du banc pour une IRF donnée (TAC de 16,4 ns, période de 12,5 ns)."""
instrument_banc(irf::Vector{Float64}, pas::Real; t₀ = 0.0, fenetre = FENETRE_BANC) =
    InstrumentParams(n_bins = length(irf), T_rep = 12.5, pas = Float64(pas), irf = irf, t₀ = Float64(t₀), fenetre = fenetre)

"""
Retard du détecteur jaune sur le bleu tiré des deux IRF : le décalage (sous-canal,
interpolation linéaire) qui superpose au mieux les fronts montants normalisés
(10 à 90 % du pic), insensible à la queue lente qui diffère entre détecteurs.
"""
function retard_irf(irf_B::Vector{Float64}, irf_J::Vector{Float64}, pas::Real)
    front(c) = begin
        k = argmax(c)
        m = c[k]
        i = k
        while i > 1 && c[i-1] < c[i] && c[i-1] > 0.1m
            i -= 1
        end
        # instant où le front franchit 50 % du pic, interpolé
        j = findfirst(t -> c[t] >= 0.5m, i:k)
        t = (i:k)[j]
        t == 1 && return Float64(t)
        t - 1 + (0.5m - c[t-1]) / (c[t] - c[t-1])
    end
    return (front(irf_J) - front(irf_B)) * pas
end

"""Largeur à mi-hauteur (ns) d'une IRF échantillonnée, interpolation linéaire."""
function largeur_mi_hauteur(c::Vector{Float64}, pas::Real)
    k = argmax(c)
    h = c[k] / 2
    i = k
    while i > 1 && c[i] > h
        i -= 1
    end
    j = k
    while j < length(c) && c[j] > h
        j += 1
    end
    a = i + (h - c[i]) / (c[i+1] - c[i])
    b = j - 1 + (c[j-1] - h) / (c[j-1] - c[j])
    return (b - a) * pas
end

# -----------------------------------------------------------------------------
# Lectures des données de banc : rapport, forme conjointe, bi-exponentiel bleu
# -----------------------------------------------------------------------------

"""
Formes affines en x des deux canaux, tirées de l'ajustement conjoint (`ajuster_bleu_jaune`) :
μ_B = A·(u_B + x·v_B) + b_B, μ_J = A·(u_J + x·v_J) + b_J sur la fenêtre. Les fonds sont
ceux de l'ajustement, remis à l'échelle des photons de l'histogramme lu.
"""
struct FormesBanc
    uB::Vector{Float64}
    vB::Vector{Float64}
    uJ::Vector{Float64}
    vJ::Vector{Float64}
    fenetre::UnitRange{Int}
end

"""
    lire_forme(F, nB, nJ, bB, bJ; mode) -> x̂ (NaN si échec)

`mode` : `:commune` (5d : une échelle pour les deux canaux, le rapport J/B informe),
`:separees` (5a : une échelle par canal, la forme seule), `:bleu` (bi-exponentiel à τ
figés sur le bleu seul). Fonds connus (b par canal). Maximum de vraisemblance de Poisson.
"""
function lire_forme(F::FormesBanc, nB, nJ, bB::Real, bJ::Real; mode::Symbol = :commune)
    W = F.fenetre
    b, y = nB[W], nJ[W]
    m = length(W)
    if mode === :bleu
        n = b
        modele!(μ, J, θ) = begin
            x, lA = θ
            A = exp(lA)
            @inbounds for i in 1:m
                g = F.uB[i] + x * F.vB[i]
                μ[i] = A * g + bB
                J[i, 1] = A * F.vB[i]
                J[i, 2] = A * g
            end
        end
        A0 = max(sum(b) - m * bB, 1.0) / sum(F.uB .+ 0.3 .* F.vB)
        r = lm_poisson(modele!, [0.3, log(A0)], [-0.5, -Inf], [1.5, Inf], n; maxiter = 200, tol_dec = 1e-9)
        return r.converge ? r.θ[1] : NaN
    end
    n = vcat(b, y)
    if mode === :commune
        modele_c!(μ, J, θ) = begin
            x, lA = θ
            A = exp(lA)
            @inbounds for i in 1:m
                gB = F.uB[i] + x * F.vB[i]
                gJ = F.uJ[i] + x * F.vJ[i]
                μ[i] = A * gB + bB
                μ[m+i] = A * gJ + bJ
                J[i, 1] = A * F.vB[i]
                J[m+i, 1] = A * F.vJ[i]
                J[i, 2] = A * gB
                J[m+i, 2] = A * gJ
            end
        end
        A0 = max(sum(b) + sum(y) - m * (bB + bJ), 1.0) / sum(F.uB .+ F.uJ .+ 0.3 .* (F.vB .+ F.vJ))
        r = lm_poisson(modele_c!, [0.3, log(A0)], [-0.5, -Inf], [1.5, Inf], n; maxiter = 200, tol_dec = 1e-9)
        return r.converge ? r.θ[1] : NaN
    end
    modele2!(μ, J, θ) = begin
        x, lA, lAJ = θ
        A, AJ = exp(lA), exp(lAJ)
        fill!(J, 0.0)
        @inbounds for i in 1:m
            gB = F.uB[i] + x * F.vB[i]
            gJ = F.uJ[i] + x * F.vJ[i]
            μ[i] = A * gB + bB
            μ[m+i] = AJ * gJ + bJ
            J[i, 1] = A * F.vB[i]
            J[m+i, 1] = AJ * F.vJ[i]
            J[i, 2] = A * gB
            J[m+i, 3] = AJ * gJ
        end
    end
    A0 = max(sum(b) - m * bB, 1.0) / sum(F.uB .+ 0.3 .* F.vB)
    AJ0 = max(sum(y) - m * bJ, 1.0) / sum(F.uJ .+ 0.3 .* F.vJ)
    r = lm_poisson(modele2!, [0.3, log(A0), log(AJ0)], [-0.5, -Inf, -Inf], [1.5, Inf, Inf], n; maxiter = 200, tol_dec = 1e-9)
    return r.converge ? r.θ[1] : NaN
end

"""
    ajuster_mono_banc(inst_u, h; pas, canaux, fenetre, tentatives) -> τ

Ton `vec_to_lifetime` sur un histogramme du banc : la période tient dans `canaux` = 195
canaux de 64 ps (convolution circulaire sur 12,48 ns), et la porte couvre `fenetre`.
L'instance doit avoir reçu l'IRF du banc tronquée aux mêmes 195 canaux.
"""
function ajuster_mono_banc(inst_u, h::Vector{Float64}; pas = 0.06400390634153652, canaux = 195, fenetre = FENETRE_BANC, tentatives = 5)
    x = h[1:canaux]
    bas = 100 * first(fenetre) / canaux                 # porte de ton code : canaux bas…haut − 1
    haut = 100 * min(last(fenetre), canaux) / canaux
    for _ in 1:tentatives
        p, _ = Base.invokelatest(inst_u.vec_to_lifetime, x; guess = [3.0, 0.0, 1e-6], bin_size = pas, histogram_resolution = canaux,
            tac_low_cut = bas, tac_high_cut = haut, first_fit = true, hist_size_threshold = 0)
        isnan(p[1]) || return p[1]
    end
    return NaN
end

"""Discriminabilité entre deux échantillons d'une lecture : |Δmoyenne| / σ moyen."""
discriminabilite(a, b) = (a = filter(isfinite, a); b = filter(isfinite, b);
    abs(mean(b) - mean(a)) / sqrt((var(a) + var(b)) / 2))

"""
    formes_banc(r, inst_B, inst_J) -> FormesBanc

Formes affines en x des deux canaux au point ajusté `r` (`ajuster_bleu_jaune`), avec
l'IRF de chaque détecteur : la base des lectures 5d, 5a et du bi-exponentiel bleu.
"""
function formes_banc(r, inst_B::InstrumentParams, inst_J::InstrumentParams)
    wsB, wsJ = EspaceTravail(inst_B), EspaceTravail(inst_J)
    n = inst_B.n_bins
    W = inst_B.fenetre
    θ = r.theta
    τD0, τDA, t0, ρD, ρA, τA, ε, dt = θ[1:8]
    τr = r.τ_r
    B = [zeros(n) for _ in 1:8]
    base!(B[1], wsB, τD0, t0); base!(B[2], wsB, τDA, t0); base!(B[3], wsB, τA, t0)
    base!(B[4], wsJ, τD0, t0 + dt); base!(B[5], wsJ, τDA, t0 + dt); base!(B[6], wsJ, τA, t0 + dt)
    base!(B[7], wsJ, τr, t0 + dt); base!(B[8], wsB, τr, t0)
    κ = (1 / τr - 1 / τD0) / (1 / τr - 1 / τA)
    c1 = ρA * τD0 / τA
    sensB = κ .* (B[3] .- B[8]) .+ ε .* B[3]
    sensY = κ .* (B[6] .- B[7]) .+ ε .* B[6]
    uB = B[1]
    vB = B[2] .- B[1] .+ r.formes.rBA * c1 .* sensB
    uJ = ρD .* B[4]
    vJ = ρD .* (B[5] .- B[4]) .+ c1 .* sensY
    return FormesBanc(uB[W], vB[W], uJ[W], vJ[W], W)
end

"""Amincissement binomial des deux canaux d'un histogramme de banc, même probabilité p."""
function amincir_banc!(bB, bJ, h::HistoBanc, p::Real, rng::AbstractRNG)
    amincir!(bB, h.bleu, p, rng)
    amincir!(bJ, h.jaune, p, rng)
    return nothing
end

"""
    ajuster_multiexp(histos, inst; τ_init, partages, montee) -> NamedTuple

Somme d'exponentielles libres (convoluées avec l'IRF de `inst`), sans modèle de
capteur : les k durées de vie sont communes aux histogrammes si `partages`, les
amplitudes (positives, ou de signe libre pour `montee`) et le fond sont propres à
chacun ; t₀ commun. Pour caractériser des déclins que M2 ne décrit pas.
"""
function ajuster_multiexp(histos::Vector{<:AbstractVector}, inst::InstrumentParams; τ_init = (0.4, 1.5, 3.5), montee::Bool = false,
                          t0_init = 0.0, maxiter = 500)
    ws = EspaceTravail(inst)
    Wf = collect(inst.fenetre)
    m = length(Wf)
    H = length(histos)
    k = length(τ_init)
    # θ = (τ_1..k, t₀, a_{h,1..k}, b_h)
    ia(h, j) = k + 1 + (h - 1) * (k + 1) + j
    ib(h) = k + 1 + (h - 1) * (k + 1) + k + 1
    P = k + 1 + H * (k + 1)
    Bs = [zeros(ws.n) for _ in 1:k]
    dBτ = [zeros(ws.n) for _ in 1:k]
    dBt = [zeros(ws.n) for _ in 1:k]
    function modele!(μ, J, θ)
        for j in 1:k
            bases_derivees!(Bs[j], dBτ[j], dBt[j], ws, θ[j], θ[k+1])
        end
        fill!(J, 0.0)
        for h in 1:H
            r0 = (h - 1) * m
            @inbounds for (i, c) in enumerate(Wf)
                v = θ[ib(h)]
                dt = 0.0
                for j in 1:k
                    a = θ[ia(h, j)]
                    v += a * Bs[j][c]
                    J[r0+i, j] += a * dBτ[j][c]
                    dt += a * dBt[j][c]
                    J[r0+i, ia(h, j)] = Bs[j][c]
                end
                μ[r0+i] = v
                J[r0+i, k+1] = dt
                J[r0+i, ib(h)] = 1.0
            end
        end
    end
    n = vcat([h[Wf] for h in histos]...)
    θ0 = zeros(P)
    θ0[1:k] .= τ_init
    θ0[k+1] = t0_init
    for (h, hist) in enumerate(histos)
        base!(Bs[1], ws, τ_init[end], t0_init)
        θ0[ib(h)] = max(0.5 * minimum(hist[Wf]), 0.0)
        tot = sum(hist[Wf]) - m * θ0[ib(h)]
        for j in 1:k
            base!(Bs[j], ws, τ_init[j], t0_init)
            θ0[ia(h, j)] = tot / k / sum(Bs[j][Wf])
        end
    end
    lo = vcat(fill(0.01, k), -1.5, repeat(vcat(fill(montee ? -Inf : 0.0, k), 0.0), H))
    hi = vcat(fill(12.0, k), 1.5, repeat(vcat(fill(Inf, k), Inf), H))
    r = lm_poisson(modele!, θ0, lo, hi, n; maxiter, tol_dec = 1e-9)
    θ = r.θ
    libres = findall(lo .< hi)
    cov = fill(NaN, P, P)
    try
        cov[libres, libres] = inv(Symmetric(r.fisher[libres, libres]))
    catch
    end
    σ = sqrt.(max.(diag(cov), 0))
    χ2 = [sum((histos[h][Wf] .- r.μ[(h-1)*m+1:h*m]) .^ 2 ./ r.μ[(h-1)*m+1:h*m]) / (m - (k + 1) - (k + 1) / H) for h in 1:H]
    amps = [[θ[ia(h, j)] for j in 1:k] for h in 1:H]
    # fraction d'intensité de chaque composante (sur la fenêtre)
    fr = [[θ[ia(h, j)] * sum(Bs[j][Wf]) for j in 1:k] ./ sum(θ[ia(h, j)] * sum(Bs[j][Wf]) for j in 1:k) for h in 1:H]
    return (τ = θ[1:k], σ_τ = σ[1:k], t₀ = θ[k+1], amplitudes = amps, fractions_intensite = fr, fonds = [θ[ib(h)] for h in 1:H],
        chi2 = χ2, converge = r.converge, μ = [r.μ[(h-1)*m+1:h*m] for h in 1:H])
end

"""IRF réduite à son cœur : la queue au-delà de `apres` ns après le pic est retirée."""
function irf_coeur(irf::Vector{Float64}, pas::Real; apres = 0.5)
    k = argmax(irf)
    c = copy(irf)
    c[k+ceil(Int, apres / pas)+1:end] .= 0.0
    return c ./ sum(c)
end

"""
    gabarits_banc(rB, rJ) -> FormesBanc

Lectures de forme sans modèle de capteur : chaque état (0 et 131 mM) est décrit par
son déclin mesuré, lissé par l'ajustement multi-exponentiel (`rB` bleu, `rJ` jaune,
fond retiré). Les deux canaux d'un état sont normalisés par le même total bleu, si bien
que le mélange μ = A·[(1 − f)·état₀ + f·état₁] garde le rapport J/B de chaque état
(lecture à échelle commune) ; f = 0 et 1 aux deux états.
"""
function gabarits_banc(rB, rJ, fenetre = FENETRE_BANC)
    T = map(1:2) do s
        b = rB.μ[s] .- rB.fonds[s]
        j = rJ.μ[s] .- rJ.fonds[s]
        nb = sum(b)
        (b ./ nb, j ./ nb)
    end
    return FormesBanc(T[1][1], T[2][1] .- T[1][1], T[1][2], T[2][2] .- T[1][2], fenetre)
end

"""Coordonnée phasor du bleu le long de l'axe état₀ → état₁ (0 et 1 aux deux états)."""
function lecture_phasor(h, P0::Complex, P1::Complex, pas::Real; fenetre = FENETRE_BANC, T = 12.5)
    ω = 2π / T
    z = sum(h[k] * cis(ω * (k - 0.5) * pas) for k in fenetre) / sum(h[k] for k in fenetre)
    d = P1 - P0
    return real((z - P0) * conj(d)) / abs2(d)
end
phasor_brut(h, pas::Real; fenetre = FENETRE_BANC, T = 12.5) =
    sum(h[k] * cis(2π / T * (k - 0.5) * pas) for k in fenetre) / sum(h[k] for k in fenetre)

"""
    lire_biexp_libre(nB, ws, fenetre, τ_init, t₀, b) -> fraction d'amplitude courte (NaN si échec)

Estimateur 3 du simulateur sur le bleu du banc : deux exponentielles à τ libres,
amplitudes libres, fond connu `b`, t₀ figé ; mêmes règles d'échec (τ sur une borne,
composantes confondues).
"""
function lire_biexp_libre(nB, ws::EspaceTravail, fenetre, τ_init::NTuple{2,Float64}, t₀::Real, b::Real; bornes = (0.01, 20.0))
    W = fenetre
    n = nB[W]
    sum(n) > 0 || return NaN
    B1, d1, B2, d2, poub = zeros(ws.n), zeros(ws.n), zeros(ws.n), zeros(ws.n), zeros(ws.n)
    modele!(μ, J, θ) = begin
        bases_derivees!(B1, d1, poub, ws, θ[1], Float64(t₀))
        bases_derivees!(B2, d2, poub, ws, θ[2], Float64(t₀))
        @inbounds for (i, k) in enumerate(W)
            μ[i] = θ[3] * B1[k] + θ[4] * B2[k] + b
            J[i, 1] = θ[3] * d1[k]
            J[i, 2] = θ[4] * d2[k]
            J[i, 3] = B1[k]
            J[i, 4] = B2[k]
        end
    end
    base!(B1, ws, τ_init[1], Float64(t₀))
    base!(B2, ws, τ_init[2], Float64(t₀))
    tot = max(sum(n) - length(W) * b, 1.0)
    θ0 = [τ_init[1], τ_init[2], 0.5tot / sum(B1[W]), 0.5tot / sum(B2[W])]
    r = lm_poisson(modele!, θ0, [bornes[1], bornes[1], 0.0, 0.0], [bornes[2], bornes[2], Inf, Inf], n; maxiter = 200, tol_dec = 1e-8)
    r.converge || return NaN
    any(r.borne[1:2]) && return NaN
    τ1, τ2, A1, A2 = r.θ
    abs(τ1 - τ2) / max(τ1, τ2) < 0.02 && return NaN
    return (τ1 < τ2 ? A1 : A2) / (A1 + A2)
end
