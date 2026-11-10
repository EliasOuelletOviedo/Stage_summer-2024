# Jaune.jl — accueillir le canal jaune mesuré (changement 3).
#
# Mesure prévue : trois concentrations (0, ≈ 30 mM, saturant), en solution, à
# temps d'acquisition fixé, les deux canaux enregistrés ensemble. Le canal bleu
# donne, par l'ajustement global (Ingest.jl), l'amplitude A_c du donneur et la
# fraction x_c à chaque concentration. Le canal jaune s'écrit alors, sans autre
# inconnue que quatre paramètres et un fond par histogramme :
#
#   μ_J,c = A_c·{ ρ_D·[x_c·B(τ_DA) + (1 − x_c)·B(τ_D0)]
#               + ρ_A·(τ_D0/τ_A)·x_c·[κ·(B(τ_A) − B(τ_DA)) + ε·B(τ_A)] } + b_c
#
# avec κ = k_t/(k_DA − k_A), k_t = 1/τ_DA − 1/τ_D0, et B la base convoluée
# (amplitude 1 à t = 0) au t₀ du bleu plus Δt (retard du détecteur jaune).
#   ρ_D = η_YD/η_BD : fuite du donneur dans le jaune, relative au bleu ;
#   ρ_A = η_YA·Φ_A/(η_BD·Φ_D) : émission sensibilisée relative au donneur ;
#   τ_A : durée de vie de l'accepteur ; ε : excitation directe (bornée à [0 ; 0,15]).
# Au Cl⁻ saturant, x ≈ 0,05 : le jaune est presque pure fuite et sa constante de
# temps doit valoir τ_D0, pas τ_A. Les gains relatifs des détecteurs sont absorbés
# dans ρ_D et ρ_A, qui sont donc directement les η effectifs dont le simulateur a
# besoin ; Φ_A effectif en découle avec η_YA des spectres.

"""Résultat de l'ajustement du canal jaune."""
struct ResultatJaune
    ρ_D::Float64
    ρ_A::Float64
    τ_A::Float64
    ε::Float64
    Δt::Float64
    b::Vector{Float64}
    σ::Dict{Symbol,Float64}
    chi2::Float64
    converge::Bool
    μ::Vector{Vector{Float64}}
end

"""
    ajuster_jaune(hY, gB, inst; τ_A_init, Δt_libre, Δt, ε_max) -> ResultatJaune

`hY` : histogrammes jaunes (un par concentration), `gB` : l'ajustement global du
canal bleu des mêmes acquisitions (ResultatGlobal, même ordre). `inst` : IRF du
canal jaune (la bleue par défaut si elle n'est pas mesurée séparément).
"""
function ajuster_jaune(hY::Vector{<:AbstractVector}, gB::ResultatGlobal, inst::InstrumentParams;
                       τ_A_init = 3.0, Δt_libre::Bool = false, Δt::Real = 0.0, ε_max = 0.15, maxiter = 500)
    ws = EspaceTravail(inst)
    W = collect(inst.fenetre)
    m = length(W)
    H = length(hY)
    τD0, τDA = gB.τ_D0, gB.τ_DA
    kt = 1 / τDA - 1 / τD0
    nf = ws.n
    BD0, dBD0, tD0 = zeros(nf), zeros(nf), zeros(nf)
    BDA, dBDA, tDA = zeros(nf), zeros(nf), zeros(nf)
    BA, dBA, tA = zeros(nf), zeros(nf), zeros(nf)
    # θ = (ρ_D, ρ_A, τ_A, ε, Δt, b_1..H)
    P = 5 + H
    function modele!(μ, J, θ)
        ρD, ρA, τA, ε, dt = θ[1], θ[2], θ[3], θ[4], θ[5]
        t = gB.t₀ + dt
        bases_derivees!(BD0, dBD0, tD0, ws, τD0, t)
        bases_derivees!(BDA, dBDA, tDA, ws, τDA, t)
        bases_derivees!(BA, dBA, tA, ws, τA, t)
        kA = 1 / τA
        κ = kt / (1 / τDA - kA)
        dκ_dτA = -kt / (1 / τDA - kA)^2 * kA^2          # ∂κ/∂τ_A = ∂κ/∂k_A · (−k_A²)
        fill!(J, 0.0)
        for h in 1:H
            A, x = gB.A[h], gB.x[h]
            c = ρA * (τD0 / τA) * x
            dc_dτA = -c / τA
            r0 = (h - 1) * m
            @inbounds for (i, k) in enumerate(W)
                fuite = x * BDA[k] + (1 - x) * BD0[k]
                sens = κ * (BA[k] - BDA[k]) + ε * BA[k]
                μ[r0+i] = A * (ρD * fuite + c * sens) + θ[5+h]
                J[r0+i, 1] = A * fuite
                J[r0+i, 2] = A * (τD0 / τA) * x * sens
                J[r0+i, 3] = A * (dc_dτA * sens + c * (dκ_dτA * (BA[k] - BDA[k]) + (κ + ε) * dBA[k]))
                J[r0+i, 4] = A * c * BA[k]
                J[r0+i, 5] = A * (ρD * (x * tDA[k] + (1 - x) * tD0[k]) + c * (κ * (tA[k] - tDA[k]) + ε * tA[k]))
                J[r0+i, 5+h] = 1.0
            end
        end
        return nothing
    end
    n = vcat([y[W] for y in hY]...)
    # départ : fuite depuis l'histogramme le plus saturé (x minimal), sensibilisée depuis le plus actif
    hs, ha = argmin(gB.x), argmax(gB.x)
    base!(BD0, ws, τD0, gB.t₀)
    base!(BDA, ws, τDA, gB.t₀)
    fuite_s = gB.A[hs] .* (gB.x[hs] .* BDA[W] .+ (1 - gB.x[hs]) .* BD0[W])
    ρD0 = 0.9 * sum(hY[hs][W]) / sum(fuite_s)
    base!(BA, ws, τ_A_init, gB.t₀)
    fuite_a = gB.A[ha] .* (gB.x[ha] .* BDA[W] .+ (1 - gB.x[ha]) .* BD0[W])
    ρA0 = max(sum(hY[ha][W]) - ρD0 * sum(fuite_a), 0.1 * sum(hY[ha][W])) / (gB.A[ha] * (τD0 / τ_A_init) * gB.x[ha] * sum(BA[W]))
    θ0 = vcat([ρD0, ρA0, τ_A_init, 0.03, Δt], [max(0.02 * minimum(y[W]), 0.0) for y in hY])
    lo = vcat([0.0, 0.0, 0.5, 0.0, Δt_libre ? -0.5 : Δt], zeros(H))
    hi = vcat([10.0, 50.0, 10.0, ε_max, Δt_libre ? 0.5 : Δt], fill(Inf, H))
    r = lm_poisson(modele!, θ0, lo, hi, n; maxiter, tol_dec = 1e-9)
    libres = findall(lo .< hi)
    cov = fill(NaN, P, P)
    try
        cov[libres, libres] = inv(Symmetric(r.fisher[libres, libres]))
    catch
    end
    σ = sqrt.(max.(diag(cov), 0))
    μs = [r.μ[(h-1)*m+1:h*m] for h in 1:H]
    χ2 = sum(sum((y[W] .- μh) .^ 2 ./ μh) for (y, μh) in zip(hY, μs)) / (H * m - length(libres))
    return ResultatJaune(r.θ[1], r.θ[2], r.θ[3], r.θ[4], r.θ[5], r.θ[6:end],
        Dict(:ρ_D => σ[1], :ρ_A => σ[2], :τ_A => σ[3], :ε => σ[4], :Δt => σ[5]), χ2, r.converge, μs)
end

"""
    ajuster_bleu_jaune(hB, hY, inst, m; inst_jaune, τ_fixes, τ_init, Δt_libre, montee_libre, ε_max) -> NamedTuple

Ajustement CONJOINT des deux canaux : x_c, A_c, fonds bleu et jaune par histogramme,
(ρ_D, ρ_A, τ_A, ε, Δt) communs, et τ_D0, τ_DA, t₀ (figés à `τ_init` si `τ_fixes`).
C'est la procédure recommandée : en deux temps (`ajuster_jaune` après
`ajuster_global`), les erreurs du bleu sur x et A — fond dégénéré avec le plancher du
déclin long — passent dans ρ_D sans être comptées (dispersion 12 fois la borne dans
l'aller-retour). Le bleu inclut la fuite de l'accepteur (rapport η_BA/η_YA des
spectres). `inst_jaune` : IRF du détecteur jaune (par défaut celle du bleu) ; Δt est
alors le retard qui reste au-delà de celui que portent les deux IRF. `montee_libre` :
la montée de l'émission sensibilisée a sa propre constante τ_r au lieu de τ_DA (M2
lent prédit τ_r = τ_DA). Jacobien par différences centrées.
"""
function ajuster_bleu_jaune(hB::Vector{<:AbstractVector}, hY::Vector{<:AbstractVector}, inst::InstrumentParams, m::Modele;
                            inst_jaune::InstrumentParams = inst, τ_fixes::Bool = true, τ_init = (m.p.τ_D0, m.p.τ_DA),
                            Δt_libre::Bool = true, montee_libre::Bool = false, ε_max = 0.15, maxiter = 300)
    H = length(hB)
    W = collect(inst.fenetre)
    nW = length(W)
    ws = EspaceTravail(inst)
    wsJ = EspaceTravail(inst_jaune)
    histos = [Histogramme("", "", NaN, 0, Vector{Float64}(h)) for h in hB]
    kw = τ_fixes ? (tau_D0_fixe = τ_init[1], tau_DA_fixe = τ_init[2]) : NamedTuple()
    gB = ajuster_global(histos, inst; τ_D0_init = τ_init[1], τ_DA_init = τ_init[2], kw...)
    rj = ajuster_jaune(hY, gB, inst_jaune; Δt_libre, ε_max)
    rBA = m.o.η[1, 2] / m.o.η[2, 2]
    # θ = (τ_D0, τ_DA, t₀, ρ_D, ρ_A, τ_A, ε, Δt, τ_r, x_1..H, A_1..H, bB_1..H, bY_1..H)
    ix(h) = 9 + h
    iA(h) = 9 + H + h
    ibB(h) = 9 + 2H + h
    ibY(h) = 9 + 3H + h
    P = 9 + 4H
    B = [zeros(ws.n) for _ in 1:8]
    function mu!(μ, θ)
        τD0, τDA, t0, ρD, ρA, τA, ε, dt, τr0 = θ[1:9]
        τr = montee_libre ? τr0 : τDA
        base!(B[1], ws, τD0, t0); base!(B[2], ws, τDA, t0); base!(B[3], ws, τA, t0)
        base!(B[4], wsJ, τD0, t0 + dt); base!(B[5], wsJ, τDA, t0 + dt); base!(B[6], wsJ, τA, t0 + dt)
        if montee_libre
            base!(B[7], wsJ, τr, t0 + dt); base!(B[8], ws, τr, t0)
        end
        Bmj, Bmb = montee_libre ? (B[7], B[8]) : (B[5], B[2])
        κ = (1 / τr - 1 / τD0) / (1 / τr - 1 / τA)
        for h in 1:H
            x, A = θ[ix(h)], θ[iA(h)]
            c = ρA * (τD0 / τA) * x
            rB, rY = (h - 1) * nW, (H + h - 1) * nW
            @inbounds for (i, k) in enumerate(W)
                sensB = κ * (B[3][k] - Bmb[k]) + ε * B[3][k]
                μ[rB+i] = A * (x * B[2][k] + (1 - x) * B[1][k] + rBA * c * sensB) + θ[ibB(h)]
                sensY = κ * (B[6][k] - Bmj[k]) + ε * B[6][k]
                μ[rY+i] = A * (ρD * (x * B[5][k] + (1 - x) * B[4][k]) + c * sensY) + θ[ibY(h)]
            end
        end
        return μ
    end
    θ0 = vcat([gB.τ_D0, gB.τ_DA, gB.t₀, rj.ρ_D, rj.ρ_A, rj.τ_A, rj.ε, rj.Δt, gB.τ_DA], gB.x, gB.A, gB.b, rj.b)
    lo = vcat([τ_fixes ? τ_init[1] : 0.2, τ_fixes ? τ_init[2] : 0.01, -1.5, 0.0, 0.0, 0.5, 0.0, Δt_libre ? -0.5 : rj.Δt,
            montee_libre ? 0.01 : gB.τ_DA], fill(-0.5, H), zeros(H), zeros(H), zeros(H))
    hi = vcat([τ_fixes ? τ_init[1] : 12.0, τ_fixes ? τ_init[2] : 6.0, 1.5, 50.0, 200.0, 10.0, ε_max, Δt_libre ? 0.5 : rj.Δt,
            montee_libre ? 3.0 : gB.τ_DA], fill(1.5, H), fill(Inf, H), fill(Inf, H), fill(Inf, H))
    libres = findall(lo .< hi)
    μp, μm = zeros(2H * nW), zeros(2H * nW)
    function modele!(μ, J, θ)
        mu!(μ, θ)
        fill!(J, 0.0)
        for j in libres
            δ = 1e-6 * max(abs(θ[j]), 1e-2)
            tp, tm = copy(θ), copy(θ)
            tp[j] += δ
            tm[j] -= δ
            mu!(μp, tp)
            mu!(μm, tm)
            @views J[:, j] .= (μp .- μm) ./ (2δ)
        end
        return nothing
    end
    n = vcat([h[W] for h in hB]..., [h[W] for h in hY]...)
    r = lm_poisson(modele!, clamp.(θ0, lo, hi), lo, hi, n; maxiter, tol_dec = 1e-9)
    cov = fill(NaN, P, P)
    try
        cov[libres, libres] = inv(Symmetric(r.fisher[libres, libres]))
    catch
    end
    σ = sqrt.(max.(diag(cov), 0))
    χ2(rg) = sum((n[rg] .- r.μ[rg]) .^ 2 ./ r.μ[rg]) / length(rg)
    θ = r.θ
    μs = r.μ
    return (ρ_D = θ[4], ρ_A = θ[5], τ_A = θ[6], ε = θ[7], Δt = θ[8], τ_D0 = θ[1], τ_DA = θ[2], t₀ = θ[3],
        τ_r = montee_libre ? θ[9] : θ[2],
        x = θ[ix(1):ix(H)], A = θ[iA(1):iA(H)], b_bleu = θ[ibB(1):ibB(H)], b_jaune = θ[ibY(1):ibY(H)],
        σ = Dict(:ρ_D => σ[4], :ρ_A => σ[5], :τ_A => σ[6], :ε => σ[7], :Δt => σ[8], :τ_D0 => σ[1], :τ_DA => σ[2],
            :τ_r => montee_libre ? σ[9] : σ[2], :x => σ[ix(1):ix(H)]),
        chi2_bleu = χ2(1:H*nW), chi2_jaune = χ2(H*nW+1:2H*nW), converge = r.converge, cov = cov,
        μ_bleu = [μs[(h-1)*nW+1:h*nW] for h in 1:H], μ_jaune = [μs[(H+h-1)*nW+1:(H+h)*nW] for h in 1:H],
        formes = (B = deepcopy(B), κ = (1 / (montee_libre ? θ[9] : θ[2]) - 1 / θ[1]) / (1 / (montee_libre ? θ[9] : θ[2]) - 1 / θ[6]), rBA = rBA),
        theta = θ)
end

"""Valeurs que le canal jaune devrait rendre pour un monde simulé (vérité du test d'aller-retour)."""
verite_jaune(m::Modele) = (ρ_D = m.o.η[2, 1] / m.o.η[1, 1], ρ_A = m.o.η[2, 2] * phi_A(m.p) / (m.o.η[1, 1] * m.p.Φ_D),
    τ_A = m.p.τ_A, ε = m.o.ε, Δt = m.inst.Δt_jaune)

"""
    parametres_mesures(rj, m; date) -> String

Le bloc à recopier dans parametres.toml après la mesure : η_YD, τ_A, Φ_A effectif
et ε au statut « mesuré », avec date et incertitude (1 σ de l'ajustement).
"""
function parametres_mesures(rj, m::Modele; date::AbstractString, source::AbstractString = "canal jaune, step12")
    ηBD, ηYA, ΦD = m.o.η[1, 1], m.o.η[2, 2], m.p.Φ_D
    ΦAe = rj.ρ_A * ηBD * ΦD / ηYA
    σΦ = rj.σ[:ρ_A] * ηBD * ΦD / ηYA
    lignes = [
        "[capteur.tau_A]", @sprintf("valeur = %.4f", rj.τ_A), "statut = \"mesuré\"", "date = \"$date\"",
        @sprintf("incertitude = %.4f", rj.σ[:τ_A]), "source = \"$source\"", "",
        "[capteur.Phi_A_eff]", @sprintf("valeur = %.4f", ΦAe), "statut = \"mesuré\"", "date = \"$date\"",
        @sprintf("incertitude = %.4f", σΦ), "source = \"$source ; ρ_A × η_BD × Φ_D / η_YA (η_YA des spectres)\"", "",
        "# [optique]", @sprintf("# eta_YD = { valeur = %.5f, statut = \"mesuré\", date = \"%s\", incertitude = %.5f, source = \"%s ; ρ_D × η_BD\" }",
            rj.ρ_D * ηBD, date, rj.σ[:ρ_D] * ηBD, source),
        @sprintf("# epsilon = { valeur = %.4f, statut = \"mesuré\", date = \"%s\", incertitude = %.4f, source = \"%s (borné à [0 ; 0,15])\" }",
            rj.ε, date, rj.σ[:ε], source),
    ]
    return join(lignes, "\n")
end

# -----------------------------------------------------------------------------
# Prédictions que la mesure confrontera
# -----------------------------------------------------------------------------

"""Rapport J/B (totaux de la fenêtre) du modèle à la concentration C."""
function ratio_JB(C::Real, m::Modele, ws::EspaceTravail)
    μB, μY = formes(C, m, ws)
    return sum(μY) / sum(μB)
end

"""
    fisher_numerique(f, θ; h) -> Matrix

Information de Fisher de Poisson d'un modèle μ = f(θ) (vecteur de comptages),
par différences centrées : I = Σ (∂μ/∂θᵢ)(∂μ/∂θⱼ)/μ.
"""
function fisher_numerique(f, θ::AbstractVector; h = 1e-5)
    μ = f(θ)
    D = zeros(length(μ), length(θ))
    for i in eachindex(θ)
        δ = h * max(abs(θ[i]), 1e-2)
        tp, tm = copy(θ), copy(θ)
        tp[i] += δ
        tm[i] -= δ
        D[:, i] = (f(tp) .- f(tm)) ./ (2δ)
    end
    ok = μ .> 0
    return D[ok, :]' * (D[ok, :] ./ μ[ok])
end

"""
    sigma_montee(τ_DA, m, N_J, ws; C, ε_libre, Δt_libre) -> σ(τ_DA)

Borne de Cramér-Rao sur la constante de montée de l'émission sensibilisée, tirée
du seul canal jaune (N_J photons à la concentration C) : amplitude de la fuite du
donneur, amplitude de la sensibilisée, τ_A et, si `ε_libre`, l'amplitude de
l'excitation directe sont libres ; t₀ est connu par le bleu. `Δt_libre` : le retard
du détecteur jaune est inconnu (pas de mesure au diffuseur) ; une montée et un
retard décalent tous deux le front montant. Sous M2, la montée a la constante
τ_DA (la même que le déclin court du donneur qui fuit dans le jaune).
"""
function sigma_montee(τ_DA::Real, m::Modele, N_J::Real, ws::EspaceTravail; C::Real = 0.0, ε_libre::Bool = true, Δt_libre::Bool = false)
    W = m.inst.fenetre
    x = fraction_active(C, m.p)
    τD0 = m.p.τ_D0
    function μ(θ)
        τda, aD, aS, τA, aE, dt = θ
        s = aS * (1 / τda - 1 / τD0) / (1 / τda - 1 / τA)
        es = ExpSum([aD * x, aD * (1 - x), -s, s + aE], [τda, τD0, τda, τA])
        return forward(es, remplacer(m.inst; t₀ = m.inst.t₀ + dt), ws)[W]
    end
    kA = 1 / m.p.τ_A
    aD = m.o.η[2, 1] * m.p.Φ_D / τD0
    aS = m.o.η[2, 2] * phi_A(m.p) * kA * x
    aE = m.o.η[2, 2] * m.o.ε * x * phi_A(m.p) * kA
    θ = [Float64(τ_DA), aD, aS, m.p.τ_A, aE, 0.0]
    S = N_J / sum(μ(θ))
    libres = [1, 2, 3, 4]
    ε_libre && push!(libres, 5)
    Δt_libre && push!(libres, 6)
    function μl(t)
        θp = copy(θ)
        θp[libres] .= t
        return S .* μ(θp)
    end
    I = fisher_numerique(μl, θ[libres])
    return sqrt(inv(Symmetric(I))[1, 1])
end

# -----------------------------------------------------------------------------
# η_YD tiré des images déjà acquises (sortie de l'ingestion, changement 3)
# -----------------------------------------------------------------------------

"""B/J : l'observable des images TIFF (C2/C1 : étiquettes inversées)."""
struct RatioBJ <: Estimateur end
nom(::RatioBJ) = "B/J"
observable(::RatioBJ, nB, nY, ctx) = (sum(nB) / sum(nY), :ok)

"""Monde dont la fuite η_YD est multipliée par f (le reste inchangé)."""
avec_fuite(m::Modele, f::Real) = remplacer(m; o = remplacer(m.o; η = SMatrix{2,2}([m.o.η[1, 1] m.o.η[1, 2]; f*m.o.η[2, 1] m.o.η[2, 2]])))

"""
    courbe_fuite(m, ctx; facteurs) -> DataFrame

Pour chaque facteur f sur η_YD : K apparent de B/J (Hill ajustée à la courbe sans
bruit, comme pour la figure 5) et plage dynamique (plateau haut / plateau bas de
cette Hill). Ce sont les deux nombres que l'analyse TIFF a mesurés.
"""
function courbe_fuite(m::Modele, ctx::Contexte; facteurs = 0.1:0.05:2.0)
    out = DataFrame()
    for f in facteurs
        r = K_apparent(RatioBJ(), avec_fuite(m, f), ctx)
        push!(out, (facteur = f, eta_YD = f * m.o.η[2, 1], K_app = r.K_app, n_app = r.n_app, plage = r.haut / r.bas, bas = r.bas, haut = r.haut))
    end
    return out
end

"""
    inverser_fuite(courbe, colonne, cible) -> facteur (NaN si hors de la courbe)

Interpolation linéaire du facteur f qui redonne la valeur mesurée `cible`
(colonne `:K_app` ou `:plage`), sur une courbe supposée monotone par morceaux ;
renvoie la première solution.
"""
function inverser_fuite(courbe::DataFrame, colonne::Symbol, cible::Real)
    f, v = courbe.facteur, courbe[!, colonne]
    for i in 1:length(f)-1
        if (v[i] - cible) * (v[i+1] - cible) <= 0 && v[i] != v[i+1]
            return f[i] + (cible - v[i]) / (v[i+1] - v[i]) * (f[i+1] - f[i])
        end
    end
    return NaN
end
