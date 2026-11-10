# Charger.jl — étape 1 : les cellules de la titration et les deux contrôles.
#
# Une cellule = une longueur d'onde × une concentration : 10 acquisitions de 1 s, deux
# canaux. Les 10 acquisitions restent séparées en mémoire ; la somme se calcule à la
# volée. L'inverse ne se fait pas, et c'est la série des 10 qui porte les contrôles et
# l'analyse du blanchiment.

const LONGUEURS_ONDE = (750, 800, 850)
const FENETRE = 14:195            # une période laser moins les canaux coupés par le TAC
const BLEU, JAUNE = 1, 2

"""
Une cellule de la grille. `acq[canal, canal TCSPC, acquisition]` (2 × 256 × 10) ;
`irf[canal, ·]` normalisée à somme 1 (une par détecteur) ; `fond[canal, ·]` : photons
attendus par canal TCSPC pour UNE acquisition, sans échantillon ; `pas_ns` : pas du TAC.
"""
struct Cellule
    λ::Int
    cl::Float64
    acq::Array{Int,3}
    irf::Matrix{Float64}
    fond::Matrix{Float64}
    duree_s::Float64
    pas_ns::Float64
end

n_acq(c::Cellule) = size(c.acq, 3)
"""Somme des acquisitions d'un canal (toutes, ou celles de `quelles`)."""
somme(c::Cellule, canal::Int; quelles = 1:n_acq(c)) = Float64.(vec(sum(view(c.acq, canal, :, quelles); dims = 2)))
acquisition(c::Cellule, canal::Int, i::Int) = Float64.(c.acq[canal, :, i])
"""Fond attendu d'un canal pour `t_s` secondes d'acquisition."""
fond(c::Cellule, canal::Int, t_s::Real) = c.fond[canal, :] .* (t_s / c.duree_s)
photons(h::AbstractVector) = sum(@view h[FENETRE])

"""
Sous-dossiers de concentration d'une titration : nom = [Cl⁻] en mM. Les noms qui ne
sont pas des nombres (`unknown`, `tampon`) sont écartés.
"""
function dossiers_concentration(dossier::AbstractString)
    out = Pair{Float64,String}[]
    for d in readdir(dossier)
        isdir(joinpath(dossier, d)) || continue
        c = tryparse(Float64, d)
        c === nothing || push!(out, c => joinpath(dossier, d))
    end
    return sort(out; by = first)
end

fichiers_acquisition(rep, λ) = sort(filter(f -> startswith(f, "$(λ)_") && endswith(f, ".sdt"), readdir(rep)))

"""
    fond_detecteurs(dossier, λ, irf_chemin, duree_s) -> (fond[2, 256], source)

Le fond « tampon seul » : si un dossier `tampon/` existe à côté des concentrations,
la moyenne de ses acquisitions à cette longueur d'onde (forme temporelle comprise).
Sinon, le niveau plat mesuré avant le pic de l'IRF de chaque détecteur (comptes
d'obscurité et lumière parasite), ramené à la durée d'une acquisition.
"""
function fond_detecteurs(dossier::AbstractString, λ::Int, irf_chemin::AbstractString, duree_s::Real)
    f = zeros(2, 256)
    rep = joinpath(dossier, "tampon")
    if isdir(rep) && !isempty(fichiers_acquisition(rep, λ))
        fs = fichiers_acquisition(rep, λ)
        for fi in fs
            d, _, dur, _ = S.lire_sdt_voies(joinpath(rep, fi))
            f .+= d[1:2, :] .* (duree_s / dur)
        end
        f ./= length(fs)
        f[:, setdiff(1:256, FENETRE)] .= 0.0
        return f, "tampon seul ($(length(fs)) acquisitions)"
    end
    _, _, dur_irf, _ = S.lire_sdt_voies(irf_chemin)
    for canal in (BLEU, JAUNE)
        _, _, niveau = S.irf_banc(irf_chemin, canal)
        f[canal, FENETRE] .= niveau * duree_s / dur_irf
    end
    return f, "niveau avant le pic des IRF (pas de tampon seul)"
end

"""IRF des deux détecteurs à une longueur d'onde (normalisées) et pas du TAC."""
function irf_detecteurs(dossier_irf::AbstractString, λ::Int)
    chemin = joinpath(dossier_irf, "IRF_$(λ)_16X.sdt")
    irf = zeros(2, 256)
    pas = NaN
    for canal in (BLEU, JAUNE)
        v, pas, _ = S.irf_banc(chemin, canal)
        irf[canal, :] .= v
    end
    return irf, pas, chemin
end

"""Une cellule : les acquisitions `<λ>_cNN.sdt` d'un dossier, gardées séparées (`cl` = [Cl⁻], NaN si inconnue)."""
function lire_cellule(rep::AbstractString, λ::Int, cl::Real, dossier::AbstractString, dossier_irf::AbstractString)
    irf, pas, irf_chemin = irf_detecteurs(dossier_irf, λ)
    fs = fichiers_acquisition(rep, λ)
    isempty(fs) && return nothing, ""
    acq = zeros(Int, 2, 256, length(fs))
    durees = Float64[]
    for (i, f) in enumerate(fs)
        d, p, dur, arret = S.lire_sdt_voies(joinpath(rep, f))
        arret || @warn "$(joinpath(rep, f)) : acquisition sans arrêt sur temps"
        abs(p - pas) < 1e-9 || error("$(f) : pas du TAC $(p) ns, IRF à $(pas) ns")
        acq[:, :, i] .= round.(Int, d[1:2, :])
        push!(durees, dur)
    end
    length(unique(durees)) == 1 || @warn "$(rep) $(λ) nm : durées d'acquisition différentes $(unique(durees))"
    f, src = fond_detecteurs(dossier, λ, irf_chemin, durees[1])
    return Cellule(λ, Float64(cl), acq, irf, f, durees[1], pas), src
end

"""
    charger(dossier; λs, concentrations, dossier_irf) -> (cellules, sources_du_fond)

`dossier` : la titration (`Data/2026-10-16/SCLM` ou `KIRIN`), un sous-dossier par concentration,
des fichiers `<λ>_cNN.sdt` à deux voies utiles (1 bleu, 2 jaune). Les IRF
`IRF_<λ>_16X.sdt` sont dans `dossier_irf` (le dossier parent par défaut).
"""
function charger(dossier::AbstractString; λs = LONGUEURS_ONDE, concentrations = nothing,
                 dossier_irf::AbstractString = dirname(normpath(dossier)))
    cellules = Cellule[]
    sources = Dict{Int,String}()
    for λ in λs, (cl, rep) in dossiers_concentration(dossier)
        concentrations === nothing || cl in concentrations || continue
        c, src = lire_cellule(rep, λ, cl, dossier, dossier_irf)
        c === nothing && continue
        sources[λ] = src
        push!(cellules, c)
    end
    return cellules, sources
end

"""L'échantillon inconnu (`<dossier>/unknown` par défaut) : une cellule par longueur d'onde, [Cl⁻] = NaN."""
function charger_inconnue(dossier::AbstractString; nom::AbstractString = "unknown", λs = LONGUEURS_ONDE,
                          dossier_irf::AbstractString = dirname(normpath(dossier)))
    out = Cellule[]
    for λ in λs
        c, _ = lire_cellule(joinpath(dossier, nom), λ, NaN, dossier, dossier_irf)
        c === nothing || push!(out, c)
    end
    return out
end

"""Pente des moindres carrés de y contre 1…n, son erreur-type et t = pente/erreur."""
function tendance(y::AbstractVector)
    n = length(y)
    x = collect(1.0:n)
    xm, ym = mean(x), mean(y)
    sxx = sum(abs2, x .- xm)
    b = sum((x .- xm) .* (y .- ym)) / sxx
    res = y .- (ym .+ b .* (x .- xm))
    se = sqrt(sum(abs2, res) / (n - 2) / sxx)
    return (pente = b, se = se, t = se > 0 ? b / se : 0.0)
end

"""
    controler(c) -> NamedTuple

Les deux contrôles avant de faire confiance à une cellule :
- la somme est-elle légitime ? tendance linéaire sur les 10 acquisitions de
  l'intensité bleue, de l'intensité jaune et de J/B ; drapeau si |pente| > 2σ ;
- le bruit est-il poissonien ? dispersion de J/B d'une acquisition à l'autre contre
  la prédiction de Poisson ; drapeau au-delà de 2 × Poisson.
"""
function controler(c::Cellule)
    n = n_acq(c)
    NB = [photons(acquisition(c, BLEU, i)) for i in 1:n]
    NJ = [photons(acquisition(c, JAUNE, i)) for i in 1:n]
    JB = NJ ./ NB
    σP_JB = mean(JB) * sqrt(mean(1 ./ NB .+ 1 ./ NJ))
    tB, tJ, tR = tendance(NB), tendance(NJ), tendance(JB)
    derives = [nom for (nom, t) in (("bleu", tB), ("jaune", tJ), ("J/B", tR)) if abs(t.t) > 2]
    return (lambda_nm = c.λ, cl_mM = c.cl, n_acq = n, duree_s = c.duree_s,
        N_bleu = mean(NB), N_jaune = mean(NJ), JB = mean(JB),
        sigma_N_bleu_sur_poisson = std(NB) / sqrt(mean(NB)), sigma_N_jaune_sur_poisson = std(NJ) / sqrt(mean(NJ)),
        sigma_JB_sur_poisson = std(JB) / σP_JB,
        pente_bleu_pct = 100 * tB.pente / mean(NB), t_bleu = tB.t,
        pente_jaune_pct = 100 * tJ.pente / mean(NJ), t_jaune = tJ.t,
        pente_JB_pct = 100 * tR.pente / mean(JB), t_JB = tR.t,
        drapeau_poisson = std(JB) / σP_JB > 2, drapeau_tendance = !isempty(derives), tendances = join(derives, " "))
end

"""Une ligne par acquisition (photons et J/B), pour les figures de tendance."""
function table_acquisitions(cellules::Vector{Cellule})
    df = DataFrame()
    for c in cellules, i in 1:n_acq(c)
        nb, nj = photons(acquisition(c, BLEU, i)), photons(acquisition(c, JAUNE, i))
        push!(df, (lambda_nm = c.λ, cl_mM = c.cl, acquisition = i, N_bleu = nb, N_jaune = nj, JB = nj / nb))
    end
    return df
end
