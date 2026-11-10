# Amincir.jl — étape 5 : amincissement des sommes et Monte-Carlo.
#
# Pour chaque cellule, la somme des 10 acquisitions est amincie par tirage binomial vers
# sept budgets de photons bleus. Le canal jaune est aminci avec la même probabilité p
# (pas vers un compte cible : J/B est l'objet de l'étude). Le jumeau tire ses photons
# (Poisson) du modèle lissé de la même cellule, à la même échelle.
#
# Les tirages amincis ne sont pas indépendants de leur source : la variance mesurée entre
# tirages vaut ≈ (1 − p) fois la variance vraie, d'où σ_vraie = σ_tirages/√(1 − p),
# appliquée systématiquement (Agreger.jl). Le budget le plus haut reste sous 20 % des
# photons de la source.

const BUDGETS = (1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5)
const P_MAX = 0.2

"""400 tirages par point sous 10³ photons, 200 au-dessus (σ connu à 3,5 % et 5 %)."""
n_tirages_defaut(N::Real) = N < 1e3 ? 400 : 200

p_amincissement(cel::Cellule, N::Real) = N / photons(somme(cel, BLEU))
"""Correction de l'amincissement : σ_vraie = σ_tirages / √(1 − p)."""
sigma_vraie(σ::Real, p::Real) = σ / sqrt(1 - p)

function splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    z = x
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

const CODE_SOURCE = Dict(:reel => 1, :jumeau => 2)
const CODE_SCENARIO = Dict(:aucun => 0, :A => 1, :B => 2, :C => 3)

"""
Graine d'un tirage, dérivée de (λ, concentration, budget, indice du tirage, source,
scénario, taux) et non d'un compteur global : deux exécutions donnent exactement les
mêmes tirages, et un tirage se rejoue seul (`rejouer`).
"""
function graine_tirage(λ::Integer, cl::Real, budget::Real, j::Integer; source::Symbol = :reel, scenario::Symbol = :aucun,
                       taux::Real = 0.0)
    h = 0x2c0f_fee1_5eed_0016
    for v in (λ, round(Int, 100cl), round(Int, budget), j, CODE_SOURCE[source], CODE_SCENARIO[scenario], round(Int, 1e4 * taux))
        h = splitmix64(h ⊻ reinterpret(UInt64, Int64(v)))
    end
    return h
end

"""Tirage de base : la somme réelle amincie (binomial), ou le jumeau (Poisson du lissage)."""
function tirage!(bB, bJ, hB, hJ, μB, μJ, p::Real, rng::AbstractRNG, source::Symbol)
    if source === :reel
        S.amincir!(bB, hB, p, rng)
        S.amincir!(bJ, hJ, p, rng)
    else
        @inbounds for k in eachindex(bB)
            bB[k] = μB[k] > 0 ? Float64(rand(rng, Poisson(p * μB[k]))) : 0.0
            bJ[k] = μJ[k] > 0 ? Float64(rand(rng, Poisson(p * μJ[k]))) : 0.0
        end
    end
    return nothing
end

"""Index d'une cellule dans son modèle global."""
indice(mg::ModeleGlobal, cl::Real) = findfirst(==(cl), mg.cl)

"""
    monte_carlo(cellules, modeles, contextes; budgets, sources, n_tirages, scenarios, lot) -> DataFrame

Une ligne par tirage : budget, p, photons tirés, observables brutes des neuf lectures
(`o_<lecture>`) et codes d'échec (`e_<lecture>`). `scenarios` : couples (scénario, taux) ;
(:aucun, 0.0) pour les tirages intacts. Un tirage blanchi part du tirage intact de même
graine, puis le modifie (Blanchir.jl) : la différence de lecture ne porte que le
blanchiment (nombres aléatoires communs). Parallélisé sur les tirages.
"""
function monte_carlo(cellules::Vector{Cellule}, modeles::Dict{Int,ModeleGlobal}, contextes::Dict{Int,Contexte};
                     budgets = BUDGETS, sources = (:reel, :jumeau), n_tirages = n_tirages_defaut,
                     scenarios = [(:aucun, 0.0)], lot::Int = 25, titre::AbstractString = "Monte-Carlo")
    taches = Tuple{Cellule,Float64,Symbol,Symbol,Float64,UnitRange{Int}}[]
    for cel in cellules, N in budgets
        p_amincissement(cel, N) <= P_MAX || continue
        nt = n_tirages(N)
        for src in sources, (sc, β) in scenarios
            (sc !== :aucun && src !== :reel) && continue
            for a in 1:lot:nt
                push!(taches, (cel, Float64(N), src, sc, β, a:min(a + lot - 1, nt)))
            end
        end
    end
    # modèles blanchis : facteurs d'injection par cellule, scénario et taux
    injections = Dict{Tuple{Int,Float64,Symbol,Float64},Any}()
    for cel in cellules, (sc, β) in scenarios
        sc === :aucun && continue
        mg = modeles[cel.λ]
        injections[(cel.λ, cel.cl, sc, β)] = facteurs_injection(mg, contextes[cel.λ], indice(mg, cel.cl), sc, β)
    end
    sorties = Vector{Vector{NamedTuple}}(undef, length(taches))
    suivant = Threads.Atomic{Int}(1)
    faits = Threads.Atomic{Int}(0)
    nf = length(first(values(contextes)).ws)
    t_debut = time()
    @sync for t in 1:nf
        Threads.@spawn begin
            bB, bJ = zeros(256), zeros(256)
            while (i = Threads.atomic_add!(suivant, 1)) <= length(taches)
                cel, N, src, sc, β, js = taches[i]
                mg, ctx = modeles[cel.λ], contextes[cel.λ]
                ic = indice(mg, cel.cl)
                hB, hJ = somme(cel, BLEU), somme(cel, JAUNE)
                p = p_amincissement(cel, N)
                t_s = p * mg.duree_somme_s[ic]
                inj = sc === :aucun ? nothing : injections[(cel.λ, cel.cl, sc, β)]
                out = NamedTuple[]
                for j in js
                    rng = StableRNG(graine_tirage(cel.λ, cel.cl, N, j; source = src))
                    tirage!(bB, bJ, hB, hJ, mg.μB[ic], mg.μJ[ic], p, rng, src)
                    if inj !== nothing
                        rng2 = StableRNG(graine_tirage(cel.λ, cel.cl, N, j; source = src, scenario = sc, taux = β))
                        injecter!(bB, bJ, inj, p, rng2)
                    end
                    o, e = lire(ctx, bB, bJ, t_s, t)
                    push!(out, (source = String(src), lambda_nm = cel.λ, cl_mM = cel.cl, budget_photons = N, p = p,
                        scenario = String(sc), taux = β, tirage = j, N_bleu = sum(@view bB[FENETRE]), N_jaune = sum(@view bJ[FENETRE]),
                        (Symbol("o_" * k) => o[q] for (q, k) in enumerate(LECTURES))...,
                        (Symbol("e_" * k) => e[q] for (q, k) in enumerate(LECTURES))...))
                end
                sorties[i] = out
                f = Threads.atomic_add!(faits, 1) + 1
                if f % max(1, length(taches) ÷ 10) == 0
                    @printf("  %s : %3.0f %% (%d tâches, %.0f s)\n", titre, 100f / length(taches), length(taches), time() - t_debut)
                end
            end
        end
    end
    return DataFrame(reduce(vcat, sorties))
end

"""
    rejouer(cel, mg, ctx, N, j; source, scenario, taux) -> NamedTuple

Un tirage seul, reconstruit depuis sa graine : ses histogrammes et ses neuf lectures,
pour déboguer un échec sans relancer le Monte-Carlo.
"""
function rejouer(cel::Cellule, mg::ModeleGlobal, ctx::Contexte, N::Real, j::Integer; source::Symbol = :reel,
                 scenario::Symbol = :aucun, taux::Real = 0.0)
    ic = indice(mg, cel.cl)
    hB, hJ = somme(cel, BLEU), somme(cel, JAUNE)
    p = p_amincissement(cel, N)
    bB, bJ = zeros(256), zeros(256)
    rng = StableRNG(graine_tirage(cel.λ, cel.cl, N, j; source))
    tirage!(bB, bJ, hB, hJ, mg.μB[ic], mg.μJ[ic], p, rng, source)
    if scenario !== :aucun
        inj = facteurs_injection(mg, ctx, ic, scenario, taux)
        injecter!(bB, bJ, inj, p, StableRNG(graine_tirage(cel.λ, cel.cl, N, j; source, scenario, taux)))
    end
    o, e = lire(ctx, bB, bJ, p * mg.duree_somme_s[ic], 1)
    return (bleu = bB, jaune = bJ, p = p, lectures = Dict(k => (o[q], e[q]) for (q, k) in enumerate(LECTURES)))
end
