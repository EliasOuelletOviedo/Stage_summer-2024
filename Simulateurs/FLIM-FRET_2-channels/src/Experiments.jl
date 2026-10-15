# Experiments.jl — balayages Monte-Carlo, parallélisation, fichier plat.
#
# Une cellule = (C vrai, budget N, normalisation, nuisance). Les réalisations
# d'une cellule sont réparties entre les fils par un compteur atomique ; chaque
# fil possède son espace FFT et son instance de ton code ; le germe de chaque
# réalisation ne dépend que de (cellule, indice), donc les résultats ne
# dépendent ni du nombre de fils ni de l'ordre d'exécution.

const GRILLE_C = [0.0, 1.0, 2.0, 5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 60.0, 100.0]
const GRILLE_N = [1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5]
repetitions(N) = N < 1e3 ? 2000 : 1000

"""Une case de la grille."""
Base.@kwdef struct Cellule
    id::Int
    C::Float64
    N::Float64
    norm::Symbol
    R::Int
    nuisance::String = "aucune"
    valeur::Float64 = NaN
end

"""
Un jeu d'estimation : un mode de calibration (« nominale », « exacte »…) et la
liste (clé, estimateur, calibration). Un même objet estimateur peut figurer
dans plusieurs jeux : son observable n'est alors calculée qu'une fois.
"""
struct Jeu
    mode::String
    entrees::Vector{Tuple{String,Estimateur,Calibration}}
end

"""Contextes de calcul, un par fil (créés en séquence : FFTW et compilation)."""
function creer_contextes(m::Modele, irf_mat; n = Threads.nthreads(), avec_utilisateur = true)
    return [Contexte(EspaceTravail(m.inst), avec_utilisateur ? charger_instance_utilisateur(irf_mat) : nothing) for _ in 1:n]
end

"""
    simuler_cellule(cell, m_vrai, jeux, ctxs; crlb=true) -> Vector{NamedTuple}

Une ligne par (jeu, estimateur) : biais, σ, RMSE (en mM, sur les réalisations
réussies), versions robustes, bornes de Cramér-Rao, taux d'échec et de
sortie de plage.
"""
function simuler_cellule(cell::Cellule, m_vrai::Modele, jeux::Vector{Jeu}, ctxs::Vector{Contexte}; crlb::Bool = true)
    ws0 = ctxs[1].ws
    μB, μY = attendus(cell.C, m_vrai, cell.N, cell.norm, ws0)
    uniques = Estimateur[]
    for j in jeux, (_, e, _) in j.entrees
        any(u -> u === e, uniques) || push!(uniques, e)
    end
    U, R = length(uniques), cell.R
    obs = fill(NaN, U, R)
    stat = fill(:echec, U, R)
    suivant = Threads.Atomic{Int}(1)
    taches = map(ctxs) do ctx
        Threads.@spawn begin
            nB, nY = zeros(length(μB)), zeros(length(μY))
            while true
                r = Threads.atomic_add!(suivant, 1)
                r > R && break
                rng = StableRNG(graine(cell.id, r))
                tirer_deux!(nB, nY, μB, μY, rng)
                for (u, e) in enumerate(uniques)
                    o, st = try
                        observable(e, nB, nY, ctx)
                    catch err
                        err isa InterruptException && rethrow()
                        (NaN, :echec)
                    end
                    obs[u, r] = o
                    stat[u, r] = st
                end
            end
        end
    end
    foreach(wait, taches)

    x_vrai = fraction_active(cell.C, m_vrai.p)
    bornes = Dict{Symbol,Float64}()
    function borne(def)
        crlb || return NaN
        get!(bornes, def) do
            try
                crlb_mM(def, cell.C, m_vrai, cell.N, cell.norm, ws0)
            catch
                NaN
            end
        end
    end
    lignes = NamedTuple[]
    for j in jeux, (cle, e, cal) in j.entrees
        u = findfirst(v -> v === e, uniques)
        ok = [stat[u, r] !== :echec for r in 1:R]
        Ĉ = Float64[]
        hp = Bool[]
        for r in 1:R
            ok[r] || continue
            c, h = inverser(cal, obs[u, r])
            push!(Ĉ, c)
            push!(hp, h)
        end
        nok = length(Ĉ)
        err = Ĉ .- cell.C
        σ = nok > 1 ? std(err) : NaN
        crlb_abs = borne(:absolue)
        crlb_p = borne(get(BORNE_PROPRE, cle, :absolue))
        o_ok = obs[u, ok]
        push!(lignes, (
            C_vrai = cell.C, x_vrai = x_vrai, estimateur = cle, nom = get(NOMS, cle, nom(e)), N_photons = cell.N,
            normalisation = String(cell.norm), N_bleu_attendu = sum(μB), N_jaune_attendu = sum(μY),
            parametre_nuisance = cell.nuisance, valeur_nuisance = cell.valeur, mode_calibration = j.mode,
            n_realisations = R, n_reussites = nok,
            biais_mM = nok > 0 ? mean(err) : NaN, sigma_mM = σ,
            RMSE_mM = nok > 0 ? sqrt(mean(abs2, err)) : NaN,
            biais_median_mM = nok > 0 ? median(err) : NaN,
            sigma_robuste_mM = nok > 1 ? 1.4826 * median(abs.(err .- median(err))) : NaN,
            CRLB_mM = crlb_abs, efficacite = crlb_abs^2 / σ^2,
            CRLB_propre_mM = crlb_p, efficacite_propre = crlb_p^2 / σ^2,
            taux_echec = 1 - nok / R, taux_hors_plage = nok > 0 ? mean(hp) : NaN,
            obs_moyenne = nok > 0 ? mean(o_ok) : NaN, obs_sigma = nok > 1 ? std(o_ok) : NaN,
        ))
    end
    return lignes
end

"""
    executer(cellules, m_vrai_de, jeux_de, ctxs; fichier) -> DataFrame

`m_vrai_de(cell)` et `jeux_de(cell)` donnent le modèle vrai et les jeux
d'estimation de chaque cellule (ils changent dans les balayages de nuisance).
Le fichier plat est réécrit après chaque cellule (reprise possible).
"""
function executer(cellules::Vector{Cellule}, m_vrai_de, jeux_de, ctxs::Vector{Contexte};
                  fichier::AbstractString, journal = stdout, deja = DataFrame())
    df = deja
    t0 = time()
    for (i, cell) in enumerate(cellules)
        t = @elapsed lignes = simuler_cellule(cell, m_vrai_de(cell), jeux_de(cell), ctxs)
        df = vcat(df, DataFrame(lignes); cols = :union)
        mkpath(dirname(fichier))
        CSV.write(fichier, df)
        écoulé = time() - t0
        @printf(journal, "[%4d/%d] C=%6.1f N=%8.0f %s %s=%-8s %6.1f s  (total %.1f min, reste ≈ %.1f min)\n",
            i, length(cellules), cell.C, cell.N, cell.norm, cell.nuisance, string(cell.valeur), t,
            écoulé / 60, écoulé / i * (length(cellules) - i) / 60)
        flush(journal)
    end
    return df
end
