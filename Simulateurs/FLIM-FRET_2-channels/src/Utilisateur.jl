# Utilisateur.jl — l'estimateur 2 appelle TON ajustement mono-exponentiel
# (`vec_to_lifetime` de lifetime_analysis.jl), pas une réimplémentation.
#
# lifetime_analysis.jl garde son état dans des globales (RUNTIME, tampons FFT,
# caches) : il n'est pas sûr en parallèle, par conception. Plutôt que de le
# modifier, on le charge dans plusieurs modules anonymes, un par fil
# d'exécution : chaque instance a ses propres globales. Le fichier est lu tel
# quel ; seuls deux noms qu'il attend de FLIMApp sont fournis (un module
# FLIMCore vide, utilisé uniquement dans des signatures de fonctions d'import
# d'IRF jamais appelées ici, et DEFAULT_HISTOGRAM_RESOLUTION, comme le fait
# analyse_SCLM_2.jl).

"""Copie de lifetime_analysis.jl placée à côté du plan (même cœur d'ajustement que tests/)."""
const CHEMIN_LIFETIME_ANALYSIS = normpath(joinpath(@__DIR__, "..", "lifetime_analysis.jl"))

"""Pas du TAC de tes .sdt (50,034 ns / gain 4 / 256) : la valeur par défaut de `vec_to_lifetime`."""
const PAS_TAC_UTILISATEUR = 0.04886091184430619

struct InstanceUtilisateur
    mod::Module
    vec_to_lifetime::Function
    conv_irf_data::Function
    get_irf_for_channels::Function
    get_gating_function::Function
    get_x_data::Function
end

"""
    charger_instance_utilisateur(irf_mat; chemin)

Charge lifetime_analysis.jl dans un module neuf, y installe l'IRF exactement
comme analyse_SCLM_2.jl (matrice [t_ns comptages] de `load_irf_from_sdt`), puis
le préchauffe (sinon le premier ajustement dépasse le time_limit d'Optim).
"""
function charger_instance_utilisateur(irf_mat::Matrix{Float64}; chemin::AbstractString = CHEMIN_LIFETIME_ANALYSIS)
    m = Module(gensym(:lifetime_analysis))
    Core.eval(m, Expr(:toplevel, :(module FLIMCore
        struct Reglages end
    end)))
    Core.eval(m, :(const DEFAULT_HISTOGRAM_RESOLUTION = 256))
    Core.eval(m, :(const LASER_PULSE_PERIOD = 12.5))
    Base.include(m, chemin)
    Core.eval(m, :(__init__()))
    rt = Core.eval(m, :(RUNTIME[]))
    rt.irf = irf_mat
    rt.irf_bin_size = Base.invokelatest(Core.eval(m, :compute_irf_bin_size), irf_mat)
    rt.tcspc_window_size = round(irf_mat[end, 1] + irf_mat[2, 1], sigdigits = 4)
    Base.invokelatest(Core.eval(m, :warmup_lifetime_fitting!))
    f(s) = Core.eval(m, s)
    return InstanceUtilisateur(m, f(:vec_to_lifetime), f(:conv_irf_data), f(:get_irf_for_channels),
        f(:get_gating_function), f(:get_x_data))
end

"""
    ajuster_mono_utilisateur(inst, h; tentatives=5) -> (params, essais)

`vec_to_lifetime(h; guess=[3.0, 0.0, 1e-6], first_fit=true)` comme
`fit_decay` d'analyse_SCLM_2.jl : jusqu'à 5 essais si le résultat est NaN.
`hist_size_threshold = 0` : sans cela ton code refuse tout histogramme de
moins de 500 photons (NaN), ce qui rendrait 10² et 3·10² ininterprétables.
`params = [τ, décalage IRF (canaux), fraction de fond]`. `fixes = (décalage,
fond)` : les fige par ton argument `fixed_parameters` (seul τ reste libre).
"""
function ajuster_mono_utilisateur(inst::InstanceUtilisateur, h::Vector{Float64}; tentatives::Int = 5,
                                  fixes::Union{Nothing,NTuple{2,Float64}} = nothing)
    guess = fixes === nothing ? [3.0, 0.0, 1e-6] : [3.0, fixes[1], fixes[2]]
    fixed = fixes === nothing ? fill(NaN, 3) : [NaN, fixes[1], fixes[2]]
    for k in 1:tentatives
        p, _ = Base.invokelatest(inst.vec_to_lifetime, h; guess = copy(guess), fixed_parameters = copy(fixed),
            histogram_resolution = 256, first_fit = true, hist_size_threshold = 0)
        isnan(p[1]) || return (p, k)
    end
    return (fill(NaN, 3), tentatives)
end

"""Bornes de τ dans `mle_reconvolution_fit` : [2·pas de l'IRF, 2·période]."""
bornes_tau_utilisateur(irf_mat) = (2 * (irf_mat[2, 1] - irf_mat[1, 1]), 2 * 12.5)

"""
    chi2_utilisateur(inst, params, h)

χ²_r de Pearson de ton ajustement mono-exponentiel, recalculé comme
`reconstruct_model_counts` + `fit_reduced_chi_squared` de tests/lifetime_analysis.jl
(mêmes fenêtres, même normalisation), avec les fonctions de ton fichier.
"""
function chi2_utilisateur(inst::InstanceUtilisateur, params::Vector{Float64}, h::Vector{Float64};
                          tac_low_cut = 5.0980392, tac_high_cut = 94.901962)
    n = length(h)
    x_data = Base.invokelatest(inst.get_x_data, n, PAS_TAC_UTILISATEUR)
    y = copy(h)
    replace!(v -> v <= 0 ? 1e-256 : v, y)
    irf_local = Base.invokelatest(inst.get_irf_for_channels, n)
    lo = max(round(Int, tac_low_cut / 100 * n), 1)
    hs = clamp(round(Int, tac_high_cut / 100 * n), 1, n + 1)
    hi_idx = n - hs
    gating = Base.invokelatest(inst.get_gating_function, n, lo, hs)
    N = float(sum(y))
    n_actifs = max(1, n - lo - hi_idx)
    fond = params[3]
    conv = Base.invokelatest(inst.conv_irf_data, x_data, (params[1], params[2], params[3]), irf_local;
        histogram_resolution = n, number_of_previous_pulses = 5, laser_pulse_period = 12.5,
        tcspc_low_cut_index = lo, tcspc_high_cut_index = hi_idx)
    modele = ((1 - fond) .* conv .+ fond / n_actifs) .* N
    replace!(v -> v <= 0 ? 1e-256 : v, modele)
    r2 = gating .* (y .- modele) .^ 2
    return sum(r2 ./ modele) / max(count(==(1), gating) - length(params), 1)
end

"""
Réserve d'instances pour les boucles parallèles : chaque tâche en prend une,
s'en sert, la rend. Créées en séquence (planification FFTW, compilation).
"""
function reserve_utilisateur(irf_mat::Matrix{Float64}, n::Int = Threads.nthreads(); chemin = CHEMIN_LIFETIME_ANALYSIS)
    ch = Channel{InstanceUtilisateur}(n)
    for _ in 1:n
        put!(ch, charger_instance_utilisateur(irf_mat; chemin))
    end
    return ch
end
