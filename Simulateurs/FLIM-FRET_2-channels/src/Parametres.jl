# Parametres.jl — parametres.toml (+ sorties de step0 et step1) → Modele.
#
# Chaque paramètre a un statut : « ajusté » (tiré de données par le simulateur),
# « supposé » (catalogue, spectres, hypothèse) ou « mesuré » (au banc, avec date
# et incertitude). Une valeur "titration" vient de l'ajustement global de la
# titration (sorties/titration/parametres_<protéine>.toml) ; "spectres" de
# sorties/spectres/eta.toml ; "fpbase" du rendement quantique FPbase. Toute
# valeur numérique écrite à la main l'emporte : c'est ainsi qu'une mesure
# remplace une valeur supposée.

const DOSSIER_SIMULATEUR = normpath(joinpath(@__DIR__, ".."))
const FICHIER_PARAMETRES = joinpath(DOSSIER_SIMULATEUR, "parametres.toml")
dossier_sorties(sous...) = joinpath(DOSSIER_SIMULATEUR, "sorties", sous...)

"""Paramètres de la titration d'une protéine (écrits par step1_titration.jl)."""
function lire_titration_toml(proteine::AbstractString)
    f = dossier_sorties("titration", "parametres_$(proteine).toml")
    isfile(f) || error("$f absent : lancer d'abord step1_titration.jl")
    return TOML.parsefile(f)
end

const STATUTS = ("ajusté", "supposé", "mesuré")

"""Valeur d'une entrée de parametres.toml : un nombre/mot-clé, ou une table `{ valeur = …, statut = … }`."""
valeur_cfg(e) = e isa AbstractDict ? e["valeur"] : e

"""
Provenance d'une entrée : (statut, source, date, incertitude). Une entrée sans
table est « supposée ». Un statut « mesuré » sans date ni incertitude est refusé :
dans six mois, il faut savoir d'où vient chaque nombre d'une figure.
"""
function provenance(nom::AbstractString, e)
    e isa AbstractDict || return (statut = "supposé", source = "", date = "", incertitude = nothing)
    st = get(e, "statut", "supposé")
    st in STATUTS || error("parametres.toml, $nom : statut « $st » inconnu (attendu : $(join(STATUTS, ", ")))")
    if st == "mesuré"
        (haskey(e, "date") && haskey(e, "incertitude")) ||
            error("parametres.toml, $nom : un paramètre « mesuré » doit porter sa date et son incertitude")
    end
    return (statut = st, source = get(e, "source", ""), date = string(get(e, "date", "")), incertitude = get(e, "incertitude", nothing))
end

"""
    charger_configuration(; fichier) -> (modele, irf_mat, config)

`modele` : le monde nominal (vrai par défaut) ; `irf_mat` : l'IRF telle que le
code de l'utilisateur la reçoit ([t_ns comptages], médiane soustraite) ;
`config` : le TOML lu, augmenté de `config["provenance"]` (statut de chaque
paramètre, voir `provenance`).
"""
function charger_configuration(; fichier::AbstractString = FICHIER_PARAMETRES)
    cfg = TOML.parsefile(fichier)
    cap, opt, ins = cfg["capteur"], cfg["optique"], cfg["instrument"]
    exp_ = get(cfg, "experience", Dict{String,Any}())
    prot = cap["proteine"]
    tit = nothing
    prov = Dict{String,Any}()
    function val(d, k; defaut = nothing)
        haskey(d, k) || return defaut
        prov[k] = provenance(k, d[k])
        v = valeur_cfg(d[k])
        if v == "titration"
            tit === nothing && (tit = lire_titration_toml(prot))
            return Float64(tit[k])
        end
        return v
    end
    η, eta = lire_eta(dossier_sorties("spectres", "eta.toml"))
    donneur = get(cap, "donneur", "mcerulean3")
    if donneur != eta["donneur"]
        e = eta["eta_donneur_$(donneur)"]
        η = SMatrix{2,2}([e["eta_BD"] η[1, 2]; e["eta_YD"] η[2, 2]])
    end
    ΦD = val(cap, "Phi_D")
    Φ_D = ΦD == "fpbase" ? Float64(eta["proprietes"][donneur]["qy"]) : Float64(ΦD)
    ΦA = val(cap, "Phi_A")
    Φ_A = ΦA == "fpbase" ? Float64(eta["proprietes"][eta["accepteur"]]["qy"]) : Float64(ΦA)
    ΦAe = val(cap, "Phi_A_eff"; defaut = "catalogue")
    Φ_A_eff = ΦAe == "catalogue" ? NaN : Float64(ΦAe)
    nombre(v) = Float64(v)
    p = SensorParams(τ_D0 = val(cap, "tau_D0"), τ_DA = val(cap, "tau_DA"), τ_A = nombre(val(cap, "tau_A")),
        x₀ = val(cap, "x0"), x_sat = val(cap, "x_sat"), K_D = val(cap, "K_D"), n = val(cap, "n"),
        Φ_D = Φ_D, Φ_A = Φ_A, Φ_A_eff = Φ_A_eff, mecanisme = Symbol(cap["mecanisme"]),
        β_A = nombre(val(exp_, "beta_A"; defaut = 0.0)), β_D = nombre(val(exp_, "beta_D"; defaut = 0.0)))
    ev = val(opt, "eta")
    if ev != "spectres"
        η = SMatrix{2,2}([ev[1][1] ev[1][2]; ev[2][1] ev[2][2]])
    end
    eYD = val(opt, "eta_YD"; defaut = "spectres")
    eYD == "spectres" || (η = SMatrix{2,2}([η[1, 1] η[1, 2]; Float64(eYD) η[2, 2]]))
    o = OpticsParams(η = η, ε = nombre(val(opt, "epsilon")), εΦ_dark = nombre(val(opt, "epsilonPhi_dark")))
    irf_mat = charger_irf_utilisateur(val(ins, "irf"))
    f = ins["fenetre"]
    inst = InstrumentParams(n_bins = ins["n_bins"], T_rep = Float64(ins["T_rep"]), irf = irf_normalisee(irf_mat),
        t₀ = nombre(val(ins, "t0")), fenetre = f[1]:f[2], Δt_jaune = nombre(val(ins, "dt_jaune"; defaut = 0.0)))
    cfg["provenance"] = prov
    return Modele(p, o, inst), irf_mat, cfg
end

"""Tableau lisible des statuts (journal des étapes, README)."""
function decrire_provenance(cfg)
    lignes = String[]
    for (k, v) in sort(collect(cfg["provenance"]); by = first)
        inc = v.incertitude === nothing ? "" : " ± $(v.incertitude)"
        dat = isempty(v.date) ? "" : " ($(v.date))"
        push!(lignes, @sprintf("  %-16s %-8s%s%s", k, v.statut, dat, inc))
    end
    return join(lignes, "\n")
end

"""Résumé lisible d'un modèle (journal des étapes)."""
function decrire(m::Modele)
    p, o = m.p, m.o
    return @sprintf("τ_D0=%.3f τ_DA=%.3f τ_A=%.2f ns  x₀=%.3f x_sat=%.3f K_D=%.1f mM n=%.2f  Φ_D=%.2f Φ_A=%.2f (effectif %.2f)  η=[%.3f %.4f; %.3f %.3f] ε=%.3f εΦ_dark=%.3f  Δt_jaune=%.0f ps  β_A=%.2f β_D=%.2f  %s",
        p.τ_D0, p.τ_DA, p.τ_A, p.x₀, p.x_sat, p.K_D, p.n, p.Φ_D, p.Φ_A, phi_A(p), o.η[1, 1], o.η[1, 2], o.η[2, 1], o.η[2, 2],
        o.ε, o.εΦ_dark, 1000 * m.inst.Δt_jaune, p.β_A, p.β_D, p.mecanisme)
end
