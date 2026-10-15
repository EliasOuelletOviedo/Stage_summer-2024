# Parametres.jl — parametres.toml (+ sorties de step0 et step1) → Modele.
#
# Une valeur écrite "titration" vient de l'ajustement global de ta titration
# (sorties/titration/parametres_<protéine>.toml) ; "spectres" vient de
# sorties/spectres/eta.toml ; "fpbase" du rendement quantique FPbase du
# donneur choisi. Toute valeur numérique écrite à la main l'emporte : c'est
# ainsi qu'une mesure remplace une valeur inférée.

const DOSSIER_SIMULATEUR = normpath(joinpath(@__DIR__, ".."))
const FICHIER_PARAMETRES = joinpath(DOSSIER_SIMULATEUR, "parametres.toml")
dossier_sorties(sous...) = joinpath(DOSSIER_SIMULATEUR, "sorties", sous...)

"""Paramètres de la titration d'une protéine (écrits par step1_titration.jl)."""
function lire_titration_toml(proteine::AbstractString)
    f = dossier_sorties("titration", "parametres_$(proteine).toml")
    isfile(f) || error("$f absent : lancer d'abord step1_titration.jl")
    return TOML.parsefile(f)
end

"""
    charger_configuration(; fichier) -> (modele, irf_mat, config)

`modele` : le monde nominal (vrai par défaut) ; `irf_mat` : l'IRF telle que ton
code la reçoit ([t_ns comptages], médiane soustraite) ; `config` : le TOML lu.
"""
function charger_configuration(; fichier::AbstractString = FICHIER_PARAMETRES)
    cfg = TOML.parsefile(fichier)
    cap, opt, ins = cfg["capteur"], cfg["optique"], cfg["instrument"]
    prot = cap["proteine"]
    tit = nothing
    val(d, k) = begin
        v = d[k]
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
    Φ_D = cap["Phi_D"] == "fpbase" ? Float64(eta["proprietes"][donneur]["qy"]) : Float64(cap["Phi_D"])
    Φ_A = cap["Phi_A"] == "fpbase" ? Float64(eta["proprietes"][eta["accepteur"]]["qy"]) : Float64(cap["Phi_A"])
    p = SensorParams(τ_D0 = val(cap, "tau_D0"), τ_DA = val(cap, "tau_DA"), τ_A = Float64(cap["tau_A"]),
        x₀ = val(cap, "x0"), x_sat = val(cap, "x_sat"), K_D = val(cap, "K_D"), n = val(cap, "n"),
        Φ_D = Φ_D, Φ_A = Φ_A, mecanisme = Symbol(cap["mecanisme"]))
    if opt["eta"] != "spectres"
        e = opt["eta"]
        η = SMatrix{2,2}([e[1][1] e[1][2]; e[2][1] e[2][2]])
    end
    o = OpticsParams(η = η, ε = Float64(opt["epsilon"]), εΦ_dark = Float64(opt["epsilonPhi_dark"]))
    irf_mat = charger_irf_utilisateur(ins["irf"])
    f = ins["fenetre"]
    inst = InstrumentParams(n_bins = ins["n_bins"], T_rep = Float64(ins["T_rep"]), irf = irf_normalisee(irf_mat),
        t₀ = Float64(ins["t0"]), fenetre = f[1]:f[2])
    return Modele(p, o, inst), irf_mat, cfg
end

"""Résumé lisible d'un modèle (journal des étapes)."""
function decrire(m::Modele)
    p, o = m.p, m.o
    return @sprintf("τ_D0=%.3f τ_DA=%.3f τ_A=%.2f ns  x₀=%.3f x_sat=%.3f K_D=%.1f mM n=%.2f  Φ_D=%.2f Φ_A=%.2f  η=[%.3f %.4f; %.3f %.3f] ε=%.3f εΦ_dark=%.3f  %s",
        p.τ_D0, p.τ_DA, p.τ_A, p.x₀, p.x_sat, p.K_D, p.n, p.Φ_D, p.Φ_A, o.η[1, 1], o.η[1, 2], o.η[2, 1], o.η[2, 2],
        o.ε, o.εΦ_dark, p.mecanisme)
end
