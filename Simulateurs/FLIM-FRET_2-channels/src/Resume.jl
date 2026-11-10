# Resume.jl — tableaux chiffrés (Markdown) tirés des fichiers plats, pour le
# README et pour relire les résultats sans ouvrir les CSV.

module Resume

using DataFrames, CSV, Printf, Statistics

sel(df; kw...) = df[reduce(.&, [df[!, k] .== v for (k, v) in kw]; init = trues(nrow(df))), :]
const ORDRE = ["1", "5d", "5a", "5c", "4", "3", "6", "2F", "2", "2H"]
const NOMS = Dict("1" => "1 ratio J/B", "2" => "2 mono (ton code)", "2H" => "2H mono + Hill (ton pipeline)",
    "2F" => "2F mono, décalage et fond figés",
    "3" => "3 biexp τ libres", "4" => "4 biexp τ figés", "5a" => "5a conjoint τ figés", "5b" => "5b conjoint τ faux",
    "5c" => "5c conjoint τ_A libre", "5d" => "5d conjoint échelle commune", "6" => "6 phasor")
const FENETRE = [5.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0]

fmt(x; d = 2) = isfinite(x) ? replace(@sprintf("%.*f", d, x), "." => ",") : "—"
fmtp(x) = isfinite(x) ? replace(@sprintf("%.0f %%", 100x), "." => ",") : "—"
budgets(N) = Dict(1e2 => "10²", 3e2 => "3·10²", 1e3 => "10³", 3e3 => "3·10³", 1e4 => "10⁴", 3e4 => "3·10⁴", 1e5 => "10⁵")[N]

function tableau(io, d, colonne, Ns; f = fmt, borne = nothing, cles = ORDRE)
    println(io, "| Estimateur | ", join(budgets.(Ns), " | "), " |")
    println(io, "| --- |", repeat(" ---: |", length(Ns)))
    for k in cles
        v = [begin
            s = sel(d; estimateur = k, N_photons = N)
            nrow(s) == 0 ? NaN : s[1, colonne]
        end for N in Ns]
        all(isnan, v) && continue
        println(io, "| ", NOMS[k], " | ", join(f.(v), " | "), " |")
    end
    if borne !== nothing
        v = [sel(d; estimateur = "1", N_photons = N)[1, borne] for N in Ns]
        println(io, "| *borne de Cramér-Rao (toute l'information)* | ", join(fmt.(v), " | "), " |")
    end
    println(io)
end

function ecrire(df, nu, chemin)
    Ns = [1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5]
    open(chemin, "w") do io
        println(io, "# Résumé chiffré\n")
        println(io, "Généré par `step6_figures.jl` à partir de `resultats_principaux.csv` et `nuisances.csv`. ",
            "Normalisation A (temps d'acquisition fixé ; N = photons bleus attendus à 20 mM), calibration exacte sauf mention. ",
            "Statistiques sur les réalisations réussies, estimations ramenées à [0, 150] mM.\n")
        d = sel(df; C_vrai = 20.0, normalisation = "A", mode_calibration = "exacte")
        println(io, "## À 20 mM : σ[Cl⁻] (mM)\n")
        tableau(io, d, :sigma_mM, Ns; borne = :CRLB_mM)
        println(io, "## À 20 mM : biais moyen (mM)\n")
        tableau(io, d, :biais_mM, Ns)
        println(io, "## À 20 mM : taux d'échec\n")
        tableau(io, d, :taux_echec, Ns; f = fmtp)
        println(io, "## À 20 mM : efficacité par rapport à la borne propre de chaque estimateur\n")
        tableau(io, d, :efficacite_propre, Ns)
        # RMSE moyenne sur la fenêtre de clamping
        println(io, "## Fenêtre de clamping (5–30 mM) : RMSE moyenne (mM), par normalisation\n")
        for norm in ("A", "B", "C")
            println(io, "**Normalisation $norm**\n")
            dn = sel(df; normalisation = norm, mode_calibration = "exacte")
            dn = dn[in.(dn.C_vrai, Ref(FENETRE)), :]
            g = combine(groupby(dn, [:estimateur, :N_photons]), :RMSE_mM => mean => :RMSE_mM, :taux_echec => maximum => :taux_echec)
            tableau(io, g, :RMSE_mM, Ns)
        end
        if "variance_parametrique_mM2" in names(df)
            println(io, "## Incertitude des paramètres ajustés (étape 11), normalisation A\n")
            println(io, "σ photonique (bruit de Poisson), σ paramétrique (200 vecteurs tirés des covariances de l'ajustement global et de ",
                "l'isotherme, lus sans bruit) et σ total = √(somme des variances), en mM.\n")
            for N in (1e3, 1e4)
                println(io, "**N = $(budgets(N))**\n")
                println(io, "| Estimateur | ", join(["$(Int(c)) mM : phot. / param. / total" for c in (5.0, 20.0, 40.0)], " | "), " |")
                println(io, "| --- |", repeat(" ---: |", 3))
                for k in ORDRE
                    v = map((5.0, 20.0, 40.0)) do c
                        s_ = sel(df; estimateur = k, N_photons = N, normalisation = "A", mode_calibration = "exacte", C_vrai = c)
                        nrow(s_) == 0 ? "—" : string(fmt(s_.sigma_mM[1]), " / ", fmt(sqrt(s_.variance_parametrique_mM2[1])), " / ", fmt(s_.sigma_total_mM[1]))
                    end
                    println(io, "| ", NOMS[k], " | ", join(v, " | "), " |")
                end
                println(io)
            end
        end
        if nu !== nothing && any(nu.parametre_nuisance .== "reference")
            println(io, "## Biais ajouté par les deux nuisances exclues jusqu'ici (étape 5), à 20 mM et 10⁴ photons\n")
            println(io, "Biais en calibration nominale moins biais du monde nominal avec les mêmes germes (balayage « reference »), en mM.\n")
            r = sel(nu; parametre_nuisance = "reference", C_vrai = 20.0, N_photons = 1e4)
            b0 = Dict(row.estimateur => row.biais_mM for row in eachrow(r))
            for p in ("decalage_jaune_ps", "blanchiment_accepteur", "blanchiment_accepteur_sombre", "blanchiment_donneur")
                dp = sel(nu; parametre_nuisance = p, C_vrai = 20.0, N_photons = 1e4, mode_calibration = "nominale")
                nrow(dp) == 0 && continue
                vals = sort(unique(dp.valeur_nuisance))
                println(io, "**$p**\n")
                println(io, "| Estimateur | ", join(fmt.(vals; d = 2), " | "), " |")
                println(io, "| --- |", repeat(" ---: |", length(vals)))
                for k in ORDRE
                    v = [begin
                        s_ = sel(dp; estimateur = k, valeur_nuisance = x)
                        nrow(s_) == 0 ? NaN : s_[1, :biais_mM] - get(b0, k, NaN)
                    end for x in vals]
                    all(isnan, v) && continue
                    println(io, "| ", NOMS[k], " | ", join(fmt.(v), " | "), " |")
                end
                println(io)
            end
        end
        if nu !== nothing && nrow(nu) > 0
            println(io, "## Nuisances, à 20 mM et 10⁴ photons\n")
            println(io, "Biais (mM) en calibration nominale (estimateurs et calibrations construits sur les valeurs nominales, données ",
                "issues du monde perturbé), puis σ (mM) en calibration exacte.\n")
            for p in unique(nu.parametre_nuisance)
                dp = sel(nu; parametre_nuisance = p, C_vrai = 20.0, N_photons = 1e4)
                vals = sort(unique(dp.valeur_nuisance))
                for (mode, col, lab) in (("nominale", :biais_mM, "biais"), ("exacte", :sigma_mM, "σ"))
                    dm = sel(dp; mode_calibration = mode)
                    nrow(dm) == 0 && continue
                    println(io, "**$p — $lab ($mode)**\n")
                    println(io, "| Estimateur | ", join(fmt.(vals; d = 2), " | "), " |")
                    println(io, "| --- |", repeat(" ---: |", length(vals)))
                    for k in vcat(ORDRE, "5b")
                        v = [begin
                            s = sel(dm; estimateur = k, valeur_nuisance = x)
                            nrow(s) == 0 ? NaN : s[1, col]
                        end for x in vals]
                        all(isnan, v) && continue
                        println(io, "| ", NOMS[k], " | ", join(fmt.(v), " | "), " |")
                    end
                    println(io)
                end
            end
        end
    end
    return chemin
end

end
