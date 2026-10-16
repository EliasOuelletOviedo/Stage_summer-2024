# Figures.jl — les figures du simulateur (CairoMakie), chargé seulement par
# les étapes qui tracent.
#
# Charte : surface #fcfcfb comme les autres analyses SCLM, grilles en filets
# pleins discrets, une couleur ET une forme de marqueur fixes par estimateur
# dans toutes les figures (la couleur suit l'estimateur, jamais son rang), texte
# dans les encres, jamais dans la couleur d'une série. Palette catégorielle à
# 8 emplacements validée (CVD ΔE ≥ 9,1 entre voisins) ; trois teintes sont sous
# 3:1 de contraste, d'où les marqueurs distincts et la légende systématique.

module Figures

using CairoMakie, DataFrames, CSV, Printf, Statistics, TOML

hexa(s) = Makie.to_color(s)

const SURFACE = hexa("#fcfcfb")
const ENCRE = hexa("#0b0b0b")
const ENCRE2 = hexa("#52514e")
const MUET = hexa("#898781")
const GRILLE = hexa("#e1e0d9")
const AXE = hexa("#c3c2b7")
const ZONE = hexa("#f0efec")

const ORDRE = ["1", "2", "3", "4", "5a", "5d", "6", "5c"]
const COULEUR = Dict(zip(ORDRE, hexa.(["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"])))
COULEUR["2H"] = COULEUR["2"]
COULEUR["2F"] = COULEUR["2"]
COULEUR["5b"] = COULEUR["5a"]
const MARQUEUR = Dict("1" => :circle, "2" => :rect, "2H" => :rect, "2F" => :rect, "3" => :diamond, "4" => :utriangle,
    "5a" => :dtriangle, "5b" => :dtriangle, "5d" => :pentagon, "6" => :hexagon, "5c" => :ltriangle)
const NOM_COURT = Dict("1" => "1 ratio J/B", "2" => "2 mono (ton code)", "2H" => "2H mono + Hill (ton pipeline)",
    "2F" => "2F mono, décalage et fond figés",
    "3" => "3 biexp τ libres", "4" => "4 biexp τ figés", "5a" => "5a conjoint τ figés", "5b" => "5b conjoint τ faux",
    "5c" => "5c conjoint τ_A libre", "5d" => "5d conjoint échelle commune", "6" => "6 phasor")

const NORMALISATION = Dict(
    "A" => "Normalisation A — temps d'acquisition fixé (N = photons bleus attendus à 20 mM ; le canal jaune s'y ajoute)",
    "B" => "Normalisation B — photons totaux fixés (N = N_bleu + N_jaune)",
    "C" => "Normalisation C — photons du canal bleu fixés (N = N_bleu)")

function theme_sim()
    return Theme(
        fontsize = 13, backgroundcolor = SURFACE, textcolor = ENCRE,
        Axis = (backgroundcolor = SURFACE, xgridcolor = GRILLE, ygridcolor = GRILLE, xgridwidth = 1, ygridwidth = 1,
            leftspinecolor = AXE, bottomspinecolor = AXE, topspinevisible = false, rightspinevisible = false,
            xtickcolor = AXE, ytickcolor = AXE, xticklabelcolor = ENCRE2, yticklabelcolor = ENCRE2,
            xlabelcolor = ENCRE2, ylabelcolor = ENCRE2, titlecolor = ENCRE, titlealign = :left, titlesize = 14,
            subtitlecolor = ENCRE2, xminorticksvisible = false, yminorticksvisible = false),
        Legend = (framevisible = false, labelcolor = ENCRE2, backgroundcolor = SURFACE, titlecolor = ENCRE2),
        Lines = (linewidth = 2, joinstyle = :round, linecap = :round),
        Scatter = (markersize = 10, strokewidth = 1.5, strokecolor = SURFACE),
    )
end

"""Coupe un texte en lignes d'au plus `n` caractères (le retour automatique de Makie décale le bloc)."""
function envelopper(s, n)
    # espaces insécables devant la ponctuation double : jamais de « ; » en début de ligne
    s = replace(s, " ;" => "\u00a0;", " :" => "\u00a0:", " ?" => "\u00a0?", " !" => "\u00a0!", "« " => "«\u00a0", " »" => "\u00a0»")
    lignes = String[]
    cur = ""
    for m in split(s, ' ')
        if isempty(cur)
            cur = m
        elseif length(cur) + 1 + length(m) <= n
            cur *= " " * m
        else
            push!(lignes, cur)
            cur = m
        end
    end
    isempty(cur) || push!(lignes, cur)
    return join(lignes, "\n")
end

"""Titre et sous-titre de figure (en haut à gauche, encres)."""
function entete!(fig, titre, sous_titre; largeur = 1)
    g = GridLayout(fig[0, 1:largeur]; tellwidth = false)
    Label(g[1, 1], titre; fontsize = 16, font = :bold, color = ENCRE, halign = :left, tellwidth = false)
    Label(g[2, 1], envelopper(sous_titre, floor(Int, size(fig.scene)[1] / 6.6)); fontsize = 12, color = ENCRE2,
        halign = :left, justification = :left, tellwidth = false)
    rowgap!(g, 1, 4)
end

"""Graduations log lisibles : 1, 2, 5, 10, 20, 50… (pas de 10^3.5)."""
function graduations_log(lo, hi)
    v = Float64[]
    for e in floor(Int, log10(lo)):ceil(Int, log10(hi)), m in (1, 2, 5)
        x = m * 10.0^e
        lo <= x <= hi && push!(v, x)
    end
    return (v, [x >= 1e4 ? replace(@sprintf("%.0f", x), r"(\d)(?=(\d{3})+$)" => s"\1 ") : x >= 1 ? @sprintf("%.0f", x) : replace(@sprintf("%g", x), "." => ",") for x in v])
end

"""
Une série par estimateur : ligne de 2 px et marqueurs cerclés de la couleur du
fond. `creux` : points dont le taux d'échec dépasse 5 % (marqueur vide) — leur
σ ne porte que sur les ajustements réussis.
"""
function serie!(ax, x, y, cle; style = :solid, label = NOM_COURT[cle], creux = falses(length(x)))
    ok = isfinite.(x) .& isfinite.(y)
    lines!(ax, x[ok], y[ok]; color = COULEUR[cle], linestyle = style, label = label)
    p = ok .& .!creux
    scatter!(ax, x[p], y[p]; color = COULEUR[cle], marker = MARQUEUR[cle], label = label)
    v = ok .& creux
    any(v) && scatter!(ax, x[v], y[v]; color = SURFACE, strokecolor = COULEUR[cle], strokewidth = 1.8, marker = MARQUEUR[cle])
end

echecs(s) = s.taux_echec .> 0.05
"""Variantes de ton code : même teinte, trait différent (2 plein, 2F pointillé, 2H tireté)."""
style_de(k) = k == "2H" ? :dash : k == "2F" ? :dot : :solid
const AVEC_2F = vcat(ORDRE, "2F")
const NOTE_CREUX = " Marqueur vide : plus de 5 % d'ajustements en échec (statistique sur les seuls réussis)."

sel(df; kw...) = df[reduce(.&, [df[!, k] .== v for (k, v) in kw]; init = trues(nrow(df))), :]

"""Bande grise de la fenêtre de clamping (5–35 mM)."""
fenetre_clamping!(ax) = vspan!(ax, 5, 35; color = (ZONE, 1.0))


"""
Bande d'incertitude systématique (changement 2) : enveloppe min–max, sur les mondes
réajustés (IRF ×0,9 et ×1,1, t₀ ± 3,5 ps, t₀ par histogramme) et le monde nominal,
de la colonne `col` selon `xcol`, pour l'estimateur `k`.
"""
function bande_mondes!(ax, mondes, nominal, k, xcol, col; filtre...)
    d = sel(mondes; estimateur = k, filtre...)
    nrow(d) == 0 && return nothing
    n0 = sel(nominal; estimateur = k, filtre...)
    tout = vcat(select(d, xcol, col), select(n0, xcol, col))
    g = combine(groupby(tout, xcol), col => (v -> minimum(filter(isfinite, v); init = Inf)) => :lo,
        col => (v -> maximum(filter(isfinite, v); init = -Inf)) => :hi)
    sort!(g, xcol)
    ok = isfinite.(g.lo) .& isfinite.(g.hi)
    band!(ax, g[ok, xcol], g.lo[ok], g.hi[ok]; color = (COULEUR[k], 0.16))
    return nothing
end
const CLES_BANDES = ["1", "5a", "4", "2H"]
const NOTE_BANDES = " Bandes teintées (1, 5a, 4, 2H) : enveloppe des mondes réajustés avec l'IRF ×0,9 et ×1,1, t₀ ± 3,5 ps et t₀ libre par histogramme (étape 11)."

# -----------------------------------------------------------------------------
# Figure 1 — σ[Cl⁻] selon [Cl⁻]
# -----------------------------------------------------------------------------
function figure1(df, chemin; N = 1e4, norm = "A", cles = AVEC_2F, mondes = nothing)
    fig = Figure(size = (980, 600))
    d = sel(df; N_photons = N, normalisation = norm, mode_calibration = "exacte")
    ax = Axis(fig[1, 1]; xlabel = "[Cl⁻] vraie (mM)", ylabel = "σ[Cl⁻] (mM, échelle log)", yscale = log10, yticks = graduations_log(0.1, 1000),
        title = "Lequel gagne, et où sur l'isotherme ?")
    fenetre_clamping!(ax)
    if mondes !== nothing
        for k in CLES_BANDES
            bande_mondes!(ax, sel(mondes; N_photons = N, normalisation = norm), d, k, :C_vrai, :sigma_mM)
        end
    end
    for k in cles
        s = sort(sel(d; estimateur = k), :C_vrai)
        nrow(s) == 0 && continue
        serie!(ax, s.C_vrai, s.sigma_mM, k; creux = echecs(s), style = style_de(k))
    end
    b = sort(sel(d; estimateur = "1"), :C_vrai)
    ok = isfinite.(b.CRLB_mM)
    lines!(ax, b.C_vrai[ok], b.CRLB_mM[ok]; color = ENCRE, linewidth = 1.5, label = "borne de Cramér-Rao (toute l'information)")
    text!(ax, 36, 1.0; text = "fenêtre de clamping", color = MUET, fontsize = 11, align = (:left, :bottom))
    i60 = findfirst(==(60.0), b.C_vrai)
    i60 === nothing || text!(ax, 61, 0.78 * b.CRLB_mM[i60]; text = "1 ratio et 5d : sur la borne", color = ENCRE2, fontsize = 11, align = (:left, :top))
    Legend(fig[1, 2], ax; merge = true, unique = true)
    entete!(fig, @sprintf("Fig. 1 — σ[Cl⁻] selon la concentration, %s photons", fmtN(N)),
        NORMALISATION[norm] * ". Calibration exacte de chaque estimateur ; σ sur les réalisations réussies, estimations ramenées à [0, 150] mM. " *
        "La borne est infinie à 0 mM (pente nulle de l'isotherme, n > 1)." * NOTE_CREUX * (mondes === nothing ? "" : NOTE_BANDES); largeur = 2)
    save(chemin, fig; px_per_unit = 2)
end

"""Budget en notation compacte : 10³, 3·10³…"""
function fmtN(N)
    e = floor(Int, log10(N) + 1e-9)
    m = round(Int, N / 10.0^e)
    exposant = map(c -> "⁰¹²³⁴⁵⁶⁷⁸⁹"[nextind("⁰¹²³⁴⁵⁶⁷⁸⁹", 0, c - '0' + 1)], string(e))
    return m == 1 ? "10" * exposant : "$(m)·10" * exposant
end

# -----------------------------------------------------------------------------
# Figure 2 — σ[Cl⁻] selon le budget (A, B, C)
# -----------------------------------------------------------------------------
function figure2(df, chemin; C = 20.0, cles = AVEC_2F, mondes = nothing)
    fig = Figure(size = (1350, 560))
    axes = Axis[]
    for (j, norm) in enumerate(("A", "B", "C"))
        d = sel(df; C_vrai = C, normalisation = norm, mode_calibration = "exacte")
        ax = Axis(fig[1, j]; xscale = log10, yscale = log10, xlabel = "budget N (photons)", ylabel = j == 1 ? "σ[Cl⁻] (mM)" : "",
            yticks = graduations_log(0.1, 1000),
            title = "Normalisation " * norm, xticks = ([1e2, 1e3, 1e4, 1e5], ["10²", "10³", "10⁴", "10⁵"]))
        if mondes !== nothing
            for k in CLES_BANDES
                bande_mondes!(ax, sel(mondes; C_vrai = C, normalisation = norm), d, k, :N_photons, :sigma_mM)
            end
        end
        for k in cles
            s = sort(sel(d; estimateur = k), :N_photons)
            nrow(s) == 0 && continue
            serie!(ax, s.N_photons, s.sigma_mM, k; creux = echecs(s), style = style_de(k))
        end
        b = sort(sel(d; estimateur = "1"), :N_photons)
        lines!(ax, b.N_photons, b.CRLB_mM; color = ENCRE, linewidth = 1.5, label = "borne de Cramér-Rao (toute l'information)")
        # pente −½ de référence
        y0 = 0.6 * b.CRLB_mM[end]
        lines!(ax, [1e3, 1e5], [y0 * 10, y0]; color = MUET, linewidth = 1)
        text!(ax, 1e3, y0 * 10; text = "pente −½", color = MUET, fontsize = 11, align = (:right, :center), offset = (-4, 0))
        push!(axes, ax)
    end
    linkyaxes!(axes...)
    Legend(fig[1, 4], axes[1]; merge = true, unique = true)
    entete!(fig, "Fig. 2 — σ[Cl⁻] selon le budget de photons, à $(Int(C)) mM",
        "Lequel tient à bas photons ? Une pente de −½ est le régime de Poisson ; un décrochage signale une rupture. " *
        "A : temps d'acquisition fixé (N = photons bleus à 20 mM) ; B : N_bleu + N_jaune = N ; C : N_bleu = N." * NOTE_CREUX *
        (mondes === nothing ? "" : NOTE_BANDES); largeur = 4)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 3 — biais selon le budget
# -----------------------------------------------------------------------------
function figure3(df, chemin; C = 20.0, norm = "A", cles = vcat(ORDRE, "2F", "2H"))
    fig = Figure(size = (980, 600))
    d = sel(df; C_vrai = C, normalisation = norm, mode_calibration = "exacte")
    ax = Axis(fig[1, 1]; xscale = log10, yscale = Makie.Symlog10(1.0), xlabel = "budget N (photons)",
        ylabel = "biais moyen sur [Cl⁻] (mM, symlog)", title = "Lesquels développent un biais quand les photons manquent ?",
        xticks = ([1e2, 1e3, 1e4, 1e5], ["10²", "10³", "10⁴", "10⁵"]), yticks = [-100, -10, -1, 0, 1, 10, 100])
    hlines!(ax, [0.0]; color = AXE, linewidth = 1)
    for k in cles
        s = sort(sel(d; estimateur = k), :N_photons)
        nrow(s) == 0 && continue
        serie!(ax, s.N_photons, s.biais_mM, k; style = style_de(k), creux = echecs(s))
    end
    Legend(fig[1, 2], ax; merge = true, unique = true)
    entete!(fig, "Fig. 3 — biais selon le budget de photons, à $(Int(C)) mM",
        NORMALISATION[norm] * ". Calibration exacte (sans bruit, sur le modèle vrai), sauf 2H : Hill ajustée sur ta grille de titration, comme ton pipeline. " *
        "Les estimations sont ramenées à [0, 150] mM : un estimateur qui sature contre 150 mM porte un biais positif." * NOTE_CREUX; largeur = 2)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 4 — biais selon la fuite spectrale η_YD
# -----------------------------------------------------------------------------
function figure4(nu, chemin; C = 20.0, cles = ORDRE, η_nominal = NaN)
    fig = Figure(size = (1250, 820))
    d = sel(nu; parametre_nuisance = "eta_YD_relatif", C_vrai = C)
    axes = Axis[]
    for (i, (mode, col, lab)) in enumerate((("nominale", :biais_mM, "biais (mM)"), ("exacte", :sigma_mM, "σ[Cl⁻] (mM, log)")))
        for (j, N) in enumerate((1e3, 1e4))
            s0 = sel(d; mode_calibration = mode, N_photons = N)
            ax = Axis(fig[i, j]; xlabel = "η_YD vrai / η_YD nominal", ylabel = j == 1 ? lab : "",
                yticks = col == :sigma_mM ? graduations_log(0.1, 1000) : [-100, -10, -1, 0, 1, 10, 100],
                yscale = col == :sigma_mM ? log10 : Makie.Symlog10(1.0),
                title = (mode == "nominale" ? "Calibration nominale — biais" : "Calibration exacte — précision") * ", N = " * fmtN(N))
            col == :biais_mM && hlines!(ax, [0.0]; color = AXE, linewidth = 1)
            vlines!(ax, [1.0]; color = AXE, linewidth = 1)
            for k in cles
                s = sort(sel(s0; estimateur = k), :valeur_nuisance)
                nrow(s) == 0 && continue
                serie!(ax, s.valeur_nuisance, s[!, col], k; creux = echecs(s))
            end
            push!(axes, ax)
        end
    end
    linkyaxes!(axes[1], axes[2])
    linkyaxes!(axes[3], axes[4])
    Legend(fig[1:2, 3], axes[1]; merge = true, unique = true)
    entete!(fig, "Fig. 4 — Fuite du donneur dans le canal jaune : le ratio s'effondre-t-il là où l'ajustement conjoint tient ?",
        "Normalisation A, $(Int(C)) mM. Les données viennent d'un monde où η_YD vaut x fois sa valeur nominale" *
        (isnan(η_nominal) ? "" : @sprintf(" (%.3f)", η_nominal)) * ". En haut : estimateurs et calibrations construits sur la valeur nominale " *
        "(calibration fausse) → biais. En bas : construits sur la vraie valeur → effet de la fuite sur la précision seule."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 4b — les autres nuisances (biais, calibration nominale)
# -----------------------------------------------------------------------------
function figure4b(nu, chemin; C = 20.0, N = 1e4, cles = ORDRE)
    panneaux = [("gain_jaune_relatif", "gain du canal jaune vrai / calibré"), ("epsilon", "excitation directe ε"),
        ("epsilonPhi_dark", "ε_dark·Φ_dark (forme protonée)"), ("Phi_A_relatif", "Φ_A effectif vrai / supposé"),
        ("tau_D0_relatif", "τ_D0 vrai / τ_D0 calibré"), ("tau_DA_relatif", "τ_DA vrai / τ_DA calibré"), ("tau_A_relatif", "τ_A vrai / τ_A supposé")]
    fig = Figure(size = (1750, 840))
    axes = Axis[]
    for (i, (p, xl)) in enumerate(panneaux)
        r, c = fldmod1(i, 4)
        s0 = sel(nu; parametre_nuisance = p, C_vrai = C, N_photons = N, mode_calibration = "nominale")
        nrow(s0) == 0 && continue
        ax = Axis(fig[r, c]; xlabel = xl, ylabel = c == 1 ? "biais (mM, symlog)" : "", yscale = Makie.Symlog10(1.0), title = xl,
            yticks = [-100, -10, -1, 0, 1, 10, 100])
        hlines!(ax, [0.0]; color = AXE, linewidth = 1)
        for k in vcat(cles, "5b")
            s = sort(sel(s0; estimateur = k), :valeur_nuisance)
            nrow(s) == 0 && continue
            serie!(ax, s.valeur_nuisance, s.biais_mM, k; creux = echecs(s))
        end
        push!(axes, ax)
    end
    isempty(axes) && return nothing
    linkyaxes!(axes...)
    Legend(fig[2, 4], axes[end]; merge = true, unique = true, tellheight = false, tellwidth = false)
    entete!(fig, "Fig. 4b — Biais sous les autres nuisances, calibration nominale",
        @sprintf("Normalisation A, %d mM, N = %s. Données issues du monde perturbé, estimateurs et calibrations construits sur les valeurs nominales ", Int(C), fmtN(N)) *
        "(gain 1 ; ε = 0,03 ; ε_dark·Φ_dark = 0 ; Φ_A = 0,60 ; τ de la titration ; τ_A = 3,0 ns). Dans les panneaux τ, l'estimateur conjoint est le 5b du plan (τ figés faux). " *
        "Le gain du jaune multiplie tout le canal jaune (dérive de détecteur) : il ne change pas sa forme ; Φ_A ne change que l'émission sensibilisée." * NOTE_CREUX; largeur = 4)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 8 — capteurs hypothétiques : σ selon τ_DA
# -----------------------------------------------------------------------------
function figure8(nu, df, chemin; C = 20.0, τ_DA_nominal = NaN, cles = ORDRE)
    fig = Figure(size = (1300, 600))
    axes = Axis[]
    for (j, N) in enumerate((1e3, 1e4))
        d = sel(nu; parametre_nuisance = "tau_DA_capteur_ns", C_vrai = C, N_photons = N, mode_calibration = "exacte")
        d0 = sel(df; C_vrai = C, N_photons = N, normalisation = "A", mode_calibration = "exacte")   # τ_DA de la titration
        ax = Axis(fig[1, j]; xlabel = "τ_DA du capteur (ns)", ylabel = j == 1 ? "σ[Cl⁻] (mM, log)" : "", yscale = log10,
            yticks = graduations_log(0.1, 1000), title = "N = " * fmtN(N))
        for k in cles
            s = sort(vcat(sel(d; estimateur = k), sel(d0; estimateur = k)[:, names(d)]; cols = :union), :valeur_nuisance)
            s.valeur_nuisance[isnan.(s.valeur_nuisance)] .= τ_DA_nominal
            s = sort(s, :valeur_nuisance)
            nrow(s) == 0 && continue
            serie!(ax, s.valeur_nuisance, s.sigma_mM, k; creux = echecs(s))
        end
        vlines!(ax, [τ_DA_nominal]; color = AXE, linewidth = 1)
        text!(ax, τ_DA_nominal, 0.12; text = " ta titration", color = MUET, fontsize = 11, align = (:left, :bottom))
        push!(axes, ax)
    end
    linkyaxes!(axes...)
    Legend(fig[1, 3], axes[1]; merge = true, unique = true)
    entete!(fig, "Fig. 8 — Si le capteur transférait moins : σ[Cl⁻] selon τ_DA, à $(Int(C)) mM",
        "Capteurs hypothétiques identiques à mCer3_Topaz (isotherme, τ_D0, τ_A, optique) sauf τ_DA, donc E = 1 − τ_DA/τ_D0 ; " *
        "calibration exacte, normalisation A. Le plan supposait τ_DA ≈ 2,1 ns ; ta titration donne 0,15 ns." * NOTE_CREUX; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 5 — K apparents
# -----------------------------------------------------------------------------
function figure5(kapp, courbes, chemin; K_vrai = NaN, mondes_K = nothing)
    fig = Figure(size = mondes_K === nothing ? (1300, 600) : (1750, 600))
    # (a) courbes de calibration normalisées : 0 à [Cl⁻] = 0, 1 à 150 mM
    ax1 = Axis(fig[1, 1]; xlabel = "[Cl⁻] (mM)", ylabel = "observable normalisée", title = "Courbes de calibration sans bruit")
    # courbes identiques une fois normalisées : regroupées sous une étiquette
    for (k, lab) in (("5a", "3, 4, 5a, 5d : lisent x (K vrai)"), ("1", "1 ratio J/B et 6 phasor (pondérés en intensité)"),
                     ("2", "2 mono (ton code)"))
        o = courbes[!, k]
        y = (o .- o[1]) ./ (o[end] - o[1])
        lines!(ax1, courbes.C_mM, y; color = COULEUR[k], label = lab)
    end
    vlines!(ax1, [K_vrai]; color = ENCRE, linewidth = 1)
    text!(ax1, K_vrai, 0.98; text = @sprintf(" K vrai = %.1f mM", K_vrai), color = ENCRE2, fontsize = 11, align = (:left, :top))
    axislegend(ax1; position = :rb)
    # (b) barres K_app / K_vrai
    ordre = ["5a", "4", "3", "5d", "6", "1", "tau_intensite", "2", "2H", "1inv"]
    lab = Dict("5a" => "5a conjoint", "4" => "4 biexp τ figés", "3" => "3 biexp τ libres", "5d" => "5d échelle commune",
        "6" => "6 phasor", "1" => "1 ratio J/B", "tau_intensite" => "⟨τ⟩ intensité (analytique)", "2" => "2 mono (ton code)",
        "2H" => "2H Hill sur ta grille", "1inv" => "ratio inversé B/J")
    rows = [k for k in ordre if any(kapp.estimateur .== k)]
    v = [kapp.rapport_K[findfirst(==(k), kapp.estimateur)] for k in rows]
    y = collect(1:length(rows))
    ax2 = Axis(fig[1, 2]; xlabel = "K apparent / K vrai", yticks = (y, [lab[k] for k in rows]), title = "Déformation de la constante de dissociation",
        ygridvisible = false)
    barplot!(ax2, y, v; direction = :x, color = COULEUR["1"], gap = 0.45, strokewidth = 0)
    vlines!(ax2, [1.0]; color = ENCRE, linewidth = 1)
    for (yi, vi) in zip(y, v)
        text!(ax2, vi, yi; text = virg(@sprintf(" %.2f", vi)), color = ENCRE2, fontsize = 11, align = (:left, :center))
    end
    # tes valeurs mesurées, en face de l'estimateur correspondant (au-dessus de la barre)
    for (mk, k) in (("mesure_tau", "2"), ("mesure_fret", "1inv"))
        i = findfirst(==(mk), kapp.estimateur)
        j = findfirst(==(k), rows)
        (i === nothing || j === nothing) && continue
        scatter!(ax2, [kapp.rapport_K[i]], [y[j] + 0.36]; color = ENCRE, marker = :diamond, markersize = 10)
        text!(ax2, kapp.rapport_K[i], y[j] + 0.36; text = virg(@sprintf("  mesuré %.2f", kapp.rapport_K[i])), color = ENCRE, fontsize = 10, align = (:left, :center))
    end
    if mondes_K !== nothing
        # plage des mondes réajustés, sous la barre
        for (yi, k) in zip(y, rows)
            d = mondes_K[mondes_K.estimateur.==k, :]
            nrow(d) == 0 && continue
            v0 = kapp.rapport_K[findfirst(==(k), kapp.estimateur)]
            lo, hi = min(minimum(d.rapport_K), v0), max(maximum(d.rapport_K), v0)
            rangebars!(ax2, [yi - 0.36], [lo], [hi]; direction = :x, color = MUET, linewidth = 2, whiskerwidth = 5)
        end
        # (c) K apparents absolus : prédits (point nominal, plage des mondes) et mesurés
        ax3 = Axis(fig[1, 3]; xlabel = "K apparent (mM)", yticks = (1:3, ["τ mono + Hill\n(grille de titration)", "J/B", "B/J (ratio TIFF C2/C1)"]),
            title = "K apparents prédits contre mesurés", ygridvisible = false)
        for (yi, k, mk) in ((1, "2H", "mesure_tau"), (2, "1", ""), (3, "1inv", "mesure_fret"))
            d = mondes_K[mondes_K.estimateur.==k, :]
            i0 = findfirst(==(k), kapp.estimateur)
            i0 === nothing && continue
            v0 = kapp.K_app[i0]
            lo, hi = min(minimum(d.K_app), v0), max(maximum(d.K_app), v0)
            rangebars!(ax3, [yi], [lo], [hi]; direction = :x, color = ENCRE2, linewidth = 2, whiskerwidth = 7)
            scatter!(ax3, [v0], [yi]; color = COULEUR["1"], markersize = 11)
            text!(ax3, hi, yi; text = virg(@sprintf("  %.1f (%.1f–%.1f)", v0, lo, hi)), align = (:left, :bottom), fontsize = 10, color = ENCRE2)
            j = findfirst(==(mk), kapp.estimateur)
            if j !== nothing
                scatter!(ax3, [kapp.K_app[j]], [yi - 0.25]; color = ENCRE, marker = :diamond, markersize = 11)
                text!(ax3, kapp.K_app[j], yi - 0.25; text = virg(@sprintf("  mesuré %.1f", kapp.K_app[j])), align = (:left, :center), fontsize = 10, color = ENCRE)
            end
        end
        xlims!(ax3, 10, 60)
        ylims!(ax3, 0.4, 3.6)
    end
    xlims!(ax2, 0, 2.35)
    colsize!(fig.layout, 1, Relative(mondes_K === nothing ? 0.45 : 0.36))
    entete!(fig, "Fig. 5 — K apparent / K vrai par estimateur",
        @sprintf("K vrai = %.1f mM (isotherme sur les x de l'ajustement global). Hill à 4 paramètres ajustée à la courbe sans bruit de chaque estimateur sur 0–50 K. ", K_vrai) *
        "Losanges : tes ajustements de Hill expérimentaux (tests/csv_exports/hill_fit_summary.csv), en face de l'estimateur qui leur correspond." *
        (mondes_K === nothing ? "" : " Barres : plage des mondes réajustés (IRF ×0,9 et ×1,1, t₀ ± 3,5 ps, t₀ par histogramme), chacun avec son propre K vrai."); largeur = mondes_K === nothing ? 2 : 3)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 6 — efficacité selon le budget
# -----------------------------------------------------------------------------
function figure6(df, chemin; C = 20.0, norm = "A", cles = AVEC_2F)
    fig = Figure(size = (1300, 580))
    d = sel(df; C_vrai = C, normalisation = norm, mode_calibration = "exacte")
    xt = ([1e2, 1e3, 1e4, 1e5], ["10²", "10³", "10⁴", "10⁵"])
    ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10, xticks = xt, xlabel = "budget N (photons)",
        ylabel = "efficacité (log)", title = "Par rapport à toute l'information (borne absolue)", yticks = graduations_log(1e-4, 10))
    ax2 = Axis(fig[1, 2]; xscale = log10, xticks = xt, xlabel = "budget N (photons)", ylabel = "efficacité",
        title = "Par rapport à l'information des données qu'il consomme")
    for ax in (ax1, ax2)
        hlines!(ax, [1.0]; color = ENCRE, linewidth = 1)
    end
    for k in cles
        s = sort(sel(d; estimateur = k), :N_photons)
        nrow(s) == 0 && continue
        serie!(ax1, s.N_photons, s.efficacite, k; creux = echecs(s), style = style_de(k))
        serie!(ax2, s.N_photons, min.(s.efficacite_propre, 1.45), k; creux = echecs(s), style = style_de(k))
    end
    ylims!(ax2, 0, 1.5)
    Legend(fig[1, 3], ax1; merge = true, unique = true)
    entete!(fig, "Fig. 6 — Efficacité σ²_CRLB / σ²_MC selon le budget, à $(Int(C)) mM",
        NORMALISATION[norm] * ". 1 = au bout de ce qui est possible ; 0,3 = 70 % de l'information laissée sur la table. À gauche, borne commune " *
        "(deux canaux, rapport J/B connu) ; à droite, borne propre (bleu seul, totaux seuls, ou deux canaux à échelles séparées selon l'estimateur). " *
        "Au-dessus de 1 : estimateur biaisé ou tronqué à [0, 150] mM (valeurs > 1,45 portées à 1,45)." * NOTE_CREUX; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Titration — le test de M2 lent sur tes données
# -----------------------------------------------------------------------------

# Rampe bleue ordinale (palette de référence, pas 250 → 650) : une teinte par concentration.
const RAMPE = hexa.(["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281"])
const CONC = ["0mM", "2.5mM", "5mM", "10mM", "20mM", "40mM", "80mM", "146mM", "F"]
const CONC_LAB = ["0", "2,5", "5", "10", "20", "40", "80", "146", "F"]
position_C(C) = C == 0 ? log10(0.8) : isinf(C) ? log10(400.0) : log10(C)
const TICKS_C = (position_C.([0, 2.5, 5, 10, 20, 40, 80, 146, Inf]), CONC_LAB)

function figure_titration(prot, dossier, chemin)
    t = CSV.read(joinpath(dossier, "titration_$(prot).csv"), DataFrame)
    a = CSV.read(joinpath(dossier, "ajustement_global_$(prot).csv"), DataFrame)
    p = TOML.parsefile(joinpath(dossier, "parametres_$(prot).toml"))
    t.pos = position_C.(t.C_mM)
    t.ic = [findfirst(==(replace(c, " " => "")), CONC) for c in t.concentration]
    fig = Figure(size = (1450, 900))

    # (a) déclins 0 mM et F, modèle global
    ax1 = Axis(fig[1, 1]; yscale = log10, xlabel = "temps (ns)", ylabel = "photons par canal", title = "Ajustement global, réplicat c1",
        yticks = graduations_log(100, 1e5))
    for (c, ic) in (("0mM_c1", 1), ("F_c1", 9))
        scatter!(ax1, a.t_ns, max.(a[!, "donnees_" * c], 0.5); color = (RAMPE[ic], 0.55), markersize = 4, strokewidth = 0)
        lines!(ax1, a.t_ns, a[!, "modele_" * c]; color = RAMPE[ic], linewidth = 2, label = c == "0mM_c1" ? "0 mM" : "F (saturant)")
    end
    axislegend(ax1; position = :rt)
    # (b) résidus normalisés
    ax2 = Axis(fig[2, 1]; xlabel = "temps (ns)", ylabel = "(n − μ)/√μ", title = "Résidus normalisés")
    hlines!(ax2, [0.0]; color = AXE, linewidth = 1)
    for (c, ic) in (("0mM_c1", 1), ("F_c1", 9))
        r = (a[!, "donnees_" * c] .- a[!, "modele_" * c]) ./ sqrt.(a[!, "modele_" * c])
        lines!(ax2, a.t_ns, r; color = RAMPE[ic], linewidth = 1)
    end
    linkxaxes!(ax1, ax2)

    # (c) x selon C et isothermes
    ax3 = Axis(fig[1, 2]; xticks = TICKS_C, xlabel = "[Cl⁻] (mM)", ylabel = "fraction active x", title = "Fraction d'amplitude courte x")
    Cf = exp10.(range(log10(0.8), log10(400); length = 300))
    iso(C, x0, xs, K, n) = xs + (x0 - xs) * K^n / (K^n + C^n)
    lines!(ax3, log10.(Cf), iso.(Cf, p["x0"], p["x_sat"], p["K_D"], p["n"]); color = ENCRE, linewidth = 1.5,
        label = @sprintf("n libre : K = %.1f mM, n = %.2f", p["K_D"], p["n"]))
    q = p["isotherme_n1"]
    lines!(ax3, log10.(Cf), iso.(Cf, q["x0"], q["x_sat"], q["K_D"], 1.0); color = MUET, linewidth = 1.5,
        label = @sprintf("n = 1 : K = %.1f mM", q["K_D"]))
    scatter!(ax3, t.pos, t.x; color = RAMPE[t.ic], markersize = 9)
    axislegend(ax3; position = :rt)

    # (d) χ² mono observé et prédit par M2 lent
    ax4 = Axis(fig[2, 2]; xticks = TICKS_C, yscale = log10, xlabel = "[Cl⁻] (mM)", ylabel = "χ²_r mono-exponentiel",
        yticks = graduations_log(1, 100),
        title = "Ton χ²_r mono, point par point, et ce que M2 lent prédit")
    scatter!(ax4, t.pos .- 0.03, t.chi2_mono; color = ENCRE, marker = :circle, markersize = 8, label = "observé (ton code)")
    scatter!(ax4, t.pos .+ 0.03, t.chi2_mono_predit_M2; color = SURFACE, strokecolor = ENCRE, strokewidth = 1.5,
        marker = :circle, markersize = 8, label = "prédit par M2 lent")
    axislegend(ax4; position = :rt)

    # (e) le test d'intensité : à temps fixé, l'amplitude du donneur à t = 0 reste constante
    #     et les comptages nets montent comme (1 − E·x)/(1 − E·x₀)
    a0 = mean(t.intensite_relative[t.C_mM.==0])
    temps_fixe = "arret_sur_temps" in names(t) && all(t.arret_sur_temps)
    ti = "test_intensite" in keys(p) ? p["test_intensite"] : nothing
    ax5 = Axis(fig[1, 3]; xticks = TICKS_C, xlabel = "[Cl⁻] (mM)", ylabel = "relatif à 0 mM",
        title = ti === nothing ? "Amplitude du donneur à t = 0" :
                virg(@sprintf("Intensité à temps fixé : ×%.2f mesuré, ×%.2f prédit", ti["mesure_sat_sur_0"], ti["predit_sat_sur_0"])))
    hlines!(ax5, [1.0]; color = AXE, linewidth = 1)
    scatter!(ax5, t.pos, t.intensite_relative ./ a0; color = SURFACE, strokecolor = RAMPE[t.ic], strokewidth = 1.8, markersize = 9,
        label = "amplitude du donneur à t = 0 (M2 : constante)")
    if "intensite_nette_sur_0mM" in names(t)
        o = sortperm(t.pos)
        lines!(ax5, t.pos[o], t.intensite_predite_sur_0mM[o]; color = ENCRE, linewidth = 1.5, label = "comptages prédits (1 − E·x)/(1 − E·x₀)")
        scatter!(ax5, t.pos, t.intensite_nette_sur_0mM; color = RAMPE[t.ic], markersize = 9, label = "comptages nets mesurés")
    end
    axislegend(ax5; position = :lt, labelsize = 10)
    ylims!(ax5, 0.6, 2.0)
    temps_fixe && text!(ax5, 0.98, 0.04; space = :relative, align = (:right, :bottom), fontsize = 10, color = MUET,
        text = virg(@sprintf("acquisitions à temps fixé (%g s, arrêt sur temps)", first(t.duree_acquisition_s))))

    # (f) phasors
    ax6 = Axis(fig[2, 3]; xlabel = "g", ylabel = "s", title = "Phasors (fond retiré, IRF corrigée) et segment M2", aspect = DataAspect())
    pD0, pDA = p["phasor_tau_D0"], p["phasor_tau_DA"]
    gmin, gmax = extrema(t.phasor_g)
    smin, smax = extrema(t.phasor_s)
    pad = 0.6 * max(gmax - gmin, smax - smin)
    d = [pDA[1] - pD0[1], pDA[2] - pD0[2]]
    λs = range(-2, 2; length = 3)
    lines!(ax6, pD0[1] .+ λs .* d[1], pD0[2] .+ λs .* d[2]; color = ENCRE, linewidth = 1, label = "droite τ_D0 → τ_DA")
    scatter!(ax6, t.phasor_g, t.phasor_s; color = RAMPE[t.ic], markersize = 9)
    limits!(ax6, gmin - pad, gmax + pad, smin - pad, smax + pad)
    axislegend(ax6; position = :lb)

    elems = [MarkerElement(color = RAMPE[i], marker = :circle, markersize = 9) for i in 1:9]
    Legend(fig[1:2, 4], elems, CONC_LAB .* [" mM", " mM", " mM", " mM", " mM", " mM", " mM", " mM", ""]; title = "[Cl⁻]")
    entete!(fig, "Titration $prot — M2 lent sur tes données",
        @sprintf("Ajustement global : τ_D0 = %.3f ns et τ_DA = %.3f ns liés sur les %d histogrammes, x, intensité et fond libres par histogramme ; χ²_r = %.2f (τ libres par histogramme : %.2f). ",
            p["tau_D0"], p["tau_DA"], nrow(t), p["chi2_r_global"], p["chi2_r_libres_moyen"]) *
        "Canal bleu (470/40), IRF irf_20260710_16x_sclm, fenêtre 14–242."; largeur = 4)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Figure 7 — cible mobile
# -----------------------------------------------------------------------------
function figure_cible_mobile(df, chemin; N = 1e4)
    fig = Figure(size = (1450, 900))
    o = df[(df.boucle.=="ouverte").&(df.N.==N), :]
    for (j, k) in enumerate(("1", "5a", "2H"))
        ax = Axis(fig[1, j]; xlabel = "temps (s)", ylabel = j == 1 ? "[Cl⁻] (mM)" : "", title = "Boucle ouverte — " * NOM_COURT[k])
        scatter!(ax, o.t, o[!, "est_"*k]; color = (COULEUR[k], 0.35), markersize = 4, strokewidth = 0, label = "estimations à 10 Hz")
        lines!(ax, o.t, o.C_vrai; color = ENCRE, linewidth = 1.5, label = "[Cl⁻] vraie")
        k == "1" && lines!(ax, o.t, o.kalman_ratio; color = COULEUR["1"], linewidth = 2, label = "Kalman, R = borne de Cramér-Rao")
        ylims!(ax, -5, 155)
        axislegend(ax; position = :lt)
    end
    f = df[(df.boucle.=="fermée").&(df.N.==N), :]
    ax = Axis(fig[2, 1:2]; xlabel = "temps (s)", ylabel = "[Cl⁻] vraie (mM)", title = "Boucle fermée — consigne 20 mM, chaque courbe commandée par un estimateur")
    hlines!(ax, [20.0]; color = ENCRE, linewidth = 1, label = "consigne")
    for k in ("1", "5d", "4", "5a", "2H")
        d = f[isfinite.(f[!, "est_"*k]), :]
        nrow(d) == 0 && continue
        lines!(ax, d.t, d.C_vrai; color = COULEUR[k], linewidth = 1.5, label = NOM_COURT[k])
    end
    axislegend(ax; position = :rb)
    # erreur statique (120 dernières secondes)
    ks = ["1", "5d", "4", "5a", "2", "2H"]
    ax2 = Axis(fig[2, 3]; xticks = (1:length(ks), ks), ylabel = "erreur statique (mM)", title = "Moyenne de [Cl⁻] vraie − consigne, t ≥ 120 s",
        xgridvisible = false)
    hlines!(ax2, [0.0]; color = AXE, linewidth = 1)
    for (iN, (NN, col)) in enumerate(((1e3, RAMPE[2]), (1e4, RAMPE[8])))
        v = Float64[]
        for k in ks
            d = df[(df.boucle.=="fermée").&(df.N.==NN).&isfinite.(df[!, "est_"*k]).&(df.t.>=120), :]
            push!(v, nrow(d) == 0 ? NaN : mean(d.C_vrai) - 20.0)
        end
        barplot!(ax2, collect(1:length(ks)) .+ (iN - 1.5) * 0.36, v; width = 0.32, color = col, strokewidth = 0,
            label = "N = " * fmtN(NN))
    end
    axislegend(ax2; position = :lb)
    entete!(fig, "Fig. 7 — Cible mobile : échelon de 5 à 25 mM (τ = 20 s) à 10 Hz, puis boucle fermée",
        "Normalisation A (flux de photons fixe ; N = " * fmtN(N) * " photons bleus par échantillon à 20 mM). Boucle fermée : intégrateur PI sur l'estimation, " *
        "actionneur ≥ 0 qui fait entrer du Cl⁻, relaxation vers 5 mM en 20 s. Un intégrateur annule l'erreur sur l'ESTIMATION : " *
        "l'erreur statique sur la vraie [Cl⁻] est l'opposé du biais de l'estimateur à la consigne."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Carte — quel estimateur minimise la RMSE
# -----------------------------------------------------------------------------
const FORME_SEULE = ["2", "2F", "2H", "3", "4", "5a", "5c", "6"]   # n'utilisent pas le rapport J/B

"""
Rangée 1 : estimateur de plus petite RMSE (taux d'échec ≤ seuil) ; ex æquo
quand d'autres sont à moins de `tolerance` (≈ 2 erreurs-types de la RMSE avec
1 000 réalisations). Rangée 2 : le meilleur des estimateurs qui n'utilisent
que la forme des déclins, et combien de fois sa RMSE dépasse la meilleure.
"""
function carte(df, chemin; seuil_echec = 0.05, tolerance = 0.05, cles = vcat(ORDRE, "2F", "2H"))
    fig = Figure(size = (1450, 1080))
    Cs = sort(unique(df.C_vrai))
    Ns = sort(unique(df.N_photons))
    gagnants = Dict{Tuple{String,Float64,Float64},String}()
    rampe = cgrad(RAMPE, scale = :linear)
    for (j, norm) in enumerate(("A", "B", "C"))
        d = sel(df; normalisation = norm, mode_calibration = "exacte")
        axs = [Axis(fig[i, j]; xticks = (1:length(Ns), [fmtN(N) for N in Ns]), yticks = (1:length(Cs), [@sprintf("%g", C) for C in Cs]),
            xlabel = "budget N (photons)", ylabel = j == 1 ? "[Cl⁻] vraie (mM)" : "",
            title = (i == 1 ? "Meilleure RMSE — " : "Forme seule : meilleur et surcoût — ") * "normalisation " * norm,
            xgridvisible = false, ygridvisible = false) for i in 1:2]
        for (iN, N) in enumerate(Ns), (iC, C) in enumerate(Cs)
            s = sel(d; N_photons = N, C_vrai = C)
            s = s[(s.taux_echec.<=seuil_echec).&in.(s.estimateur, Ref(cles)), :]
            nrow(s) == 0 && continue
            rmin = minimum(s.RMSE_mM)
            ex = sort(s.estimateur[s.RMSE_mM.<=(1 + tolerance)*rmin]; by = k -> findfirst(==(k), cles))
            rect = Rect(iN - 0.47, iC - 0.47, 0.94, 0.94)
            if length(ex) == 1
                k = ex[1]
                gagnants[(norm, N, C)] = k
                poly!(axs[1], rect; color = COULEUR[k])
                lum = 0.2126 * COULEUR[k].r + 0.7152 * COULEUR[k].g + 0.0722 * COULEUR[k].b
                text!(axs[1], iN, iC; text = k, align = (:center, :center), fontsize = 10, color = lum > 0.5 ? ENCRE : SURFACE)
            else
                gagnants[(norm, N, C)] = join(ex, "≈")
                poly!(axs[1], rect; color = ZONE)
                text!(axs[1], iN, iC; text = join(ex, "≈"), align = (:center, :center), fontsize = 9, color = ENCRE)
            end
            f = s[in.(s.estimateur, Ref(FORME_SEULE)), :]
            nrow(f) == 0 && continue
            i = argmin(f.RMSE_mM)
            r = f.RMSE_mM[i] / rmin
            c = rampe[clamp(log10(r) / 2, 0, 1)]
            poly!(axs[2], rect; color = c)
            lum = 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b
            text!(axs[2], iN, iC; text = @sprintf("%s\n×%.0f", f.estimateur[i], r), align = (:center, :center), fontsize = 8,
                color = lum > 0.45 ? ENCRE : SURFACE)
        end
        for ax in axs
            xlims!(ax, 0.5, length(Ns) + 0.5)
            ylims!(ax, 0.5, length(Cs) + 0.5)
        end
    end
    presents = Set(k for v in values(gagnants) for k in split(v, "≈") if !occursin("≈", v))
    elems = vcat([PolyElement(color = COULEUR[k]) for k in cles if k in presents], [PolyElement(color = ZONE)])
    labs = vcat([NOM_COURT[k] for k in cles if k in presents], [@sprintf("ex æquo (à %d %% près)", round(Int, 100tolerance))])
    Legend(fig[1, 4], elems, labs; title = "plus petite RMSE", tellheight = false)
    Colorbar(fig[2, 4]; colormap = rampe, limits = (0, 2), ticks = ([0, 1, 2], ["×1", "×10", "×100"]),
        label = "RMSE forme seule / meilleure RMSE", height = Relative(0.7), tellheight = false)
    entete!(fig, "Carte — quel estimateur minimise la RMSE sur [Cl⁻], et ce que coûte de renoncer au rapport d'intensités",
        @sprintf("Parmi les estimateurs dont le taux d'échec est ≤ %d %% ; calibration exacte ; estimations ramenées à [0, 150] mM. ", round(Int, 100seuil_echec)) *
        "En bas, « forme seule » = estimateurs qui n'utilisent pas le rapport J/B (2, 2F, 2H, 3, 4, 5a, 5c, 6), donc insensibles à un photoblanchiment de l'accepteur " *
        "ou à une dérive de gain entre détecteurs. A : temps fixé (N = photons bleus à 20 mM) ; B : photons totaux ; C : photons bleus."; largeur = 4)
    save(chemin, fig; px_per_unit = 2)
    return gagnants
end


# -----------------------------------------------------------------------------
# Batterie H1–H5 — le module d'ingestion retrouve-t-il la vérité ? (changement 1)
# -----------------------------------------------------------------------------

# Familles d'hypothèses : teinte catégorielle (ordre fixe de la palette) et forme fixes.
const FAMILLE = [
    ("H1", "H1 modèle actuel (contrôle)", "#2a78d6", :circle),
    ("H2", "H2/H2b accepteur sombre, complexe statique", "#eb6834", :rect),
    ("H3b", "H3b donneur nu (4,1 ns)", "#1baf7a", :diamond),
    ("H4", "H4/H5 τ_DA = 1,0 ns (IRF juste ou fausse)", "#eda100", :utriangle),
    ("IRF", "H1 ajusté avec une IRF fausse", "#e87ba4", :dtriangle),
    ("H3a", "H3a état FRET intermédiaire", "#008300", :pentagon),
    ("gigue", "gigue de t₀ entre échantillons", "#4a3aa7", :hexagon),
]
function famille(jeu)
    occursin("gigue", jeu) && !startswith(jeu, "H3a") && return "gigue"
    startswith(jeu, "H3a") && return "H3a"
    startswith(jeu, "H3b") && return "H3b"
    startswith(jeu, "H2") && return "H2"
    (startswith(jeu, "H4") || startswith(jeu, "H5")) && return "H4"
    occursin("IRF", jeu) && return "IRF"
    return "H1"
end
couleur_famille(f) = hexa(FAMILLE[findfirst(t -> t[1] == f, FAMILLE)][3])
"""Virgule décimale dans un texte de figure (0.15 → 0,15)."""
virg(s) = replace(string(s), r"(\d)\.(\d)" => s"\1,\2")
marqueur_famille(f) = FAMILLE[findfirst(t -> t[1] == f, FAMILLE)][4]

function figure_batterie(dossier, chemin)
    bat = CSV.read(joinpath(dossier, "batterie.csv"), DataFrame; types = Dict(:jeu => String, :verdict => String))
    reel = CSV.read(joinpath(dossier, "reel_perturbations.csv"), DataFrame)
    prof = CSV.read(joinpath(dossier, "profil_tau_DA.csv"), DataFrame)
    t0 = CSV.read(joinpath(dossier, "t0_par_histogramme.csv"), DataFrame; types = Dict(:jeu => String, :dossier => String))
    res = TOML.parsefile(joinpath(dossier, "batterie_resume.toml"))
    χr = res["chi2_r_reel"]
    moy = combine(groupby(bat[bat.jeu.!="réel", :], :jeu, sort = false), :chi2_r => mean => :chi2_r, :tau_DA => mean => :tau_DA,
        :tau_DA_vrai => first => :tau_DA_vrai, :dtau_mono_0 => mean => :dt0, :dtau_mono_sat => mean => :dts,
        :phasor_mediane => mean => :ph, :biais_ratio_20mM => mean => :b20, :s_irf_ajustement => first => :s)
    moy.fam = famille.(moy.jeu)
    sr = res["signatures_reelles"]
    fig = Figure(size = (1500, 1000))

    # (a) profil de χ²ᵣ en τ_DA sur les données réelles
    ax1 = Axis(fig[1, 1]; xscale = log10, xticks = ([0.05, 0.1, 0.2, 0.5, 1.0, 2.0], ["0,05", "0,1", "0,2", "0,5", "1", "2"]),
        xlabel = "τ_DA imposé (ns)", ylabel = "χ²ᵣ des 27 histogrammes réels", title = "τ_DA tient : les données rejettent τ_DA ≥ 0,3 ns")
    for (sv, col, st, lab) in ((1.0, ENCRE, :solid, "IRF mesurée"), (0.9, MUET, :dash, "IRF × 0,9"), (1.1, MUET, :dot, "IRF × 1,1"))
        d = sort(prof[prof.s_irf.==sv, :], :tau_DA)
        lines!(ax1, d.tau_DA, d.chi2_r; color = col, linestyle = st, linewidth = sv == 1 ? 2 : 1.5, label = lab)
        sv == 1 && scatter!(ax1, d.tau_DA, d.chi2_r; color = col, markersize = 7)
    end
    d1 = prof[(prof.s_irf.==1.0).&(prof.tau_DA.==1.0), :]
    dmin = prof[prof.s_irf.==1.0, :]
    nrow(d1) == 1 && text!(ax1, 0.03, 0.62; space = :relative, align = (:left, :top), fontsize = 11, color = ENCRE2,
        text = virg(@sprintf("τ_DA = 1 ns (hypothèse H4) : χ²ᵣ = %.1f\ncontre %.2f au minimum (τ_DA = %.2f ns)", d1.chi2_r[1],
            minimum(dmin.chi2_r), dmin.tau_DA[argmin(dmin.chi2_r)])))
    axislegend(ax1; position = :lt)

    # (b) τ_DA retrouvé contre vrai : H1, H4 et leurs ajustements avec une IRF fausse
    lignes_b = ["H1", "H1 (IRF ×0.90)", "H1 (IRF ×1.10)", "H4", "H5 (IRF ×0,9)", "H5 (IRF ×1,1)"]
    etiq_b = ["H1", "H1, IRF × 0,9", "H1, IRF × 1,1", "H4", "H5, IRF × 0,9", "H5, IRF × 1,1"]
    ax2 = Axis(fig[1, 2]; xscale = log10, xticks = ([0.1, 0.2, 0.5, 1.0], ["0,1", "0,2", "0,5", "1"]),
        yticks = (1:length(lignes_b), etiq_b), yreversed = true, xlabel = "τ_DA (ns)",
        title = "Une IRF fausse de 10 % ne fabrique pas E = 0,96")
    for (i, j) in enumerate(lignes_b)
        r = moy[moy.jeu.==j, :]
        nrow(r) == 1 || continue
        f = r.fam[1]
        scatter!(ax2, [r.tau_DA_vrai[1]], [i]; color = SURFACE, strokecolor = ENCRE2, strokewidth = 1.5, marker = :circle, markersize = 12)
        scatter!(ax2, [r.tau_DA[1]], [i]; color = couleur_famille(f), marker = marqueur_famille(f), markersize = 12)
        gauche = r.tau_DA[1] < r.tau_DA_vrai[1]
        text!(ax2, r.tau_DA[1], i; text = virg(@sprintf("%.3f", r.tau_DA[1])), align = (gauche ? :right : :left, :center),
            offset = (gauche ? -10 : 10, 0), fontsize = 10, color = ENCRE2)
    end
    xlims!(ax2, 0.07, 2.2)
    ylims!(ax2, length(lignes_b) + 0.6, 0.4)
    axislegend(ax2, [MarkerElement(color = SURFACE, strokecolor = ENCRE2, strokewidth = 1.5, marker = :circle, markersize = 12),
            MarkerElement(color = ENCRE2, marker = :circle, markersize = 12)], ["vrai", "retrouvé"]; position = :rt)

    # (c) ce que vaut χ²ᵣ = 1,96
    jeux_c = filter(j -> !occursin("imposé", j), moy.jeu)
    d = moy[in.(moy.jeu, Ref(jeux_c)), :]
    sort!(d, :chi2_r)
    nom_jeu(j) = j == "H3a" ? "H3a τ=1,5 20 %" : j == "H3b" ? "H3b 10 %" : virg(j)
    ax3 = Axis(fig[1, 3]; yticks = (1:nrow(d), nom_jeu.(d.jeu)), xlabel = "χ²ᵣ de l'ajustement global à deux états", ygridvisible = false,
        title = "Ce que vaut χ²ᵣ = 1,96", yticklabelsize = 10)
    σ = res["sigma_chi2_poisson"]
    vspan!(ax3, 1 - 3σ, 1 + 3σ; color = (ZONE, 1.0))
    vlines!(ax3, [χr]; color = ENCRE, linewidth = 1.5)
    text!(ax3, χr, 0.5; text = virg(@sprintf("données réelles %.2f", χr)), align = (:right, :bottom), offset = (-4, 0), fontsize = 10, color = ENCRE)
    text!(ax3, 1 + 3σ, 0.5; text = "bruit de\nPoisson", align = (:left, :bottom), offset = (4, 0), fontsize = 9, color = MUET)
    for (i, r) in enumerate(eachrow(d))
        scatter!(ax3, [r.chi2_r], [i]; color = couleur_famille(r.fam), marker = marqueur_famille(r.fam), markersize = 10)
    end
    ylims!(ax3, 0.4, nrow(d) + 0.6)
    xlims!(ax3, 0.9, max(2.1, maximum(d.chi2_r) + 0.1))

    # (d) signatures : τ mono et phasors
    ax4 = Axis(fig[2, 1]; xlabel = "τ mono, données − modèle ajusté, à 0 mM (ns)", ylabel = "même écart à saturation (ns)",
        title = "Signature : un état intermédiaire la reproduit, une IRF fausse non")
    hlines!(ax4, [0.0]; color = AXE, linewidth = 1)
    vlines!(ax4, [0.0]; color = AXE, linewidth = 1)
    for r in eachrow(moy)
        occursin("imposé", r.jeu) && continue
        scatter!(ax4, [r.dt0], [r.dts]; color = couleur_famille(r.fam), marker = marqueur_famille(r.fam), markersize = 10)
    end
    for (j, lab, al, off) in (("H3a", "H3a τ=1,5 ns 20 %", (:left, :center), (8, 0)),
                              ("H3a τ=2.5 40 %", "H3a τ=2,5 ns 40 %", (:left, :bottom), (6, 4)),
                              ("H3a 30 % + gigue", "H3a 30 % + gigue de t₀", (:center, :top), (0, -9)),
                              ("H1 (IRF ×1.10)", "IRF × 1,1", (:right, :center), (-8, 0)))
        r = moy[moy.jeu.==j, :]
        nrow(r) == 1 && text!(ax4, r.dt0[1], r.dts[1]; text = lab, align = al, offset = off, fontsize = 10, color = ENCRE2)
    end
    scatter!(ax4, [sr["dtau_mono_0"]], [sr["dtau_mono_sat"]]; color = ENCRE, marker = :star5, markersize = 18)
    text!(ax4, sr["dtau_mono_0"], sr["dtau_mono_sat"]; text = "données réelles", align = (:left, :bottom), offset = (8, 4), fontsize = 11, color = ENCRE)
    xlims!(ax4, -0.13, 0.035)

    # (e) x₀ sous perturbation de l'IRF et de t₀ (données réelles)
    rr = reel[.!occursin.("canal", reel.perturbation), :]
    ax5 = Axis(fig[2, 2]; yticks = (1:nrow(rr), virg.(rr.perturbation)), yreversed = true, xlabel = "Δx₀ (données réelles réajustées)",
        title = "x₀ sous perturbation de l'IRF et de t₀", yticklabelsize = 10, ygridvisible = false)
    vlines!(ax5, [0.0]; color = AXE, linewidth = 1)
    for (i, r) in enumerate(eachrow(rr))
        rej = r.chi2_r > χr + 0.01
        scatter!(ax5, [r.delta_x0], [i]; color = rej ? SURFACE : ENCRE, strokecolor = ENCRE, strokewidth = 1.5, markersize = 10)
        text!(ax5, r.delta_x0, i; text = virg(@sprintf("χ²ᵣ %.2f", r.chi2_r)), align = (:left, :center), offset = (9, 0), fontsize = 9, color = MUET)
    end
    ylims!(ax5, nrow(rr) + 0.6, 0.4)
    xlims!(ax5, min(-0.08, minimum(rr.delta_x0) - 0.02), max(0.12, maximum(rr.delta_x0) + 0.06))
    rc = reel[occursin.("canal", reel.perturbation), :]
    nrow(rc) > 0 && text!(ax5, 0.0, nrow(rr) + 0.5; text = @sprintf("t₀ décalé d'un canal : χ²ᵣ %.0f–%.0f, hors cadre", minimum(rc.chi2_r), maximum(rc.chi2_r)),
        align = (:center, :bottom), fontsize = 9, color = MUET)

    # (f) t₀ libre par histogramme
    ax6 = Axis(fig[2, 3]; xticks = TICKS_C, xlabel = "[Cl⁻] (mM)", ylabel = "t₀ − moyenne (ps)", title = "t₀ libre par histogramme : gigue sans tendance")
    hlines!(ax6, [0.0]; color = AXE, linewidth = 1)
    for (j, col, dx, lab) in (("H1 (contrôle)", MUET, 0.04, "H1 simulé (contrôle)"), ("réel", ENCRE, -0.04, "données réelles"))
        d = t0[t0.jeu.==j, :]
        y = 1000 .* (d.t0_ns .- mean(d.t0_ns))
        scatter!(ax6, position_C.(d.C_mM) .+ dx, y; color = col, markersize = 8, label = lab)
    end
    q = res["t0_par_histogramme"]["reel"]
    text!(ax6, 0.02, 0.98; text = virg(@sprintf("ρ de Spearman %+.2f (p ≈ %.2f)\nentre échantillons %.1f ps, entre réplicats %.1f ps", q["rho_spearman"],
        q["p_permutation"], q["dispersion_inter_ps"], q["dispersion_intra_ps"])), space = :relative, align = (:left, :top), fontsize = 10, color = ENCRE2)
    axislegend(ax6; position = :rb)

    elems = [MarkerElement(color = hexa(c), marker = mk, markersize = 11) for (_, _, c, mk) in FAMILLE]
    Legend(fig[3, 1:3], vcat(elems, [MarkerElement(color = ENCRE, marker = :star5, markersize = 15)]),
        vcat([l for (_, l, _, _) in FAMILLE], ["données réelles"]); orientation = :horizontal, nbanks = 2, tellwidth = false)
    entete!(fig, "Batterie H1–H5 — l'ingestion du canal bleu face à des titrations dont la vérité est connue",
        @sprintf("Mêmes IRF, photons, fonds et concentrations que la titration du 2026-07-20 ; moyenne de 3 tirages pour les jeux principaux. ") *
        "H2 (transfert vers un accepteur non émissif) est indiscernable du modèle actuel par le canal bleu : seul le canal jaune peut le trancher. " *
        "Un état FRET intermédiaire (H3a) reproduit l'écart du τ mono et la distance des phasors ; une IRF fausse ne les reproduit pas."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end


# -----------------------------------------------------------------------------
# Facteur 8 — mono + Hill contre bi-exponentiel à τ figés, sur les données réelles
# -----------------------------------------------------------------------------
const STYLE_SOURCE = Dict("reel" => (:solid, true, "données réelles amincies"), "jumeau" => (:dash, false, "jumeau (tirages du modèle ajusté)"),
    "simulateur" => (:dot, nothing, "simulateur (étape 4, normalisation C)"))

function figure_facteur8(res, chemin; C = 20.0)
    fig = Figure(size = (1350, 960))
    Ns = sort(unique(res.N))
    serie_src!(ax, x, y, k, src; positif = true) = begin
        st, plein, _ = STYLE_SOURCE[src]
        ok = isfinite.(y) .& (positif ? (y .> 0) : trues(length(y)))
        lines!(ax, x[ok], y[ok]; color = COULEUR[k], linestyle = st, linewidth = src == "simulateur" ? 1.5 : 2)
        plein === nothing && return
        scatter!(ax, x[ok], y[ok]; color = plein ? COULEUR[k] : SURFACE, strokecolor = COULEUR[k], strokewidth = 1.8, marker = MARQUEUR[k])
    end
    lab = Dict("4R" => "4", "2H" => "2H")
    # (a) σ selon N à C
    ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10, xticks = (Ns, fmtN.(Ns)), yticks = graduations_log(1, 100),
        xlabel = "photons par histogramme", ylabel = "σ[Cl⁻] (mM)", title = virg(@sprintf("σ[Cl⁻] à %g mM : le facteur se matérialise sur les données", C)))
    for k in ("4R", "2H"), src in ("simulateur", "jumeau", "reel")
        d = sort(res[(res.source.==src).&(res.C_vrai.==C), :], :N)
        serie_src!(ax1, d.N, d[!, "sigma_"*k], lab[k], src)
    end
    # (b) facteur selon N, trois concentrations
    ax2 = Axis(fig[1, 2]; xscale = log10, xticks = (Ns, fmtN.(Ns)), xlabel = "photons par histogramme",
        ylabel = "σ(2H) / σ(4)", title = "Le facteur dépend du budget : il vient de la saturation du mono")
    hlines!(ax2, [1.0]; color = AXE, linewidth = 1)
    for (c, col) in ((5.0, RAMPE[2]), (20.0, RAMPE[5]), (40.0, RAMPE[8])), src in ("simulateur", "jumeau", "reel")
        d = sort(res[(res.source.==src).&(res.C_vrai.==c), :], :N)
        st, plein, _ = STYLE_SOURCE[src]
        ok = isfinite.(d.facteur)
        lines!(ax2, d.N[ok], d.facteur[ok]; color = col, linestyle = st, linewidth = src == "simulateur" ? 1.5 : 2)
        plein === nothing || scatter!(ax2, d.N[ok], d.facteur[ok]; color = plein ? col : SURFACE, strokecolor = col, strokewidth = 1.8)
    end
    axislegend(ax2, [LineElement(color = c, linewidth = 2) for c in (RAMPE[2], RAMPE[5], RAMPE[8])], ["5 mM", "20 mM", "40 mM"];
        position = :lt, title = "[Cl⁻]")
    # (c) réel contre jumeau, toutes les cellules
    ax3 = Axis(fig[2, 1]; xscale = log10, yscale = log10, xticks = graduations_log(1, 100), yticks = graduations_log(1, 100),
        xlabel = "σ[Cl⁻] du jumeau (mM)", ylabel = "σ[Cl⁻] des données réelles (mM)", title = "Réel contre jumeau : 4R sur la diagonale, mono mieux que prévu",
        aspect = DataAspect())
    lines!(ax3, [1.0, 100.0], [1.0, 100.0]; color = AXE, linewidth = 1)
    for k in ("4R", "2H")
        r = res[res.source.=="reel", :]
        j = res[res.source.=="jumeau", :]
        cle = string.(r.C_vrai, "_", r.N)
        jj = Dict(string.(j.C_vrai, "_", j.N) .=> j[!, "sigma_"*k])
        x = [get(jj, c, NaN) for c in cle]
        y = r[!, "sigma_"*k]
        ok = isfinite.(x) .& isfinite.(y) .& (x .> 0) .& (y .> 0) .& (r.hors_plage_2H .< 0.5)
        scatter!(ax3, x[ok], y[ok]; color = COULEUR[lab[k]], marker = MARQUEUR[lab[k]], markersize = 10)
    end
    limits!(ax3, 1.5, 100, 1.5, 100)
    text!(ax3, 0.04, 0.96; space = :relative, align = (:left, :top), fontsize = 10, color = ENCRE2,
        text = "cellules où moins de la moitié des lectures\nmono + Hill saturent à 150 mM")
    # (d) taux de saturation du pipeline mono + Hill
    ax4 = Axis(fig[2, 2]; xscale = log10, xticks = (Ns, fmtN.(Ns)), xlabel = "photons par histogramme",
        ylabel = "lectures ramenées à 150 mM (%)", title = virg(@sprintf("Mono + Hill : estimations saturées à %g mM", C)))
    for k in ("2H", "4R"), src in ("simulateur", "jumeau", "reel")
        d = sort(res[(res.source.==src).&(res.C_vrai.==C), :], :N)
        serie_src!(ax4, d.N, 100 .* d[!, "hors_plage_"*k], lab[k], src; positif = false)
    end
    ylims!(ax4, -3, 100)
    elems = [[LineElement(color = COULEUR["4"], linestyle = :solid), MarkerElement(color = COULEUR["4"], marker = MARQUEUR["4"], markersize = 10)],
        [LineElement(color = COULEUR["2H"], linestyle = :solid), MarkerElement(color = COULEUR["2H"], marker = MARQUEUR["2H"], markersize = 10)],
        [LineElement(color = ENCRE2, linestyle = :solid), MarkerElement(color = ENCRE2, marker = :circle, markersize = 10)],
        [LineElement(color = ENCRE2, linestyle = :dash), MarkerElement(color = SURFACE, strokecolor = ENCRE2, strokewidth = 1.8, marker = :circle, markersize = 10)],
        [LineElement(color = ENCRE2, linestyle = :dot)]]
    Legend(fig[1:2, 3], elems, ["4R bi-exponentiel, τ figés, fond libre", "2H mono + Hill (pipeline actuel)", STYLE_SOURCE["reel"][3],
        STYLE_SOURCE["jumeau"][3], STYLE_SOURCE["simulateur"][3]]; tellheight = false)
    entete!(fig, "Fig. 10 — Le « facteur 8 » sur les données réelles",
        "Les 27 histogrammes de la titration (≈ 10⁶ photons) amincis par tirage binomial, 200 tirages par histogramme et par budget ; σ = écart-type des " *
        "estimations d'un même histogramme (bruit de photons seul), moyenné sur les réplicats. Chaque lecture est calibrée sur les histogrammes complets " *
        "(isotherme pour 4R, Hill pour 2H, comme au banc). Le jumeau refait tout sur des tirages du modèle ajusté."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end


# -----------------------------------------------------------------------------
# Changement 5 — les deux biais exclus : décalage inter-canaux et blanchiment
# -----------------------------------------------------------------------------

"""
Tolérance sur Δt : plus petit |Δt| (de part et d'autre de 0) où |biais| atteint
`seuil` mM, par interpolation linéaire depuis (0, 0) ; Inf si jamais atteint.
"""
function tolerance_decalage(dts, biais; seuil = 1.0)
    tol = Inf
    for signe in (-1, 1)
        idx = sortperm(abs.(dts))
        x, y = [0.0], [0.0]
        for i in idx
            sign(dts[i]) == signe || continue
            push!(x, abs(dts[i]))
            push!(y, abs(biais[i]))
        end
        for j in 2:length(x)
            if y[j] >= seuil
                t = x[j-1] + (seuil - y[j-1]) / (y[j] - y[j-1]) * (x[j] - x[j-1])
                tol = min(tol, t)
                break
            end
        end
    end
    return tol
end

"""
Biais AJOUTÉ par une nuisance : biais(valeur) − biais du monde nominal avec les
mêmes germes (balayage « reference » de l'étape 5), par estimateur, [Cl⁻] et N ;
la référence entre au point `x0` (Δt = 0, β = 0). Le biais propre à N fini (mono
saturé, phasor…) n'est pas causé par la nuisance et se retranche.
"""
function biais_ajoute(nu, param; x0 = 0.0)
    d = sel(nu; parametre_nuisance = param, mode_calibration = "nominale")
    r = sel(nu; parametre_nuisance = "reference", mode_calibration = "nominale")
    nrow(r) == 0 && return nothing
    ref = Dict((row.estimateur, row.C_vrai, row.N_photons) => row.biais_mM for row in eachrow(r))
    out = DataFrame(estimateur = String[], C_vrai = Float64[], N_photons = Float64[], valeur = Float64[], delta = Float64[], echec = Float64[])
    for row in eachrow(d)
        b0 = get(ref, (row.estimateur, row.C_vrai, row.N_photons), NaN)
        push!(out, (row.estimateur, row.C_vrai, row.N_photons, row.valeur_nuisance, row.biais_mM - b0, row.taux_echec))
    end
    for row in eachrow(r)
        push!(out, (row.estimateur, row.C_vrai, row.N_photons, x0, 0.0, row.taux_echec))
    end
    return out
end

function figure_decalage(nu, chemin; N = 1e4, cles = ORDRE, seuil = 1.0)
    d = biais_ajoute(nu, "decalage_jaune_ps")
    (d === nothing || nrow(d) == 0) && return nothing
    d = d[d.N_photons.==N, :]
    fig = Figure(size = (1450, 640))
    ax1 = Axis(fig[1, 1]; xlabel = "retard du détecteur jaune sur le bleu, Δt (ps)", ylabel = "biais ajouté à 20 mM (mM, symlog)",
        yscale = Makie.Symlog10(1.0), yticks = [-100, -10, -1, 0, 1, 10, 100], title = "Biais ajouté par Δt, calibration à Δt = 0")
    hspan!(ax1, -seuil, seuil; color = (ZONE, 1.0))
    hlines!(ax1, [0.0]; color = AXE, linewidth = 1)
    for k in vcat(cles, "2H")
        s_ = sort(d[(d.estimateur.==k).&(d.C_vrai.==20.0), :], :valeur)
        nrow(s_) == 0 && continue
        serie!(ax1, s_.valeur, s_.delta, k; creux = s_.echec .> 0.05, style = style_de(k))
    end
    text!(ax1, 195, -seuil; text = virg(@sprintf("± %.0f mM", seuil)), align = (:right, :top), offset = (0, -3), fontsize = 10, color = MUET)
    ks = [k for k in vcat(cles, "2H") if any(d.estimateur .== k)]
    ax2 = Axis(fig[1, 2]; xscale = log10, xticks = ([5, 10, 20, 50, 100, 200], ["5", "10", "20", "50", "100", "200"]),
        yticks = (1:length(ks), [NOM_COURT[k] for k in ks]), yreversed = true, ygridvisible = false,
        xlabel = virg(@sprintf("|Δt| tolérable pour un biais ajouté ≤ %.0f mM (ps)", seuil)), title = "Tolérance de calibration du banc")
    for (i, k) in enumerate(ks), (Cv, col, dy) in ((5.0, RAMPE[2], -0.2), (20.0, RAMPE[5], 0.0), (40.0, RAMPE[8], 0.2))
        s_ = d[(d.estimateur.==k).&(d.C_vrai.==Cv).&(d.valeur.!=0), :]
        nrow(s_) == 0 && continue
        t = tolerance_decalage(s_.valeur, s_.delta; seuil)
        if isfinite(t)
            scatter!(ax2, [max(t, 4.0)], [i + dy]; color = col, marker = MARQUEUR[k], markersize = 11)
        else
            scatter!(ax2, [200.0], [i + dy]; color = col, marker = :rtriangle, markersize = 11)
        end
    end
    xlims!(ax2, 4, 260)
    axislegend(ax2, [MarkerElement(color = c, marker = :circle, markersize = 10) for c in (RAMPE[2], RAMPE[5], RAMPE[8])],
        ["5 mM", "20 mM", "40 mM"]; position = :rb, title = "[Cl⁻]")
    Legend(fig[1, 3], ax1; merge = true, unique = true)
    entete!(fig, "Fig. 13 — Décalage temporel entre les deux détecteurs : la tolérance à exiger du banc",
        @sprintf("Normalisation A, N = %s, données d'un monde où le détecteur jaune retarde de Δt ; estimateurs calibrés à Δt = 0. ", fmtN(N)) *
        "Biais ajouté = biais à Δt moins biais à Δt = 0, mêmes germes. Le ratio n'utilise que les totaux et les estimateurs à un canal ignorent le jaune : " *
        "insensibles (triangle à droite : tolérance > 200 ps). Les ajustements de forme à deux canaux (5a, 5c) lisent Δt comme de l'amplitude."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

function figure_blanchiment(nu, chemin; C = 20.0, cles = ORDRE)
    fig = Figure(size = (1750, 640))
    axes = Axis[]
    for (j, (p, xl, tt)) in enumerate((("blanchiment_accepteur", "accepteurs blanchis en fin d'acquisition, β_A", "Accepteur blanchi non absorbant (donneur nu)"),
                                      ("blanchiment_accepteur_sombre", "accepteurs blanchis en fin d'acquisition, β_A", "Accepteur blanchi en extincteur sombre"),
                                      ("blanchiment_donneur", "donneurs blanchis en fin d'acquisition, β_D", "Donneur blanchi")))
        d = biais_ajoute(nu, p)
        (d === nothing || nrow(d) == 0) && continue
        ax = Axis(fig[1, j]; xlabel = xl, ylabel = j == 1 ? "biais ajouté à 20 mM (mM, symlog)" : "", yscale = Makie.Symlog10(1.0),
            yticks = [-100, -10, -1, 0, 1, 10, 100], title = tt)
        hlines!(ax, [0.0]; color = AXE, linewidth = 1)
        for k in vcat(cles, "2H")
            s_ = sort(d[(d.estimateur.==k).&(d.C_vrai.==C).&(d.N_photons.==1e4), :], :valeur)
            nrow(s_) == 0 && continue
            serie!(ax, s_.valeur, s_.delta, k; creux = s_.echec .> 0.05, style = style_de(k))
        end
        push!(axes, ax)
    end
    isempty(axes) && return nothing
    linkyaxes!(axes...)
    Legend(fig[1, length(axes)+1], axes[1]; merge = true, unique = true)
    entete!(fig, "Fig. 14 — Photoblanchiment pendant l'acquisition, calibration sans blanchiment",
        virg(@sprintf("%g mM, N = 10⁴, normalisation A (temps d'acquisition du capteur intact) ; 1 000 réalisations, erreur-type ≈ 2 mM pour 2H (σ ≈ 60 mM), < 0,3 mM pour les autres. ", C)) *
        "Décroissance exponentielle des molécules intactes, donneur et accepteur indépendamment ; β = fraction blanchie en fin d'acquisition. " *
        "Biais ajouté = biais moins biais sans blanchiment, mêmes germes. Le produit du blanchiment de l'accepteur décide qui résiste : " *
        "s'il n'absorbe plus, le donneur se déséteint et la forme bouge autant que le ratio ; s'il reste extincteur, seul le jaune change."; largeur = 4)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Changement 4 — fusion ratio + forme (extension de la figure 7)
# -----------------------------------------------------------------------------
const COULEUR_FILTRE = Dict("ratio" => COULEUR["1"], "forme" => COULEUR["5a"], "fusion" => ENCRE)
const NOM_FILTRE = Dict("ratio" => "ratio seul (gain supposé fixe)", "forme" => "forme seule (5a)", "fusion" => "fusion, état ([Cl⁻], log g)")

function figure_fusion(res, traj, saut, chemin)
    fig = Figure(size = (1500, 960))
    ax1 = Axis(fig[1, 1:2]; xlabel = "temps (min)", ylabel = "[Cl⁻] (mM)",
        title = "Une réalisation : dérive de gain τ_g = 5 min, 10⁴ photons bleus par échantillon à 10 Hz")
    lines!(ax1, traj.t ./ 60, traj.C_vrai; color = (MUET, 0.8), linewidth = 4, label = "[Cl⁻] vraie")
    for f in ("forme", "ratio", "fusion")
        lines!(ax1, traj.t ./ 60, traj[!, f]; color = COULEUR_FILTRE[f], linewidth = f == "fusion" ? 1.8 : 1.4, label = NOM_FILTRE[f])
    end
    ax2 = Axis(fig[1, 3]; xlabel = "temps (min)", ylabel = "gain du jaune relatif au bleu", title = "Le gain, estimé par la fusion")
    lines!(ax2, traj.t ./ 60, traj.g_vrai; color = (MUET, 0.8), linewidth = 4, label = "vrai")
    lines!(ax2, traj.t ./ 60, traj.g_fusion; color = ENCRE, linewidth = 1.5, label = "estimé")
    axislegend(ax2; position = :lt)
    caps = unique(res.capteur)
    axes = Axis[]
    for (j, cap) in enumerate(caps)
        d = res[res.capteur.==cap, :]
        τs = sort(unique(d.tau_g_s))
        ax = Axis(fig[2, j]; xscale = log10, yscale = log10, xticks = (τs, [t < 60 ? @sprintf("%.0f s", t) : @sprintf("%.0f min", t / 60) for t in τs]),
            yticks = graduations_log(0.1, 100), xlabel = "échelle de temps de la dérive, τ_g", ylabel = j == 1 ? "RMSE sur [Cl⁻] (mM)" : "",
            title = virg(cap))
        for f in ("ratio", "forme", "fusion"), (N, st) in ((1e4, :solid), (1e3, :dash))
            g = combine(groupby(d[(d.filtre.==f).&(d.N.==N), :], :tau_g_s), :RMSE_mM => mean => :r)
            sort!(g, :tau_g_s)
            lines!(ax, g.tau_g_s, g.r; color = COULEUR_FILTRE[f], linestyle = st, linewidth = 2)
            scatter!(ax, g.tau_g_s, g.r; color = st == :solid ? COULEUR_FILTRE[f] : SURFACE, strokecolor = COULEUR_FILTRE[f], strokewidth = 1.5, markersize = 8)
        end
        j == 1 && axislegend(ax, [LineElement(color = ENCRE2, linestyle = :solid), LineElement(color = ENCRE2, linestyle = :dash)],
            ["10⁴ photons bleus par échantillon", "10³ photons"]; position = :rt, labelsize = 10)
        push!(axes, ax)
    end
    linkyaxes!(axes...)
    ax5 = Axis(fig[2, 3]; xlabel = "temps après le saut (s)", ylabel = "biais sur [Cl⁻] (mM)", title = "Saut de gain de 5 % à 20 mM, 10⁴ photons")
    hlines!(ax5, [0.0]; color = AXE, linewidth = 1)
    retours = String[]
    for (cap, st) in zip(caps, (:solid, :dash))
        d = saut[(saut.capteur.==cap).&(saut.N.==1e4), :]
        t = d.t .- 60.0
        lines!(ax5, t, d.biais_ratio_seul_mM; color = COULEUR_FILTRE["ratio"], linestyle = st, linewidth = 1.5)
        lines!(ax5, t, d.biais_fusion_mM; color = ENCRE, linestyle = st, linewidth = 2)
        τr = d.temps_retour_s[1]
        isfinite(τr) && push!(retours, virg(@sprintf("%.0f s (%s)", τr, startswith(cap, "mCer3") ? "mCer3_Topaz" : "E = 0,43")))
    end
    xlims!(ax5, -20, 180)
    axislegend(ax5, [LineElement(color = ENCRE2, linestyle = :solid), LineElement(color = ENCRE2, linestyle = :dash)],
        ["mCer3_Topaz (E = 0,96)", "capteur hypothétique (E = 0,43)"]; position = :rc, labelsize = 10)
    isempty(retours) || text!(ax5, 0.97, 0.62; space = :relative, align = (:right, :top), fontsize = 10, color = ENCRE2,
        text = "fusion : retour à 1/e du pic en\n" * join(retours, ", "))
    elems = vcat([[LineElement(color = (MUET, 0.8), linewidth = 4)]], [[LineElement(color = COULEUR_FILTRE[f], linewidth = 2)] for f in ("ratio", "forme", "fusion")])
    Legend(fig[3, 1:3], elems, vcat(["[Cl⁻] vraie / gain vrai"], [NOM_FILTRE[f] for f in ("ratio", "forme", "fusion")]);
        orientation = :horizontal, tellwidth = false)
    entete!(fig, "Fig. 7b — Fusion du ratio et de la forme : le gain des détecteurs dans l'état du filtre de Kalman",
        "[Cl⁻] = 20 + 8·sin(2πt/5 min), 10 Hz, 30 min ; gain du jaune en marche aléatoire sur log g, 5 % d'écart-type après τ_g. RMSE sur la seconde moitié, " *
        "moyenne de 6 réalisations. Le filtre connaît la statistique de la dérive (q_ℓ = 0,05²/τ_g par seconde) et la variance de Fisher de chaque lecture."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

# -----------------------------------------------------------------------------
# Changement 3 — le canal jaune : ce qu'on en tire déjà, ce que la mesure tranchera
# -----------------------------------------------------------------------------
function figure_canal_jaune(dossier, chemin; K_mes = NaN, plage_mes = NaN)
    courbes = CSV.read(joinpath(dossier, "eta_YD_courbes.csv"), DataFrame; types = Dict(:monde => String))
    inv_ = CSV.read(joinpath(dossier, "eta_YD_inversion.csv"), DataFrame; types = Dict(:monde => String))
    pred = CSV.read(joinpath(dossier, "predictions_JB.csv"), DataFrame; types = Dict(:scenario => String))
    mt = CSV.read(joinpath(dossier, "montee.csv"), DataFrame)
    ar = CSV.read(joinpath(dossier, "aller_retour.csv"), DataFrame; types = Dict(:protocole => String, :parametre => String))
    fig = Figure(size = (1500, 960))
    # (a), (b) inversions de η_YD
    for (j, (col, cible, yl, tt)) in enumerate(((:K_app, K_mes, "K apparent de B/J (mM)", "η_YD par le K apparent du ratio TIFF"),
                                                (:plage, plage_mes, "plage dynamique de B/J (×)", "η_YD par la plage dynamique du ratio TIFF")))
        ax = Axis(fig[1, j]; xlabel = "η_YD / valeur des spectres", ylabel = yl, title = tt)
        for mo in unique(courbes.monde)
            d = courbes[courbes.monde.==mo, :]
            lines!(ax, d.facteur, d[!, col]; color = mo == "nominal" ? ENCRE : MUET, linewidth = mo == "nominal" ? 2 : 1.2,
                linestyle = mo == "nominal" ? :solid : :dash, label = mo == "nominal" ? "modèle nominal" : "mondes IRF / t₀")
        end
        hlines!(ax, [cible]; color = COULEUR["1"], linewidth = 1.5, label = "mesuré (TIFF)")
        f = inv_[inv_.monde.=="nominal", col == :K_app ? :facteur_par_K : :facteur_par_plage][1]
        vlines!(ax, [1.0]; color = AXE, linewidth = 1)
        isfinite(f) && scatter!(ax, [f], [cible]; color = COULEUR["1"], markersize = 11)
        lo, hi = extrema(filter(isfinite, inv_[!, col == :K_app ? :facteur_par_K : :facteur_par_plage]))
        isfinite(f) && text!(ax, f, cible; text = virg(@sprintf("η_YD × %.2f\n(mondes IRF/t₀ : %.2f–%.2f)", f, lo, hi)), align = (:left, :bottom), offset = (8, 8), fontsize = 11, color = ENCRE2)
        axislegend(ax; position = col == :K_app ? :rt : :rt, merge = true, unique = true)
    end
    # (c) J/B prédit selon [Cl⁻]
    ax3 = Axis(fig[1, 3]; xticks = TICKS_C, xlabel = "[Cl⁻] (mM)", ylabel = "J/B prédit (gains égaux)", title = "Ce que la mesure confrontera : J/B et sa plage")
    garde(C) = C == 0 || isinf(C) || 0.8 <= C <= 300
    sup = pred[startswith.(pred.scenario, "supposés").&garde.(pred.C_mM), :]
    env = combine(groupby(sup, :C_mM), :JB => minimum => :lo, :JB => maximum => :hi)
    sort!(env, :C_mM)
    band!(ax3, position_C.(env.C_mM), env.lo, env.hi; color = (ZONE, 1.0))
    text!(ax3, position_C(1.0), env.hi[findfirst(>=(1.0), env.C_mM)]; text = "paramètres supposés :\nη_YD ±50 %, Φ_A ±30 %, τ_A ±20 %, ε",
        align = (:left, :top), offset = (4, -4), fontsize = 10, color = MUET)
    for (sc, col, st, lab) in (("nominal", ENCRE, :solid, "modèle nominal"), ("H2 accepteur non émissif", hexa("#eb6834"), :solid, "H2 : accepteur non émissif"),
                               ("H3a τ=1,5 ns 20 %", hexa("#008300"), :dash, "H3a : état intermédiaire"))
        d = sort(pred[(pred.scenario.==sc).&garde.(pred.C_mM), :], :C_mM)
        lines!(ax3, position_C.(d.C_mM), d.JB; color = col, linestyle = st, linewidth = 2, label = lab)
    end
    axislegend(ax3; position = :rt, labelsize = 11)
    # (d) montée observable ?
    ax4 = Axis(fig[2, 1]; xscale = log10, yscale = log10, xticks = ([0.05, 0.1, 0.2, 0.5, 1, 2], ["0,05", "0,1", "0,2", "0,5", "1", "2"]),
        yticks = graduations_log(0.001, 10), xlabel = "τ_DA vrai (ns)", ylabel = "σ(τ_DA) tiré du jaune seul (ns)", title = "La montée est-elle observable ? (0 mM)")
    τg = [0.04, 2.5]
    band!(ax4, τg, τg ./ 3, fill(10.0, 2); color = (ZONE, 1.0))
    lines!(ax4, τg, τg ./ 3; color = MUET, linewidth = 1)
    text!(ax4, 1.2, 3.5; text = "σ > τ_DA/3 :\nmontée non résolue", align = (:right, :top), fontsize = 10, color = MUET)
    for (N, col) in ((1e5, RAMPE[3]), (1e6, RAMPE[8])), (dl, st) in ((false, :solid), (true, :dash))
        d = sort(mt[(mt.N_jaune.==N).&(mt.dt_libre.==dl), :], :tau_DA)
        lab = (N == 1e5 ? "10⁵" : "10⁶") * " photons jaunes, Δt " * (dl ? "inconnu" : "mesuré au diffuseur")
        lines!(ax4, d.tau_DA, d.sigma_tau_DA; color = col, linewidth = 2, linestyle = st, label = lab)
        scatter!(ax4, d.tau_DA, d.sigma_tau_DA; color = dl ? SURFACE : col, strokecolor = col, strokewidth = 1.5, markersize = 8)
    end
    vlines!(ax4, [0.1515]; color = ENCRE, linewidth = 1)
    text!(ax4, 0.1515, 0.0032; text = " titration", align = (:left, :bottom), fontsize = 10, color = ENCRE2)
    ylims!(ax4, 0.003, 5)
    axislegend(ax4, [[LineElement(color = RAMPE[3], linewidth = 2)], [LineElement(color = RAMPE[8], linewidth = 2)],
            [LineElement(color = ENCRE2, linestyle = :solid)], [LineElement(color = ENCRE2, linestyle = :dash)]],
        ["10⁵ photons jaunes", "10⁶ photons jaunes", "Δt mesuré au diffuseur", "Δt inconnu"]; position = :lt, labelsize = 10)
    # (e) aller-retour
    params = ["ρ_D", "ρ_A", "τ_A", "ε", "Δt"]
    ax5 = Axis(fig[2, 2]; yticks = (1:length(params), params), yreversed = true, xlabel = "écart normé (estimé − vrai)/σ de Fisher",
        title = "Aller-retour sur 20 bancs simulés", ygridvisible = false)
    vspan!(ax5, -1, 1; color = (ZONE, 1.0))
    vlines!(ax5, [0.0]; color = AXE, linewidth = 1)
    for (pr, col, dy) in (("deux temps", MUET, 0.15), ("conjoint", ENCRE, -0.15))
        for (i, q) in enumerate(params)
            d = ar[(ar.protocole.==pr).&(ar.parametre.==q), :]
            nrow(d) == 0 && continue
            μ, σ = mean(d.ecart_norme), std(d.ecart_norme)
            lines!(ax5, [clamp(μ - σ, -14, 14), clamp(μ + σ, -14, 14)], [i + dy, i + dy]; color = col, linewidth = 2)
            scatter!(ax5, [clamp(μ, -14, 14)], [i + dy]; color = col, markersize = 9, label = pr == "conjoint" ? "conjoint (recommandé)" : "en deux temps")
        end
    end
    xlims!(ax5, -15, 15)
    axislegend(ax5; position = :rc, merge = true, unique = true)
    # (f) le protocole, en clair
    conj = ar[ar.protocole.=="conjoint", :]
    σrel(q) = mean(conj.sigma[conj.parametre.==q]) / mean(conj.vrai[conj.parametre.==q])
    txt = "Protocole : 3 concentrations (0, ≈ 30 mM, saturant),\ntemps d'acquisition fixé, les deux canaux ensemble.\n\n" *
        "Ajustement conjoint, τ_D0 et τ_DA de la titration :\n" *
        virg(@sprintf("  ρ_D = η_YD/η_BD   à ± %.1f %%\n", 100σrel("ρ_D"))) *
        virg(@sprintf("  ρ_A (→ Φ_A effectif)   à ± %.1f %%\n", 100σrel("ρ_A"))) *
        virg(@sprintf("  τ_A   à ± %.3f ns\n", mean(conj.sigma[conj.parametre.=="τ_A"]))) *
        virg(@sprintf("  Δt du détecteur jaune   à ± %.1f ps\n", 1000mean(conj.sigma[conj.parametre.=="Δt"]))) *
        virg(@sprintf("  ε   à ± %.3f seulement : borné\n\n", mean(conj.sigma[conj.parametre.=="ε"]))) *
        "Tests : la plage de J/B (sans gain) tranche\nentre M2 et H2 ; à saturation, le jaune doit\ndécroître en τ_D0, pas en τ_A."
    Label(fig[2, 3], txt; halign = :left, valign = :top, justification = :left, fontsize = 12, color = ENCRE2, tellwidth = false, tellheight = false)
    entete!(fig, "Fig. 12 — Le canal jaune : ce que les images déjà acquises en disent, et ce qu'une demi-journée de banc tranchera",
        "En haut à gauche : η_YD est le seul paramètre purement supposé du K apparent et de la plage du ratio ; les deux inversions, indépendantes, " *
        "le placent sous la valeur des spectres. Les images TIFF ne sont pas photoniques : un décalage d'intensité déplacerait les deux nombres."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end


# -----------------------------------------------------------------------------
# Figure 11 — incertitude des paramètres : variance paramétrique et classement
# -----------------------------------------------------------------------------
function figure11(df, nu, chemin; N = 1e4, norm = "A", C = 20.0)
    fig = Figure(size = (1600, 900))
    d = sel(df; N_photons = N, normalisation = norm, mode_calibration = "exacte")
    ax1 = Axis(fig[1, 1:2]; yscale = log10, yticks = graduations_log(0.1, 100), xlabel = "[Cl⁻] vraie (mM)", ylabel = "σ[Cl⁻] (mM)",
        title = "Bruit de photons (trait plein) et incertitude des paramètres ajustés (tireté)")
    for k in ("1", "5a", "4")
        s_ = sort(sel(d; estimateur = k), :C_vrai)
        ok = s_.C_vrai .> 0
        lines!(ax1, s_.C_vrai[ok], s_.sigma_mM[ok]; color = COULEUR[k], linewidth = 2, label = NOM_COURT[k] * " : photons")
        scatter!(ax1, s_.C_vrai[ok], s_.sigma_mM[ok]; color = COULEUR[k], marker = MARQUEUR[k])
        lines!(ax1, s_.C_vrai[ok], sqrt.(s_.variance_parametrique_mM2[ok]); color = COULEUR[k], linewidth = 2, linestyle = :dash)
    end
    text!(ax1, 0.98, 0.04; space = :relative, align = (:right, :bottom), fontsize = 10, color = ENCRE2,
        text = "tireté : 200 vecteurs tirés des covariances (τ de l'ajustement global ×χ²ᵣ, isotherme),\nlus sans bruit par chaque estimateur calibré sur le point estimé")
    axislegend(ax1; position = :lt)
    panneaux = (("tau_A_relatif", "τ_A vrai / supposé"), ("Phi_A_relatif", "Φ_A effectif vrai / supposé"), ("eta_YD_relatif", "η_YD vrai / supposé"))
    cles = ["1", "5d", "5a", "5c", "4", "2H"]
    axes = Axis[]
    for (j, (p, xl)) in enumerate(panneaux), (i, mode) in enumerate(("exacte", "nominale"))
        s0 = sel(nu; parametre_nuisance = p, C_vrai = C, N_photons = N, mode_calibration = mode)
        nrow(s0) == 0 && continue
        ax = Axis(fig[i, j+2]; xlabel = xl, ylabel = j == 1 ? (mode == "exacte" ? "RMSE (mM), calibration exacte" : "RMSE (mM), calibration supposée") : "",
            yscale = log10, yticks = graduations_log(0.5, 100), title = mode == "exacte" ? "Précision si la valeur était connue" : "Erreur si l'on garde la valeur supposée")
        vlines!(ax, [1.0]; color = AXE, linewidth = 1)
        for k in cles
            kk = mode == "nominale" && k == "5a" && startswith(p, "tau") ? "5b" : k
            s_ = sort(sel(s0; estimateur = kk), :valeur_nuisance)
            nrow(s_) == 0 && continue
            r0 = sel(df; N_photons = N, normalisation = norm, mode_calibration = "exacte", C_vrai = C, estimateur = k)
            x = vcat(s_.valeur_nuisance, nrow(r0) == 1 ? 1.0 : Float64[])
            y = vcat(s_.RMSE_mM, nrow(r0) == 1 ? r0.RMSE_mM[1] : Float64[])
            o = sortperm(x)
            serie!(ax, x[o], y[o], k; style = style_de(k))
        end
        push!(axes, ax)
    end
    length(axes) > 1 && linkyaxes!(axes...)
    Legend(fig[2, 1:2], axes[1]; merge = true, unique = true, nbanks = 2, tellheight = false)
    entete!(fig, "Fig. 11 — Incertitude des paramètres d'entrée : ce qu'elle ajoute à σ, et si le classement y survit",
        virg(@sprintf("Normalisation %s, N = %s. À gauche, selon [Cl⁻] : la variance paramétrique s'ajoute à la variance photonique (colonnes séparées dans resultats_principaux.csv). ", norm, fmtN(N))) *
        virg(@sprintf("À droite, à %g mM : τ_A ±20 %%, Φ_A ±30 %%, η_YD ±50 %%. En haut, estimateurs construits sur la vraie valeur ; en bas, sur la valeur supposée.", C)); largeur = 5)
    save(chemin, fig; px_per_unit = 2)
end


# -----------------------------------------------------------------------------
# Banc à deux canaux du 2026-10-16 (autre version de SCLM, traitée seule)
# -----------------------------------------------------------------------------

"""Lecture tolérante d'un CSV : en-têtes et valeurs éventuellement alignés par un éditeur (espaces)."""
function lire_csv(f)
    df = CSV.read(f, DataFrame; stripwhitespace = true)
    rename!(df, [strip(n) for n in names(df)])
    return df
end

const TECH_BANC = ["1", "5d", "5a", "6", "2J", "4", "2", "3"]
const NOM_BANC = Dict("1" => "1 rapport J/B", "2" => "2 mono, bleu (ton code)", "2J" => "2J mono, jaune (ton code)",
    "3" => "3 bi-exp. τ libres (bleu)", "4" => "4 bi-exp. τ figés (bleu)", "5a" => "5a forme 2 canaux, échelles séparées",
    "5d" => "5d forme 2 canaux, échelle commune", "6" => "6 phasor (bleu)")
couleur_banc(k) = k == "2J" ? COULEUR["2"] : COULEUR[k]
marqueur_banc(k) = k == "2J" ? :star5 : MARQUEUR[k]
style_banc(k) = k == "2J" ? :dash : :solid

function figure_banc_declins(dossier, dossier_donnees, chemin; λ = 800)
    mexp = lire_csv(joinpath(dossier, "sommes_trois_exponentielles.csv"))
    irf = lire_csv(joinpath(dossier, "irf.csv"))
    acq = lire_csv(joinpath(dossier, "acquisitions.csv"))
    ph = lire_csv(joinpath(dossier, "sommes_phasors.csv"))
    S_ = Main.SimFLIMFRET
    hs = S_.lire_banc(dossier_donnees)
    pas = 0.06400390634153652
    W = S_.FENETRE_BANC
    fig = Figure(size = (1500, 960))
    cols = Dict(0.0 => RAMPE[2], 131.0 => RAMPE[8])
    for (j, (voie, nomv)) in enumerate(((1, "bleu"), (2, "jaune")))
        ax = Axis(fig[1, j]; yscale = log10, yticks = graduations_log(1, 1e6), xlabel = "temps depuis le pic de l'IRF (ns)",
            ylabel = j == 1 ? "photons par canal (somme des 10 acquisitions)" : "",
            title = virg(@sprintf("Canal %s, %d nm : sommes et trois exponentielles", nomv, λ)))
        irfv, _, _ = S_.irf_banc(joinpath(dossier_donnees, "IRF_$(λ)_16X.sdt"), voie)
        inst = S_.instrument_banc(irfv, pas)
        H = [sum(voie == 1 ? h.bleu : h.jaune for h in hs if h.C == C && h.λ == λ) for C in (0.0, 131.0)]
        r = S_.ajuster_multiexp(H, inst; τ_init = (0.1, 0.9, 3.0), montee = voie == 2)
        t = ((W) .- argmax(irfv)) .* pas
        for (s, C) in enumerate((0.0, 131.0))
            scatter!(ax, t, max.(H[s][W], 0.5); color = (cols[C], 0.45), markersize = 4, strokewidth = 0)
            lines!(ax, t, r.μ[s]; color = cols[C], linewidth = 2, label = virg(@sprintf("%g mM", C)))
        end
        lines!(ax, t, max.(irfv[W] .* maximum(H[1][W]) ./ maximum(irfv), 0.5); color = MUET, linewidth = 1, label = "IRF (même pic)")
        ylims!(ax, 10, 3 * maximum(maximum.(H)))
        axislegend(ax; position = :rt)
    end
    # fractions d'intensité des trois composantes
    ax3 = Axis(fig[1, 3]; xticks = (1:4, ["bleu 0 mM", "bleu 131 mM", "jaune 0 mM", "jaune 131 mM"]), ylabel = "fraction de l'intensité",
        title = virg(@sprintf("Trois composantes, %d nm", λ)), xgridvisible = false)
    d = mexp[mexp.lambda_nm.==λ, :]
    vals = [d[(d.canal.==c).&(d.C_mM.==C), [:fraction_1, :fraction_2, :fraction_3]] for (c, C) in (("bleu", 0.0), ("bleu", 131.0), ("jaune", 0.0), ("jaune", 131.0))]
    τs = [d[d.canal.==c, [:tau_1, :tau_2, :tau_3]][1, :] for c in ("bleu", "jaune")]
    tons = (RAMPE[2], RAMPE[5], RAMPE[9])
    for (i, v) in enumerate(vals)
        base = 0.0
        for c in 1:3
            f = max(v[1, c], 0.0)
            poly!(ax3, Rect(i - 0.35, base, 0.7, f); color = tons[c], strokecolor = SURFACE, strokewidth = 1)
            f > 0.06 && text!(ax3, i, base + f / 2; text = virg(@sprintf("%.2f ns", (i <= 2 ? τs[1] : τs[2])[c])), align = (:center, :center),
                fontsize = 10, color = c == 1 ? ENCRE : SURFACE)
            base += f
        end
    end
    ylims!(ax3, 0, 1.02)
    # IRF des deux détecteurs
    ax4 = Axis(fig[2, 1]; yscale = log10, yticks = graduations_log(1e-4, 1), xlabel = "temps depuis le pic (ns)", ylabel = "IRF normalisée au pic",
        title = "IRF des deux détecteurs : queue lente, retard du jaune")
    for (λx, st) in ((750, :dot), (800, :solid), (850, :dash)), (voie, col) in ((1, RAMPE[5]), (2, COULEUR["4"]))
        irfv, _, _ = S_.irf_banc(joinpath(dossier_donnees, "IRF_$(λx)_16X.sdt"), voie)
        k = argmax(irfv)
        lines!(ax4, ((W) .- k) .* pas, max.(irfv[W] ./ irfv[k], 1e-5); color = col, linestyle = st, linewidth = λx == 800 ? 2 : 1.2,
            label = (voie == 1 ? "bleu " : "jaune ") * string(λx) * " nm")
    end
    xlims!(ax4, -1, 8)
    ylims!(ax4, 1e-4, 1.5)
    r8 = irf[(irf.lambda_nm.==800), :]
    text!(ax4, 0.97, 0.97; space = :relative, align = (:right, :top), fontsize = 11, color = ENCRE2,
        text = virg(@sprintf("800 nm : largeur %.0f ps (bleu), %.0f ps (jaune)\nqueue au-delà de 1 ns : %.0f %% (bleu), %.0f %% (jaune)\nretard du jaune : %+.0f ps",
            r8.largeur_mi_hauteur_ps[1], r8.largeur_mi_hauteur_ps[2], 100r8.queue_apres_1ns[1], 100r8.queue_apres_1ns[2], r8.retard_sur_bleu_ps[2])))
    axislegend(ax4; position = :rc, labelsize = 10, nbanks = 2)
    # stabilité d'une acquisition à l'autre
    g5 = fig[2, 2] = GridLayout()
    ax5 = Axis(g5[1, 1]; xlabel = "acquisition (1 s chacune)", ylabel = "écart à la moyenne de l'échantillon (%)",
        title = "J/B d'une acquisition à l'autre, contre les intensités")
    hspan!(ax5, -0.4, 0.4; color = (ZONE, 1.0))
    for (C, col) in ((0.0, cols[0.0]), (131.0, cols[131.0]))
        g = [h for h in hs if h.C == C && h.λ == λ]
        sort!(g; by = h -> h.cycle)
        JB = [sum(h.jaune[W]) / sum(h.bleu[W]) for h in g]
        NB = [sum(h.bleu[W]) for h in g]
        lines!(ax5, 1:length(g), 100 .* (JB ./ mean(JB) .- 1); color = col, linewidth = 2, label = virg(@sprintf("J/B, %g mM", C)))
        scatter!(ax5, 1:length(g), 100 .* (JB ./ mean(JB) .- 1); color = col, markersize = 8)
        lines!(ax5, 1:length(g), 100 .* (NB ./ mean(NB) .- 1); color = col, linewidth = 1.2, linestyle = :dash, label = virg(@sprintf("photons bleus, %g mM", C)))
    end
    a8 = acq[acq.lambda_nm.==λ, :]
    text!(ax5, 0.03, 0.04; space = :relative, align = (:left, :bottom), fontsize = 10, color = ENCRE2,
        text = virg(@sprintf("σ(J/B) = %.1f et %.1f × Poisson ; σ(photons bleus) = %.1f et %.1f × Poisson",
            a8.sigma_JB[1] / a8.sigma_poisson_JB[1], a8.sigma_JB[2] / a8.sigma_poisson_JB[2],
            a8.sigma_N_bleu[1] / a8.sigma_poisson_bleu[1], a8.sigma_N_bleu[2] / a8.sigma_poisson_bleu[2])))
    Legend(g5[2, 1], ax5; orientation = :horizontal, nbanks = 2, labelsize = 10, tellwidth = false)
    rowgap!(g5, 4)
    # phasors
    ax6 = Axis(fig[2, 3]; xlabel = "g", ylabel = "s", title = "Phasors des sommes (IRF de chaque détecteur corrigée)", aspect = DataAspect())
    θc = range(0, π; length = 200)
    lines!(ax6, 0.5 .+ 0.5 .* cos.(θc), 0.5 .* sin.(θc); color = AXE, linewidth = 1)
    for (c, mk) in (("bleu", :circle), ("jaune", :diamond)), λx in (750, 800, 850)
        d = ph[(ph.canal.==c).&(ph.lambda_nm.==λx), :]
        lines!(ax6, d.g, d.s; color = MUET, linewidth = 1)
        for r in eachrow(d)
            scatter!(ax6, [r.g], [r.s]; color = cols[r.C_mM], marker = mk, markersize = λx == 800 ? 12 : 8)
        end
    end
    for (c, x, y) in (("bleu", 0.73, 0.33), ("jaune", 0.25, 0.40))
        text!(ax6, x, y; text = c, fontsize = 11, color = ENCRE2, align = (:center, :top))
    end
    limits!(ax6, 0, 1, 0, 0.55)
    entete!(fig, "Fig. 15 — Banc du 2026-10-16 (autre version de SCLM) : sommes des 10 acquisitions",
        "Deux échantillons (0 et 131 mM de Cl⁻), excitation à 750, 800 et 850 nm, 10 acquisitions de 1 s, deux canaux ensemble (voie 1 bleu, voie 2 jaune), " *
        "TAC de 16,4 ns : canaux de 64 ps, fenêtre 14–195. Les durées de vie sont communes aux deux échantillons ; les fractions d'intensité changent."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

function figure_banc_stabilite(dossier, chemin; λ = 800)
    st = lire_csv(joinpath(dossier, "stabilite.csv"))
    st.technique = string.(st.technique)
    fig = Figure(size = (1550, 1000))
    Ns = sort(unique(st.N_bleu))
    tN = (Ns, fmtN.(Ns))
    d = st[st.lambda_nm.==λ, :]
    ax1 = Axis(fig[1, 1:2]; xscale = log10, yscale = log10, xticks = tN, yticks = graduations_log(0.001, 10),
        xlabel = "photons bleus par lecture", ylabel = "σ, en fraction de l'écart 0 → 131 mM",
        title = virg(@sprintf("Bruit de chaque lecture, %d nm : sommes réelles amincies (plein), jumeau (tireté), borne de Fisher (fin)", λ)))
    for k in TECH_BANC
        r = sort(d[(d.source.=="reel").&(d.technique.==k), :], :N_bleu)
        j = sort(d[(d.source.=="jumeau").&(d.technique.==k), :], :N_bleu)
        f = sort(d[(d.source.=="Fisher").&(d.technique.==k), :], :N_bleu)
        lines!(ax1, r.N_bleu, r.sigma; color = couleur_banc(k), linestyle = style_banc(k), linewidth = 2.2, label = NOM_BANC[k])
        scatter!(ax1, r.N_bleu, r.sigma; color = r.echec .> 0.05, colormap = [couleur_banc(k), SURFACE], strokecolor = couleur_banc(k),
            strokewidth = 1.5, marker = marqueur_banc(k), markersize = 10)
        lines!(ax1, j.N_bleu, j.sigma; color = (couleur_banc(k), 0.55), linestyle = :dash, linewidth = 1.2)
        nrow(f) > 0 && lines!(ax1, f.N_bleu, f.sigma; color = ENCRE, linewidth = 0.8)
    end
    lines!(ax1, [3e3, 1e5], [1.6, 1.6 / sqrt(1e5 / 3e3)]; color = MUET, linewidth = 1)
    text!(ax1, 3e3, 1.6; text = "pente −½", align = (:right, :center), offset = (-4, 0), fontsize = 10, color = MUET)
    Legend(fig[1, 3], ax1; merge = true, unique = true, title = "technique")
    # biais
    ax2 = Axis(fig[2, 1]; xscale = log10, xticks = tN, yscale = Makie.Symlog10(0.01),
        yticks = ([-1, -0.1, -0.01, 0, 0.01, 0.1, 1], ["−1", "−0,1", "−0,01", "0", "0,01", "0,1", "1"]),
        xlabel = "photons bleus par lecture", ylabel = "biais moyen (fraction de l'écart)", title = "Biais quand les photons manquent")
    hspan!(ax2, -0.02, 0.02; color = (ZONE, 1.0))
    hlines!(ax2, [0.0]; color = AXE, linewidth = 1)
    for k in TECH_BANC
        r = sort(d[(d.source.=="reel").&(d.technique.==k), :], :N_bleu)
        b = (r.biais_0mM .+ r.biais_131mM) ./ 2
        lines!(ax2, r.N_bleu, b; color = couleur_banc(k), linestyle = style_banc(k), linewidth = 2)
        scatter!(ax2, r.N_bleu, b; color = couleur_banc(k), marker = marqueur_banc(k), markersize = 9)
    end
    # échecs
    ax3 = Axis(fig[2, 2]; xscale = log10, xticks = tN, xlabel = "photons bleus par lecture", ylabel = "lectures en échec (%)",
        title = "Ajustements qui échouent")
    for k in TECH_BANC
        r = sort(d[(d.source.=="reel").&(d.technique.==k), :], :N_bleu)
        lines!(ax3, r.N_bleu, 100 .* r.echec; color = couleur_banc(k), linestyle = style_banc(k), linewidth = 2)
        scatter!(ax3, r.N_bleu, 100 .* r.echec; color = couleur_banc(k), marker = marqueur_banc(k), markersize = 9)
    end
    ylims!(ax3, -1, max(5, 100 * maximum(filter(isfinite, d.echec)) + 2))
    # discriminabilité à 10⁴ photons, trois longueurs d'onde
    ax4 = Axis(fig[2, 3]; xscale = log10, xticks = ([0.5, 1, 3, 10, 30, 100], ["0,5", "1", "3", "10", "30", "100"]), yticks = (1:length(TECH_BANC), [NOM_BANC[k] for k in TECH_BANC]),
        yreversed = true, xlabel = "discriminabilité d' à 10⁴ photons bleus", title = "Selon la longueur d'onde (réel)", ygridvisible = false,
        yticklabelsize = 11)
    for (i, k) in enumerate(TECH_BANC), (λx, col, dy) in ((750, RAMPE[2], -0.22), (800, RAMPE[5], 0.0), (850, RAMPE[9], 0.22))
        r = st[(st.source.=="reel").&(st.lambda_nm.==λx).&(st.N_bleu.==1e4).&(st.technique.==k), :]
        nrow(r) == 1 && scatter!(ax4, [max(r.d_prime[1], 0.4)], [i + dy]; color = col, marker = marqueur_banc(k), markersize = 11)
    end
    xlims!(ax4, 0.35, 200)
    axislegend(ax4, [MarkerElement(color = c, marker = :circle, markersize = 10) for c in (RAMPE[2], RAMPE[5], RAMPE[9])],
        ["750 nm", "800 nm", "850 nm"]; position = :rb)
    entete!(fig, "Fig. 16 — Banc du 2026-10-16 : quelle lecture reste stable quand les photons manquent ?",
        "Sommes des 10 acquisitions amincies par tirage binomial (400 tirages par point) ; chaque lecture est ramenée sur l'écart entre les deux échantillons " *
        "(0 à 0 mM, 1 à 131 mM, valeurs des sommes à pleins photons). Marqueur vide : plus de 5 % d'échecs. Le jumeau tire ses photons du lissage à trois " *
        "exponentielles ; la borne de Fisher suppose les fonds connus."; largeur = 3)
    save(chemin, fig; px_per_unit = 2)
end

end
