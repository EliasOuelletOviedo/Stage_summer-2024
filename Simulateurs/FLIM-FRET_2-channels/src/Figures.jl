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

# -----------------------------------------------------------------------------
# Figure 1 — σ[Cl⁻] selon [Cl⁻]
# -----------------------------------------------------------------------------
function figure1(df, chemin; N = 1e4, norm = "A", cles = AVEC_2F)
    fig = Figure(size = (980, 600))
    d = sel(df; N_photons = N, normalisation = norm, mode_calibration = "exacte")
    ax = Axis(fig[1, 1]; xlabel = "[Cl⁻] vraie (mM)", ylabel = "σ[Cl⁻] (mM, échelle log)", yscale = log10, yticks = graduations_log(0.1, 1000),
        title = "Lequel gagne, et où sur l'isotherme ?")
    fenetre_clamping!(ax)
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
        "La borne est infinie à 0 mM (pente nulle de l'isotherme, n > 1)." * NOTE_CREUX; largeur = 2)
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
function figure2(df, chemin; C = 20.0, cles = AVEC_2F)
    fig = Figure(size = (1350, 560))
    axes = Axis[]
    for (j, norm) in enumerate(("A", "B", "C"))
        d = sel(df; C_vrai = C, normalisation = norm, mode_calibration = "exacte")
        ax = Axis(fig[1, j]; xscale = log10, yscale = log10, xlabel = "budget N (photons)", ylabel = j == 1 ? "σ[Cl⁻] (mM)" : "",
            yticks = graduations_log(0.1, 1000),
            title = "Normalisation " * norm, xticks = ([1e2, 1e3, 1e4, 1e5], ["10²", "10³", "10⁴", "10⁵"]))
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
        "A : temps d'acquisition fixé (N = photons bleus à 20 mM) ; B : N_bleu + N_jaune = N ; C : N_bleu = N." * NOTE_CREUX; largeur = 4)
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
        ("epsilonPhi_dark", "ε_dark·Φ_dark (forme protonée)"), ("tau_D0_relatif", "τ_D0 vrai / τ_D0 calibré"),
        ("tau_DA_relatif", "τ_DA vrai / τ_DA calibré"), ("tau_A_relatif", "τ_A vrai / τ_A supposé")]
    fig = Figure(size = (1500, 840))
    axes = Axis[]
    for (i, (p, xl)) in enumerate(panneaux)
        r, c = fldmod1(i, 3)
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
    Legend(fig[1:2, 4], axes[end]; merge = true, unique = true)
    entete!(fig, "Fig. 4b — Biais sous les autres nuisances, calibration nominale",
        @sprintf("Normalisation A, %d mM, N = %s. Données issues du monde perturbé, estimateurs et calibrations construits sur les valeurs nominales ", Int(C), fmtN(N)) *
        "(gain 1 ; ε = 0,03 ; ε_dark·Φ_dark = 0 ; τ de la titration ; τ_A = 3,0 ns). Dans les panneaux τ, l'estimateur conjoint est le 5b du plan (τ figés faux). " *
        "Le gain du jaune multiplie tout le canal jaune (dérive de détecteur) : il ne change pas sa forme." * NOTE_CREUX; largeur = 4)
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
function figure5(kapp, courbes, chemin; K_vrai = NaN)
    fig = Figure(size = (1300, 600))
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
        text!(ax2, vi, yi; text = @sprintf(" %.2f", vi), color = ENCRE2, fontsize = 11, align = (:left, :center))
    end
    # tes valeurs mesurées, en face de l'estimateur correspondant
    for (mk, k, txt) in (("mesure_tau", "2", "mesuré (ton Hill sur τ)"), ("mesure_fret", "1inv", "mesuré (ton ratio TIFF C2/C1)"))
        i = findfirst(==(mk), kapp.estimateur)
        j = findfirst(==(k), rows)
        (i === nothing || j === nothing) && continue
        scatter!(ax2, [kapp.rapport_K[i]], [y[j] + 0.32]; color = ENCRE, marker = :diamond, markersize = 11)
        text!(ax2, kapp.rapport_K[i], y[j] + 0.32; text = @sprintf("  %s : %.2f", txt, kapp.rapport_K[i]), color = ENCRE, fontsize = 11, align = (:left, :center))
    end
    xlims!(ax2, 0, 2.35)
    colsize!(fig.layout, 1, Relative(0.45))
    entete!(fig, "Fig. 5 — K apparent / K vrai par estimateur",
        @sprintf("K vrai = %.1f mM (isotherme sur les x de l'ajustement global). Hill à 4 paramètres ajustée à la courbe sans bruit de chaque estimateur sur 0–50 K. ", K_vrai) *
        "Losanges : tes ajustements de Hill expérimentaux (tests/csv_exports/hill_fit_summary.csv), en face de l'estimateur qui leur correspond."; largeur = 2)
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

    # (e) intensité relative (validation croisée gratuite)
    ax5 = Axis(fig[1, 3]; xticks = TICKS_C, xlabel = "[Cl⁻] (mM)", ylabel = "A / moyenne",
        title = "Intensité du donneur à t = 0 (doit être constante)")
    hspan!(ax5, 0.9, 1.1; color = (ZONE, 1.0))
    hlines!(ax5, [1.0]; color = AXE, linewidth = 1)
    scatter!(ax5, t.pos, t.intensite_relative; color = RAMPE[t.ic], markersize = 9)
    ylims!(ax5, 0.6, 1.5)

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

end
