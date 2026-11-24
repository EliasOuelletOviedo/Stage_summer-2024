# Étape 20 — biais sur [Cl⁻] sous autofluorescence (figure A.2), SCLM (Cerulean + Topaz), 850 nm.
#
# L'argument à démontrer : un ajustement à durées de vie caractérisées peut séparer
# l'autofluorescence du signal du capteur, alors qu'une mesure ratiométrique la confond
# avec lui. On l'éprouve plutôt que de l'affirmer :
#
#   1. Le monde : le capteur de la titration du 2026-10-16 à 850 nm, tel que l'ajustement
#      global de l'étape 15 le décrit (bleu : τ_DA = τ₂, τ_D = τ₃ ; jaune : τ_J, τ_A ;
#      diffusion et fond), IRF mesurées des deux voies. Entre 0 et 131 mM, mélange
#      statique des deux états mesurés : fraction g(c) = c(131 + K)/(131(c + K)), K ajusté
#      sur les 8 concentrations.
#   2. L'autofluorescence (valeurs typiques de la littérature, à deux photons) :
#        NAD(P)H dans le bleu : libre 0,4 ns (75 % de l'amplitude), lié 2,0 ns ;
#        flavines (FAD) dans le jaune : liées, éteintes, 0,3 ns (75 %), libres 2,6 ns.
#      Intensité absolue constante (indépendante du Cl⁻), exprimée en fraction p des photons
#      du canal à 20 mM.
#   3. Quatre lectures, chacune calibrée sur le monde sans autofluorescence, puis appliquée
#      au monde contaminé ; biais = [Cl⁻] lue − [Cl⁻] vraie :
#        - rapport J/B brut ;
#        - τ du donneur, mono-exponentielle sur le bleu seul (fond uniforme libre) ;
#        - 2 canaux, τ caractérisés, échelle commune (le 5d de l'étape 15 : x et S libres),
#          sans terme d'autofluorescence (témoin : sépare l'effet du terme de celui des 2 canaux) ;
#        - le même, plus un terme d'autofluorescence par canal à τ_auto caractérisés
#          (amplitudes libres, de signe libre).
#      IRF, t₀, diffusion et fond du détecteur connus pour toutes les lectures.
#   4. Robustesse : autofluorescence réelle différente de celle supposée par le modèle
#      (fractions et τ décalés dans les plages publiées).
#   5. Le prix : écart-type de la [Cl⁻] lue à N photons bleus (Monte-Carlo).
#
#   julia --project=. step20_autofluorescence.jl [--niveau 0.05] [--N 1e4] [--tirages 300]
#
# Sorties (sorties/autofluorescence_2026-10-16_SCLM_850nm/) : autofluorescence.png / .pdf,
# biais_selon_cl.csv, biais_selon_niveau.csv, precision.csv, monde.toml

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
include(joinpath(@__DIR__, "src", "deux_canaux", "DeuxCanaux.jl"))
using .SimFLIMFRET, .DeuxCanaux
const S = SimFLIMFRET
const D = DeuxCanaux
using Printf, Statistics, Random, StableRNGs, CSV, DataFrames

args = copy(ARGS)
valeur(nom, defaut) = (i = findfirst(==(nom), args); i === nothing ? defaut : args[i+1])
p_ref = parse(Float64, valeur("--niveau", "0.05"))
N_prec = parse(Float64, valeur("--N", "1e4"))
ntir = parse(Int, valeur("--tirages", "600"))
λ = 850
donnees = joinpath(@__DIR__, "Data", "2026-10-16", "SCLM")
dossier = joinpath(@__DIR__, "sorties", "autofluorescence_2026-10-16_SCLM_$(λ)nm")
mkpath(dossier)

# -----------------------------------------------------------------------------
# 1. Le monde : ajustement global (étape 15, modèle libre) et mélange à deux états
# -----------------------------------------------------------------------------
cellules, _ = D.charger(donnees; λs = (λ,))
lie = D.ajuster_global(cellules; lien = true)
mg = D.ajuster_global(cellules; lien = false, τ_init = (lie.τ2, lie.τ3, lie.τA, 1.0), t0_init = (lie.t0B, lie.t0J))
c1 = cellules[1]
Δ = c1.pas_ns
W = collect(D.FENETRE)
nb = size(c1.irf, 2)
instB = S.instrument_banc(c1.irf[D.BLEU, :], Δ; t₀ = mg.t0B, fenetre = D.FENETRE)
instJ = S.instrument_banc(c1.irf[D.JAUNE, :], Δ; t₀ = mg.t0J, fenetre = D.FENETRE)
wsB, wsJ = S.EspaceTravail(instB), S.EspaceTravail(instJ)
function base(ws, τ, t₀)
    b = zeros(nb)
    S.base!(b, ws, τ, t₀)
    return S.masquer!(b, D.FENETRE)
end
B2, B3, Bd = base(wsB, mg.τ2, mg.t0B), base(wsB, mg.τ3, mg.t0B), base(wsB, D.TAU_DIFFUSION, mg.t0B)
JJ, JA, Jd = base(wsJ, mg.τJ, mg.t0J), base(wsJ, mg.τA, mg.t0J), base(wsJ, D.TAU_DIFFUSION, mg.t0J)
cl = mg.cl
i0, i1 = findfirst(==(0.0), cl), argmax(cl)
cmax = cl[i1]
FB(i) = mg.aB[i, 2] .* B2 .+ mg.aB[i, 3] .* B3          # fluorescence du capteur, bleu
FJ(i) = mg.aJ[i, 2] .* JJ .+ mg.aJ[i, 3] .* JA          # fluorescence du capteur, jaune
DB, DJ = mean(mg.aB[:, 1]) .* Bd, mean(mg.aJ[:, 1]) .* Jd   # diffusion : indépendante du Cl⁻
fB, fJ = S.masquer!(mean(mg.fondB), D.FENETRE), S.masquer!(mean(mg.fondJ), D.FENETRE)
FB0, FB1, FJ0, FJ1 = FB(i0), FB(i1), FJ(i0), FJ(i1)

# fraction g de l'état « 131 mM » qui reproduit au mieux chaque concentration (moindres carrés pondérés)
gmes = map(eachindex(cl)) do i
    yb, yj = FB(i), FJ(i)
    wb, wj = 1 ./ (yb .+ DB .+ fB), 1 ./ (yj .+ DJ .+ fJ)
    eb, ej = FB1 .- FB0, FJ1 .- FJ0
    num = sum(wb[W] .* (yb[W] .- FB0[W]) .* eb[W]) + sum(wj[W] .* (yj[W] .- FJ0[W]) .* ej[W])
    den = sum(wb[W] .* eb[W] .^ 2) + sum(wj[W] .* ej[W] .^ 2)
    num / den
end
gK(c, K) = c * (cmax + K) / (cmax * (c + K))
Ks = exp.(range(log(1.0), log(500.0), length = 5000))
K = Ks[argmin([sum(abs2, gK.(cl, k) .- gmes) for k in Ks])]
@printf("SCLM %d nm : τ_DA = %.3f, τ_D = %.3f, τ_J = %.3f, τ_A = %.3f ns ; mélange à deux états, K = %.1f mM\n",
    λ, mg.τ2, mg.τ3, mg.τJ, mg.τA, K)
for (c, g) in zip(cl, gmes)
    @printf("   %6g mM : g mesuré %.3f, isotherme %.3f\n", c, g, gK(c, K))
end

function monde(c)
    g = gK(c, K)
    return (FB0 .* (1 - g) .+ FB1 .* g .+ DB .+ fB, FJ0 .* (1 - g) .+ FJ1 .* g .+ DJ .+ fJ)
end
totB20, totJ20 = sum(monde(20.0)[1][W]), sum(monde(20.0)[2][W])

# -----------------------------------------------------------------------------
# 2. Autofluorescence
# -----------------------------------------------------------------------------
forme(α, τ, inst, ws) = (m = S.masquer!(S.forward(S.ExpSum(α, τ), inst, ws), D.FENETRE); m ./ sum(m[W]))
const NADH = (α = [0.75, 0.25], τ = [0.4, 2.0])          # supposée par le modèle (valeurs typiques)
const FAD = (α = [0.75, 0.25], τ = [0.3, 2.6])
const NADH_AUTRE = (α = [0.60, 0.40], τ = [0.5, 2.5])    # réelle, différente de la supposée
const FAD_AUTRE = (α = [0.60, 0.40], τ = [0.4, 2.9])
nadh, fad = forme(NADH..., instB, wsB), forme(FAD..., instJ, wsJ)
nadh_autre, fad_autre = forme(NADH_AUTRE..., instB, wsB), forme(FAD_AUTRE..., instJ, wsJ)

"""Monde à [Cl⁻] = c, plus pN (bleu) et pF (jaune) d'autofluorescence ; `autre` : formes différentes de celles du modèle."""
function contamine(c, pN, pF; autre = false)
    yb, yj = monde(c)
    return (yb .+ pN * totB20 .* (autre ? nadh_autre : nadh), yj .+ pF * totJ20 .* (autre ? fad_autre : fad))
end

# -----------------------------------------------------------------------------
# 3. Les lectures
# -----------------------------------------------------------------------------
lire_ratio(yb, yj) = sum(yj[W]) / sum(yb[W])

const Bt, dBτ, dBt = zeros(nb), zeros(nb), zeros(nb)
const B_DEPART = base(wsB, 2.0, mg.t0B)
function lire_mono(yb, _)
    yW = yb[W]
    function modele!(μ, J, θ)
        S.bases_derivees!(Bt, dBτ, dBt, wsB, θ[1], mg.t0B)
        @inbounds for (i, k) in enumerate(W)
            μ[i] = θ[2] * Bt[k] + θ[3]
            J[i, 1] = θ[2] * dBτ[k]
            J[i, 2] = Bt[k]
            J[i, 3] = 1.0
        end
    end
    θ0 = [2.0, 0.9 * sum(yW) / sum(B_DEPART[W]), 0.02 * sum(yW) / length(W)]
    r = S.lm_poisson(modele!, θ0, [0.05, 0.0, 0.0], [20.0, Inf, Inf], yW; maxiter = 400, tol_dec = 1e-12)
    return r.θ[1]
end

const uB, vB, uJ, vJ = FB0[W], (FB1.-FB0)[W], FJ0[W], (FJ1.-FJ0)[W]
const kB, kJ = (DB.+fB)[W], (DJ.+fJ)[W]
const nN, nF = nadh[W], fad[W]
"""2 canaux, τ caractérisés, échelle commune : x (= g) et S libres ; `auto` : + amplitudes NADH (bleu) et FAD (jaune)."""
function lire_deux(yb, yj; auto::Bool, echelle = 1.0)
    y = vcat(yb[W], yj[W])
    m = length(W)
    function modele!(μ, J, θ)
        x, s = θ[1], θ[2]
        fill!(J, 0.0)
        @inbounds for i in 1:m
            μ[i] = s * (uB[i] + x * vB[i]) + echelle * kB[i]
            μ[m+i] = s * (uJ[i] + x * vJ[i]) + echelle * kJ[i]
            J[i, 1], J[i, 2] = s * vB[i], uB[i] + x * vB[i]
            J[m+i, 1], J[m+i, 2] = s * vJ[i], uJ[i] + x * vJ[i]
            if auto
                μ[i] += θ[3] * nN[i]
                μ[m+i] += θ[4] * nF[i]
                J[i, 3], J[m+i, 4] = nN[i], nF[i]
            end
        end
        μ .= max.(μ, 1e-12)
    end
    θ0 = auto ? [0.5, echelle, 0.0, 0.0] : [0.5, echelle]
    lo = auto ? [-3.0, 0.0, -Inf, -Inf] : [-3.0, 0.0]
    hi = auto ? [4.0, Inf, Inf, Inf] : [4.0, Inf]
    r = S.lm_poisson(modele!, θ0, lo, hi, y; maxiter = 400, tol_dec = 1e-12)
    return r.θ[1]
end

LECTURES = [
    ("ratio", "rapport J/B (brut)", (yb, yj; e = 1.0) -> lire_ratio(yb, yj)),
    ("mono", "τ du donneur, mono-exponentielle (bleu seul)", (yb, yj; e = 1.0) -> lire_mono(yb, yj)),
    ("deux", "2 canaux, τ caractérisés, sans terme d'autofluorescence", (yb, yj; e = 1.0) -> lire_deux(yb, yj; auto = false, echelle = e)),
    ("deux_auto", "2 canaux, τ caractérisés + τ_auto", (yb, yj; e = 1.0) -> lire_deux(yb, yj; auto = true, echelle = e)),
]

# calibration de chaque lecture sur le monde propre : grille dense, g de −0,25 à g(400 mM)
c_bas = -0.25 * cmax * K / (1.25cmax + K)
cgrille = collect(range(c_bas, 400.0, length = 800))
cal = Dict{String,Vector{Float64}}()
for (cle, _, f) in LECTURES
    cal[cle] = [f(monde(c)...) for c in cgrille]
    d = diff(cal[cle])
    (all(d .> 0) || all(d .< 0)) || @warn "calibration non monotone : $cle"
end
function inverser(cle, X)
    xs, cs = cal[cle], cgrille
    isfinite(X) || return NaN
    if xs[end] < xs[1]
        xs, cs = reverse(xs), reverse(cs)
    end
    (X < xs[1] || X > xs[end]) && return NaN               # hors de la calibration
    j = clamp(searchsortedlast(xs, X), 1, length(xs) - 1)
    return cs[j] + (X - xs[j]) * (cs[j+1] - cs[j]) / (xs[j+1] - xs[j])
end
lire_c(cle, f, yb, yj; e = 1.0) = inverser(cle, f(yb, yj; e = e))

# contrôle : sans autofluorescence, chaque lecture relit la vraie [Cl⁻]
for (cle, _, f) in LECTURES
    e = maximum(abs(lire_c(cle, f, monde(c)...) - c) for c in (0.0, 5.0, 20.0, 80.0))
    @printf("contrôle sans autofluorescence, %-10s : écart max %.1e mM\n", cle, e)
end

# fraction des photons de l'autofluorescence autour du pic (±0,25 ns) : comparaison avec la composante prompte des tranches
pic = argmax(monde(20.0)[1])
fen_pic = filter(k -> abs(k - pic) * Δ <= 0.25, W)
part_pic(pN) = (a = pN * totB20 .* nadh; sum(a[fen_pic]) / (sum(monde(20.0)[1][fen_pic]) + sum(a[fen_pic])))
@printf("NADH à %.0f %% des photons bleus : %.0f %% des photons à ±0,25 ns du pic\n", 100p_ref, 100part_pic(p_ref))

# -----------------------------------------------------------------------------
# 4. Biais
# -----------------------------------------------------------------------------
SCENARIOS = [("NADH", "NADH dans le bleu", (p) -> (p, 0.0)), ("FAD", "flavines (FAD) dans le jaune", (p) -> (0.0, p)),
             ("les_deux", "NADH et FAD", (p) -> (p, p))]
ccourbe = collect(0.0:2.5:80.0)
lignes = DataFrame(scenario = String[], formes = String[], niveau = Float64[], cl_vraie = Float64[], lecture = String[], cl_lue = Float64[], biais = Float64[])
for (sc, _, niv) in SCENARIOS, autre in (false, true), c in ccourbe
    yb, yj = contamine(c, niv(p_ref)...; autre)
    for (cle, _, f) in LECTURES
        (autre && cle != "deux_auto") && continue
        cr = lire_c(cle, f, yb, yj)
        push!(lignes, (sc, autre ? "différentes du modèle" : "celles du modèle", p_ref, c, cle, cr, cr - c))
    end
end
CSV.write(joinpath(dossier, "biais_selon_cl.csv"), lignes)
niveaux = [0.01, 0.02, 0.05, 0.10, 0.20]
lv = DataFrame(scenario = String[], niveau = Float64[], cl_vraie = Float64[], lecture = String[], biais = Float64[])
for (sc, _, niv) in SCENARIOS, p in niveaux, c in (5.0, 20.0)
    yb, yj = contamine(c, niv(p)...)
    for (cle, _, f) in LECTURES
        push!(lv, (sc, p, c, cle, lire_c(cle, f, yb, yj) - c))
    end
end
CSV.write(joinpath(dossier, "biais_selon_niveau.csv"), lv)
println("\nBiais (mM) à $(round(Int, 100p_ref)) % d'autofluorescence :")
for (sc, _, _) in SCENARIOS, c in (5.0, 20.0)
    l = lignes[(lignes.scenario .== sc) .& (lignes.cl_vraie .== c), :]
    @printf("  %-9s %4.0f mM : ", sc, c)
    for r in eachrow(l)
        @printf("%s%s %+.2f   ", r.lecture, r.formes == "celles du modèle" ? "" : " (formes autres)", r.biais)
    end
    println()
end

# -----------------------------------------------------------------------------
# 5. Précision à N photons bleus (Monte-Carlo, autofluorescence NADH + FAD au niveau de référence)
# -----------------------------------------------------------------------------
k_N = N_prec / totB20
cprec = [0.0, 2.5, 5.0, 10.0, 20.0, 40.0, 80.0]
rng = StableRNG(2020)
prec = DataFrame(cl_vraie = Float64[], lecture = String[], N_bleu = Float64[], ecart_type_robuste = Float64[], biais_median = Float64[], hors_echelle = Int[])
for c in cprec
    yb, yj = contamine(c, p_ref, p_ref)
    μb, μj = k_N .* yb, k_N .* yj
    nbt, njt = zeros(nb), zeros(nb)
    lus = Dict(cle => Float64[] for (cle, _, _) in LECTURES)
    for _ in 1:ntir
        S.tirer!(nbt, μb, rng)
        S.tirer!(njt, μj, rng)
        for (cle, _, f) in LECTURES
            push!(lus[cle], lire_c(cle, f, nbt, njt; e = k_N))
        end
    end
    for (cle, _, _) in LECTURES
        v = lus[cle]
        ok = filter(isfinite, v)
        σ = length(ok) > 10 ? (quantile(ok, 0.75) - quantile(ok, 0.25)) / 1.349 : NaN
        push!(prec, (c, cle, N_prec, σ, length(ok) > 0 ? median(ok) - c : NaN, length(v) - length(ok)))
    end
    @printf("précision, %5.1f mM : %s\n", c, join([@sprintf("%s σ %.1f", r.lecture, r.ecart_type_robuste) for r in eachrow(prec[prec.cl_vraie .== c, :])], " · "))
end
CSV.write(joinpath(dossier, "precision.csv"), prec)
open(joinpath(dossier, "monde.toml"), "w") do io
    println(io, "# Monde de l'étape 20 (SCLM, $(λ) nm, titration du 2026-10-16)")
    @printf(io, "tau_DA = %.4f\ntau_D = %.4f\ntau_J = %.4f\ntau_A = %.4f\nK_melange_mM = %.2f\n", mg.τ2, mg.τ3, mg.τJ, mg.τA, K)
    println(io, "NADH = { alpha = $(NADH.α), tau = $(NADH.τ) }\nFAD = { alpha = $(FAD.α), tau = $(FAD.τ) }")
    println(io, "NADH_autre = { alpha = $(NADH_AUTRE.α), tau = $(NADH_AUTRE.τ) }\nFAD_autre = { alpha = $(FAD_AUTRE.α), tau = $(FAD_AUTRE.τ) }")
    @printf(io, "niveau_reference = %.3f\nN_bleu_precision = %.0f\n", p_ref, N_prec)
end

# -----------------------------------------------------------------------------
# 6. Figure
# -----------------------------------------------------------------------------
include(joinpath(@__DIR__, "src", "Figures.jl"))
using CairoMakie
const F = Figures
set_theme!(F.theme_sim())
virgs(v) = [F.virg(@sprintf("%g", x)) for x in v]
STYLE = Dict("ratio" => (F.COULEUR["1"], :solid, 2.0), "mono" => (F.COULEUR["2"], :solid, 2.0),
             "deux" => (F.COULEUR["5d"], :dash, 1.8), "deux_auto" => (F.COULEUR["6"], :solid, 3.0))
NOM = Dict(cle => nom for (cle, nom, _) in LECTURES)
pc = round(Int, 100p_ref)
fig = Figure(size = (1300, 1040))
Label(fig[0, 1:2], "Biais sur [Cl⁻] sous autofluorescence · SCLM (Cerulean + Topaz) · $(λ) nm · monde de la titration du 2026-10-16";
    halign = :left, fontsize = 16, color = F.ENCRE)

# (a) les formes
tc = ((0:nb-1) .+ 0.5) .* Δ
ga = fig[1, 1] = GridLayout()
for (j, (titre, s0, s1, a, sous)) in enumerate((
        ("(a) Formes des déclins · canal bleu", FB0, FB1, nadh,
         F.virg(@sprintf("capteur : τ_DA = %.2f, τ_D = %.2f ns · NADH : %.0f %% à %.1f + %.0f %% à %.1f ns", mg.τ2, mg.τ3, 100NADH.α[1], NADH.τ[1], 100NADH.α[2], NADH.τ[2]))),
        ("canal jaune", FJ0, FJ1, fad,
         F.virg(@sprintf("capteur : τ_J = %.2f, τ_A = %.2f ns · FAD : %.0f %% à %.1f + %.0f %% à %.1f ns", mg.τJ, mg.τA, 100FAD.α[1], FAD.τ[1], 100FAD.α[2], FAD.τ[2])))))
    ax = Axis(ga[j, 1]; yscale = log10, ylabel = "normalisé (max = 1)", xtickformat = virgs, title = titre, subtitle = sous,
        xlabel = j == 2 ? "temps du TAC (ns)" : "")
    for (v, col, st) in ((s0, F.ENCRE, :solid), (s1, F.ENCRE2, :dash), (a, F.COULEUR["6"], :solid))
        ok = v[W] .> 1e-4 * maximum(v[W])
        lines!(ax, tc[W][ok], v[W][ok] ./ maximum(v[W]); color = col, linestyle = st, linewidth = 1.6)
    end
    ylims!(ax, 1e-3, 1.5)
    ax.yticks = ([1e-3, 1e-2, 1e-1, 1], ["0,001", "0,01", "0,1", "1"])
    j == 1 && hidexdecorations!(ax; grid = false, ticks = false)
end
Legend(ga[3, 1], [LineElement(color = F.ENCRE, linestyle = :solid), LineElement(color = F.ENCRE2, linestyle = :dash),
        LineElement(color = F.COULEUR["6"], linestyle = :solid)],
    ["capteur, 0 mM", "capteur, 131 mM", "autofluorescence (NADH en haut, FAD en bas)"];
    orientation = :horizontal, framevisible = false, labelsize = 11, tellheight = true)
rowgap!(ga, 4)

# (b) et (c) biais contre [Cl⁻]
function panneau_biais!(pos, sc, titre)
    ax = Axis(pos; xlabel = "[Cl⁻] vraie (mM)", ylabel = "biais : [Cl⁻] lue − vraie (mM)", xtickformat = virgs, ytickformat = virgs,
        title = titre, subtitle = "$(pc) % des photons du canal à 20 mM, intensité indépendante du Cl⁻")
    hlines!(ax, [0.0]; color = F.MUET, linewidth = 1)
    for (cle, _, _) in LECTURES
        l = lignes[(lignes.scenario .== sc) .& (lignes.lecture .== cle) .& (lignes.formes .== "celles du modèle"), :]
        col, st, lw = STYLE[cle]
        p = lines!(ax, l.cl_vraie, l.biais; color = col, linestyle = st, linewidth = lw, label = NOM[cle])
        cle == "mono" && (translate!(p, 0, 0, 10); p.linestyle = :dash; p.linewidth = 2.6)      # au-dessus du trait épais à 0 (flavines : biais nul)
    end
    l = lignes[(lignes.scenario .== sc) .& (lignes.formes .== "différentes du modèle"), :]
    lines!(ax, l.cl_vraie, l.biais; color = F.COULEUR["6"], linestyle = :dot, linewidth = 2.2,
        label = "2 canaux + τ_auto, autofluorescence réelle ≠ supposée")
    xlims!(ax, 0, 80)
    return ax
end
axb = panneau_biais!(fig[2, 1], "NADH", "(c) NADH dans le bleu")
axc = panneau_biais!(fig[2, 2], "FAD", "(d) Flavines (FAD) dans le jaune")

# (d) précision
axd = Axis(fig[1, 2]; xlabel = "[Cl⁻] vraie (mM)", ylabel = "écart-type de la [Cl⁻] lue (mM)", yscale = log10,
    xtickformat = virgs, title = F.virg(@sprintf("(b) Le prix : précision à %s photons bleus", N_prec == 1e4 ? "10⁴" : @sprintf("%.0f", N_prec))),
    subtitle = "NADH et FAD à $(pc) % ; écart-type robuste de $(ntir) tirages de Poisson")
for (cle, _, _) in LECTURES
    l = prec[prec.lecture .== cle, :]
    col, st, lw = STYLE[cle]
    ok = isfinite.(l.ecart_type_robuste)
    scatterlines!(axd, l.cl_vraie[ok], l.ecart_type_robuste[ok]; color = col, linestyle = st, linewidth = lw, markersize = 7)
end
yt = [0.5, 1, 2, 5, 10, 20, 50, 100, 200]
axd.yticks = (yt, virgs(yt))
xlims!(axd, -2, 82)

Legend(fig[3, 1:2], axb; orientation = :horizontal, nbanks = 2, framevisible = false, labelsize = 12)
note = "Lectures calibrées sur le monde sans autofluorescence, puis appliquées au monde contaminé ; IRF, t₀, diffusion et fond du détecteur connus pour toutes.\n" *
       "Autofluorescence (valeurs typiques à deux photons) : NAD(P)H libre 0,3–0,4 ns, liée 1,9–5,7 ns ; FAD libre 2,3–2,9 ns, liée et éteinte, courte. " *
       "À 850 nm, le NADH est peu excité (efficace à 710–780 nm) : les flavines dominent.\n" *
       "Pointillé : autofluorescence réelle différente de la supposée (NADH 60 % à 0,5 + 40 % à 2,5 ns ; FAD 60 % à 0,4 + 40 % à 2,9 ns). " *
       @sprintf("Biais proportionnels au niveau ; 5 %% des photons bleus = %.0f %% des photons à ±0,25 ns du pic.\n", 100part_pic(p_ref)) *
       @sprintf("Isotherme à n = 1, K = %.0f mM (8 points) : sous 10 mM, la titration est 2 à 3 fois moins sensible que l'isotherme, ", K) *
       "donc les mM y sont sous-estimés d'autant, pour toutes les lectures à la fois (leur classement n'en dépend pas)."
Label(fig[4, 1:2], F.virg(note); halign = :left, fontsize = 11, color = F.ENCRE2, justification = :left)
save(joinpath(dossier, "autofluorescence.png"), fig; px_per_unit = 2)
save(joinpath(dossier, "autofluorescence.pdf"), fig)
println("-> ", joinpath(dossier, "autofluorescence.png"))
