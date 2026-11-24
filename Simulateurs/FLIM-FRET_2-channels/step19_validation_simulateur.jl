# Étape 19 — validation du simulateur sur la titration SCLM (régulier) du 2026-10-16, 850 nm.
#
# Réponse à « comment sais-tu que ton simulateur est juste ? », canal bleu :
#
# (a) Aller-retour : on injecte une fraction d'amplitude FRET α dans le déclin bleu
#       d(t) = α·e^{−t/τ₂} + (1 − α)·e^{−t/τ₃}   (+ diffusion et fond, comme dans les mesures)
#     avec le moteur du simulateur (S.forward : noyau périodisé, convolution par l'IRF
#     mesurée, intégration par canal), on tire du Poisson pour N photons bleus, et on
#     récupère α par maximum de vraisemblance, τ et IRF connus. α récupéré contre α
#     injecté, pour N = 10³, 10⁴, 10⁵ : la diagonale, et la dispersion qui rétrécit en 1/√N.
#
# (b) Réel contre simulé : la somme réelle à 20 mM (10 acquisitions), et une décroissance
#     simulée avec les mêmes τ, la même IRF, le même t₀, les mêmes amplitudes (ajustement
#     global de l'étape 15 sur les 8 concentrations) et exactement le même N (tirage
#     multinomial). Résidus de Pearson des deux contre le même modèle : si le simulateur
#     reproduit la mesure, ils ont la même allure et le même χ²ᵣ ; l'écart de χ²ᵣ mesure
#     ce que le simulateur ne capte pas.
#
#   julia --project=. step19_validation_simulateur.jl [--cl 20] [--tirages 500]
#
# Sorties (sorties/validation_simulateur_2026-10-16_SCLM_850nm/) :
#   validation_simulateur.png / .pdf, aller_retour.csv, reel_simule.csv

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
include(joinpath(@__DIR__, "src", "deux_canaux", "DeuxCanaux.jl"))
using .SimFLIMFRET, .DeuxCanaux
const S = SimFLIMFRET
const D = DeuxCanaux
using Printf, Statistics, Random, StableRNGs, Distributions, CSV, DataFrames

args = copy(ARGS)
valeur(nom, defaut) = (i = findfirst(==(nom), args); i === nothing ? defaut : args[i+1])
cl_b = parse(Float64, valeur("--cl", "20"))
ntir = parse(Int, valeur("--tirages", "500"))
λ = 850
donnees = joinpath(@__DIR__, "Data", "2026-10-16", "SCLM")
dossier = joinpath(@__DIR__, "sorties", "validation_simulateur_2026-10-16_SCLM_$(λ)nm")
mkpath(dossier)

# -----------------------------------------------------------------------------
# 1. Paramètres tirés des mesures : ajustement global (modèle libre, étape 15)
# -----------------------------------------------------------------------------
cellules, _ = D.charger(donnees; λs = (λ,))
lie = D.ajuster_global(cellules; lien = true)
mg = D.ajuster_global(cellules; lien = false, τ_init = (lie.τ2, lie.τ3, lie.τA, 1.0), t0_init = (lie.t0B, lie.t0J))
ic = findfirst(c -> c.cl == cl_b, cellules)
c = cellules[ic]
W = collect(D.FENETRE)
Δ = c.pas_ns
@printf("SCLM %d nm : τ₂ = %.3f ns, τ₃ = %.3f ns, t₀ bleu = %+.1f ps, χ²ᵣ global %.2f\n", λ, mg.τ2, mg.τ3, 1000mg.t0B, mg.χ2)

# -----------------------------------------------------------------------------
# 2. Le monde simulé : moteur du simulateur, IRF mesurée du bleu
# -----------------------------------------------------------------------------
inst = S.instrument_banc(c.irf[D.BLEU, :], Δ; t₀ = mg.t0B, fenetre = D.FENETRE)
ws = S.EspaceTravail(inst)
base(τ) = (b = zeros(inst.n_bins); S.base!(b, ws, τ, inst.t₀); b)
Bd, B2, B3 = base(D.TAU_DIFFUSION), base(mg.τ2), base(mg.τ3)
fond = mg.fondB[ic]
aB = mg.aB[ic, :]
# déclin du capteur par le simulateur (S.forward), plus diffusion et fond tels qu'ajustés
μ_sim = S.forward(S.ExpSum([aB[2], aB[3]], [mg.τ2, mg.τ3]), inst, ws) .+ aB[1] .* Bd .+ fond
S.masquer!(μ_sim, D.FENETRE)
ecart_moteur = maximum(abs.(μ_sim[W] .- mg.μB[ic][W])) / maximum(mg.μB[ic][W])
@printf("Moteur du simulateur contre modèle de l'ajustement : écart relatif max %.1e\n", ecart_moteur)

# -----------------------------------------------------------------------------
# 3. (b) réel contre simulé, même N
# -----------------------------------------------------------------------------
y_reel = D.somme(c, D.BLEU)
N = round(Int, sum(y_reel[W]))
rng = StableRNG(19)
y_sim = zeros(length(y_reel))
y_sim[W] .= rand(rng, Multinomial(N, μ_sim[W] ./ sum(μ_sim[W])))
μN = μ_sim .* (N / sum(μ_sim[W]))
pearson(y) = (y[W] .- μN[W]) ./ sqrt.(μN[W])
r_reel, r_sim = pearson(y_reel), pearson(y_sim)
χ2r(r) = sum(abs2, r) / (length(W) - 1)
@printf("%g mM, N = %d photons bleus : χ²ᵣ réel %.2f, simulé %.2f\n", cl_b, N, χ2r(r_reel), χ2r(r_sim))

# Le χ²ᵣ d'un seul tirage simulé fluctue (écart-type ≈ √(2/181) ≈ 0,11) : moyenne sur 200 tirages.
# Et aux N des simulations : le réel aminci (tirage binomial canal par canal) contre le simulé au même N.
p_mod = μ_sim[W] ./ sum(μ_sim[W])
function khi2(y)
    Nt = sum(y[W])
    m = p_mod .* Nt
    return sum(abs2, (y[W] .- m) ./ sqrt.(m)) / (length(W) - 1)
end
rng2 = StableRNG(2019)
function khi2_moyen(Nt; reel::Bool, n = 200)
    v = Float64[]
    y = zeros(length(y_reel))
    for _ in 1:n
        if reel
            y[W] .= [rand(rng2, Binomial(round(Int, k), Nt / N)) for k in y_reel[W]]
        else
            y[W] .= rand(rng2, Multinomial(Nt, p_mod))
        end
        push!(v, khi2(y))
    end
    return mean(v)
end
Nk = (1_000, 10_000, 100_000)
tab_k = DataFrame(N = [Nk..., N], khi2_reel = [[khi2_moyen(n; reel = true) for n in Nk]..., khi2(y_reel)],
                  khi2_simule = [khi2_moyen(n; reel = false) for n in (Nk..., N)])
CSV.write(joinpath(dossier, "khi2_selon_N.csv"), tab_k)
for l in eachrow(tab_k)
    @printf("N = %8d : χ²ᵣ moyen réel (aminci) %.3f, simulé %.3f\n", l.N, l.khi2_reel, l.khi2_simule)
end
tc = ((0:inst.n_bins-1) .+ 0.5) .* Δ
CSV.write(joinpath(dossier, "reel_simule.csv"),
    DataFrame(t_ns = tc[W], reel = y_reel[W], simule = y_sim[W], modele = μN[W], pearson_reel = r_reel, pearson_simule = r_sim))

# -----------------------------------------------------------------------------
# 4. (a) aller-retour sur α
# -----------------------------------------------------------------------------
# Proportions de diffusion et de fond fixées à celles de la mesure ; N = photons bleus de la fenêtre.
# Récupération : Poisson, τ₂, τ₃, t₀ et IRF connus ; amplitudes (diffusion, τ₂, τ₃) libres ; fond connu.
Iw(b) = sum(@view b[W])
part_d = aB[1] * Iw(Bd) / Iw(μ_sim)
part_f = Iw(fond) / Iw(μ_sim)
αs = collect(0.05:0.10:0.95)
Ns = (1_000, 10_000, 100_000)
function attendu(α, Nt)
    d = α .* B2 .+ (1 - α) .* B3
    d .*= (1 - part_d - part_f) / Iw(d)
    μ = d .+ part_d .* Bd ./ Iw(Bd) .+ part_f .* fond ./ Iw(fond)
    S.masquer!(μ, D.FENETRE)
    return μ .* Nt
end
function recuperer(y, Nt)
    f = fond[W] .* (part_f * Nt / Iw(fond))
    X = hcat(Bd[W], B2[W], B3[W])
    yW = y[W]
    function modele!(μ, J, θ)
        μ .= X * θ .+ f
        J .= X
    end
    θ0 = max.(X \ max.(yW .- f, 0.0), 1e-6 .* Nt)
    r = S.lm_poisson(modele!, θ0, zeros(3), fill(Inf, 3), yW; maxiter = 200, tol_dec = 1e-10)
    return r.θ[2] / (r.θ[2] + r.θ[3])
end
lignes = DataFrame(N = Int[], alpha_injecte = Float64[], moyenne = Float64[], ecart_type = Float64[],
                   q025 = Float64[], q975 = Float64[], biais = Float64[])
rng = StableRNG(1919)
for Nt in Ns, α in αs
    μ = attendu(α, Nt)
    est = Float64[]
    y = zeros(length(μ))
    for _ in 1:ntir
        S.tirer!(y, μ, rng)
        push!(est, recuperer(y, Nt))
    end
    push!(lignes, (Nt, α, mean(est), std(est), quantile(est, 0.025), quantile(est, 0.975), mean(est) - α))
end
CSV.write(joinpath(dossier, "aller_retour.csv"), lignes)
for Nt in Ns
    l = lignes[lignes.N .== Nt, :]
    @printf("N = %6d : |biais| max %.4f, écart-type %.3f–%.3f\n", Nt, maximum(abs.(l.biais)), minimum(l.ecart_type), maximum(l.ecart_type))
end

# -----------------------------------------------------------------------------
# 5. Figure
# -----------------------------------------------------------------------------
include(joinpath(@__DIR__, "src", "Figures.jl"))
using CairoMakie
const F = Figures
set_theme!(F.theme_sim())
virgs(v) = [F.virg(@sprintf("%g", x)) for x in v]
milliers(n) = replace(string(n), r"(?<=\d)(?=(\d{3})+$)" => " ")
puiss(n) = n == 1_000 ? "10³" : n == 10_000 ? "10⁴" : n == 100_000 ? "10⁵" : milliers(n)
fig = Figure(size = (1250, 620))
Label(fig[0, 1:2], "Validation du simulateur · SCLM (Cerulean + Topaz) · $(λ) nm · canal bleu (IRF mesurée IRF_$(λ)_16X)";
    halign = :left, fontsize = 16, color = F.ENCRE)
coulN = Dict(1_000 => F.COULEUR["2"], 10_000 => F.COULEUR["1"], 100_000 => F.COULEUR["3"])

ga = fig[1, 1] = GridLayout()
axa = Axis(ga[1, 1]; xlabel = "α injecté (fraction d'amplitude FRET du déclin bleu)", ylabel = "α récupéré",
    xticks = 0:0.2:1, yticks = 0:0.2:1, xtickformat = virgs, ytickformat = virgs,
    title = "(a) Aller-retour : α injecté contre α récupéré",
    subtitle = F.virg(@sprintf("τ₂ = %.3f ns, τ₃ = %.3f ns, IRF et t₀ connus ; moyenne ± écart-type, %d tirages", mg.τ2, mg.τ3, ntir)))
lines!(axa, [0, 1], [0, 1]; color = F.MUET, linestyle = :dash, linewidth = 1.2, label = "diagonale")
dx = Dict(1_000 => -0.012, 10_000 => 0.0, 100_000 => 0.012)
for Nt in Ns
    l = lignes[lignes.N .== Nt, :]
    errorbars!(axa, l.alpha_injecte .+ dx[Nt], l.moyenne, l.ecart_type; color = coulN[Nt], whiskerwidth = 6)
    scatter!(axa, l.alpha_injecte .+ dx[Nt], l.moyenne; color = coulN[Nt], markersize = 8,
        label = "N = $(puiss(Nt)) photons")
end
limits!(axa, 0, 1, -0.05, 1.05)
axislegend(axa; position = :lt, framevisible = false)

gb = fig[1, 2] = GridLayout()
yt = filter(v -> minimum(μN[W]) / 1.5 <= v <= 1.5maximum(μN[W]), [200, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000])
axb = Axis(gb[1, 1]; yscale = log10, ylabel = "photons / canal", yticks = (yt, milliers.(yt)),
    title = F.virg("(b) Réel contre simulé : $(F.virg(@sprintf("%g", cl_b))) mM, N = $(milliers(N)) photons bleus"),
    subtitle = "mêmes τ, IRF, t₀, amplitudes et N ; trait : modèle commun")
ok = y_reel[W] .> 0
lines!(axb, tc[W], μN[W]; color = (F.MUET, 0.9), linewidth = 1.0, label = "modèle")
oks = y_sim[W] .> 0
scatter!(axb, tc[W][oks], y_sim[W][oks]; color = (F.COULEUR["2"], 0.45), markersize = 9, marker = :rect, label = "simulé")
scatter!(axb, tc[W][ok], y_reel[W][ok]; color = F.ENCRE, markersize = 4, label = "réel (somme de 10 × 1 s)")
axislegend(axb; position = :rt, framevisible = false)
axr = Axis(gb[2, 1]; xlabel = "temps du TAC (ns)", ylabel = "résidu de Pearson", xtickformat = virgs, ytickformat = virgs,
    subtitle = F.virg(@sprintf("χ²ᵣ réel %.2f · χ²ᵣ simulé %.2f (ce tirage ; moyenne de 200 tirages %.2f)",
        χ2r(r_reel), χ2r(r_sim), tab_k.khi2_simule[end])))
hlines!(axr, [0.0]; color = F.MUET, linewidth = 1)
scatter!(axr, tc[W], r_sim; color = (F.COULEUR["2"], 0.6), markersize = 6, marker = :rect)
scatter!(axr, tc[W], r_reel; color = F.ENCRE, markersize = 4)
linkxaxes!(axb, axr)
hidexdecorations!(axb; grid = false)
rowsize!(gb, 2, Relative(0.32))
rowgap!(gb, 6)
note = "(a) points décalés horizontalement selon N. (b) Aux N des simulations, χ²ᵣ moyen (200 tirages) du réel aminci / du simulé contre le même modèle : " *
    join([F.virg(@sprintf("%s photons %.2f / %.2f", puiss(l.N), l.khi2_reel, l.khi2_simule)) for l in eachrow(tab_k) if l.N != N], " · ")
Label(fig[2, 1:2], note; halign = :left, fontsize = 12, color = F.ENCRE2)
save(joinpath(dossier, "validation_simulateur.png"), fig; px_per_unit = 2)
save(joinpath(dossier, "validation_simulateur.pdf"), fig)
println("-> ", joinpath(dossier, "validation_simulateur.png"))
