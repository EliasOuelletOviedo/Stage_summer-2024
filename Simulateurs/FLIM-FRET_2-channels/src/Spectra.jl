# Spectra.jl — spectres FPbase, matrice de recouvrement spectral η, cache CSV.
# Exécuté une fois (step0_spectres.jl) ; les autres étapes lisent le cache.

url_proteine(slug) = "https://www.fpbase.org/api/proteins/$slug/?format=json"
url_emission(slug) = "https://www.fpbase.org/api/proteins/$slug/?include_fields=states__em_spectrum&format=json"
const FPBASE_GRAPHQL = "https://www.fpbase.org/graphql/"

function _telecharger_json(url::AbstractString; corps = nothing)
    io = IOBuffer()
    if corps === nothing
        Downloads.download(url, io; timeout = 60)
    else
        Downloads.request(url; method = "POST", output = io, input = IOBuffer(corps),
            headers = ["Content-Type" => "application/json"], timeout = 60)
    end
    return JSON.parse(String(take!(io)))
end

"""Rendement quantique, durée de vie et maxima d'une protéine FPbase (état par défaut)."""
function proprietes_fpbase(slug::AbstractString)
    d = _telecharger_json(url_proteine(slug))
    s = d["states"][1]
    return Dict{String,Any}("nom" => d["name"], "slug" => slug, "qy" => s["qy"],
        "lifetime" => something(s["lifetime"], NaN), "em_max" => s["em_max"], "ex_max" => s["ex_max"])
end

"""Spectre d'émission FPbase : paires (λ en nm, valeur normalisée au maximum)."""
function emission_fpbase(slug::AbstractString)
    d = _telecharger_json(url_emission(slug))
    sp = d["states"][1]["em_spectrum"]
    return Float64[p[1] for p in sp], Float64[p[2] for p in sp]
end

"""Courbe de transmission d'un filtre FPbase (identifiant de spectre, GraphQL)."""
function filtre_fpbase(id::Integer)
    q = JSON.json(Dict("query" => "{ spectrum(id: $id) { owner { name } data } }"))
    d = _telecharger_json(FPBASE_GRAPHQL; corps = q)["data"]["spectrum"]
    return d["owner"]["name"], Float64[p[1] for p in d["data"]], Float64[p[2] for p in d["data"]]
end

function ecrire_spectre(chemin, λ, v; entete = "lambda_nm,valeur")
    mkpath(dirname(chemin))
    open(chemin, "w") do io
        println(io, entete)
        for (a, b) in zip(λ, v)
            println(io, a, ",", b)
        end
    end
end

function lire_spectre(chemin)
    lignes = readlines(chemin)[2:end]
    λ = [parse(Float64, split(l, ",")[1]) for l in lignes]
    v = [parse(Float64, split(l, ",")[2]) for l in lignes]
    return λ, v
end

"""Interpolation linéaire de (λ, v) sur une grille ; 0 hors du domaine mesuré."""
function interpoler(λ::Vector{Float64}, v::Vector{Float64}, grille::AbstractVector)
    out = zeros(length(grille))
    for (i, g) in enumerate(grille)
        (g < λ[1] || g > λ[end]) && continue
        j = searchsortedlast(λ, g)
        j >= length(λ) && (out[i] = v[end]; continue)
        w = (g - λ[j]) / (λ[j+1] - λ[j])
        out[i] = (1 - w) * v[j] + w * v[j+1]
    end
    return out
end

"""Filtre rectangulaire idéal (transmission 1 dans [a, b])."""
rectangle(a, b) = (λ -> a <= λ <= b ? 1.0 : 0.0)

"""
    calculer_eta(em_D, em_A, T_bleu, T_jaune; pas=0.1) -> SMatrix{2,2}

η[canal, émetteur] = ∫ S_e(λ)·T_c(λ) dλ / ∫ S_e(λ) dλ : fraction des photons
émis par chaque émetteur qui atteint chaque canal (détection plate supposée).
`em_*` : (λ, valeur) ; `T_*` : (λ, transmission) ou une fonction λ -> T.
Les spectres FPbase sont traités comme des spectres en photons.
"""
function calculer_eta(em_D, em_A, T_bleu, T_jaune; pas = 0.1)
    grille = 300.0:pas:800.0
    S = [interpoler(em_D..., grille), interpoler(em_A..., grille)]
    T = [t isa Function ? t.(grille) : clamp.(interpoler(t..., grille), 0, 1) for t in (T_bleu, T_jaune)]
    η = zeros(2, 2)
    for c in 1:2, e in 1:2
        η[c, e] = sum(S[e] .* T[c]) / sum(S[e])
    end
    return SMatrix{2,2}(η)
end

"""
    fuite_donneur_jaune(η, x, p, o)

Part du canal jaune due à la fuite du donneur (η_YD·I_D)/I_jaune, à comparer
au contrôle « ≈ 24 % » du plan.
"""
function fuite_donneur_jaune(x, p::SensorParams, o::OpticsParams)
    _, esY = expsums(x, p, o)
    I_D = integrale_donneur(x, p)
    return o.η[2, 1] * I_D / integrale(esY)
end

"""
    preparer_spectres(dossier; donneur, accepteur, filtre_bleu, filtre_jaune)

Télécharge (ou relit du cache) les spectres et propriétés, écrit `eta.toml`
avec la provenance de chaque valeur. Renvoie le dictionnaire écrit.
"""
function preparer_spectres(dossier::AbstractString; donneur = "mcerulean3", accepteur = "topaz",
                           filtre_bleu = 4887, filtre_jaune = 7399, autres_donneurs = ["cerulean"])
    mkpath(dossier)
    cache(nom) = joinpath(dossier, nom * ".csv")
    function emission(slug)
        f = cache("emission_" * slug)
        isfile(f) && return lire_spectre(f)
        λ, v = emission_fpbase(slug)
        ecrire_spectre(f, λ, v)
        return λ, v
    end
    function filtre(id)
        f = cache("filtre_$(id)")
        if isfile(f)
            nom = strip(readline(f)[3:end])
            l = readlines(f)[3:end]
            return nom, [parse(Float64, split(x, ",")[1]) for x in l], [parse(Float64, split(x, ",")[2]) for x in l]
        end
        nom, λ, v = filtre_fpbase(id)
        mkpath(dirname(f))
        open(f, "w") do io
            println(io, "# ", nom)
            println(io, "lambda_nm,transmission")
            foreach(((a, b),) -> println(io, a, ",", b), zip(λ, v))
        end
        return nom, λ, v
    end
    props = Dict(slug => proprietes_fpbase(slug) for slug in unique([donneur, accepteur, autres_donneurs...]))
    emD, emA = emission(donneur), emission(accepteur)
    nomB, λB, TB = filtre(filtre_bleu)
    nomY, λY, TY = filtre(filtre_jaune)
    η = calculer_eta(emD, emA, (λB, TB), (λY, TY))
    η_plan = calculer_eta(emD, emA, rectangle(460, 500), rectangle(520, 560))
    η_autres = Dict(s => calculer_eta(emission(s), emA, (λB, TB), (λY, TY)) for s in autres_donneurs)

    res = Dict{String,Any}(
        "source" => "FPbase (https://www.fpbase.org), téléchargé le $(Dates.today())",
        "donneur" => donneur, "accepteur" => accepteur,
        "filtre_bleu" => "$nomB (spectre FPbase $filtre_bleu)", "filtre_jaune" => "$nomY (spectre FPbase $filtre_jaune)",
        "hypotheses" => "détection plate ; spectres FPbase traités comme des spectres en photons",
        "eta_BD" => η[1, 1], "eta_BA" => η[1, 2], "eta_YD" => η[2, 1], "eta_YA" => η[2, 2],
        "eta_plan_rectangles_460_500_520_560" => Dict("eta_BD" => η_plan[1, 1], "eta_BA" => η_plan[1, 2],
                                                      "eta_YD" => η_plan[2, 1], "eta_YA" => η_plan[2, 2]),
        "proprietes" => props,
    )
    for (s, e) in η_autres
        res["eta_donneur_$(s)"] = Dict("eta_BD" => e[1, 1], "eta_YD" => e[2, 1])
    end
    open(io -> TOML.print(io, res; sorted = true), joinpath(dossier, "eta.toml"), "w")
    return res
end

"""Matrice η lue dans `eta.toml` (écrit par `preparer_spectres`)."""
function lire_eta(chemin::AbstractString)
    d = TOML.parsefile(chemin)
    return SMatrix{2,2}([d["eta_BD"] d["eta_BA"]; d["eta_YD"] d["eta_YA"]]), d
end
