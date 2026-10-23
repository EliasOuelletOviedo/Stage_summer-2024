# Étape 0 — spectres FPbase et matrice de recouvrement spectral η.
#
# Télécharge une fois les spectres d'émission (mCerulean3, Topaz, Cerulean)
# et les transmissions des filtres des deux canaux FLIM (Chroma D470/40m pour
# le bleu, Chroma ET530/30m pour le jaune), les met en cache dans
# sorties/spectres/ et écrit sorties/spectres/eta.toml. Relancer ne
# retélécharge rien : effacer le cache pour forcer.
#
#   julia --project=. step0_spectres.jl

include(joinpath(@__DIR__, "src", "SimFLIMFRET.jl"))
using .SimFLIMFRET
const S = SimFLIMFRET
using Printf

dossier = joinpath(@__DIR__, "sorties", "spectres")
res = S.preparer_spectres(dossier)

println("Filtres : bleu = ", res["filtre_bleu"], " ; jaune = ", res["filtre_jaune"])
@printf("η (vrais filtres)        BD = %.4f  BA = %.4f  YD = %.4f  YA = %.4f\n",
    res["eta_BD"], res["eta_BA"], res["eta_YD"], res["eta_YA"])
r = res["eta_plan_rectangles_460_500_520_560"]
@printf("η (rectangles du plan)   BD = %.4f  BA = %.4f  YD = %.4f  YA = %.4f\n",
    r["eta_BD"], r["eta_BA"], r["eta_YD"], r["eta_YA"])
for (k, v) in res["proprietes"]
    @printf("%-10s Φ = %.2f  τ = %s ns  em_max = %d nm\n", v["nom"], v["qy"], string(v["lifetime"]), v["em_max"])
end
