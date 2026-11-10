"""
    DeuxCanaux

Analyse FLIM à deux canaux d'une titration de banc (document « Analyse FLIM à deux
canaux : modifications du code », 2026-10-16) : sept modules, une responsabilité chacun,
une sortie persistée après chaque étape. Utilise les briques du simulateur (convolution
avec l'IRF, Levenberg-Marquardt de Poisson, ton `vec_to_lifetime`, lecteur .sdt).
"""
module DeuxCanaux

using ..SimFLIMFRET
const S = SimFLIMFRET
using LinearAlgebra, Statistics, Random, Printf, CSV, DataFrames, StableRNGs, Distributions

include("Charger.jl")
include("Global.jl")
include("Lectures.jl")
include("Calibrer.jl")
include("Amincir.jl")
include("Blanchir.jl")
include("Agreger.jl")

end
