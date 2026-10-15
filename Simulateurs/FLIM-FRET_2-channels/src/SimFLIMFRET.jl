"""
    SimFLIMFRET

Simulateur FLIM-FRET à deux canaux (mCer3 → mTopaz, modèle M2 lent) : erreur
sur [Cl⁻] de six estimateurs en fonction du budget de photons et de la
concentration vraie, borne de Cramér-Rao comprise. Voir le README du dossier
et le plan de conception.
"""
module SimFLIMFRET

using LinearAlgebra, Statistics, Random, Printf, Dates, TOML, Downloads
using FFTW, StaticArrays, Distributions, StableRNGs, ForwardDiff
using CSV, DataFrames, JSON

include("Types.jl")
include("Photophysics.jl")
include("Instrument.jl")
include("Spectra.jl")
include("Generate.jl")
include("Fisher.jl")
include("Ajustement.jl")
include("Utilisateur.jl")
include("Ingest.jl")
include("Parametres.jl")
include("Estimators.jl")
include("Calibrate.jl")
include("Experiments.jl")

end
