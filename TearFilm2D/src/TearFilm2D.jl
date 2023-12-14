module TearFilm2D

using Parameters
using DifferentialEquations
using LinearAlgebra
using ProgressMeter
using Krylov, LinearSolve

export PhysicalConstants
include("constants.jl")

export Grid1D, Grid1DSym, radial_solve
include("solvers-1D.jl")

export twodim_solve
include("solvers-2D.jl")

end # module TearFilm2D
