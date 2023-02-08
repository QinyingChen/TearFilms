using Test
using LinearAlgebra

include("spectral.jl")

@testset "fderiv" begin
    x = trig(44)[1]
    u = @. exp(sin(2x) + cos(x))
    du = @. (2cos(2x) - sin(x)) * u
    ddu = @. (2cos(2x) - sin(x)) * du + (-4sin(2x) - cos(x)) * u
    @test fderiv(u) ≈ du
    @test fderiv(u, 2) ≈ ddu

    F, Finv, mult = plan_fderiv(44)
    @test fderiv(u, F, Finv, mult) ≈ du
    F, Finv, mult = plan_fderiv(44, 2)
    @test fderiv(u, F, Finv, mult) ≈ ddu
end


