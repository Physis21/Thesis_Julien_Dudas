using Plots
using LinearAlgebra
using QuantumOptics



include("Qreservoir.jl")
figpath = "C:/Users/julie/Downloads/"

ψtest = coherentstate(basis_a, 1) ⊗ coherentstate(basis_b, 1)
test_ρ = tensor( ψtest, dagger( ψtest ))
# println(test_ρ)
test_result = Qmeasure_shot(test_ρ)
println(test_result)
