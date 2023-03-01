using Plots
using LinearAlgebra



include("Qreservoir.jl")
figpath = "C:/Users/julie/Downloads/"

test_ρ = tensor( (fockstate(basis_a,1)+fockstate(basis_a,2))⊗fockstate(basis_b,1) / sqrt(2), dagger((fockstate(basis_a,1)+fockstate(basis_a,2))⊗fockstate(basis_b,1) / sqrt(2)))
# print(test_ρ)
test_result = Qmeasure_shot(test_ρ)
print(test_result)