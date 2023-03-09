using Plots
using LinearAlgebra
using QuantumOptics



include("Qreservoir.jl")
figpath = "C:/Users/julie/Downloads/"

ψtest = coherentstate(basis_a, 1) ⊗ coherentstate(basis_b, 2)
# ψtest = fockstate(basis_a, 0) ⊗ fockstate(basis_b, 2)

test_ρ = tensor( ψtest, dagger( ψtest ))
println("printing test_ρ matrix , size is ", size(test_ρ))
for i1 in 0:meas_max
    for i2 in 0:meas_max
        print(round(real(expect(test_ρ,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))), digits = 5), " ")
    end
    println(" ")
end

# println(test_ρ)
# test_result, possible_shot_results = Qmeasure_shot(test_ρ, 3)
# for i in 1:length(test_result)
#     println(test_result[i])
# end
println("Qmeasure_shot test")
Qmeasure_shot_error(test_ρ, (Ndim+1)^2, (Ndim+1)^4, 1)
@time Qmeasure_shot_error_theory(test_ρ, (Ndim+1)^2, 1)
