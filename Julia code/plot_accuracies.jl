using Plots
using Printf

neuron_numbers = [4, 9, 16, 25]
neuron_numbers2 = 2*[4, 9, 16, 25]
accuracies_classical_static = [50, 59, 80.5, 92]
accuracies_classical_dynamic = [61, 85, 96.5, 99.5]
# accuracies_quantum = [68.59, 93.7, 100.0, 100.0] #eA = eB = 2e6
# accuracies_quantum = [52.51, 100, 96.98, 100.0] #eA = 5e6 , eB = 2e6
# accuracies_quantum = [70.92, 83.45, 94.49, 100.0] #eA = 9e5 , eB = 5e5
# accuracies_quantum = [73.18, 94.49, 94.49, 99.74] #eA = 9e5, eB = 9e5
# accuracies_quantum = [84.21, 94.73, 94.73, 100] #eA = 1.2e6 , eB = 1.2e6
# accuracies_quantum = [68.17, 81.20, 77.94, 100] #50ns normal kappa
# accuracies_quantum = [75.93, 99.74, 99.74, 99.74] #sampling = 2
# accuracies_quantum = [74.69, 83.71, 89.72, 87.47] #100ns half kappa
# accuracies_quantum = [88.24, 94.12, 94.12, 100] # dataset of 400
# accuracies_quantum = [77.69, 80.20, 94.73, 100] # half sine offset = 0.5
# accuracies_quantum = [78.19, 89.22, 99.74, 99.74] # half sine offset = 0.6

pre_text = " "

accuracies = hcat(accuracies_classical_static, accuracies_classical_dynamic, accuracies_quantum)
neuron_numbers_hcat = hcat(neuron_numbers,neuron_numbers,neuron_numbers2)
# ! change file name manually
println("start plotting")
gr()
p = plot(neuron_numbers,
accuracies,
size=(700,400),
title = "halved kappas",
margin = 5Plots.mm,
shape = :circle,
tickfontsize = 10,legendfontsize=12,fmt=:pdf,
label=["classical static" "classical dynamic" "quantum basis states"],linewidth = 3,
legend=:bottomright)
xlabel!(p,"Number of neurons")
ylabel!("Accuracy")
display(p)

figname = string(pre_text, "accuracies_eA=1.2e6_eB=1.2e6.pdf")
savefig(figname)
