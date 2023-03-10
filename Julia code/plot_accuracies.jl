using Plots
using Printf

neuron_numbers = [4, 9, 16, 25]
accuracies_classical_static = [50, 59, 80.5, 92]
accuracies_classical_dynamic = [61, 85, 96.5, 99.5]
# accuracies_quantum = [68.59, 93.7, 100.0, 100.0]
# accuracies_quantum = [52.51, 100, 96.98, 100.0]
# accuracies_quantum = [70.92, 83.45, 94.49, 100.0]
# accuracies_quantum = [73.18, 94.49, 94.49, 99.74]
# accuracies_quantum = [84.21, 94.73, 94.73, 100]
accuracies_quantum = [68.17, 81.20, 77.94, 100]

accuracies = hcat(accuracies_classical_static, accuracies_classical_dynamic, accuracies_quantum)

# ! change file name manually
println("start plotting")
gr()
p = plot(neuron_numbers,
accuracies,
size=(700,400),
margin = 5Plots.mm,
shape = :circle,
tickfontsize = 10,legendfontsize=12,fmt=:pdf,
label=["classical static" "classical dynamic" "quantum basis states"],linewidth = 3,
legend=:bottomright)
xlabel!(p,"Number of neurons")
ylabel!("Accuracy")
display(p)

figname = string("accuracies_time_interval=50ns_eA=1.2e6_eB=1.2e6_plot.pdf")
savefig(figname)
