using Plots
using Printf

neuron_numbers = [4, 9, 16, 25]
accuracies_classical_static = [50, 59, 80.5, 92]
accuracies_classical_dynamic = [61, 85, 96.5, 99.5]
accuracies_quantum = [68.59, 93.7, 100.0, 100.0]

accuracies = hcat(accuracies_classical_static, accuracies_classical_dynamic, accuracies_quantum)


println("start plotting")
gr()
p = plot(neuron_numbers,
accuracies,
size=(700,400),
margin = 5Plots.mm,
shape = :circle,
tickfontsize = 10,legendfontsize=12,fmt=:pdf,
label=["classical static" "classical dynamic" "quantum basis states"],linewidth = 3,
legend=:right)
xlabel!(p,"Number of neurons")
ylabel!("Accuracy")
display(p)

figname = string("accuracies_plot.pdf")
savefig(figname)
