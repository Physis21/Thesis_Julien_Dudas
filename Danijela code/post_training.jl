using Plots
using LinearAlgebra
include("sin_square_functions.jl")
using JLD2


figpath = "/Users/markovic/Documents/VSC Julia/sin_square/"


function compute_accuracy(predicted_test, test_y)
    correct0 = 0
    error = 0
    for m in 3:length(predicted_test)
        correct0 += ((predicted_test[m]>0.5) == test_y[m])
        error += (predicted_test[m] - test_y[m])^2

    end
    accuracy=correct0/size(test_y)[1]
    rmse = sqrt(error/size(test_y)[1])
    
    println("accuracy= " , accuracy)
    print("rmse=", rmse)
    return(accuracy)
end
wA = 10
wB=9.5
kappaA=17.0
kappaB=21.0
g = 700e6
meas_max=2
sampling = 1
gr()

figure_title = string("wA =", wA/1e9, "GHz, wB = ", wB/1e9, "GHz, g =", g/1e6, "MHz,
κA =", kappaA/1e6, "MHz, κB =", kappaB/1e6, "MHz, sampling=", sampling, "_meas_max=", meas_max )
filename = string("sin_square_features_coupling=", g, "MHz_kappas=", kappaA,"_", kappaB, "MHz_mesmax=", meas_max , "_sampling=", sampling, ".jld2")

figname = string("sin_square_coupling=", g, "MHz_kappas=", kappaA,"_", kappaB, "MHz_mesmax=", meas_max , "_sampling=", sampling, ".pdf")

time_plot = load(filename, "time_plot")
X = load(filename, "X")
X_test = load(filename, "X_test")
Y = load(filename, "Y")
Y_test = load(filename, "Y_test")

W = pinv(X) * Y
   
predicted_test = X_test * W

compute_accuracy(predicted_test, Y_test)

final_plot=hcat(predicted_test, Y_test)

index_min = 1
index_max = 81

index_min =25
index_max = 128

gr()
p = plot(1e6*time_plot[index_min-24:index_max-24],
final_plot[index_min:index_max,:],
size=(1132,400),
margin = 5Plots.mm,
tickfontsize = 14,legendfontsize=12,fmt=:pdf,
label=["prediction" "data" "target"],linewidth = 3,
legend=:right)
xlabel!(p,"Time (us)")
display(p)
savefig("performance_9_neurons.pdf")






