using Plots
using LinearAlgebra
include("sin_square_task.jl")



figpath = "C:/Users/julie/Downloads/"


function compute_accuracy(predicted_test, test_y)
    correct0 = 0
    error = 0
    for m in 3:length(predicted_test)
        correct0 += ((predicted_test[m]>0.5) == test_y[m])
        error += (predicted_test[m] - test_y[m])^2

    end
    accuracy=correct0/(size(test_y)[1] - 2)
    rmse = sqrt(error/(size(test_y)[1] - 2))
    
    println("accuracy= " , accuracy)
    println("rmse=", rmse)
    return(accuracy)
end
wA = 10
wB= 9
κA= 17e6
κB= 21e6
g = 700e6
meas_max=5
sampling = 1
gr()

figure_title = string("wA =", wA, "GHz, wB = ", wB, "GHz, g =", g/1e6, "MHz,
κA =", κA/1e6, "MHz, κB =", κB/1e6, "MHz, sampling=", sampling, "_meas_max=", meas_max )
filename = string("sin_square_features_eA=9e5_eB=5e5_coupling=", g/1e6, "MHz_kappas=", κA/1e6,"_", κB/1e6, "MHz_mesmax=", meas_max , "_sampling=", sampling, ".jld2")
# filename = string("sin_square_features_coupling=", g/1e6, "MHz_kappas=", κA/1e6,"_", κB/1e6, "MHz_mesmax=", meas_max , "_sampling=", sampling, ".jld2")

figname = string("sin_square_coupling=", g/1e6, "MHz_κs=", κA/1e6,"_", κB/1e6, "MHz_mesmax=", meas_max , "_sampling=", sampling, ".pdf")

time_plot = load(filename, "time_plot")
X = load(filename, "X")
X_test = load(filename, "X_test")
Y = load(filename, "Y")
Y_test = load(filename, "Y_test")

target = load(filename, "target")
# target is practically equal to Y_test for high enough neuron number
# target = [(y>0.5) for y in Y_test]
# W = pinv(X) * Y

# now that they are loaded, can change meas_max with same dataset to truncate higher level probs
meas_max_new = 1

X_new = zeros(size(X)[1], sampling*(meas_max_new+1)^2 )
X_test_new = zeros(size(X_test)[1], sampling*(meas_max_new+1)^2 )
for i1 in 1:size(X)[1]
    for i2 in 1:(sampling*(meas_max_new+1)^2)
        X_new[i1,i2] = X[i1,i2]
        X_test_new[i1,i2] = X_test[i1,i2]
    end

end

W = pinv(X_new) * Y

predicted_test = X_test_new * W

accuracy = compute_accuracy(predicted_test, target)

final_plot=hcat(predicted_test, Y_test, target)

index_min = 1
index_max = 81

index_min =25
index_max = 128

println("start plotting")
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

figname = string("sin_square_eA=9e5_eB=5e5_coupling=", g/1e6, "MHz_κs=", κA/1e6,"_", κB/1e6, "MHz_mesmaxnew=", meas_max_new , "_sampling=", sampling, "_accuracy=", accuracy, ".pdf")
savefig(figname)






