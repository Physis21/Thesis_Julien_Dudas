using Plots
using LinearAlgebra
using QuantumOptics
using Printf
using BenchmarkTools
using JLD2
include("sin_square_functions.jl")

#PARAMETERS

meas_max = 2
wA = 10.0e9
wB = 9.5e9
g = 700.0e6

kappaA = 17.0e6
kappaB = 21.0e6
eA = 20e5
eB = 20e5

training_x,training_y,test_x,test_y = createDataset(100, 8)

# FIGURE DATASET
gr()
p = plot(training_x[1:75],
seriestype=:scatter,markersize = 8,
size=(1132,700), tickfontsize = 12, legend = false, fmt=:pdf)
display(p)
savefig("sinesquare_dataset.pdf")


time_ar = 1:length(training_x)
plotly()
display(plot(time_ar,training_x))
time_resolution = 300
time_interval = 100  # ns
offset = 1
sampling = 1

@time time_plot, X, X_test, Y, Y_test = classification_task(wA, wB, g, kappaA, kappaB, eA, eB, training_x, training_y, test_x, test_y, time_interval,time_resolution,sampling,offset)
filename = string("sin_square_features_coupling=", g, "MHz_kappas=", kappaA/1e6, "_", kappaB/1e6, "MHz_mesmax=", meas_max, "_sampling=", sampling, ".jld2")    
save(filename, "time_plot", time_plot, "X", X, "X_test", X_test, "Y", Y, "Y_test", Y_test) 