using Plots
using LinearAlgebra



include("Qreservoir.jl")
figpath = "C:/Users/julie/Downloads/"
pre_text = "half_sine_"

# dataset with arbitrary resolution

# function createDataset(nb_exemples,resolution)
#     data = rand([0,1],nb_exemples)
#     time_function = LinRange(0.,2*pi,resolution)
#     square_ar = vcat([1 for i in 1:resolution÷2],[-1 for i in 1:(resolution - resolution÷2)])
#     sin_ar = sin.(time_function)

#     dataset = Float64[]
#     target = Float64[]
#     for i in data
#         if i == 0
#             dataset = vcat(dataset,square_ar)
#             target = vcat(target,zeros(resolution))
#         else
#             dataset = vcat(dataset,sin_ar)
#             target = vcat(target,ones(resolution))
#         end
#     end

#     training_x = dataset[1:length(dataset)÷2]
#     training_y = target[1:length(target)÷2]
#     test_x = dataset[length(dataset)÷2:end]
#     test_y = target[length(target)÷2:end]

#     return training_x,training_y,test_x,test_y
# end

# dataset with resolution = 8, and sine values exact 0, 0.7071, 1 !
# divide values by 2
function createDataset(nb_exemples)
    data = rand([0,1],nb_exemples)
    time_function = LinRange(0.,2*pi,8)
    square_ar = vcat([1 for i in 1:4],[-1 for i in 1:4])/2
    sin_ar = [0,0.7071,1,0.7071,0,-0.7071,-1,-0.7071]/2

    dataset = Float64[]
    target = Float64[]
    for i in data
        if i == 0
            dataset = vcat(dataset,square_ar)
            target = vcat(target,zeros(8))
        else
            dataset = vcat(dataset,sin_ar)
            target = vcat(target,ones(8))
        end
    end

    training_x = dataset[1:length(dataset)÷2]
    training_y = target[1:length(target)÷2]
    test_x = dataset[length(dataset)÷2:end]
    test_y = target[length(target)÷2:end]

    return training_x,training_y,test_x,test_y
end
training_x,training_y,test_x,test_y = createDataset(100)
# training_x,training_y,test_x,test_y = createDataset(25)

# gr()
# p = plot(training_x[1:75],
# seriestype=:scatter,markersize = 8,
# size=(1132,700),tickfontsize = 12,legend = false,fmt=:pdf)
# display(p)
# savefig(figpath*"sinesquare.pdf")


# time_ar = 1:length(training_x)
# plotly()
# display(plot(time_ar,training_x))
## model
time_resolution = 300
time_interval = 100  # ns , before used time_interval = 100ns


offset = 0.5 # value of 1.1 so that input doesn't become negative nor 0 when sin(x) = -1
# sampling = 6
sampling = 1

function classification_task()
    println("calculation initialized")

    global time_plot, X = Qmodel(training_x,time_interval,time_resolution,sampling,offset)
    global Y = training_y
    # Train the model by Moore-Penrose pseudoinversion.
    global W = pinv(X) * Y
    # Evaluate the model on the test set
    # We pass the latest training state in order to avoid the need for another washout
    println("training complete")
    ##model
    global time_plot, X_test = Qmodel(test_x,time_interval,time_resolution,sampling,offset)
    global Y_test = X_test * W

    println("testing complete")
    ## Compute and print the accuracy
    global correct0 = 0
    for m in 1:length(Y_test)
        global correct0 += ((Y_test[m]>0.5) == test_y[m])
    end
    global correct0=correct0/size(test_y)[1]
    global final_plot = hcat(Y_test,test_x,test_y)
    println("accuracy= " ,correct0)

    # index_min = 400
    # index_max = 600

    index_min = 100
    index_max = 200
    gr()
    p = plot(1e9*time_plot[index_min:index_max],
    final_plot[index_min:index_max,:],
    size=(1132,700),title=figure_title,
    margin = 5Plots.mm,
    tickfontsize = 12,legendfontsize=12,fmt=:pdf,
    label=["prediction" "data" "target"],linewidth = 3,
    legend=:right,titlefontsize=12)
    xlabel!(p,"Time (nanoseconds)")
    display(p)
    savefig(p,figpath*"classification.pdf")

    #save JLD2 file

    parameters = Dict([("wA", wA),
    ("wB", wB),
    ("g", g),
    ("κA", κA),
    ("κB", κB),
    ("eA", eA),
    ("eB", eB),
    ("meas_max", meas_max),
    ("sampling", sampling),
    ("time_interval", time_interval)])

    filename = string(pre_text, "sin_square_features_eA=", eA, "_eB=" , eB,"_coupling=", g/1e6, "MHz_kappas=", κA/1e6, "_", κB/1e6, "MHz_mesmax=", meas_max, "_sampling=", sampling, ".jld2")    
    #target (test_y) is useful for post processing
    save(filename, "time_plot", time_plot, "X", X, "X_test", X_test, "Y", Y, "Y_test", Y_test, "target", test_y, "parameters", parameters)

    println("end of calculation")
end

function classification_task_shots( shot_nb=(Ndim+1)^4 )
    println("calculation initialized")

    global time_plot, X = Qmodel_shots(training_x,time_interval,time_resolution,sampling,shot_nb,offset)
    global Y = training_y
    # Train the model by Moore-Penrose pseudoinversion.
    global W = pinv(X) * Y
    # Evaluate the model on the test set
    # We pass the latest training state in order to avoid the need for another washout
    println("training complete")
    ##model
    global time_plot, X_test = Qmodel_shots(test_x,time_interval,time_resolution,sampling,shot_nb,offset)
    global Y_test = X_test * W
    println("testing complete")
    ## Compute and print the accuracy
    global correct0 = 0
    for m in 1:length(Y_test)
        global correct0 += ((Y_test[m]>0.5) == test_y[m])
    end
    global correct0=correct0/size(test_y)[1]
    global final_plot = hcat(Y_test,test_x,test_y)
    println("accuracy= " ,correct0)

    # index_min = 400
    # index_max = 600

    index_min = 20
    index_max = 80
    gr()
    p = plot(1e9*time_plot[index_min:index_max],
    final_plot[index_min:index_max,:],
    size=(1132,700),title=figure_title,
    margin = 5Plots.mm,
    tickfontsize = 12,legendfontsize=12,fmt=:pdf,
    label=["prediction" "data" "target"],linewidth = 3,
    legend=:right,titlefontsize=12)
    xlabel!(p,"Time (nanoseconds)")
    display(p)
    savefig(p,figpath*"classification_shots.pdf")

    #save JLD2 file

    filename = string(pre_text, "shots_sin_square_features_eA=", eA, "_eB=" , eB,"_coupling=", g/1e6, "MHz_kappas=", κA/1e6, "_", κB/1e6, "MHz_mesmax=", meas_max, "_sampling=", sampling, ".jld2")    
    #target (test_y) is useful for post processing

    parameters = Dict([("wA", wA),
    ("wB", wB),
    ("g", g),
    ("κA", κA),
    ("κB", κB),
    ("eA", eA),
    ("eB", eB),
    ("meas_max", meas_max),
    ("sampling", sampling),
    ("time_interval", time_interval),
    ("shot_nb", shot_nb)])

    save(filename, "time_plot", time_plot, "X", X, "X_test", X_test, "Y", Y, "Y_test", Y_test, "target", test_y, "parameters", parameters)

    println("end of calculation")
end

