using Plots
using LinearAlgebra



include("Qreservoir.jl")
figpath = "C:/Users/julie/Downloads/"
function createDataset(nb_exemples,resolution)
    data = rand([0,1],nb_exemples)
    time_function = LinRange(0.,2*pi,resolution)
    square_ar = vcat([1 for i in 1:resolution÷2],[-1 for i in 1:(resolution - resolution÷2)])
    sin_ar = sin.(time_function)

    dataset = Float64[]
    target = Float64[]
    for i in data
        if i == 0
            dataset = vcat(dataset,square_ar)
            target = vcat(target,zeros(resolution))
        else
            dataset = vcat(dataset,sin_ar)
            target = vcat(target,ones(resolution))
        end
    end

    training_x = dataset[1:length(dataset)÷2]
    training_y = target[1:length(target)÷2]
    test_x = dataset[length(dataset)÷2:end]
    test_y = target[length(target)÷2:end]

    return training_x,training_y,test_x,test_y
end

training_x,training_y,test_x,test_y = createDataset(100,20)

#gr()
#p = plot(training_x[1:75],
#seriestype=:scatter,markersize = 8,
#size=(1132,700),tickfontsize = 12,legend = false,fmt=:pdf)
#display(p)
#savefig(figpath*"sinesquare.pdf")


# time_ar = 1:length(training_x)
# plotly()
# display(plot(time_ar,training_x))
## model
time_resolution = 300
time_interval = 100  # ns
offset = 1
sampling = 6

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
    global predicted_test = X_test * W

    println("testing complete")
    ## Compute and print the accuracy
    global correct0 = 0
    for m in 1:length(predicted_test)
        global correct0 += ((predicted_test[m]>0.5) == test_y[m])
    end
    global correct0=correct0/size(test_y)[1]
    global final_plot = hcat(predicted_test,test_x,test_y)
    println("accuracy= " ,correct0)

    index_min = 400
    index_max = 600

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

    println("end of calculation")
end
