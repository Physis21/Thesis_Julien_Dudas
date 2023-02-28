using Plots
using LinearAlgebra
include("Qreservoir.jl")
# watch out, DifferentialEquations also uses ⊗ which induces errors
using DifferentialEquations



function createDataset(dataset_start=100,dataset_span=1000,delay = 2)
    β0 = 0.2
    θ = 1
    n = 10
    τ = 17
    γ = 0.1
    h(p,t) = 1
    f(u,h,p,t) = -u*γ + β0*(θ^n)*h(p,t-τ)/((θ^n)+(h(p,t-τ)^n))
    u0 = 0
    tspan = (0,2*dataset_span+2*dataset_start+3)
    prob = DDEProblem(f,u0,h,tspan,saveat=1)
    dataset = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

    training_x = Float64[]
    training_y = Float64[]
    test_x = Float64[]
    test_y = Float64[]
    for i in 1:dataset_span
        push!(training_x,dataset[dataset_start-delay+i])
        push!(training_y,dataset[dataset_start+i])
        push!(test_x,dataset[dataset_span+2*dataset_start-delay+i])
        push!(test_y,dataset[dataset_span+2*dataset_start+i])
    end
    return training_x,training_y,test_x,test_y

end

training_x,training_y,test_x,test_y = createDataset(100,1000,10)

# time_ar = 1:length(training_x)
# plotly()
# display(plot(time_ar,training_x))
## model
time_resolution = 300
time_interval = 100  # ns
sampling = 6

function classification_task()
    println("calculation initialized")

    global time_plot, X = Qmodel(training_x,time_interval,time_resolution,sampling)
    global Y = training_y
    # Train the model by Moore-Penrose pseudoinversion.
    global W = pinv(X) * Y
    # Evaluate the model on the test set
    # We pass the latest training state in order to avoid the need for another washout
    println("training complete")
    ##model
    global time_plot, X_test = Qmodel(test_x,time_interval,time_resolution,sampling)
    global predicted_test = X_test * W

    println("testing complete")
    ## Compute and print the mean square root error
    global error0 = 0
    for m in 30:length(predicted_test)
        error0 += (predicted_test[m]-test_y[m])^2
    end
    error0=sqrt(error0/length(predicted_test))
    println("prediction error: " ,error0)

    # list_accuracies.append(accuracy)

    index_min = 220
    index_max = 500
    gr()
    final_plot = hcat(predicted_test,test_y)

    p = plot(time_plot[index_min:index_max],final_plot[index_min:index_max,:],
    size=(1132,700),title=figure_title,label=["prediction" "target"],
    legend = :best,titlefontsize=12)

    display(p)
    savefig(figpath*"test.pdf")
    println("end of calculation")
end
