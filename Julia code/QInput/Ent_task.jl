using Plots
using LinearAlgebra
include("QreservoirE.jl")

L = 50

function createDataset(L)
    train_r, test_r = 2*rand(Float32,L,1), 2*rand(Float32,L,1)
    train_ϕ, test_ϕ = 2*π*rand(Float32,L,1), 2*π*rand(Float32,L,1)
    train_Y, test_Y = [Int(r>1) for r in train_r], [Int(r>1) for r in test_r]
    return train_r, test_r, train_ϕ, test_ϕ, train_Y, test_Y
end
train_r, test_r, train_ϕ, test_ϕ, train_Y, test_Y = createDataset(L)

time_resolution = 300
time_interval = 100  # ns
sampling = 10

function classification_task()
    println("calculation initialized")

    global X = Qmodel(train_r,train_ϕ,time_interval,time_resolution,sampling)
    global Y = train_Y
    # Train the model by Moore-Penrose pseudoinversion.
    global W = pinv(X) * Y
    # Evaluate the model on the test set
    # We pass the latest training state in order to avoid the need for another washout
    println("training complete")
    ##model
    global X_test = Qmodel(test_r,test_ϕ,time_interval,time_resolution,sampling)
    global predicted_test = X_test * W

    println("testing complete")
    ## Compute and print the accuracy
    global deviation = 0
    global entangled, unentangled  = [], []
    for m in 1:L
        global deviation += (predicted_test[m]-test_Y[m])^2
        if predicted_test[m]>1
            push!(entangled, test_r[m]*exp(im*test_ϕ[m]))
        else
            push!(unentangled, test_r[m]*exp(im*test_ϕ[m]))
        end
    end
    global deviation = sqrt(deviation/L)
    println("error= ", deviation, "%")

    ## plot the classification accuracy
    plotly()
    x_circle = [cos(θ) for θ in 0.0:2*π/100:2*π+0.001]
    y_circle = [sin(θ) for θ in 0.0:2*π/100:2*π+0.001]
    p = plot(real(entangled),imag(entangled),seriestype = :scatter,
    size=(1132,700),title=figure_title,titlefontsize=12,aspect_ratio=:equal)
    plot!(p,real(unentangled),imag(unentangled),seriestype = :scatter, legend=:right)
    plot!(p, x_circle, y_circle)
    xlabel!(p,"Re(r)")
    ylabel!(p,"Im(r)")
    display(p)

    println("end of calculation")
end
