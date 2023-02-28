using Plots
using LinearAlgebra
using QuantumOptics
using Printf
using BenchmarkTools

Ndim = 7
basis_a = FockBasis(Ndim)
basis_b = FockBasis(Ndim)
basis = basis_a ⊗ basis_b
a = embed(basis,1,destroy(basis_a))
b = embed(basis,2,destroy(basis_b))
a_d = embed(basis,1,create(basis_a))
b_d = embed(basis,2,create(basis_b))
Na = embed(basis,1,number(basis_a))
Nb = embed(basis,2,number(basis_b))
Xa = ((a_d+a)/sqrt(2))
Pa = (im*(a_d-a)/sqrt(2))
Xb = ((b_d+b)/sqrt(2))
Pb = (im*(b_d-b)/sqrt(2))

vacuum = fockstate(basis_a,0)⊗fockstate(basis_b,0)



#some useful testing functions

T = LinRange(0,2*pi,20)
sinT = sin.(T)
square = vcat([1 for i in 1:10],[-1 for i in 1:10])
squaresin = vcat(square,sinT,square,sinT)




function Qmodel(wA, wB, g, kappaA, kappaB, eA, eB, training_x, time_interval, resolution, sampling, offset=1)
    train_length = length(training_x)
    
    eA_ar = eA.*training_x.+ offset
    eB_ar = eB.*training_x.+ offset
    
    H_0 = wA*a_d*a + wB*b_d*b + g*(a_d*b + a*b_d)
    
    c_ops = [a,b]
    
    ψ = vacuum
    reservoir_output = zeros(train_length, sampling*(meas_max+1)^2 )
    input_times = zeros(train_length, resolution)

    for i in 1:train_length
        for j in 1:resolution
            input_times[i,j] = ((i-1)*time_interval + j*(time_interval/resolution))*(10^-9)
        end
        H_drive= im*sqrt(2*kappaA)*eA_ar[i]*(a-a_d) + im*sqrt(2*kappaB)*eB_ar[i]*(b-b_d)
        H_full = H_0+ H_drive

        tout,states = timeevolution.master(input_times[i,:],ψ,H_full,c_ops,rates=[kappaA,kappaB])
        temp_ar = Float32[]
        for k in 1:sampling
            rho_mes = states[k*(resolution÷sampling)]
            normalize!(rho_mes)
            for i1 in 0:meas_max
                for i2 in 0:meas_max
                    # !!!
                    temp = expect(rho_mes,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))
                    push!(temp_ar,abs(temp)^2)
                end
            end
        end
        reservoir_output[i,:] = temp_ar
        ψ = states[end]
        if i%100 == 0
            println("100 simulations passed, i = ", i)
        end
    end
    input_times = zeros(train_length)
    for i in 1:train_length
        input_times[i] = (i*time_interval)*(10^-9)
    end
    return input_times,reservoir_output
end





function createDataset(nb_exemples,resolution)
    data = rand([0,1],nb_exemples)
    time_function = LinRange(0.,2*pi,resolution+1)
    square_ar = vcat([1 for i in 1:resolution÷2],[-1 for i in 1:(resolution - resolution÷2)])
    sin_ar = sin.(time_function[1:8])

    dataset = Float32[]
    target = Float32[]
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


function classification_task(wA, wB, g, kappaA, kappaB, eA, eB, training_x, training_y, test_x, test_y, time_interval,time_resolution,sampling,offset)
    println("calculation initialized")

    time_plot, X = Qmodel(wA, wB, g, kappaA, kappaB, eA, eB,training_x, time_interval, time_resolution, sampling, offset)
    Y = training_y
    # Train the model by Moore-Penrose pseudoinversion.
    W = pinv(X) * Y
    # Evaluate the model on the test set
    # We pass the latest training state in order to avoid the need for another washout
    println("training complete")
    ##model
    time_plot, X_test = Qmodel(wA, wB, g, kappaA, kappaB, eA, eB,test_x,time_interval,time_resolution,sampling,offset)
    predicted_test = X_test * W

    println("testing complete")
   
    return(time_plot, X, X_test, Y, test_y)
end


