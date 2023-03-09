using QuantumOptics
using Printf
using Plots
using BenchmarkTools
using Random
using JLD2
using LaTeXStrings

#usually keep Ndim = 7 for computing speed and accuracy of quantum simulations
Ndim = 8
basis_a = FockBasis(Ndim)
basis_b = FockBasis(Ndim)
basis = basis_a ⊗ basis_b
a = embed(basis,1,destroy(basis_a))
b = embed(basis,2,destroy(basis_b))
a_d = embed(basis,1,create(basis_a))
b_d = embed(basis,2,create(basis_b))
Na = embed(basis,1,number(basis_a))
Nb = embed(basis,2,number(basis_b))

vacuum = fockstate(basis_a,0)⊗fockstate(basis_b,0)

# meas_max = 7
# wA = 10e9
# wB = 9.5e9
# g = 7e8

# κA = 1e8
# κB = 1e8
# eA = 5e6
# eB = 1e6

meas_max = 5
wA = 10e9
wB = 9e9
g = 7e8

κA = 17e6
κB = 21e6
eA = 1.2e6
eB = 1.2e6

# !!!!! for these parameters and Ndim = 7, eA and eB shouldn't exceed 0.9e6
# !!!!! for these parameters and Ndim = 8, eA and eB shouldn't exceed 1.2e6
#some useful testing functions

T = LinRange(0,2*pi,20)
sinT = sin.(T)
square = vcat([1 for i in 1:10],[-1 for i in 1:10])
squaresin = vcat(square,sinT,square,sinT)


figpath = "C:/Users/julie/Downloads/"
#figure_title = "Reservoir simulation for: wA =  $(@sprintf("%.2e", wA)), wB =  $(@sprintf("%.2e", wB)), g =  $(@sprintf("%.2e", g))
#κA =  $(@sprintf("%.2e", κA)), κB =  $(@sprintf("%.2e", κB))
#eA =  $(@sprintf("%.2e", eA)), eB =  $(@sprintf("%.2e", eB))"

figure_title = string("Reservoir simulation for wA = ", wA, ", wB = ", wB, ", g = ", g, "
kA = ", κA, ", kB = ", κA,"
eA = ", eA,", eB = ", eB, ", Ndim = ", Ndim)

## Set Hamiltonian functions

# param[0]= Wa ; param[1]= Wb ; param[2]= gR
H_0(param) = param[1]*a_d*a + param[2]*b_d*b + param[3]*(a_d*b + a*b_d)

# param[0]= kappaA; param[1] = kappaB, param[2] = eA, param[3] = eB
H_drive(param) = im*sqrt(2*param[1])*param[3]*(a-a_d) + im*sqrt(2*param[2])*param[4]*(b-b_d)

param_0 = [wA,wB,g]
param_drive = [κA,κB,eA,eB]

H_full = H_0(param_0) + H_drive(param_drive)
c_ops = [a,b]

function constant_drive(show = 1,resolution = 10000, max_time = 50*(10^-9))
    tspan = LinRange(0,max_time,resolution)
    tout,states = timeevolution.master(tspan,vacuum,H_full,c_ops,rates=[κA,κB])
    Z1 = fockstate(basis_a,0)⊗fockstate(basis_b,0)⊗dagger(fockstate(basis_a,0)⊗fockstate(basis_b,0))
    Z2 = fockstate(basis_a,0)⊗fockstate(basis_b,1)⊗dagger(fockstate(basis_a,0)⊗fockstate(basis_b,1))
    #plot_means = hcat(real(expect(Z1,states)),
    #real(expect(Z2,states)))
    plot_means = hcat(real(expect(Na,states)),
    real(expect(Nb,states)))

    #plot means Na and Nb
    if show == 1
        gr()

        p = plot(tspan*1e9,plot_means,title=figure_title,size=(1132,700),
        legendfontsize=20,linewidth = 4,legend=:best,fmt=:pdf,
        margin=5Plots.mm,
        xguidefontsize=22,
        yguidefontsize=22,
        titlefontsize=20,tickfontsize = 15,label=["Na" "Nb"])
        xlabel!(p,"Time (nanoseconds)")
        ylabel!(p,"Average populations")
        display(p)
        #savefig(figpath*"constantdrive_high_g.pdf")
        savefig(figpath*"test_eA=1.2e6_eB=1.2e6_Ndim=8.pdf")
    end

end

function constant_drive_change(train=squaresin,time_interval=100,offset=1,show=1,resolution = 1000)
    #time_interval in nanoseconds
    train_length = length(train)
    eA_ar = [eA*(i+offset) for i in train]
    eB_ar = [eB*(i+offset) for i in train]
    ψ = vacuum
    reservoir_output = zeros(train_length*resolution,(meas_max+1)^2 )
    Na_mean = Float64[]
    Nb_mean = Float64[]
    input_times = zeros(train_length,resolution)
    for i in 1:train_length
        for j in 1:resolution
            input_times[i,j] = ((i-1)*time_interval + j*(time_interval/resolution))*(10^-9)
        end
        param_drive = [κA,κB,eA_ar[i],eB_ar[i]]
        H_full = H_0(param_0) + H_drive(param_drive)

        tout,states = timeevolution.master(input_times[i,:],ψ,H_full,c_ops,rates=[κA,κB])

        for k in 1:resolution
            temp_ar = Float64[]
            rho_mes = states[k]
            normalize!(rho_mes)
            push!(Na_mean,real(tr(rho_mes*Na)))
            push!(Nb_mean,real(tr(rho_mes*Nb)))
            for i1 in 0:meas_max
                for i2 in 0:meas_max
                    # !!!
                    temp = expect(rho_mes,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))
                    push!(temp_ar,real(temp))
                end
            end
            reservoir_output[(i-1)*resolution + k,:] = temp_ar
        end

        ψ = states[end]
#        if i%100 == 0
#            println("100 simulations passed, i = ", i)
#        end
    end
    input_times = zeros(train_length*resolution)
    for i in 1:train_length
        for j in 1:resolution
            input_times[(i-1)*resolution+j] = ((i-1)*time_interval+j*(time_interval/resolution))*(10^-9)
        end
    end
    train_plot = Float64[]
    for i in 1:train_length
        for j in 1:resolution
            push!(train_plot, train[i])
        end
    end
    #plot a layout
    if show ==1
        gr()
        p1 = plot(input_times,hcat(Na_mean,Nb_mean,train_plot),
        label=["Na_mean" "Nb_mean" "input"],
        fmt=:pdf,size=(1132,700),title=figure_title,titlefontsize=14,
        legendfontsize=14,margin = 5Plots.mm,legend=:right,tickfontsize = 12,
        xlabel="time",ylabel="mean populations",linewitdh = 5,
        xguidefontsize = 14, yguidefontsize =14)
        # push!(input_times,train_length*time_interval*(10^-9))
        # C = cgrad(:inferno, scale = :log)
        # p2 = heatmap(input_times,0:(meas_max+1)^2+1,transpose(reservoir_output),
        # xlabel="time",ylabel="occupation levels",c= cgrad(scale= :log))
        display(p1)
    end
    println("algorithm finished")
end

function Qmodel(train,time_interval,resolution,multiplexing,offset=0)
    train_length = length(train)
    eA_ar = [eA*(i+offset) for i in train]
    eB_ar = [eB*(i+offset) for i in train]
    ψ = vacuum
    reservoir_output = zeros(train_length,multiplexing*(meas_max+1)^2 )
    input_times = zeros(train_length,resolution)
    for i in 1:train_length
        for j in 1:resolution
            input_times[i,j] = ((i-1)*time_interval + j*(time_interval/resolution))*(10^-9)
        end
        param_drive = [κA,κB,eA_ar[i],eB_ar[i]]
        H_full = H_0(param_0) + H_drive(param_drive)

        tout,states = timeevolution.master(input_times[i,:],ψ,H_full,c_ops,rates=[κA,κB])
        temp_ar = Float64[]
        for k in 1:multiplexing
            rho_mes = states[k*(resolution÷multiplexing)]
            normalize!(rho_mes)
            for i1 in 0:meas_max
                for i2 in 0:meas_max
                    # !!!
                    temp = expect(rho_mes,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))
                    # push!(temp_ar,abs(temp)^2)   !!! previous error because temp is imaginary
                    push!(temp_ar,real(temp))
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

function Qmeasure_shot(ρ, shot_nb=1)
    #start temp_ar at 0 for later cumulative prob distribution
    result_ar = Tuple{Int64, Int64}[]
    temp_ar = [0.0]
    for i1 in 0:meas_max
        for i2 in 0:meas_max
            temp = expect(ρ, fockstate(basis_a,i1)⊗fockstate(basis_b,i2))
            push!(temp_ar,real(temp))
            push!(result_ar, (i1,i2))
            
        end
    end
    # println("creation of probability matrix")

    ## normalize!(temp_ar)
    ## printf("normalization finished")
    #we have probability distribution, now sum proba into cumulative distribution and locate a random [0,1) in the distribution to see where it lands
    temp_ar = accumulate(+, temp_ar)

    # println("probability distribution = ", temp_ar)

    # create a number of projections equal to shot_nb
    localize_ar = Float64[]
    shots = Tuple{Int64, Int64}[]
    for i in 1:shot_nb
        # !!
        localize = rand()
        push!(localize_ar,localize)
    end
    
    # println("localize rands = ", localize_ar)
    for i in 1:shot_nb
        result = (0,0)
        found = 0
        localize = localize_ar[i]
        for i in 1:(length(temp_ar)-1)
            if (localize >= temp_ar[i]) && (localize <= temp_ar[i+1])
                result = result_ar[i]
                found = 1
            end
        end
        #if no shot is measured, assume the projected energy state is beyond what we can measure, so return highest measurable by default
        if found == 1
            push!(shots, result)
        end
        
    end
    # println("result computed")
    return shots, result_ar
    
end

function Qmeasure_shot_mean(ρ, shot_nb=100, print_output = 0, heatmap_output = 0)
    
    # if all the shots are different, can't make error bars
    shots, shot_possibilities = Qmeasure_shot(ρ, shot_nb)
    shot_means = Float64[]
    # search probability of occupation for each energy level
    for shot1 in shot_possibilities
        apparition = 0
        # search frequency of apparition of each projected state
        for shot2 in shots
            if shot1 == shot2
                apparition += 1
            end
        end
        push!(shot_means,apparition/shot_nb)

    end

    # test to see the probabilities
    if print_output == 1

        for i in 1:length(shot_possibilities)
            print_output = string("state ", shot_possibilities[i], " has prob ", shot_means[i])
            println(print_output)
        end
        
    end

    if heatmap_output == 1

        println("print heatmap")
        means_histogram = reshape(shot_means, meas_max+1, meas_max+1)
        gr()
        axis_level = 0:meas_max
        h = heatmap(
        axis_level, axis_level, (x, y)->means_histogram[x+1, y+1], c=:viridis,
        nx=50, ny=50, 
        )
        display(h)

    end
    
    return shot_means
end

function Qmeasure_shot_error(ρ, shot_nb=100, N = (Ndim+1)^4, heatmap_output = 0)
    real_values = Float64[]

    for i1 in 0:meas_max
        for i2 in 0:meas_max
            push!(real_values, real(expect(test_ρ,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))) )
        end
    end

    var = zeros((meas_max+1)^2 )
    for i1 in 1:N

        shot_means = Qmeasure_shot_mean(ρ, shot_nb)
        for i2 in 1:length(shot_means)
            var[i2] += ((shot_means[i2] - real_values[i2])^2) / N
        end
        
    end

    diff = [sqrt(v) for v in var]
    
    if heatmap_output == 1

        println("print heatmap")
        diff_histogram = reshape(diff, meas_max+1, meas_max+1)

        gr()

        axis_level = 0:meas_max

        h = heatmap(
        axis_level, axis_level, (x, y)->log10.(diff_histogram[x+1, y+1]), c=:viridis,
        nx=50, ny=50, 
        xlabel="b levels", ylabel = "a levels", title = string("shot error for ",L"shot_{nb} = {%$shot_nb}", "and ", L"N = {%$N}")
        )
        

        display(h)
        savefig(figpath*"experiment.pdf")

    end
end

function Qmeasure_shot_error_theory(ρ, N = (Ndim+1)^4, heatmap_output = 0)
    real_values = Float64[]

    for i1 in 0:meas_max
        for i2 in 0:meas_max
            push!(real_values, real(expect(test_ρ,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))) )
        end
    end

    #theoretical value of error with multinomial law

    diff = [sqrt(p*(1-p)/N) for p in real_values]

    if heatmap_output == 1

        println("print heatmap")
        diff_histogram = reshape(diff, meas_max+1, meas_max+1)

        gr()

        axis_level = 0:meas_max

        h = heatmap(
        axis_level, axis_level, (x, y)->log10.(diff_histogram[x+1, y+1]), c=:viridis,
        nx=50, ny=50, 
        xlabel="b levels", ylabel = "a levels", title = string("theoretical shot error for ",L"N = {%$N}")
        )
        
        display(h)
        savefig(figpath*"theory.pdf")

    end

end

function Qmodel_shots(train,time_interval,resolution,multiplexing,offset=0, shot_nb = Ndim^4)
    train_length = length(train)
    eA_ar = [eA*(i+offset) for i in train]
    eB_ar = [eB*(i+offset) for i in train]
    ψ = vacuum
    reservoir_output = zeros(train_length,multiplexing*(meas_max+1)^2 )
    input_times = zeros(train_length,resolution)
    for i in 1:train_length
        for j in 1:resolution
            input_times[i,j] = ((i-1)*time_interval + j*(time_interval/resolution))*(10^-9)
        end
        param_drive = [κA,κB,eA_ar[i],eB_ar[i]]
        H_full = H_0(param_0) + H_drive(param_drive)

        tout,states = timeevolution.master(input_times[i,:],ψ,H_full,c_ops,rates=[κA,κB])
        temp_ar = Float64[]
        for k in 1:multiplexing
            rho_mes = states[k*(resolution÷multiplexing)]
            normalize!(rho_mes)
            # for i1 in 0:meas_max
            #     for i2 in 0:meas_max
            #         !!!
            #         temp = expect(rho_mes,fockstate(basis_a,i1)⊗fockstate(basis_b,i2))
            #         push!(temp_ar,abs(temp)^2)  #  !!! previous error because temp is imaginary
            #         push!(temp_ar,real(temp))
            #     end
            # end

            # add shot 
            temp_ar = Qmeasure_shot_mean(rho_mes, shot_nb)
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

