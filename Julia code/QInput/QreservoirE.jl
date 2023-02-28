using QuantumOptics
using Printf
using Plots


#useful quantum operators and kets

Ndim = 4
Ncd = 3

basis_a, basis_b = FockBasis(Ndim), FockBasis(Ndim)
basis_c, basis_d = FockBasis(Ncd), FockBasis(Ncd)


basis = basis_a ⊗ basis_b ⊗ basis_c ⊗ basis_d
a, b = embed(basis,1,destroy(basis_a)), embed(basis,2,destroy(basis_b))
a_d, b_d= embed(basis,1,create(basis_a)), embed(basis,2,create(basis_b))
Na, Nb = embed(basis,1,number(basis_a)), embed(basis,2,number(basis_b))
c, d = embed(basis,3,destroy(basis_c)), embed(basis,4,destroy(basis_d))
c_d, d_d= embed(basis,3,create(basis_c)), embed(basis,4,create(basis_d))
Nc, Nd = embed(basis,3,number(basis_c)), embed(basis,4,number(basis_d))

vacuum = fockstate(basis_a,0)⊗fockstate(basis_b,0)⊗fockstate(basis_c,0)⊗fockstate(basis_d,0)

meas_max = 3
wA, wB= 10e9, 9.5e9
g = 30e6
κA, κB = 17e6, 21e6
r, ϕp = 0.1,0.5*π

figure_title = "Reservoir simulation for: wA =  $(@sprintf("%.2e", wA)), wB =  $(@sprintf("%.2e", wB)), g =  $(@sprintf("%.2e", g))
κA =  $(@sprintf("%.2e", κA)), κB =  $(@sprintf("%.2e", κB)), r =  $(@sprintf("%.2e", r)), ϕp =  $(@sprintf("%.2e", ϕp))"


## Set Hamiltonian functions

# p = [wA,wB,g]
H_0(p) = p[1]*a_d*a + p[2]*b_d*b + p[3]*(a_d*b + a*b_d)

function H_input(p)
    # p = [κA,κB,r,ϕp]
    cout = cosh(p[3])*c + exp(im*p[4])*sinh(p[3])*d_d
    dout = cosh(p[3])*d + exp(im*p[4])*sinh(p[3])*c_d
    cout_d, dout_d = dagger(cout), dagger(dout)
    return 2*p[1]*(cout_d*a + a_d*cout) + 2*p[2]*(dout_d*b + b_d*dout)
end

param_0 = [wA,wB,g]
param_input = [κA,κB,r,ϕp]
H_full = H_0(param_0) + H_input(param_input)
c_ops = [a,b]


function constant_drive(show = 1,resolution = 3000, max_time = 50*(10^-9))
    println("simulation initialized")
    tspan = LinRange(0,max_time,resolution)
    ψ = fockstate(basis_a,0)⊗fockstate(basis_b,0)⊗coherentstate(basis_c,1)⊗coherentstate(basis_d,1)
    tout,states = timeevolution.master(tspan,ψ,H_full,c_ops,rates=[κA,κB])

    plot_means = hcat(real(expect(Na,states)),real(expect(Nb,states)))

    #plot means Na and Nb
    if show == 1
        plotly()
        p = plot(tspan,plot_means,title=figure_title,size=(1132,700),
        titlefontsize=10,label=["Na" "Nb"],xlabel = "time (seconds)", ylabel="populations")
        display(p)
    end
    println("end of simulation")
end

function Qmodel(r_list,ϕ_list,time_interval=100,resolution=300,sampling=10)
    L = length(r_list)

    ψ = fockstate(basis_a,0)⊗fockstate(basis_b,0)⊗coherentstate(basis_c,1)⊗coherentstate(basis_d,1)
    reservoir_output = zeros(L,sampling*(meas_max+1)^2 )
    ## training phase
    for i in 1:L
        ψ = fockstate(basis_a,0)⊗fockstate(basis_b,0)⊗coherentstate(basis_c,1)⊗coherentstate(basis_d,1)
        input_times = [j*(time_interval/resolution)*(10^-9) for j in 1:resolution]
        param_input = [κA,κB,r_list[i],ϕ_list[i]]
        H_full = H_0(param_0) + H_input(param_input)

        tout,states = timeevolution.master(input_times,ψ,H_full,c_ops,rates=[κA,κB])
        temp_ar = Float64[]
        for k in 1:sampling
            rho_mes = states[k*(resolution÷sampling)]
            normalize!(rho_mes)
            for i1 in 0:meas_max
                for i2 in 0:meas_max
                    # !!!
                    temp = expect(rho_mes,fockstate(basis_a,i1)⊗fockstate(basis_b,i2)⊗fockstate(basis_c,0)⊗fockstate(basis_d,0))
                    push!(temp_ar,abs(temp)^2)
                end
            end
        end
        reservoir_output[i,:] = temp_ar
        if i%10 == 0
            println("10 simulations passed, i = ", i)
        end
    end

    return reservoir_output
end
