using QuantumOptics
using Plots

N=10
wA = 10e9
κA = 17e6
eA = 1e6
basis = FockBasis(N)
a,a_d,Na = destroy(basis),create(basis),number(basis)
ϕ0 = fockstate(basis,1)
H = wA*Na + sqrt(2*κA)*im*eA*(a - a_d)
c_ops = [sqrt(κA)*a]
tspan = LinRange(0,50e-9,3000)
L = -im*spre(H) + im*spost(H) + κA *(spre(a)*spost(a_d) - 0.5*spre(Na) - 0.5*spost(Na))
@time tout,states1 = timeevolution.master(tspan,ϕ0,L)
@time tout,states2 = timeevolution.master(tspan,ϕ0,H,c_ops)
plot_means = hcat(real(expect(Na,states1)),real(expect(Na,states2)))

plotly()
p = plot(tout,plot_means,title="test",size=(500,350),
titlefontsize=12,label=["liouvillian" "linblad"],legend=:topright)
display(p)
