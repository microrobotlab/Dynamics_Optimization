using DrWatson
@quickactivate "active-brownian-particles"
using Optim
using Plots


# test function
f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

X = LinRange(-5, 5, 1000)
Y = LinRange(-5 ,5, 1000)
Z = f.(Iterators.product(X,Y))
heatmap(X,Y,Z, c=:thermal, title="Nelder-Mead")

x0 = [0., 0.]
opt = Optim.Options(store_trace = true, show_trace=true, extended_trace=true)
res = optimize(f, x0, LBFGS(), opt)
x_trace = Optim.x_trace(res)

# animated plot
anim = @animate for i in eachindex(x_trace)
    contourf(X,Y,Z, c=:thermal, title="LBFGS, i=$i")
    plot!([x_trace[j][1] for j=1:i], [x_trace[j][2] for j=1:i], label=nothing, c=:red, marker=:star)
end

gif(anim, fps=2)