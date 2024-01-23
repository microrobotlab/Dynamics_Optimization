# OPTIMIZATION TESTS WITH VIZUALIZATION OF THE ITERATIONS

using DrWatson
using Optim
using Plots


# test function
f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2

X = LinRange(-5, 5, 1000)
Y = LinRange(-5 ,5, 1000)
Z = f.(Iterators.product(X,Y))

# initial condition for optimization
x0 = [0., 0.]
# optimization options
opt = Optim.Options(store_trace = true, show_trace=true, extended_trace=true)
res = optimize(f, x0, LBFGS(), opt) # use LBFGS() algorithm here
x_trace = Optim.x_trace(res)        # variables values explored

# animated plot
anim = @animate for i in eachindex(x_trace)
    # plot function as heatmap with contours
    contourf(X,Y,Z, c=:thermal, title="LBFGS, i=$i")
    # to see optimization path (red lines between points explored)
    plot!([x_trace[j][1] for j=1:i], [x_trace[j][2] for j=1:i], label=nothing, c=:red, marker=:star)
end

gif(anim, fps=2)