using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using Optim
using Random

include("objective_functions.jl")
include("utils_optimization.jl")


# ----- SIMULATOR HYPERPARAMETERS AND FIXED PARAMETERS
# N, M are the number of vertical, horizontal divisions of the space for parallel computing; 
# :auto mode infer maximum values for it. Those maximum values depend on L and R,
#  so it is useful to have an auto mode since R will vary through optimization process.
p =(
    wall_condition = "periodic",
    collision_correction = true,
    nb_runs = 10,
    N = 4, 
    M = 4    
)

Nt = 1000
L = 100.

# ----- OPTIMIZATION PARAMETERS
objective_function = mean_pf

# Initial conditions
u0 = [
    0.1, # packing_fraction
    1.5, # R
    5. # v
]

# Boundaries
lb = [ 
    0.2, # packing_fraction
    1.5, # R
    5. # v
] # low
ub = [ 
    0.2, # packing_fraction
    1.5, # R
    5. # v
] # up


# function F(u, p)
#     packing_fraction, R, v = u
#     # wall_condition, nb_runs, N, M = p
#     lb = [0.1, 0.55, 1.]
#     ub = [0.3, 1., 10.]
#     for i in eachindex(u)
#         if u[i] < lb[i] || u[i] > ub[i]
#             return Inf
#         end
#     end
#     return mean_pf([Nt, packing_fraction, L, R, v]; p)
#     # return packing_fraction*R - v
# end

# function F_small(u)
#     return u.x^2 + u.y
# end


# Optimization problem (OptimizationProblem instance)
# First modify the objective function, because it must be of form F(u, p)
# where u are the state variables and p other parameters (here simulator hyperparameters)
# (see https://docs.sciml.ai/Optimization/stable/API/optimization_function/)
function optf(u,p)
    packing_fraction, R, v = u
    wall_condition, collision_correction, nb_runs, N, M = p
    return objective_function((Nt=Nt, packing_fraction=packing_fraction, L=L, R=R, v=v); wall_condition=wall_condition, collision_correction=collision_correction, nb_runs=nb_runs, N=N, M=M)
end
prob = OptimizationProblem(optf, u0, p) #lb=lb, ub=ub)
# Specific parameters for the optimization methods in the different packages are passed in `solve()` as keywords arguments
sol = solve(prob, NelderMead(), store_trace=true, extended_trace=true, show_trace=true, allow_f_increases=true, iterations=20)

centroids_explored = Optim.centroid_trace(sol.original)
f_explored = Optim.f_trace(sol.original)

println(centroids_explored)

# scatter(xyz[:,1], xyz[:,2], xyz[:,3], c=colormap("Reds",length(t)), markersize=2.5)
# scatter!(first.(Optim.centroid_trace(sol)), Optim.f_trace(sol))