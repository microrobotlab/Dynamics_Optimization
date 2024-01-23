# MAIN FILE FOR THE OPTIMIZATION OF PARTICLES BEHAVIOR (w.r.t. objective functions defined in `objective_functions.jl`)
# ON THE SET OF PHYSICAL PARAMETERS OF THE SIMULATOR (see src/main_parallel.jl).


using DrWatson
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using Optim
using Random

include("objective_functions.jl")
include("utils_optimization.jl")


# ----- SIMULATOR HYPERPARAMETERS AND PARAMETERS THAT WILL REMAIN FIXED DURING OPTIMIZATION; 
# N, M are the number of vertical, horizontal divisions of the space for parallel computing (see `ABP main_parallel.jl`) 
# TODO: maximum possible values for N and M depend on L and R. It is thus useful to have an auto mode to automatically determine 
# those maximum values since R will vary during optimization (see `max_space_division()` function from `utils_simulation.jl` for an 
# implementation that could be used). Fixed values for N and M can lead to overlapping as R varies (see `hardsphere!` function in 
# ABP main_parallel.jl for more details) and make the simulation fail.

# Simulator 'hyperparameters'
p =(
    # /!\ Do not use "elliptical" for the moment because `packing_fraction_to_Np` cannot handle this shape 
    # and the number of particles might be wrong. See comment below in OPTIMIZED PARAMETERS.
    wall_condition = "periodic", # / "open" / "squared" 
    collision_correction = true,
    nb_runs = 10,
    N = 4, 
    M = 4    
)

# Physical parameters that remain fixed during optimization 
Nt = 1000
L = 100.

# ----- OPTIMIZED PARAMETERS 
# `packing_fraction`, represents a density of particles in space. It is used instead of the number of 
# particles directly. It transforms a discrete parameter into a continuous one for our code.
# See `packing_fraction_to_Np` for a precise definition of the packing fraction; this function is 
# used to convert back to the number of particles.  

# Initial conditions, determines parameters to be optimized
u0 = [
    0.1, # packing_fraction
    1.5, # particle radius R
    5.   # particle speed v
]

# # Boundaries
# lb = [ 
#     0.2, # packing_fraction
#     1.5, # particle radius R
#     5. # particle speed v
# ] # low
# ub = [ 
#     0.2, # packing_fraction
#     1.5, # particle radius R
#     5. # particle speed v
# ] # up


# ATTEMPT TO DEFINE BOX CONSTRAINTS ON PHYSICAL PARAMETERS
# SHOULD USE SCIML PACKAGE INSTEAD
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


# ----- OBJECTIVE FUNCTION (see `objective_functions.jl`) 
objective_function = mean_pf


# Optimization problem (OptimizationProblem instance)
# First modify the objective function, because it must fit SciML objective functions format
# F(u, p) where u are the state variables and p other parameters (here simulator hyperparameters)
# (see https://docs.sciml.ai/Optimization/stable/API/optimization_function/)
function optf(u,p)
    packing_fraction, R, v = u
    wall_condition, collision_correction, nb_runs, N, M = p
    return objective_function((Nt=Nt, packing_fraction=packing_fraction, L=L, R=R, v=v); wall_condition=wall_condition, collision_correction=collision_correction, nb_runs=nb_runs, N=N, M=M)
end

# Problem instance
prob = OptimizationProblem(optf, u0, p)
# Specific parameters for the optimization methods of the different packages used by SciML are passed in `solve()` as keywords arguments
# Other optimization functions can be tried
sol = solve(prob, NelderMead(), store_trace=true, extended_trace=true, show_trace=true, allow_f_increases=true, iterations=20)



# To vizualise optimization output with NelderMead method
# centroids of the successive simplexes used during optimization with NelderMead 
# (see https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
centroids_explored = Optim.centroid_trace(sol.original)
println(centroids_explored)

# f_explored = Optim.f_trace(sol.original)