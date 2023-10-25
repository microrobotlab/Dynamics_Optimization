using Optimization
using OptimizationOptimJL

include("objective_functions.jl")
include("utils.jl")


# Optimization parameters
# Nt = 100; packing_fraction = 0.1; L = 100.; R = 1.5; v = 10.
# wall_condition = "periodic"; nb_runs = 10
# N = 16; M = 16;

# Initial conditions
u0 = [0.2, 0.75, 5.]
# Boundaries
lb = [0.1, 0.55, 1.]
ub = [0.3, 1., 10.]

function F(u, p)
    packing_fraction, R, v = u
    @show packing_fraction
    @show packing_fraction_to_Np(packing_fraction, R, 100.)
    @show R
    @show v
    # wall_condition, nb_runs, N, M = p
    return mean_pf((Nt=100, Np=packing_fraction_to_Np(packing_fraction, R, 100.), L=100., R=R, v=v); wall_condition="periodic", nb_runs=3, N=16, M=16)
end

# Select parameters to be optimized
# variable_params = [:Np, :R, :v]

# Function to be optimized (OptimizationFunction instance)
optf = OptimizationFunction(F, Optimization.AutoForwardDiff())

# Optimization problem (OptimizationProblem instance)
prob = OptimizationProblem(optf, u0, lb=lb, ub=ub)
sol = solve(prob, NelderMead())