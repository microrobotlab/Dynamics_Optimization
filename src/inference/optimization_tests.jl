using Optimization
using OptimizationOptimJL

# Objective function, initial value, parameters
f(u, p) =  (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
optf = OptimizationFunction(f, Optimization.AutoForwardDiff())
u0 = zeros(2)
p = [1., 100.]

# Definition of the problem
prob = Optimization.OptimizationProblem(optf, u0, p)

# Solve 
sol = solve(prob, NelderMead())
sol = solve(prob, BFGS())