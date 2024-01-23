# OPTIMIZATION TESTS USING UNIFIED 'Optimization' PACKAGE FRAMEWORK 
# (see https://docs.sciml.ai/Optimization/stable/)

using Optimization
using OptimizationOptimJL

# Objective function 
# must have variable x and parameters p as parameters
f(u, p) =  (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
optf = OptimizationFunction(f, Optimization.AutoForwardDiff())
# Initial value
u0 = zeros(2)
# Parameters
p = [1., 100.]

# Definition of the problem
prob = Optimization.OptimizationProblem(optf, u0, p)

# Solve using different optimization methods
sol = solve(prob, NelderMead())
sol = solve(prob, BFGS())