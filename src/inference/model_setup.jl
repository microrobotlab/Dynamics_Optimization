# SciML Tools
using DifferentialEquations, OrdinaryDiffEq, ModelingToolkit, DataDrivenDiffEq, SciMLSensitivity, DataDrivenSparse
using Optimization, OptimizationOptimisers, OptimizationOptimJL

# Standard Libraries
using LinearAlgebra, Statistics

# External Libraries
using ComponentArrays, Lux, Zygote, Plots, StableRNGs
gr()


# Set a random seed for reproducible behaviour
rng = StableRNG(5)

# ----- MODELS
"""
# Random differential equation (ABP model)
function ABP!(du, u, p, t, W)
    v0, DT, DR = p
    x, y, θ = u
    du[1] = v0 * cos(θ) + sqrt(2*DT) * W[1]
    du[2] = v0 * sin(θ) + sqrt(2*DT) * W[2]
    du[3] = sqrt(2*DR) * W[3]
end

# Define the experimental parameter
tspan = (0.0, 5.0)
u0 = [0., 0., 0.]
p_ = [10., 0.5, 0.5]

# problem setting and solving
prob = RODEProblem(ABP!, u0, tspan, p_)
sol = solve(prob, RandomEM(), dt=1/10)
"""

# Deterministic model similar to ABP for inference
function deterministic_ABP!(du, u, p, t)
    # translational speed, translational / rotational agitation
    v0, DT, DR, f = p
    x, y, θ = u
    # we use sine functions to imitate diffusion
    du[1] = v0 * cos(θ) + sqrt(2*DT) * sin((1/f)*t)
    du[2] = v0 * sin(θ) + sqrt(2*DT) * sin((1/f)*t)
    du[3] = sqrt(2*DR)  * sin((1/f)*t)
end

# Define the experimental parameter
exp_duration = 10. # duration of the caption (s)
img_per_sec = 10 # number of images per second in video (gives sampling frequency)
exp_noise_magnitude = 0.5 # estimated measurement noise magnitude (μm)
v0 = 5. # speed of the particles (μm/s)
noise_trans = 10. # noise for translational part of the movement (`deterministic_ABP!` model)
noise_rot = 10. # noise for rotational part of the movement (`deterministic_ABP!` model)
f = 500. # oscillation frequency for our brownian motion imitation (`deterministic_ABP!` model) (
tspan = (0.0, exp_duration)
u0 = [0., 0., 0.]
p_ = [v0, noise_trans, noise_rot, f]

# ----- ORIGINAL PROBLEM
prob = ODEProblem(deterministic_ABP!, u0, tspan, p_)
sol = solve(prob, Vern7(), abstol = 1e-12, reltol = 1e-12, saveat = 1/img_per_sec)

X = Array(sol)
# println(X)
t = sol.t
plot(sol[1:100], label = ["x" "y" "θ"])

# Adding noise
x̄ = mean(X, dims = 2)
noise_magnitude = 5e-3
Xₙ = X .+ (noise_magnitude * x̄) .* randn(rng, eltype(X), size(X))
scatter!(t, transpose(Xₙ), label = ["x noisy" "y noisy" "θ noisy"], markersize=3.)