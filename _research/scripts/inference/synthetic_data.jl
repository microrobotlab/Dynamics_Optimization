# TEST TO GENERATE SYNTHETIC POSITIONAL DATA BASED ON DIFFERENTIAL EQUATION RESOLUTION FOR DYNAMIC PARTICLES. THE GOAL 
# WAS TO IMITATE DATA FROM REAL MICROSCOPE OBSERVATIONS, USING PARAMETERS SUCH AS THE NUMBER OF FRAMES PER SECOND.
# GENERATED DATA IS GIVEN AS INPUT FOR MODEL IN NN_inference.jl 
# USES CODE FROM DOCUMENTATION TUTORIAL https://docs.sciml.ai/Overview/dev/showcase/missing_physics/


"""
MIT License

Copyright (c) 2020 Christopher Rackauckas

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


# SciML Tools
using DifferentialEquations

# External Libraries
using Plots, StableRNGs, Statistics
gr()


# Set random seed for reproducibility 
rng = StableRNG(5)

# ----- MODEL WITH RANDOMNESS
# Random differential equation using Active Brownian Particle (ABP) model
function ABP!(du, u, p, t, W)
    v0, DT, DR = p
    x, y, θ = u
    du[1] = v0 * cos(θ) + sqrt(2*DT) * W[1]
    du[2] = v0 * sin(θ) + sqrt(2*DT) * W[2]
    du[3] = sqrt(2*DR) * W[3]
end

# Define experimental parameters
# Time span
tspan = (0.0, 5.0)
# Initial state (x₀, y₀, θ₀)
u0 = [0., 0., 0.]
# Model parameters: v0 (particle initial speed), 
# DT (translational diffusion coefficient), DR (rotational diffusion coefficient)
p_ = [10., 0.5, 0.5]

# Problem setting and generation of positional data using random differential equation solver
prob = RODEProblem(ABP!, u0, tspan, p_)
sol = solve(prob, RandomEM(), dt=1/10)



# ----- DETERMINISTIC MODEL (similar to ABP random model)
function deterministic_ABP!(du, u, p, t)
    # f 
    v0, DT, DR, f = p
    x, y, θ = u
    # we use sine functions to imitate random diffusion in ABP model
    du[1] = v0 * cos(θ) + sqrt(2*DT) * sin((1/f)*t)
    du[2] = v0 * sin(θ) + sqrt(2*DT) * sin((1/f)*t)
    du[3] = sqrt(2*DR)  * sin((1/f)*t)
end

# Define experimental parameters that would fit real experiment parameters
exp_duration = 10.          # duration of the caption (s)
img_per_sec = 10            # number of images per second in video (gives sampling frequency)
exp_noise_magnitude = 0.5   # estimated measurement noise magnitude (μm)
v0 = 5.                     # speed of the particles (μm/s)
noise_trans = 10.           # noise for translational part of the movement (`deterministic_ABP!` model)
noise_rot = 10.             # noise for rotational part of the movement (`deterministic_ABP!` model)
f = 500.                    # oscillation frequency for our brownian motion imitation (`deterministic_ABP!` model)
tspan = (0.0, exp_duration)
u0 = [0., 0., 0.]
p_ = [v0, noise_trans, noise_rot, f]    # grouped parameters similar to previous random model

# problem setting and generation of positional data using differential equation solver
prob = ODEProblem(deterministic_ABP!, u0, tspan, p_)
sol = solve(prob, Vern7(), abstol = 1e-12, reltol = 1e-12, saveat = 1/img_per_sec)

X = Array(sol)
t = sol.t

# add noise 
x̄ = mean(X, dims = 2)
noise_magnitude = 5e-3
Xₙ = X .+ (noise_magnitude * x̄) .* randn(rng, eltype(X), size(X))

plot(sol, alpha = 0.75, color = :black, label = ["True Data" nothing])
scatter!(t, transpose(Xₙ), color = :red, label = ["Noisy Data" nothing])