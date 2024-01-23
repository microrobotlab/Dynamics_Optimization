# TESTS ON DYNAMICS MODEL INFERENCE USING NEURAL NETWORK / DIFFERENTIAL EQUATION MIXED MODELS 
# SEE "Universal Differential Equations" 
# THE DATA COMES FROM synthetic_data.jl
# STRONGLY INSPIRED BY DOCUMENTATION TUTORIAL https://docs.sciml.ai/Overview/dev/showcase/missing_physics/


"""MIT License

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


using Optimization, OptimizationOptimisers, OptimizationOptimJL
using Lux, Zygote, SciMLSensitivity
using ComponentArrays

include("synthetic_data.jl")


# Activation function
rbf(x) = exp.(-(x .^ 2))
# The neural network
const U = Lux.Chain(Lux.Dense(3, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 3))
p, st = Lux.setup(rng, U) # rng defined in synthetic_data.jl
const _st = st

# Hybrid Differential Equation / Neural Network model : ITS WEIGHTS WILL BE THE PARAMETERS OF THE DIFFERENTIAL EQUATION
function hybrid_deterministic_ABP!(du, u, p, t, p_true)
    # Initial parameters for Active Brownian Particles (ABP) model
    # -> speed and translational / rotational diffusion coefficients
    v0, DT, DR = p_true
    # Network prediction, takes u and parameters p as input
    û = U(u, p, _st)[1] 
    x, y, θ = u
    # The NN will try to model the random (diffusion) part of the ABP model 
    # see https://en.wikipedia.org/wiki/Active_Brownian_particle#Equations_of_motion
    du[1] = v0 * cos(θ) + û[1]
    du[2] = v0 * sin(θ) + û[2]
    du[3] = û[3]
end

# Closure with known parameters
nn_deterministic_ABP!(du, u, p, t) = hybrid_deterministic_ABP!(du, u, p, t, p_)

# Define the problem
prob_nn = ODEProblem(nn_deterministic_ABP!, Xₙ[:, 1], tspan, p) # Xₙ defined in synthetic_data.jl

# Prediction
function predict(θ, X = Xₙ[:, 1], T = t)
    # Create differential equation problem with new parameters 
    # (they are updated while training of the neural network)
    _prob = remake(prob_nn, u0=X, tspan=(T[1], T[end]), p=θ)
    # solve new problem
    Array(solve(_prob, Vern7(), saveat=T, 
                abstol=1e-6, reltol=1e-6,
                verbose=false)) 
end

# Loss compared to original data 
function loss(θ)
    X̂ = predict(θ)
    mean(abs2, Xₙ .- X̂)
end

# Callback is the function called after each loop iteration
losses = Float64[]
callback = function (p, l)
    push!(losses, l)
    if length(losses) % 50 == 0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end

# Build optimization problem based on loss
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentVector{Float64}(p))

# Training
# FIRST WITH ADAM
res1 = Optimization.solve(optprob, ADAM(), callback=callback, maxiters=6000)
println("Training loss after $(length(losses)) iterations: $(losses[end])")

# FINISH WITH LBFGS
optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(optprob2, LBFGS(), callback=callback, maxiters=1000)
println("Final training loss after $(length(losses)) iterations: $(losses[end])")
# Best parameters
p_trained = res2.u

# Loss plot
# Plot loss evolution with ADAM optimizer 
pl_losses = plot(1:6000, losses[1:6000], yaxis=:log10, xaxis=:log10, xlabel="Iterations", ylabel="Loss", label="ADAM", color=:blue)
# Plot loss evolution with LBFGS optimizer 
plot!(6001:length(losses), losses[6001:end], yaxis=:log10, xaxis=:log10, xlabel="Iterations", ylabel="Loss", label="LBFGS", color=:red)

# Reconstruction analysis
ts = first(sol.t):(mean(diff(sol.t)) / 2):last(sol.t)
X̂ = predict(p_trained, Xₙ[:, 1], ts)
pl_trajectory = plot(ts, transpose(X̂), xlabel="t", ylabel="x(t), y(t)", color=:red, label=["UDE Approximation" nothing])
scatter!(sol.t, transpose(Xₙ), color=:black, label=["Measurements" nothing])
