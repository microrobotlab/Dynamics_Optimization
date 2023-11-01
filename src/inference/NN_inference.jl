include("model_setup.jl")


# ----- INFERENCE WITH HYBRID DEFINED / NEURAL NETWORK MODEL (`Universal Differential Equation`)
# Neural Network model : THE WEIGHTS WILL BE THE PARAMETERS OF THE ODE 
rbf(x) = exp.(-(x .^ 2))
const U = Lux.Chain(Lux.Dense(3, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 5, rbf), Lux.Dense(5, 3))
p, st = Lux.setup(rng, U)
const _st = st

# Hybrid ODE 
function hybrid_deterministic_ABP!(du, u, p, t, p_true)
    v0, DT, DR = p_true
    # Network prediction, takes u and also parameters p
    û = U(u, p, _st)[1] 
    x, y, θ = u
    # The NN will try to model the "noise part"
    du[1] = v0 * cos(θ) + û[1]
    du[2] = v0 * sin(θ) + û[2]
    du[3] = û[3]
end

# Closure with the known parameters
nn_deterministic_ABP!(du, u, p, t) = hybrid_deterministic_ABP!(du, u, p, t, p_)

# Define the problem
prob_nn = ODEProblem(nn_deterministic_ABP!, Xₙ[:, 1], tspan, p)

# Prediction
function predict(θ, X = Xₙ[:, 1], T = t)
    _prob = remake(prob_nn, u0=X, tspan=(T[1], T[end]), p=θ)
    Array(solve(_prob, Vern7(), saveat=T, 
                abstol=1e-6, reltol=1e-6,
                verbose=false)) 
end

# Loss
function loss(θ)
    X̂ = predict(θ)
    mean(abs2, Xₙ .- X̂)
end

# The callback will be called after each loop iteration
losses = Float64[]
callback = function (p, l)
    push!(losses, l)
    if length(losses) % 50 == 0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end

# Build optimization problem
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentVector{Float64}(p))

# Training
# FIRST WITH ADAM
res1 = Optimization.solve(optprob, ADAM(), callback=callback, maxiters=6000)
println("Training loss after $(length(losses)) iterations: $(losses[end])")

# FINISHING WITH LBFGS
optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(optprob2, LBFGS(), callback=callback, maxiters=1000)
println("Final training loss after $(length(losses)) iterations: $(losses[end])")
# Best parameters
p_trained = res2.u

# Loss plot
pl_losses = plot(1:6000, losses[1:6000], yaxis=:log10, xaxis=:log10, xlabel="Iterations", ylabel="Loss", label="ADAM", color=:blue)
plot!(6001:length(losses), losses[6001:end], yaxis=:log10, xaxis=:log10, xlabel="Iterations", ylabel="Loss", label="LBFGS", color=:red)

# Reconstruction analysis
ts = first(sol.t):(mean(diff(sol.t)) / 2):last(sol.t)
X̂ = predict(p_trained, Xₙ[:, 1], ts)
pl_trajectory = plot(ts, transpose(X̂), xlabel="t", ylabel="x(t), y(t)", color=:red, label=["UDE Approximation" nothing])
scatter!(sol.t, transpose(Xₙ), color=:black, label=["Measurements" nothing])
