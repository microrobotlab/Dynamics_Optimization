# PRELIMINARY TESTS ON AUTOMATIC MODEL INFERENCE USING SINDy ALGORITHM 
# (see "Discovering governing equations from data: Sparse identification of nonlinear 
# dynamical systems", Steven L. Brunton, Joshua L. Proctor, J. Nathan Kutz)
# SIMILAR TO python_sindy.jl FILE BUT USES JULIA PACKAGES DIRECTLY


using DrWatson
using DataDrivenDiffEq
using ModelingToolkit
using DataDrivenSparse
using Plots
using ProgressBars

# `srcdir` from DrWatson provides path to src/ to which we add elements provided in arguments
include(srcdir("inference", "utils_inference.jl"))
include(srcdir("ABP output.jl"))


# Timestep
δt = 1e-3
# Number of timesteps
Nt = 1000
# Time
times = collect(range(start=0, length=Nt+1, step=δt))


# ----- TESTS ON VERY SIMPLE DATA FROM DYNAMICS WITH LINEAR TIME DEPENDENCE 
# Trajectory
X_simple = hcat(-2.0times, 3.0times)'
# Automatically extract derivatives from positions (derivatives are constants here)
# (seems to smooth it)  
DX_simple_interpolated, X_simple_interpolated, times_interpolated = collocate_data(X_simple ,times)
# Problem definition
prob_simple = ContinuousDataDrivenProblem(X_simple, times, DX_simple_interpolated)

@parameters t
@variables u(t), v(t)
# features library that will be used to find dynamics
basis_simple = Basis(monomial_basis([u,v], 1), [u,v])
println(basis_simple)

# Solve
res_simple = solve(prob_simple, basis_simple, STLSQ())
println(res_simple.basis)
# Print found parameters
println("Found parameters : $(res_simple.prob.p)")



# ----- TESTS ON DATA FROM THE SIMULATOR IN WHICH WE CAN CHANGE PARTICLE DYNAMICS USING DIFFERENT MODELS IN src/step.jl.
# (note that in this example we use one particle only with `Np`=1, no boundaries with  `wall_condition`="open" 
# and no collision between particles with `collision_correction`=false)
# Generate data from simulator
sim_output = run((Nt=Nt, Np=1, L=100., R=1.5, v=5.); wall_condition="open", collision_correction=false, N=16, M=16)
# Output formatting
X = tup2matrix_states(sim_output) 
# Same as before
DX_interpolated, X_interpolated, times_interpolated = collocate_data(X ,times)
prob = ContinuousDataDrivenProblem(X, times, DX_interpolated)

@variables (u[1:3])(t) 
basis = Basis([1.; u[1]; u[2]; u[3]], u, iv=t)
# basis = Basis(Num[1.], [x, y, θ])
println(basis)

# Same as before
res = solve(prob, basis, STLSQ())
println(res.basis)
println("Found parameters : $(res.prob.p)")

# Just to observe data after applying `collocate_data()`
plot(X[3,:], label="X")
plot!(X_interpolated[3,:], label="X_interpolated")
plot!(DX_interpolated[3,:], label="DX_interpolated")


# ----- Training loop with various thresholds for Sequentially Thresholded Least Squares algorithm (STLSQ) 
# used for sparsity (described in original SINDy paper cited at the beginning of the file)
# random data (must be changed)
X = rand(2, 10)
N = 10
# Time
t = Float64.(collect(1:N))
X = Matrix(2.0t') + 0.5 * rand(N)'
# Problem definition
prob = ContinuousDataDrivenProblem(X, t)

sampler = DataProcessing(split = 0.8, shuffle = true, batchsize = 30)
# Thresholds that will be tested
λ_list = exp10.(-10:0.1:0)
# Thresholds that work for the algorithm
λ_list_ok = []

# Accuracy of the infered models (Residual Sum of Squares metric)
acc_list = []
# Complexity for the different infered models are also stored
# It corresponds to the number of elements used in the features library,
# that is the number of non-zero parameters in the model
complexity_list = []
for λ in ProgressBar(λ_list)
    try
        res = solve(prob, basis, STLSQ(λ), options = DataDrivenCommonOptions(data_processing = sampler, digits = 1))
        push!(acc_list, res.residuals)
        system = get_basis(res)
        # length to have number of non-zero parameters in the model 
        push!(complexity_list, length(get_parameter_values(system)))
        push!(λ_list_ok, λ)
    catch
        push!(acc_list, -1.)
        println("error : λ=$λ")
    end
end

plot(complexity_list, acc_list, yaxis=:log10)


# Data X corresponds to positions x, y and angle θ (see ABP main_parallel.jl)
# plot(X[1,:], label="x")
# plot!(X[2,:], label="y")
# plot!(X[3,:], label="θ")
# plot!(DX_interpolated[1,:], label="∂x")
# plot!(DX_interpolated[2,:], label="∂y")
# plot!(DX_interpolated[3,:], label="∂θ")