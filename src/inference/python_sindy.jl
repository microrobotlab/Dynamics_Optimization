# PRELIMINARY TESTS ON AUTOMATIC MODEL INFERENCE USING SINDy ALGORITHM 
# (see "Discovering governing equations from data: Sparse identification of nonlinear 
# dynamical systems", Steven L. Brunton, Joshua L. Proctor, J. Nathan Kutz)
# PYTHON IMPLEMENTATION "pysindy" PACKAGE IS USED ON DATA COMING FROM THE SIMULATOR 
# IN WHICH WE CAN CHANGE PARTICLE DYNAMICS USING DIFFERENT MODELS IN src/step.jl.
# STRONGLY INSPIRED BY THE TUTORIAL IN THE pysindy CORRESPONDING DOCUMENTATION
# (see https://pysindy.readthedocs.io/en/latest/examples/15_pysindy_lectures/example.html)

# NOTE THAT PYTHON AND PYTHON PACKAGES "pysindy" and "sklearn.metrics"
# MUST BE INSTALLED TO USE THE FOLLOWING


"""
MIT License

Copyright (c) for portions of project PySINDy are held by Markus Quade, 2019 as
part of project sparsereg. All other copyright for project PySINDy are held by
Brian de Silva and Kathleen Champion 2019.

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


using PyCall 
using DrWatson
using ProgressBars

# `srcdir` from DrWatson provides path to src/. We add to the string folder/file names provided in arguments.
include(srcdir("ABP output.jl"))
include("utils_inference.jl")
pysindy = pyimport("pysindy")
metrics = pyimport("sklearn.metrics")


# Timestep
δt = 1e-3
# Number of timesteps
Nt = 10000
# Time
t = collect(range(start=0., length=Nt+1, step=δt))
t_span = (t[1], t[end])
# Train set -> multiple trajectories with different initial conditions (but same dynamics)
n_trajectories = 3
# Generate data from simulator
sim_output_train = [run((Nt=Nt, Np=1, L=100., R=1.5, v=10.); wall_condition="open", collision_correction=false, N=16, M=16) for _=1:n_trajectories]
# (format output to fit pysindy)
X_train = [tup2matrix_states(s)' for s in sim_output_train]
# Test set 
sim_output_test = run((Nt=Nt, Np=1, L=100., R=1.5, v=10.); wall_condition="open", collision_correction=false, N=16, M=16)
X_test = tup2matrix_states(sim_output_test)' 


# ----- Training loop with various thresholds for Sequentially Thresholded Least Squares algorithm (STLSQ) 
# used for sparsity (described in original SINDy paper cited at the beginning of the file)

# Thresholds
threshold_scan = collect(range(0.0,1.0,100))
# Will store mean squared error on predicted trajectories 
mse_traj_list = []
# Will store mean squared error on predicted derivatives 
mse_deriv_list = []
models_list = []

for threshold in ProgressBar(threshold_scan)
    # Create STLSQ optimizer with given threshold 
    stlsq_opt = pysindy.STLSQ(threshold=threshold)
    # Functions library for inference : functions that will be combined
    # to find dynamics corresponding to the provided data
    feature_library = pysindy.PolynomialLibrary(1)
    # Model definition and fitting on train data
    model = pysindy.SINDy(feature_library=feature_library, optimizer=stlsq_opt, feature_names=["x", "y", "θ"], t_default=δt)
    model.fit(X_train, t=δt, quiet=true, multiple_trajectories=true)

    # Error computation one test data
    # 1. On trajectories using `simulate()`
    X_test_sim = model.simulate(X_test[1,:], t, integrator="odeint")
    # if(any(X_test_sim .> 1e4)) fill!(X_test_sim, 1e4) end
    mse_traj = metrics.mean_squared_error(X_test, X_test_sim)
    # 2. on derivatives using `score()`
    mse_deriv = model.score(X_test, t=δt, metric=metrics.mean_squared_error)
    push!(mse_traj_list, mse_traj)
    push!(mse_deriv_list, mse_deriv)
    push!(models_list, model)
end

# Complexity for the different infered models are also stored
# The complexity of a model corresponds to the number of elements used in the features library
# that is the number of non-zero parameters in the model
complexity_list = [model.complexity for model in models_list]

# Plotting
# MSE on trajectory prediction w.r.t STLSQ threshold λ
p_acc_λ = plot(threshold_scan, mse_traj_list, label="MSE (trajectory) vs. λ", xlabel="Threshold (λ)", ylabel="Mean squared error", yscale=:log10, ylims=:auto, color=:red, legend=:topright)
# MSE on derivatives prediction w.r.t. STLSQ threshold λ
p_score_λ = plot(threshold_scan, mse_deriv_list, label="MSE (derivatives) vs. λ", xlabel="Threshold (λ)", ylabel="Mean squared error", yscale=:log10, ylims=:auto, color=:red, legend=:topright)
# Model complexity w.r.t. STLSQ threshold λ
p_complexity_λ = plot(threshold_scan, complexity_list, label="Model complexity vs. λ", xlabel="Threshold (λ)", ylabel="Model complexity", color=:orange)
# Model complexity w.r.t. MSE on derivatives prediction
p_complexity_acc = scatter(complexity_list, mse_deriv_list, label="Model complexity \nvs. MSE (derivatives)",minorgrid=true, xlabel="Model complexity", ylabel="Mean squared error", yscale=:log10, ylims=:auto, color=:green, legend=:bottomright)

# Everything together
plot(p_acc_λ, p_score_λ, p_complexity_λ, p_complexity_acc, layout=layout = grid(4, 2, heights=[0.5, 0.5, 0.5, 0.5]), xguidefonthalign = :right, yguidefontvalign = :top, guidefontsize=8, tickfontsize=6, legendfontsize=5, size=(600, 400), minorgrid=false)


# To plot best models predictions
#=
# index corresponding to best model w.r.t. trajectory error 
best_traj_idx = argmin(mse_traj_list)
# corresponding complexity
best_traj_complexity = complexity_list[best_traj_idx]
# corresponding model
best_model_traj = models_list[best_traj_idx]
# same w.r.t. derivatives
best_deriv_idx = argmin(mse_deriv_list)
best_deriv_complexity = complexity_list[best_deriv_idx]
best_model_deriv = models_list[best_deriv_idx]
# simulate using the two best models above with intitial conditions of the test set for reconstruction evalutation
X_test_sim_traj = best_model_traj.simulate(X_test[1,:], t, integrator="odeint")
X_test_sim_deriv = best_model_deriv.simulate(X_test[1,:], t, integrator="odeint")
plot(X_test, label=["x" "y" "θ"], color=[:red :orange :blue])
plot!(X_test_sim_traj, labels=["x_pred_traj" "y_pred_traj" "θ_pred_traj"], color=[:red :orange :blue], linestyle=:dash)
plot!(X_test_sim_deriv, labels=["x_pred_deriv" "y_pred_deriv" "θ_pred_deriv"], color=[:red :orange :blue], linestyle=:dashdot)
=#