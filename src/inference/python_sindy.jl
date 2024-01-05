using PyCall 
using DrWatson
using ProgressBars
include(srcdir("ABP output.jl"))
include("utils_inference.jl")
pysindy = pyimport("pysindy")
metrics = pyimport("sklearn.metrics")
integrate = pyimport("scipy.integrate")

# py"""
# import numpy as np
# from scipy.integrate import solve_ivp
# from pysindy.utils import lorenz, lorenz_control, enzyme
# np.random.seed(100)

# # Initialize integrator keywords for solve_ivp to replicate the odeint defaults
# integrator_keywords = {}
# integrator_keywords['rtol'] = 1e-12
# integrator_keywords['method'] = 'LSODA'
# integrator_keywords['atol'] = 1e-12
# def get_train_test(δt):
#     t_train = np.arange(0, 10, δt)
#     x0_train = [-8, 8, 27]
#     t_train_span = (t_train[0], t_train[-1])
#     x_train = solve_ivp(
#         lorenz, t_train_span, x0_train, t_eval=t_train, **integrator_keywords
#     ).y.T

#     t_test = np.arange(0, 15, δt)
#     t_test_span = (t_test[0], t_test[-1])
#     x0_test = np.array([8, 7, 15])
#     x_test = solve_ivp(
#         lorenz, t_test_span, x0_test, t_eval=t_test, **integrator_keywords
#     ).y.T
#     return t_train, x_train, t_test, x_test
# """
# t_train, X_train, t_test, X_test = py"get_train_test"(δt)

δt = 1e-3
Nt = 10000
# Time
t = collect(range(start=0., length=Nt+1, step=δt))
t_span = (t[1], t[end])
# Train set -> multiple trajectories with di
n_trajectories = 3
sim_output_train = [run((Nt=Nt, Np=1, L=100., R=1.5, v=10.); wall_condition="open", collision_correction=false, N=16, M=16) for _=1:n_trajectories]
X_train = [tup2matrix_states(s)' for s in sim_output_train]
sim_output_test = run((Nt=Nt, Np=1, L=100., R=1.5, v=10.); wall_condition="open", collision_correction=false, N=16, M=16)
X_test = tup2matrix_states(sim_output_test)' 

# Training loop with various thresholds for STLSQ method
threshold_scan = collect(range(0.0,1.0,100))
mse_traj_list = []
mse_deriv_list = []
complexity_list = []
models_list = []
for threshold in ProgressBar(threshold_scan)
    # Optimizer with given threshold 
    stlsq_opt = pysindy.STLSQ(threshold=threshold)
    # Functions library
    feature_library = pysindy.PolynomialLibrary(1)
    # Model definition and fitting on train data
    model = pysindy.SINDy(feature_library=feature_library, optimizer=stlsq_opt, feature_names=["x", "y", "θ"], t_default=δt)
    model.fit(X_train, t=δt, quiet=true, multiple_trajectories=true)

    # Error computation one test data
    # 1. On trajectory using `simulate()`
    X_test_sim = model.simulate(X_test[1,:], t, integrator="odeint")
    # if(any(X_test_sim .> 1e4)) fill!(X_test_sim, 1e4) end
    mse_traj = metrics.mean_squared_error(X_test, X_test_sim)
    # 2. on derivatives using `score()`
    mse_deriv = model.score(X_test, t=δt, metric=metrics.mean_squared_error)
    push!(mse_traj_list, mse_traj)
    push!(mse_deriv_list, mse_deriv)
    push!(models_list, model)
end

complexity_list = [model.complexity for model in models_list]

# Plotting
# MSE on trajectory prediction w.r.t STLSQ threshold λ
p_acc_λ = plot(threshold_scan, mse_traj_list, label="MSE (trajectory) vs. λ", xlabel="Threshold (λ)", ylabel="Mean squared error", yscale=:log10, ylims=:auto, color=:red, legend=:topright)
# MSE on derivatives prediction w.r.t. STLSQ threshold λ
p_score_λ = plot(threshold_scan, mse_deriv_list, label="MSE (derivatives) vs. λ", xlabel="Threshold (λ)", ylabel="Mean squared error", yscale=:log10, ylims=:auto, color=:red, legend=:topright)
# Model complexity w.r.t. STLSQ threshold λtrue
p_complexity_λ = plot(threshold_scan, complexity_list, label="Model complexity vs. λ", xlabel="Threshold (λ)", ylabel="Model complexity", color=:orange)
# Model complexity w.r.t. MSE on derivatives prediction
p_complexity_acc = scatter(complexity_list, mse_deriv_list, label="Model complexity \nvs. MSE (derivatives)",minorgrid=true, xlabel="Model complexity", ylabel="Mean squared error", yscale=:log10, ylims=:auto, color=:green, legend=:bottomright)

# Everything together
plot(p_acc_λ, p_score_λ, p_complexity_λ, p_complexity_acc, layout=layout = grid(4, 2, heights=[0.5, 0.5, 0.5, 0.5]), xguidefonthalign = :right, yguidefontvalign = :top, guidefontsize=8, tickfontsize=6, legendfontsize=5, size=(600, 400), minorgrid=false)



# To plot best models predictions
# @show best_traj_idx = argmin(mse_traj_list)
# @show best_deriv_idx = argmin(mse_deriv_list)
# @show best_traj_complexity = complexity_list[best_traj_idx]
# @show best_deriv_complexity = complexity_list[best_deriv_idx]
# best_model_traj = models_list[best_traj_idx]
# best_model_deriv = models_list[best_deriv_idx]
# X_test_sim_traj = best_model_traj.simulate(X_test[1,:], t, integrator="odeint")
# X_test_sim_deriv = best_model_deriv.simulate(X_test[1,:], t, integrator="odeint")
# plot(X_test, label=["x" "y" "θ"], color=[:red :orange :blue])
# plot!(X_test_sim_traj, labels=["x_pred_traj" "y_pred_traj" "θ_pred_traj"], color=[:red :orange :blue], linestyle=:dash)
# plot!(X_test_sim_deriv, labels=["x_pred_deriv" "y_pred_deriv" "θ_pred_deriv"], color=[:red :orange :blue], linestyle=:dashdot)