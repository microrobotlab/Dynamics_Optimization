using DrWatson
using DataDrivenDiffEq
using ModelingToolkit
using LinearAlgebra
using DataDrivenSparse
using ProgressBars

include("utils_inference.jl")
include(srcdir("ABP output.jl"))

δt = 1e-3
Nt = 1000
times = collect(range(start=0, length=Nt+1, step=δt))


# X_simple = hcat(-2.0t, 3.0t)'
# DX_simple_interpolated, X_simple_interpolated, t_interpolated = collocate_data(X_simple ,t)
# prob = ContinuousDataDrivenProblem(X_simple, t, DX_simple_interpolated)

# @parameters t
# @variables u(t), v(t)
# basis = Basis(monomial_basis([u,v], 1), [u,v])
# println(basis)

# res = solve(prob, basis, STLSQ())
# system = res.basis
# println(system)


# sim_output = run((Nt=Nt, Np=1, L=100., R=1.5, v=5.); wall_condition="open", collision_correction=false, N=16, M=16)
# X = tup2matrix_states(sim_output) 
DX_interpolated, X_interpolated, times_interpolated = collocate_data(X ,times)
prob = ContinuousDataDrivenProblem(X, times, DX_interpolated)

@parameters t
@variables u[1:3](t)
basis = Basis([1.; u[1]; u[2]; u[3]], u, iv=t)
# basis = Basis(Num[1.], [x, y, θ])
println(basis)

res = solve(prob, basis, STLSQ())
system = res.basis
println(system)
println(res)


# using Plots
# plot(X[3,:], label="X")
# plot!(X_interpolated[3,:], label="X_interpolated")
# plot!(DX_interpolated[3,:], label="DX_interpolated")

# X = rand(2, 10)
# N = 10
# t = Float64.(collect(1:N))
# X = Matrix(2.0t') + 0.5 * rand(N)'
# prob = ContinuousDataDrivenProblem(X, t)

# sampler = DataProcessing(split = 0.8, shuffle = true, batchsize = 30)
# λ_list = exp10.(-10:0.1:0)
# λ_list_ok = []

# acc_list = []
# complexity_list = []
# for λ in ProgressBar(λ_list)
#     try
#         res = solve(prob2, basis, STLSQ(λ), options = DataDrivenCommonOptions(data_processing = sampler, digits = 1))
#         push!(acc_list, res.residuals)
#         system = get_basis(res)
#         push!(complexity_list, length(get_parameter_values(system)))
#         push!(λ_list_ok, λ)
#     catch
#         # push!(acc_list, -1.)
#         # println("error : λ=$λ")
#     end
# end

# plot(complexity_list, acc_list, yaxis=:log10)

# plot(X[1,:], label="x")
# plot!(X[2,:], label="y")
# plot!(X[3,:], label="θ")
# plot!(DX_interpolated[1,:], label="∂x")
# plot!(DX_interpolated[2,:], label="∂y")
# plot!(DX_interpolated[3,:], label="∂θ")