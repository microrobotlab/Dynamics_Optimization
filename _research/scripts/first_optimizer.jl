using DrWatson
using Optim
using Plots

include(srcdir("ABP output.jl"))
include(srcdir("ABP VOP.jl"))


# OBJECTIVE FUNCTION
function pf_simulation(parameters::NamedTuple, nb_runs; wall_condition::String)
    # call simulator
    simulation_folder_path = run(parameters, nb_runs; wall_condition=wall_condition)
    # compute mean polarization factor from outputed files
    mean_pf_vec = Array{Float64}([])
    for filename in readdir(simulation_folder_path)
        data_path = joinpath(simulation_folder_path, filename)
        # averaged over time by polarization_factor
        mean_pf = polarization_factor(data_path; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    # remove folder (because a lot will be generated over optimization)
    rm(simulation_folder_path; recursive=true)
    # return /!\ MINUS THE AVERAGE (here we minimize so need to invert)
    return - mean(mean_pf_vec)
end

# SIMULATION PARAMETERS
# wall condition
wall_condition = "periodic"
# number of runs
nb_runs = 10
# other fixed parameters 
Np = 50; Nt = 10000; L = 100.

# here we optimize w.r.t. radius R and velocity v with fixed other parameters
plot()
for i=1:2
    x0 = [1.5 + rand() * (5. - 1.5), 1. + rand() * (15. - 1.)]
    opt = Optim.Options(store_trace = true, show_trace=true, extended_trace=true, iterations=10)
    res = optimize(x -> pf_simulation((Nt=Nt, Np=Np, L=L, R=x[1], v=x[2]), nb_runs; wall_condition=wall_condition), x0, opt)
    f_trace = Optim.f_trace(res)
    plot!(- f_trace, labels=x0)
end
