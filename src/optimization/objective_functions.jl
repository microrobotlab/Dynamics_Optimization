# The objective functions defined for the optimizer

using DrWatson
@quickactivate "active-brownian-particles"
using Optim
using Plots

include(projectdir("src", "ABP output.jl"))
include(projectdir("src", "ABP VOP.jl"))


function mean_pf(parameters::NamedTuple, nb_runs; wall_condition::String)
    # call simulator and get corresponding generated folder path
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