# The objective functions defined for the optimizer
# MUST HAVE THE FORM F(u, p) WITH u the state variables and p other parameters (here those of the simulator)
# SEE https://docs.sciml.ai/Optimization/stable/API/optimization_function/

using DrWatson
@quickactivate "active-brownian-particles"
using Optim
using Plots

include(projectdir("src", "ABP output.jl"))
include(projectdir("src", "ABP VOP.jl"))


# /!\ For optimization, the number of particle Np is given by packing_fraction 
# (see utils.jl file) through the function `packing_fraction_to_Np()`
function mean_pf(parameters; wall_condition, nb_runs, N, M)
    # packing fraction to Np to fit simulator parameters
    run_parameters = (
        Nt=parameters.Nt, 
        Np=packing_fraction_to_Np(parameters.packing_fraction, parameters.R, parameters.L), 
        L=parameters.L, 
        R=parameters.R, 
        v=parameters.v
    )
    # call simulator and get corresponding generated folder path
    simulation_folder_path = run_multiple(
        run_parameters; wall_condition=wall_condition, 
        nb_runs=nb_runs, 
        save_stride=1,
        animate=false,
        N=N, M=M,
        verbose=false
    )
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