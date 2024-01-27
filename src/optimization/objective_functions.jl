# The objective functions defined for the optimizer (see `optimizer.jl`)

using DrWatson
using Optim
using Plots

# `projectdir` from DrWatson provides path to current project. We add to the string folder/file names provided in arguments.
include(projectdir("src", "ABP output.jl"))
include(projectdir("src", "ABP VOP.jl"))


# /!\ For optimization, the number of particles Np is indirectly provided by
# the packing fraction (see utils.jl file) through the function `packing_fraction_to_Np()`

"""
    mean_pf(parameters; wall_condition, collision_correction, nb_runs, N, M)

Compute polarization factor averaged over time and over `nb_runs` simulations using `polarization_factor` function.

`N` and `M` still allow for selecting the number of vertical and horizontal divisions of the space for parallel computation (see 'ABP main_parallel.jl' or 'ABP output.jl')
""" 
function mean_pf(parameters; wall_condition, collision_correction, nb_runs, N, M)
    # Packing fraction to Np to fit simulator parameters
    run_parameters = (
        Nt=parameters.Nt, 
        Np=packing_fraction_to_Np(parameters.packing_fraction, parameters.R, parameters.L), 
        L=parameters.L, 
        R=parameters.R, 
        v=parameters.v
    )
    # Call simulator and get corresponding generated folder path
    simulation_folder_path = run_multiple(
        run_parameters; 
        wall_condition=wall_condition, collision_correction=collision_correction,
        nb_runs=nb_runs, 
        save_stride=1,
        animate=false,
        N=N, M=M,
        verbose=false
    )
    # Compute mean polarization factor from outputed files
    mean_pf_vec = Array{Float64}([])
    for filename in readdir(simulation_folder_path)
        data_path = joinpath(simulation_folder_path, filename)
        # Compute time averaged polarization factor for one simulation 
        mean_pf = polarization_factor(data_path; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    # Remove folder (because a lot will be generated over optimization)
    rm(simulation_folder_path; recursive=true)
    # Mean over all simulations
    # /!\ return MINUS the result because we will minimize 
    # and we want the polarization factor to increase
    return - mean(mean_pf_vec)
end