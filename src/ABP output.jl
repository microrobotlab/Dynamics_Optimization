using DrWatson
using Plots,Distances,NaNStatistics,CSV, DataFrames

# `srcdir` from DrWatson provides path to src/. We add to the string folder/file names provided in arguments.
include(srcdir("ABP main_parallel.jl"))
include(srcdir("ABP file.jl"))
include(srcdir("ABP analysis.jl"))
include(srcdir("ABP multifolder.jl"))
include(srcdir("ABP animate.jl"))


"""
    run(
    parameters::NamedTuple; 
    wall_condition::String="periodic", collision_correction::Bool=true,
    animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

Simulator envelope with parameters given as NamedTuple.  

`wall_condition` can be 'open', 'periodic', 'squared', or 'elliptical'. 
`collision_correction` to unable / disable correction of the collistions in simulation. 
`N` and `M` for vertical (resp. horizontal) number of divisions of the space for parallel computation of the simulation.
`animate` (false by default) for animation generation based on output. Adjustable animation filename (based on experiment parameters by default) and animation stride (display one out of n steps in animation).
`verbose` to see simulation and animation progress.
Return tabular data corresponding to the simulation.
"""
function run(
    parameters::NamedTuple; 
    wall_condition::String="periodic", collision_correction::Bool=true,
    animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

    # Wall condition
    if(wall_condition in ["open", "periodic"])
        simulation_output = multiparticleE(;parameters..., wall_condition=wall_condition, collision_correction=collision_correction, N=N, M=M, verbose=verbose)
    elseif(wall_condition in ["squared", "elliptical"])
        simulation_output = multiparticleE_wall(;parameters..., wall_condition=wall_condition, collision_correction=collision_correction, N=N, M=M, verbose=verbose)
    else 
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    # Animate if required
    if(animate)
        animate!(simulation_output; parameters=parameters, wall_condition=wall_condition, collision_correction=collision_correction, animation_filename=animation_filename, stride=animation_stride, verbose=verbose)
    end

    return simulation_output
end


"""
    run_multiple(
    parameters::NamedTuple; 
    wall_condition::String="periodic", collision_correction::Bool=true,
    nb_runs::Integer=1, 
    save_stride::Integer=1, # stride for timesteps in file export
    animate::Bool=false, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

To use `run()` function in batch, with saved output in CSV file in 'data/sims/'. Number of runs and stride of the output file (save one out of n state snapshots) tunable.
Return the corresponding folderpath.
"""
function run_multiple(
    parameters::NamedTuple; 
    wall_condition::String="periodic", collision_correction::Bool=true,
    nb_runs::Integer=1, 
    save_stride::Integer=1, # stride for timesteps in file export
    animate::Bool=false, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )
    
    # Used for the name of simulator output to make it recognizable and unique
    experiment_marker = instance_marker(parameters=parameters, wall_condition=wall_condition, collision_correction=collision_correction)
    # automatic output file name generation depending on parameters 
    # `datadir` from DrWatson provides path to data directory. We add to the string folder/file names provided in arguments.
    simulation_folder_path = datadir("sims", experiment_marker)
    mkdir(simulation_folder_path)
    # for each run
    for i_run in 1:nb_runs
        filename = experiment_marker * "_" * "run_$i_run"
        file_path = joinpath(simulation_folder_path, filename)
        # the animation filename will be based corresponding CSV filename
        simulation_output = run(
            parameters; 
            wall_condition=wall_condition, collision_correction=collision_correction,
            animate=animate, animation_filename=filename*"_animation.gif", animation_stride=animation_stride, 
            N=N, M=M,
            verbose=verbose
        )
        file_store(simulation_output, parameters.Nt, file_path, save_stride)
    end
    # return folder name containing run outputs
    return simulation_folder_path
end


"""
    run_from_file(
    param_file_name::String; 
    wall_condition::String="periodic", collision_correction::Bool=true,
    nb_runs::Integer=1, 
    save_stride::Integer=1, # stride for timesteps in file export
    animate::Bool=false, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

To run `run_multiple()` with simulation parameters stored in a CSV file (`param_file_name`). Number of runs also tunable (`nb_runs`).
Return a list of folderpaths (`run_multiple()` output for each set of parameters in `param_file_name`).
"""
function run_from_file(
    param_file_name::String; 
    wall_condition::String="periodic", collision_correction::Bool=true,
    nb_runs::Integer=1, 
    save_stride::Integer=1, # stride for timesteps in file export
    animate::Bool=false, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

    # extract parameters instances from given file
    # `datadir` from DrWatson provides path to data directory. We add to the string folder/file names provided in arguments.
    parameter_instances = CSV.read(datadir("parameters", param_file_name), DataFrame)
    # to keep track of created folders
    folder_paths_list = Array{String}([])
    # for each parameter instance
    for p in eachrow(parameter_instances)
        simulation_folder_path = run_multiple(
            NamedTuple(p); 
            wall_condition=wall_condition, collision_correction=collision_correction,
            nb_runs=nb_runs, 
            save_stride=save_stride, 
            animate=animate, animation_stride=animation_stride,
            N=N, M=M,
            verbose=verbose
        )
        # add folder to the list 
        push!(folder_paths_list, simulation_folder_path)
    end
    return(folder_paths_list)
end