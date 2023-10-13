# PURPOSE: Output of ABP main 
# all codes in repository are complied in this function
# VARIABLES: Destination folder path and filename

using DrWatson
@quickactivate "active-brownian-particles"
using Plots,Distances,NaNStatistics,CSV, DataFrames

include(srcdir("ABP main_parallel.jl"))
include(srcdir("ABP file.jl"))
include(srcdir("ABP analysis.jl"))
include(srcdir("ABP multifolder.jl"))
include(srcdir("ABP animate.jl"))


"""
    run(parameters::NamedTuple; wall_condition::String="periodic", animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1)

Simulator envelope with parameters given as NamedTuple and wall condition choice. 
Set `animate` to true (false by default) for animation generation based on output, animation filename and adjustable stride (`animation_stride`).
"""
function run(
    parameters::NamedTuple; wall_condition::String="periodic",
    animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

    if(wall_condition == "periodic")
        simulation_output = multiparticleE(;parameters..., N=N, M=M, verbose=verbose)
    elseif(wall_condition in ["squared", "elliptical"])
        simulation_output = multiparticleE_wall(;parameters..., wall_condition, N=N, M=M, verbose=verbose)
    else 
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    # animate if required
    if(animate)
        animate!(simulation_output; parameters, wall_condition, animation_filename=animation_filename, stride=animation_stride, verbose=verbose)
    end

    return simulation_output
end


"""
    run_multiple(parameters::NamedTuple; wall_condition::String="periodic", 
    nb_runs::Integer=1, 
    save_stride::Integer=1,
    animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1)

Run simulator with parameters given as NamedTuple, wall condition choice and variable number of runs. 
Save output in CSV files (one per run) with adjustable timestep stride (`save_stride`) and return corresponding folder path.
Set `animate` to true (false by default) for animation generation based on output, animation filename and adjustable stride (`animation_stride`).
"""
function run_multiple(
    parameters::NamedTuple; wall_condition::String="periodic", 
    nb_runs::Integer=1, 
    save_stride::Integer=1, # stride for timesteps in file export
    animate::Bool=false, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )
    
    experiment_marker = instance_marker(parameters=parameters, wall_condition=wall_condition)
    # automatic output file name generation depending on parameters 
    simulation_folder_path = datadir("sims", experiment_marker)
    mkdir(simulation_folder_path)
    # for each run
    for i_run in 1:nb_runs
        filename = experiment_marker * "_" * "run_$i_run"
        file_path = joinpath(simulation_folder_path, filename)
        # the animation filename will be the same as the corresponding CSV file
        simulation_output = run(parameters; wall_condition=wall_condition, 
                                animate=animate, animation_filename=filename*"_animation.gif", animation_stride=animation_stride, 
                                N=N, M=M,
                                verbose=verbose)
        file_store(simulation_output, parameters.Nt, file_path, save_stride)
    end
    # return folder name containing run outputs
    return simulation_folder_path
end


"""
    run_from_file(param_file_name::String; wall_condition::String="periodic", 
    nb_runs::Integer=1, 
    save_stride::Integer=1,
    animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1)

Run simulator with file stored parameters (`param_file_name`), wall condition choice and variable number of runs. 
Save output in CSV files (one per run) with adjustable timestep stride (`save_stride`) and return corresponding folder path list (one per set of parameters).
Set `animate` to true (false by default) for animation generation based on output, animation filename and adjustable stride (`animation_stride`).
"""
function run_from_file(
    param_file_name::String; wall_condition::String="periodic", 
    nb_runs::Integer=1, 
    save_stride::Integer=1, # stride for timesteps in file export
    animate::Bool=false, animation_stride::Integer=1,
    N::Integer, M::Integer,
    verbose::Bool=true
    )

    # extract parameters instances from given file
    parameter_instances = CSV.read(datadir("parameters", param_file_name), DataFrame)
    # to keep track of created folders
    folder_paths_list = Array{String}([])
    # for each parameter instance
    for p in eachrow(parameter_instances)
        simulation_folder_path = run_multiple(NamedTuple(p); wall_condition=wall_condition, 
                                              nb_runs=nb_runs, 
                                              save_stride=save_stride, 
                                              animate=animate, animation_stride=animation_stride,
                                              N=N, M=M,
                                              verbose=verbose)
        # add folder to the list 
        push!(folder_paths_list, simulation_folder_path)
    end
    return(folder_paths_list)
end