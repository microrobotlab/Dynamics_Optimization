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
    run(parameters::NamedTuple; wall_condition::String="periodic", animate::Bool=false, animation_filename::String=nothing, animation_stride::Integer=1)

Simulator envelope with parameters given as NamedTuple and wall condition choice. 
Set `animate` to true (false by default) for animation generation based on output, animation filename and stride adjustable
"""
function run(parameters::NamedTuple; wall_condition::String="periodic", animate::Bool=false, animation_filename=nothing, animation_stride::Integer=1)
    if(wall_condition == "periodic")
        simulation_output = multiparticleE(;parameters...)
    elseif(wall_condition in ["squared", "elliptical"])
        simulation_output = multiparticleE_wall(;parameters..., wall_condition)
    else 
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    # animate if required
    if(animate)
        animate!(simulation_output; parameters, wall_condition, animation_filename=animation_filename, stride=animation_stride)
    end

    return simulation_output
end


"""
    run(parameters::NamedTuple, nb_runs::Integer; wall_condition::String="periodic", stride::Integer=1)


Run simulator with parameters given as NamedTuple, wall condition choice and variable number of runs. 
Save output in CSV files (one per run) with adjustable `stride` and return corresponding folder path.
"""
function run(parameters::NamedTuple, nb_runs::Integer; wall_condition::String="periodic", stride::Integer=1)
    experiment_marker = instance_marker(parameters, wall_condition)
    # automatic output file name generation depending on parameters 
    simulation_folder_path = datadir("sims", experiment_marker)
    mkdir(simulation_folder_path)
    # for each run
    for i_run in 1:nb_runs
        file_path = joinpath(simulation_folder_path, experiment_marker * "_" * "run_$i_run")
        simulation_output = run(parameters; wall_condition=wall_condition)
        file_store(simulation_output, parameters.Nt, file_path, stride)
    end
    # return folder name containing run outputs
    return simulation_folder_path
end


"""
    run(param_file_name::String, nb_runs::Integer; wall_condition::String="periodic", stride::Integer=1)

Run simulator with parameters provided by file (`param_file_name`), wall condition choice and variable number of runs. 
Save output in CSV files (one per run) with adjustable `stride` and return corresponding folder path list (one per set of parameters).
"""
function run(param_file_name::String, nb_runs::Integer; wall_condition::String="periodic", stride::Integer=1)
    # extract parameters instances from given file
    parameter_instances = CSV.read(datadir("parameters", param_file_name), DataFrame)
    # to keep track of created folders
    folder_paths_list = Array{String}([])
    # for each parameter instance
    for p in eachrow(parameter_instances)
        simulation_folder_path = run(NamedTuple(p), nb_runs; wall_condition=wall_condition, stride=stride)
        # add the folder corresponding to the last simulation
        push!(folder_paths_list, simulation_folder_path)
    end
    return(folder_paths_list)
end