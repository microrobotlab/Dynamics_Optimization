# CODE FOR PLOTTING PERFORMANCE MEASUREMENTS (see benchmarking_run.jl)

using DrWatson
using Plots
using JSON3
using DataFrames

# `srcdir` from DrWatson provides path to src/. We add to the string folder/file names provided in arguments.
include(srcdir("ABP file.jl"))
include("benchmarking_struct.jl")

using CSV


# VARIABLES TO BE PLOTTED
# choose a variable among simulation parameters (-> x axis of the plot)
simulator_variable = :Np 
# choose a variable from benchmarking results, must be a field of BenchmarkTools.Trial
benchmarking_variable = :times

# STATISTICAL PROCESSING THAT WILL BE APPLIED TO THE BENCHMARKING RAW OUTPUT
# minimum among execution times in seconds (divide by 1.0e9, the output times are in nanoseconds)
# good measurement because it is close to the "real" execution time without noise from external programs
# (see https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Which-estimator-should-I-use?)
F(X) = minimum(X) / 10^9 

# EXTRACT ALL RAW DATA IN ONE SINGLE DATAFRAME 
# all outputs are located at data/processed/benchmarking_outputs
# `datadir` from DrWatson provides path to data directory. We add to the string folder/file names provided in arguments.
data_filename_list = readdir(datadir("processed", "benchmarking_outputs"))
# to have full paths of all files corresponding to the outputs
data_filepath_list = datadir.("processed", "benchmarking_outputs", data_filename_list)

# data extraction in one dataframe
benchmarking_df = DataFrame()
json = ""
for data_filepath in data_filepath_list
    json_string = read(data_filepath, String)
    benchmarking_instance_list = JSON3.read(json_string, Vector{Dict})
    append!(benchmarking_df, DataFrame(benchmarking_instance_list))
end

# Extract one sub-dataframe per value of (N,M) as an ensemble (e.g. (14,15) != (14,16))
# N (resp. M) is the number of vertical (resp. horizontal) divisions of the 
# space for parallel computing as mentioned above (see `hardsphere!()` in ABP main_parallel.jl for example)
benchmarking_df_per_nb_cells = groupby(benchmarking_df, [:N, :M])
plot(legend=:outerright)

# plot 
X = []
Y = []
labels = Vector{String}([])
for g in benchmarking_df_per_nb_cells
    # select variable for x axis
    x = g[:,simulator_variable]
    push!(X, x)
    # select variable for axis and apply F
    y = [F(b[String(benchmarking_variable)]) for b in g[:, :benchmarking]]
    push!(Y, y)
    # extract N and M
    N = g[1, :N]
    M = g[1, :M]
    # for the non-parallel version of the simulator there are no N, M divisions of the space
    N == "original" && M=="original" ? label = "original" : label = "$N x $M" 
    push!(labels, label)
end


# plot results for non parallel version of the simulator
plot!(X[1], Y[1], label=labels[1], color=:black, linestyle=:dash, lw=2)
# for parallel version
plot!(X[2:end], Y[2:end], label=permutedims(labels[2:end]), line_z=(1:length(X[2:end]))', color=cgrad(:thermal), colorbar=false, lw=1)

# SAVED IN PLOTS FOLDER
# `projectdir` from DrWatson provides path to current project. We add to the string folder/file names provided in arguments.
savefig(projectdir("plots", instance_marker() * "_" * "benchmark.svg"))
