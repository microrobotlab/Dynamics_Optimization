using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")
using Plots
using JSON3
using DataFrames

include("benchmarking_struct.jl")
include(srcdir("ABP file.jl"))



using CSV



# VARIABLES TO BE PLOTTED
# choose a variable among simulation parameters (-> x-axis of the plot)
simulator_variable = :Np 
# choose a variable from benchmarking results, must be a field of BenchmarkTools.Trial
benchmarking_variable = :times

# FUNCTION THAT WILL BE APPLIED TO THE BENCHMARKING DATA
F(X) = minimum(X) / 10^9 # minimum among execution times in seconds (divide by 1.0e9, the output times are in nanoseconds)

# WHERE THE DATA IS LOCATED 
data_filename_list = readdir(datadir("processed", "benchmarking_outputs"))
data_filepath_list = datadir.("processed", "benchmarking_outputs", data_filename_list)


# data extraction as one dataframe
benchmarking_df = DataFrame()
json = ""
for data_filepath in data_filepath_list
    json_string = read(data_filepath, String)
    benchmarking_instance_list = JSON3.read(json_string, Vector{Dict})
    append!(benchmarking_df, DataFrame(benchmarking_instance_list))
end

# extract one sub-dataframe per value of (N,M) as an ensemble (e.g. (14,15) != (14,16))
# N (resp. M) is the number of vertical (resp. horizontal) slices of the space during parallelization
benchmarking_df_per_nb_cells = groupby(benchmarking_df, [:N, :M])
plot(legend=:outerright)

# test = DataFrame()

X = []
Y = []
labels = Vector{String}([])
for g in benchmarking_df_per_nb_cells
    x = g[:,simulator_variable]
    push!(X, x)
    # Y = g[:,:v ]
    # Z = g[:, :benchmarking]
    # println(benchmarking_df_per_nb_cells)
    # println(X)
    # println([F(b[String(benchmarking_variable)]) for b in g[:, :benchmarking]])
    y = [F(b[String(benchmarking_variable)]) for b in g[:, :benchmarking]]
    push!(Y, y)
    # println(Z[1].times)
    N = g[1, :N]
    M = g[1, :M]
    N == "original" && M=="original" ? label = "original" : label = "$N x $M" 
    push!(labels, label)
    # plot!(X, Y, yaxis=:log, label=label, legend=:outerright)

    # append!(test, g)
end

# CSV.write("test.csv", test)
# plot!()

# original
plot!(X[1], Y[1], label=labels[1], color=:black, linestyle=:dash, lw=2)
# with parallel
plot!(X[2:end], Y[2:end], label=permutedims(labels[2:end]), line_z=(1:length(X[2:end]))', color=cgrad(:thermal), colorbar=false, lw=1)

savefig(projectdir("plots", instance_marker() * "_" * "benchmark.svg"))
