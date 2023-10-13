using Plots
using DataFrames

include("benchmarking.jl")

data_filename = "2023-10-10_16:53:10:665_benchmarking.json"
data_folderpath = datadir("processed", "benchmarking_outputs")
data_filepath = joinpath(data_folderpath, data_filename)