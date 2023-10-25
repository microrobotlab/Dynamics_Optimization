using BenchmarkTools


# ----- THE STRUCTURE THAT WILL BE USED TO SAVE RESULTS
# For more information about BenchmarkTools setup, see "Parameters" section of the BenchmarkTools Manual at
# https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Benchmark-Parameters and definitions at
# https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Manual)
struct BenchmarkingInstance
    # simulation parameters
    Nt::Int64
    Np::Int64                   
    L::Float64                  
	R::Float64  
	v::Float64 	
    # parallelization parameters 
    # number of horizontal/vertical divisions,
    # set to "original" for original main
    N::Union{Int64, String}
    M::Union{Int64, String}
    # parameters and results of the benchmarking
    benchmarking::BenchmarkTools.Trial
end