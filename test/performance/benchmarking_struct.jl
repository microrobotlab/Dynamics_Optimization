using BenchmarkTools


# ----- STRUCTURE THAT WILL BE USED TO SAVE PERFORMANCE MEASURE RESULTS
# For more information about BenchmarkTools setup, see "Parameters" section of the BenchmarkTools Manual at
# https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Benchmark-Parameters and definitions at
# https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Manual)
struct BenchmarkingInstance
    # Simulation parameters
    Nt::Int64
    Np::Int64                   
    L::Float64                  
	R::Float64  
	v::Float64 	
    # Parallelization parameters 
    # number of horizontal/vertical divisions (see `ABP main_parallel.jl`),
    # Set to "original" to use first (non parallel) version of simulator main 
    N::Union{Int64, String}
    M::Union{Int64, String}
    # Parameters and results of the benchmarking
    benchmarking::BenchmarkTools.Trial
end