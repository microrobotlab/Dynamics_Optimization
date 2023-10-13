using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")
using Plots
using BenchmarkTools, Random
using ProgressBars
using JSON3

include(srcdir("ABP output.jl"))
# for comparison with original version
include("ABP main_original.jl")
include(srcdir("ABP file.jl"))


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


# ----- FIXED PARAMETERS
Nt = 1000
L = 100.
R = 1.5
v = 10.

wall_condition = "periodic"


# ----- PARAMETERS THAT WILL VARY DURING BENCHMARKING
Np_list = [10]
nb_div_cells_list = vcat("original", 1:16)


# FIX RANDOM SEED
Random.seed!(3)

# SETUP PARAMETERS
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 7200 # max time per sample (in seconds)
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10


# ----- LOOP
folderpath = datadir("processed", "benchmarking_outputs")
# time marker for filename (see instance_marker function at `src/ABP file.jl`)
filename = instance_marker() * "_benchmarking" * ".json" 
filepath = joinpath(folderpath, filename)

open(filepath, "w") do f
    iter = ProgressBar(enumerate(Iterators.product(Np_list, nb_div_cells_list)))
    println(f, '[')
    for (i, (Np, nb_div_cells)) in iter
        if(nb_div_cells == "original")
            set_description(iter, "BENCHMARKING (original version: Np=$Np)")
            b = @benchmarkable multiparticleE($Np, $L, $R, $v, $Nt)
        else
            set_description(iter, "BENCHMARKING (parallel version: Np=$Np, N=$nb_div_cells, M=$nb_div_cells)")
            b = @benchmarkable run((Nt=$Nt,Np=$Np,L=$L,R=$R,v=$v); wall_condition=$wall_condition, animate=false, N=$nb_div_cells, M=$nb_div_cells, verbose=false)
        end
        t = BenchmarkTools.run(b)
        JSON3.write(f, BenchmarkingInstance(Nt, Np, L, R, v, 
                                        nb_div_cells, nb_div_cells,
                                        t))
        if(i<length(iter)) 
            println(f, ',')
        end
    end
    print(f, "\n]")
end

# loop completed
println("EVALUATION COMPLETED (results saved at $filepath")