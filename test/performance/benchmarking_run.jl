# CODE FOR RUNTIME EVALUATION OF THE SIMULATOR WITH DIFFERENT CONFIGURATIONS
# USED IN PARTICULAR TO COMPARE PERFORMANCES BETWEEN NON-PARALLEL AND PARALLEL VERSIONS OF THE 
# SIMULATOR WITH DIFFERENT NUMBER OF DIVISIONS OF THE SPACE FOR PARALLEL COMPUTING (see src/ABP main_parallel.jl)
# BenchmarkTools package is used (see https://juliaci.github.io/BenchmarkTools.jl/stable/)

# OUTPUT IS SAVED IN data/processed/benchmarking_outputs FOLDER USING JSON FORMAT WITH CUSTOM STRUCTURE 
# DEFINED IN benchmarking_struct.jl FILE.

using DrWatson
using Plots
using BenchmarkTools, Random
using ProgressBars
using JSON3

# `srcdir` from DrWatson provides path to src/. We add to the string folder/file names provided in arguments.
include(srcdir("ABP file.jl"))
include(srcdir("ABP output.jl"))
# for comparison with original version
# same as `srcdir` above, with general project folder
include(projectdir("test", "ABP main_original.jl"))
include(srcdir("ABP file.jl"))
include("benchmarking_struct.jl")


# ----- FIXED PARAMETERS
Nt = 1000
L = 100.
R = 1.5
v = 10.

# /!\ IN THE FIRST VERSION OF THE SIMULATOR, BOUNDARY CONDITIONS IS NOT A PARAMETER
# LIKE IN THE NEW (PARALLEL) VERSION. IF YOU WANT TO COMPARE PARALLEL AND NON-PARALLEL
# VERSIONS OF THE SIMULATOR, set `wall_condition` AND PUT THE SAME VALUE IN "ABP main_original.jl"
# -> uncomment `wall_condition!` / `elliptical_wall_condition!` for squared / elliptical boundary 
# in `update_wall` and use `multiparticleE_wall` instead of `multiparticleE` (see "ABP main_parallel.jl")
# in the following. 
# uncomment `periodic_BC_array!` for periodic boundary in `update`
# comment every boundary condition function in update to simulate open boundaries
# Also `initABPE` must be chosen properly. ELLIPSE for elliptical boundary and SQUARE for others.
wall_condition = "periodic"


# ----- PARAMETERS THAT WILL VARY DURING BENCHMARKING
# Number of particles 
Np_list = [10, 50, 100]
# Number of divisions of the space (horizontally and vertically) for parallel execution
# (see src/ABP main_parallel.jl)
nb_div_cells_list = vcat(10:16)


# FIX RANDOM SEED
Random.seed!(3)

# SETUP PARAMETERS (see https://juliaci.github.io/BenchmarkTools.jl/stable/manual/)
# Number of executions for one runtime measurement, i.e. one sample
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1  
# Number of samples for a specific parameters configuration
# Statistical (mean, standard deviation, etc.) study is then based on the different samples obtained  
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10  
# Maximum time for one sample in seconds
# The number of samples can be reduced (keeping at least one sample) if exceeded
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 7200    


# ----- RUNTIME MEASUREMENTS
# `datadir` from DrWatson provides path to data directory. We add to the string folder/file names provided in arguments.
folderpath = datadir("processed", "benchmarking_outputs")
# time marker for filename (see `instance_marker()` function at `src/ABP file.jl`)
filename = instance_marker() * "_benchmarking" * ".json" 
filepath = joinpath(folderpath, filename)

open(filepath, "w") do f
    # To vizualise progress
    iter = ProgressBar(enumerate(Iterators.product(Np_list, nb_div_cells_list)))
    # Start by writing '[' in the output file for JSON format (to have a list)
    println(f, '[')
    for (i, (Np, nb_div_cells)) in iter
        # non parallel version of the simulator
        if(nb_div_cells == "original")
            # just for progress bar description
            set_description(iter, "BENCHMARKING (original version: Np=$Np)")
            # runable measurement of exectution time (see https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Benchmarking-basics)
            b = @benchmarkable multiparticleE($Np, $L, $R, $v, $Nt)
        else
            # progress bar description
            set_description(iter, "BENCHMARKING (parallel version: Np=$Np, N=$nb_div_cells, M=$nb_div_cells)")
            # runable measurement of exectution time (see https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Benchmarking-basics)
            b = @benchmarkable run((Nt=$Nt,Np=$Np,L=$L,R=$R,v=$v); wall_condition=$wall_condition, animate=false, N=$nb_div_cells, M=$nb_div_cells, verbose=false)
        end
        # run measurement
        t = BenchmarkTools.run(b)
        # write result in the file using custom struct (see benchmarking_struct.jl)
        JSON3.write(f, BenchmarkingInstance(Nt, Np, L, R, v, 
                                        nb_div_cells, nb_div_cells,
                                        t))
        # close list for JSON format
        if(i<length(iter)) 
            println(f, ',')
        end
    end
    print(f, "\n]")
end

# loop completed
println("EVALUATION COMPLETED (results saved at $filepath")