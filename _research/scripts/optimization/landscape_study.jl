# WE WANT TO OPTIMIZE COLLECTIVE MOTION OF THE SIMULATED PARTICLES. OPTIMIZATION IS BASED ON AVERAGED 
# POLARIZATION FACTOR see src/ABP VOP.jl AS THE OBJECTIVE FUNCTION. HERE THIS FACTOR IS MEASURED FOR 
# DIFFERENT POINTS IN PARAMETERS SPACE; THE RESULTS ARE PLOTTED IN landscape_plot.jl


using DrWatson
using DataFrames, CSV
using ProgressBars

# `projectdir` from DrWatson provides path to current project to which we add elements provided in arguments
include(projectdir("src", "ABP output.jl"))
include(projectdir("src", "ABP VOP.jl"))


# Polarization factor is averaged over time and multiple simulation runs (for one specific set of parameters)
NB_RUNS = 3
# Simulator parameters 
WALL_CONDITION = "elliptical"

# ----- VALUES THAT WILL BE TESTED IN PARAMETERS SPACE SOME VARIABLES CAN REMAIN FIXED 
# (see `ABPE2` struct in src/ABP main_parallel.jl to know simulator parameters)
parameters_spaces = Dict(
    :v => [10.],                             # particle speed 
    :R => [0.55, 1.5, 5.0],                  # particle radius
    :Nt => [10000],                          # number of timesteps
    :Np => [2],                              # number of particles
    :L => [100.0]                            # size of the space (reduced for elliptical boundary conditions)
    )   

# Generate list of all combinations of parameters from parameters_spaces
# (see https://juliadynamics.github.io/DrWatson.jl/dev/run&list/#DrWatson.dict_list)
all_params = dict_list(parameters_spaces)

# Using a progress bar
for params in ProgressBar(all_params)
    # Store time averaged polarization factor for each run (there will be NB_RUNS runs) 
    mean_pf_vec = Vector{Float64}([])
    
    # Obtain output from simulator (folder with output files from multiple runs)
    generated_folder = run_multiple((Nt=params[:Nt], Np=params[:Np], L=params[:L], R=params[:R], v=params[:v]); nb_runs=NB_RUNS, wall_condition=WALL_CONDITION)
    for filename in readdir(generated_folder)
        # to build path to each output file (corresponding to one run); `datadir` function gives path 
        # for data folder to which we can add "generated_folder/" and filename 
        # (see https://juliadynamics.github.io/DrWatson.jl/dev/project/#Navigating-a-Project-1) 
        data_path = datadir("sims", generated_folder, filename)
        # compute time averaged polarization factor for one run
        mean_pf = polarization_factor(data_path; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    
    # average over runs and save in data/procesed/saved_valuations with corresponding parameters
    CSV.write(datadir("processed", "saved_valuations.csv"), DataFrame(Nt=params[:Nt], Np=params[:Np], L=params[:L], R=params[:R], v=params[:v], nb_runs=NB_RUNS, mean_polarization_factor=mean(mean_pf_vec)), append=true) 
    # remove simulation output at each step to avoid memory explosion 
    rm(generated_folder[1]; recursive=true)
end
