# SCRIPT TO CALL `run_multiple` FUNCTIONS FROM "ABP output.jl" WITH PARAMETERS PROVIDED IN A FILE
# MUST BE IN data/parameters/ FOLDER


using DrWatson

# `srcdir` from DrWatson provides path to src/. We add to the string folder/file names provided in arguments.
include(srcdir("ABP output.jl"))


#---------------------------------------------------------------------------------------------------------------------
# SIMULATOR PARAMETERS / "HYPERPARAMETERS" (see src/ABP main_parallel.jl)

# NAME OF THE FILE CONTAINING SIMULATOR PARAMETERS; THE FILE MUST BE IN data/parameters/ FOLDER.
param_file_name = "instance1.csv"

# "hyperparameters"
wall_condition = "periodic"
collision_correction = true
# /!\ `nb_runs` RUNS FOR EACH SET OF PARAMETERS PROVIDED IN `param_file_name`
nb_runs = 1


# Parameters for CSV export in data/sims folder :
# unlike in scripts/run_params.jl, in this version simulator output 
# (particles position / orientation in time) is automatically saved in a CSV file.
# stride used for CSV file (one postition / orientation in time saved out of save_stride) 
save_stride = 1 

# Animation parameters: 
# to generate gif animation in plots/ folder from simulator output 
# /!\ if true there will be `nb_runs`*(number of sets of parameters in `param_file_name`) animations 
animate = false
# stride for the animation
animation_stride = 1
# by default `animation_filename` will be based on a marker using current time and chosen parameters, 
# similarly to CSV data filename (see `instance_marker` function in "ABP file.jl"). If `save` = true,
# `animation_filename` will be automatically based on the CSV data filename (see `run_multiple` 
# in src/ABP output.jl and `animate!` in src/ABP animate.jl)
animation_filename = nothing; 

# Parallelization parameters: (N, M) number of cell divisions (horizontally /vertically) for parallel computing
# see src/ABP main_parallel.jl
N = 16; M = 16;

# Display progress bar for simulation
verbose = true


# RUN (see "src/ABP output.jl")
run_from_file(
    param_file_name;
    wall_condition=wall_condition, collision_correction=collision_correction, 
    nb_runs=nb_runs, 
    save_stride=save_stride,
    animate=animate, animation_stride=animation_stride,
    N=N, M=M,
    verbose=verbose
);