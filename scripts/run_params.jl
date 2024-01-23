# SCRIPT TO CALL `run` and `run_multiple` FUNCTIONS FROM "ABP output.jl" WITH PROVIDED PARAMETERS 

using DrWatson

# `srcdir` from DrWatson provides path to src/ to which we add elements provided in arguments
include(srcdir("ABP output.jl"))


#---------------------------------------------------------------------------------------------------------------------
# SIMULATOR PARAMETERS / "HYPERPARAMETERS" (see src/ABP main_parallel.jl)

# Simulator parameters
Nt = 100
Np = 50
L = 100.
R = 1.5
v = 10.

# "hyperparameters"
wall_condition = "periodic"
collision_correction = true
nb_runs = 1


# Parameters for CSV export in data/sims folder:
# to save simulator output (particles position / orientation in time) in CSV file
# /!\ `save` WILL BE CONSIDERED AS true FOR `nb_runs` > 1 REGARLESS OF THE CHOSEN VALUE
save = false
# stride used for CSV file (one postition / orientation in time saved out of save_stride) 
save_stride = 1 

# Animation parameters: 
# to generate gif animation in plots/ folder from simulator output 
# /!\ if true, there will be one animation per run 
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


# RUN
# For one run only, and save set to true, call run which simply returns particles position / orientation 
# through time. For multiple runs, the outputs will be saved in CSV files even if save is set to true
# (i.e. we never directly return multiple simulator outputs through a variable, see "src/ABP output.jl")
if(nb_runs==1 && !save)
    run(
        (Nt=Nt,Np=Np,L=L,R=R,v=v); 
        wall_condition=wall_condition, collision_correction=collision_correction,
        animate=animate, animation_stride=animation_stride, 
        N=N, M=M, 
        verbose=verbose
    );
# Otherwise return the folder path containing saved runs
else
    run_multiple(
        (Nt=Nt,Np=Np,L=L,R=R,v=v); 
        wall_condition=wall_condition, collision_correction=collision_correction,  
        nb_runs=nb_runs, 
        save_stride=save_stride, 
        animate=animate, animation_stride=animation_stride,
        N=N, M=M,
        verbose=verbose
    );
end