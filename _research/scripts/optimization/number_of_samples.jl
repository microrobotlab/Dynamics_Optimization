using DrWatson
using Random

# `projectdir` from DrWatson provides path to current project. We add to the string folder/file names provided in arguments.
include(projectdir("src", "ABP VOP.jl"))


# Compute polarization factor averaged over time and multiple simulation runs 
# from data in `filename_list` (similar as `pf_simulation` function in first_optimizer.jl)
function averaged_runs_pf(filename_list::Vector{String})
    mean_pf_vec = Array{Float64}([])
    for filename in filename_list
        # time averaged polarization factor for one run
        mean_pf = polarization_factor(filename; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    # Return the mean over all results from different runs
    return mean(mean_pf_vec)
end


# FOLDER CONTAINING SIMULATION RUNS OUTPUTS
EXPERIMENT_FOLDER_NAME = "2023-10-25_18:54:37:134_periodic_L=100.0_Np=100_Nt=100_R=1.5_v=10.0"
# Extract run outputs from given folder located in data/sims/
# `datadir` function gives path to /data. We add to the string folder/file names provided in arguments.
filename_list = readdir(datadir("sims", EXPERIMENT_FOLDER_NAME))
# Complete paths (broadcast completion with "datadir.(...)")
filename_list = datadir.("sims", EXPERIMENT_FOLDER_NAME, filename_list)


# Results for full number of samples (i.e. all runs in folder)
full_samples_pf = averaged_runs_pf(filename_list)

# We compute polarization factor averaged over time and over several (NB_ITER) random subsets of simulation outputs 
# and observe error and relative error with through averaged polarization factor (ie over all outputs)  
NB_ITER = 10
NB_SAMPLES = 5

for iter=1:NB_ITER
    # random subset of run outputs
    random_sub_filename_list = shuffle(filename_list)[1:NB_SAMPLES]
    # corresponding averaged polarization factor
    sub_averaged_pf = averaged_runs_pf(random_sub_filename_list)
    # result, error and relative error
    println("result : $sub_averaged_pf / error : $(sub_averaged_pf - full_samples_pf) / relative error : $((sub_averaged_pf - full_samples_pf) / full_samples_pf)")
end