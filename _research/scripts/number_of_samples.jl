using DrWatson
@quickactivate "active-brownian-particles"
using Random

include(projectdir("src", "ABP VOP.jl"))


function averaged_runs_pf(filename_list::Vector{String})
    mean_pf_vec = Array{Float64}([])
    for filename in filename_list
        mean_pf = polarization_factor(filename; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    # return the mean over all results from different runs
    # (which are themselves means over time)
    return mean(mean_pf_vec)
end


# FOLDER CONTAINING SIMULATION DATA FOR THE PLOT
EXPERIMENT_FOLDER_NAME = "2023-08-17_17:25:55_L=100.0_Np=50_Nt=10000_R=5.0_v=10.0"
# extract run outputs from given folder
filename_list = readdir(datadir("sims", EXPERIMENT_FOLDER_NAME))
# complete path (broadcast completion with "datadir.(...)")
filename_list = datadir.("sims", EXPERIMENT_FOLDER_NAME, filename_list)

# see NB_ITER times the result for mean polarization factor with NB_SAMPLES samples 
# also gives the relative error compared to the full number of samples (i.e. all runs in folder)
NB_ITER = 10
NB_SAMPLES = 999

# result for full number of samples (i.e. all runs in folder)
full_samples_pf = averaged_runs_pf(filename_list)

for iter=1:NB_ITER
    random_sub_filename_list = shuffle(filename_list)[1:NB_SAMPLES]
    sub_averaged_pf = averaged_runs_pf(random_sub_filename_list)
    println("result : $sub_averaged_pf / error : $(sub_averaged_pf - full_samples_pf) / relative error : $((sub_averaged_pf - full_samples_pf) / full_samples_pf)")
end