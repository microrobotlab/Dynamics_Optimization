using DrWatson
@quickactivate "active-brownian-particles"
using DataFrames, CSV
using ProgressBars

include(projectdir("src", "ABP output.jl"))
include(projectdir("src", "ABP VOP.jl"))

# Output averaged (over time and runs) polarization factor varying to significant parameters : R and v 
# the results are saved in data/procesed/saved_valuations to avoid wasting long computation
NB_SAMPLES_LANDSCAPE = 70
NB_RUNS = 200
WALL_CONDITION = "elliptical"

# values to be tested
parameters_spaces = Dict(
    v => [10.],
    R => [1.5],
    Nt => [10000],
    # pf => 0.1 .+ (0.6 - 0.1)*rand(25), # packing fraction between .1 and .6 so in average 5 points per decimal 
    Np => [10, 20, 40, 50, 100, 200, 300, 400, 500, 700, 850]
    L => [100.0]
    )   

# using a progress bar
for p in ProgressBar(dict_list(parameters_spaces)) 
    # we will average on 200 runs
    mean_pf_vec = Vector{Float64}([])
    generated_folder = run((Nt=Nt, Np=Np, L=L, R=R, v=v); nb_runs=NB_RUNS, wall_condition=WALL_CONDITION)
    for filename in readdir(generated_folder[1])
        data_path = datadir("sims", generated_folder[1], filename)
        mean_pf = polarization_factor(data_path; averaged=true)
        # store time-averaged polarization factor for current run
        push!(mean_pf_vec, mean_pf)
    end
    # average over runs and save in data/procesed/saved_valuations with corresponding parameters
    CSV.write(datadir("processed", "saved_valuations.csv"), DataFrame(Nt=Nt, Np=Np, L=L, R=R, v=v, nb_runs=NB_RUNS, mean_polarization_factor=mean(mean_pf_vec)), append=true) 
    # remove simulation output at each step to avoid memory explosion 
    rm(generated_folder[1]; recursive=true)
end
