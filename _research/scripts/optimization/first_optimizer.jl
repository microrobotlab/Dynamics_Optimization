# ONE INITIAL VERSION FOR THE OPTIMIZATION OF PARTICLES BEHAVIOR ON THE SET OF PHYSICAL PARAMETERS OF THE SIMULATOR
# W.R.T. POLARIZATION FACTOR AVERAGED OVER TIME AND MULTIPLE SIMULATIONS.

# NOTE : In the following, N and M are the number of vertical and horizontal divisions of the space for parallel computing (see `ABP main_parallel.jl`) 
# TODO: maximum possible values for N and M depend on L and R. It is thus useful to have an auto mode to automatically determine 
# those maximum values since R will vary during optimization (see `max_space_division()` function from `utils_simulation.jl` for an 
# implementation that could be used). Fixed values for N and M can lead to overlapping as R varies (see `hardsphere!` function in 
# ABP main_parallel.jl for more details) and make the simulation fail.


using DrWatson
using Optim
using Plots

# `srcdir` from DrWatson provides path to src/ to which we add elements provided in arguments
include(srcdir("ABP output.jl"))
include(srcdir("ABP VOP.jl"))


# ----- OBJECTIVE FUNCTION
# Compute polarization factor averaged over time and multiple simulations (`nb_runs`)
# see `polarization_factor` function of ABP VOP.jl for polarization factor metric definition.
function pf_simulation(parameters::NamedTuple, nb_runs; wall_condition::String)
    # call simulator, `run_multiple` returns the path of the folder containing simulator outputs 
    simulation_folder_path = run_multiple(parameters; nb_runs=nb_runs, N=5, M=5, wall_condition=wall_condition)
    # compute averaged polarization factor from outputed files
    mean_pf_vec = Array{Float64}([])
    for filename in readdir(simulation_folder_path)
        # generate complete path for one file containing a simulation output (simulation folder path + simulation filename)
        data_path = joinpath(simulation_folder_path, filename)
        # averaged over time by polarization_factor (see `polarization_factor` definition in ABP VOP.jl)
        mean_pf = polarization_factor(data_path; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    # remove folder (because a lot will be generated over optimization)
    rm(simulation_folder_path; recursive=true)
    # return /!\ MINUS THE AVERAGE over multiple simulations runs (here we minimize so need to invert)
    return - mean(mean_pf_vec)
end


# ----- SIMULATION HYPERPARAMETERS AND FIXED PARAMETERS (THOSE THAT REMAIN CONSTANT DURING OPTIMIZATION) 
# Wall condition
wall_condition = "periodic"
# Number of runs
nb_runs = 5

# Other fixed parameters 
Np = 50     # number of particles 
Nt = 10000  # number of timesteps
L = 100.    # size of for the space where particles move


# ----- OPTIMIZATION
# Here we optimize w.r.t. radius R and velocity v with fixed other parameters
# initial parameters values based on valid ranges for them
x0 = [1.5 + rand() * (5. - 1.5), 1. + rand() * (15. - 1.)]
opt = Optim.Options(store_trace = true, show_trace=true, extended_trace=true, iterations=10)
# optimization based on `pf_simulation`
res = optimize(x -> pf_simulation((Nt=Nt, Np=Np, L=L, R=x[1], v=x[2]), nb_runs; wall_condition=wall_condition), x0, opt)
# values of the objective function `pf_simulation` taken during optimization
f_trace = Optim.f_trace(res)
plot(- f_trace, labels="averaged polarization factor")
