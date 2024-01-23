using DrWatson
using Plots, Indicators, Temporal, Statistics

# `projectdir` from DrWatson provides path to current project to which we add elements provided in arguments
include(projectdir("src", "ABP VOP.jl"))

# Plot polarization factor over time without time averaging (see src/optimization/ABP VOP.jl).
# Also apply moving average and add it to the plot -> we want a better observation of potential
# polarization factor trends by removing noise. `window` to choose window size of moving average
# (see `SMA` in https://github.com/dysonance/Indicators.jl/blob/master/README.md)
function plot_pf(filename::String, output_name::String, window::Integer)
    # Compute polarization factor without time averaging
    pf = polarization_factor(filename)
    # Apply moving average with provided window size
    filtered_pf = sma(pf, n=window)
    plot(pf, label="polarization factor")
    # Add moving average to the plot
    plot!(filtered_pf, label="filtered p.f.", c=:red)
    title!("polarization factor w.r.t. time, with moving average (window=$window)", titlefontsize=9)
    ylims!(0, 1.0) # polarization factor ∈ [0,1] so plot in this range to stay objective
    xlabel!("time")
    ylabel!("polarization factor")
    # `projectdir` from DrWatson provides path to current project to which we add elements provided in arguments
    savefig(projectdir("plots", output_name))    
end


# FOLDER CONTAINING SIMULATION DATA FOR THE PLOT
# Choose a folder in data/sims containing simulator output files
EXPERIMENT_FOLDER_NAME = "2023-11-01_15:55:49:445_periodic_collision-correction=true_L=100.0_Np=283_Nt=1000_R=1.5_v=5.0"
# `datadir` function gives path for data folder (provided by DrWatson package). We can add subfolder, "sims",
# EXPERIMENT_FOLDER_NAME and EXPERIMENT_FOLDER_NAME to path to generate complete paths to simulator output file
# (see https://juliadynamics.github.io/DrWatson.jl/dev/project/#Navigating-a-Project-1)
# "_run_1.csv" to select first file in the folder with multiple simulation outputs 
data_path = datadir("sims", EXPERIMENT_FOLDER_NAME, EXPERIMENT_FOLDER_NAME * "_run_1.csv")
plot_pf(data_path, EXPERIMENT_FOLDER_NAME * "_pf_plot" * ".svg", 50)


# Analysis of the polarization factor 
pf = polarization_factor(data_path)
println("Polarization factor mean : ", mean(pf))
println("Polarization factor standard deviation : ", std(pf))