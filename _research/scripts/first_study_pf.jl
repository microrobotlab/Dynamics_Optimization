using DrWatson
@quickactivate "active-brownian-particles"
using Plots, Indicators, Temporal, Statistics

include(projectdir("src", "ABP VOP.jl"))

# to visualize polarization factor over time, with moving average
function plot_pf(filename::String, output_name::String, window::Integer)
    pf = polarization_factor(filename)
    # apply moving average
    filtered_pf = sma(pf, n=window)
    plot(pf, label="polarization factor")
    plot!(filtered_pf, label="filtered p.f.", c=:red)
    title!("polarization factor w.r.t. time, with moving average (window=$window)", titlefontsize=9)
    ylims!(0, 1.0) # pf ∈ [0,1] so plot with this range to stay objective
    xlabel!("time")
    ylabel!("polarization factor")
    savefig(projectdir("plots", output_name))    
end


# FOLDER CONTAINING SIMULATION DATA FOR THE PLOT
EXPERIMENT_FOLDER_NAME = "2023-08-22_16:28:08_squared_L=100.0_Np=50_Nt=10000_R=5.0_v=10.0"
data_path = datadir("sims", EXPERIMENT_FOLDER_NAME, EXPERIMENT_FOLDER_NAME * "_run_1.csv")
plot_pf(data_path, EXPERIMENT_FOLDER_NAME * "_pf_plot" * ".svg", 50)


# analysis of the polarization factor 
pf = polarization_factor(data_path)
println("Polarization factor mean : ", mean(pf))
println("Polarization factor standard deviation : ", std(pf))