using DrWatson
using Printf
using Plots
using Distributions, StatsPlots

# `projectdir` from DrWatson provides path to current project to which we add elements provided in arguments
include(projectdir("src", "ABP VOP.jl"))


# To visualize the amount of variation (using a histogram) of time averaged polarization factor 
# (see src/optimization/ABP VOP.jl file) for different runs -> initial conditions, i.e. initial positions 
# of the particles vary between runs but simulation parameter values remain the same. 
# `folder_name` corresponds to the folder containing simulator output for different runs
function hist_plot_mean_pf(folder_name::String, output_name::String; fit::Bool=false)
    mean_pf_vec = Array{Float64}([])
    # `datadir` function gives path for data folder (provided by DrWatson package). We can add subfolder, "sims" and 
    # folder_name to path to generate complete paths to simulator output files 
    # (see https://juliadynamics.github.io/DrWatson.jl/dev/project/#Navigating-a-Project-1)
    for filename in readdir(datadir("sims", folder_name))
        data_path = datadir("sims", folder_name, filename)
        # Compute time averaged polarization factor
        mean_pf = polarization_factor(data_path; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    nb_samples = length(mean_pf_vec)
    # plot histogram of all obtained values for different runs
    histogram(mean_pf_vec, label="repartition ($nb_samples samples)", normalize=:pdf)
    final_output_name = output_name
    # Fit and add to plot a normal distribution if required
    if fit
        mean_pf_fit = fit_mle(Normal, mean_pf_vec)
        # compute mean / standard deviation to add to plot
        μ = @sprintf "%.3f" mean(mean_pf_fit)
        σ = @sprintf "%.3f" std(mean_pf_fit)
        plot!(mean_pf_fit, label="gaussian fit (μ≈$μ, σ≈$σ)", c=:red, linewidth=2.0)
        # add "_fitted" before extension in filename
        final_output_name = final_output_name[1:findlast('.', output_name)] * "_fitted" * final_output_name[findlast('.', output_name):end]
    end
    # plot
    title!("Repartition of time-averaged polarization factor for different initial conditions", titlefontsize=9)
    xlims!(0, 1.0)
    xlabel!("time-averaged polarization factor")
    ylabel!("probability")
    # `projectdir` from DrWatson provides path to current project. We add to the string folder/file names provided in arguments.
    savefig(projectdir("plots", final_output_name))  
end


# ----- FOLDER CONTAINING SIMULATION DATA FOR THE HISTOGRAM
EXPERIMENT_FOLDER_NAME = "" # choose
# Use EXPERIMENT_FOLDER_NAME + "_mean_pf_hist" in name of the output histogram
hist_plot_mean_pf(EXPERIMENT_FOLDER_NAME, EXPERIMENT_FOLDER_NAME * "_mean_pf_hist" * ".svg", fit=true)