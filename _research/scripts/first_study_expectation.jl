using DrWatson
using Printf
using Plots
using Distributions, StatsPlots

include(projectdir("src", "ABP VOP.jl"))


# to visualize the repartition (histogram) of the averaged polarization factor 
# for different runs (i.e. different initial conditions)
function hist_plot_mean_pf(folder_name::String, output_name::String; fit::Bool=false)
    mean_pf_vec = Array{Float64}([])
    for filename in readdir(datadir("sims", folder_name))
        data_path = datadir("sims", folder_name, filename)
        mean_pf = polarization_factor(data_path; averaged=true)
        push!(mean_pf_vec, mean_pf)
    end
    nb_samples = length(mean_pf_vec)
    histogram(mean_pf_vec, label="repartition ($nb_samples samples)", normalize=:pdf)
    final_output_name = output_name
    # fit and add to plot a normal distribution if required
    if fit
        mean_pf_fit = fit_mle(Normal, mean_pf_vec)
        μ = @sprintf "%.3f" mean(mean_pf_fit)
        σ = @sprintf "%.3f" std(mean_pf_fit)
        plot!(mean_pf_fit, label="gaussian fit (μ≈$μ, σ≈$σ)", c=:red, linewidth=2.0)
        # add "_fitted" before extension in filename
        final_output_name = final_output_name[1:findlast('.', output_name)] * "_fitted" * final_output_name[findlast('.', output_name):end]
    end
    title!("Repartition of time-averaged polarization factor for different initial conditions", titlefontsize=9)
    xlims!(0, 1.0)
    xlabel!("time-averaged polarization factor")
    ylabel!("probability")
    savefig(projectdir("plots", final_output_name))  
end


# FOLDER CONTAINING SIMULATION DATA FOR THE PLOT
EXPERIMENT_FOLDER_NAME = "2023-08-22_16:22:20_squared_L=100.0_Np=50_Nt=10000_R=0.55_v=10.0"
hist_plot_mean_pf(EXPERIMENT_FOLDER_NAME, EXPERIMENT_FOLDER_NAME * "_mean_pf_hist" * ".svg", fit=true)