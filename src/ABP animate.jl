using DrWatson
using Plots 
using CSV, DataFrames
using .Threads
using ProgressBars

"""
    animate!(xyθ; parameters::NamedTuple, wall_condition::String, collision_correction::Bool, animation_filename=nothing, stride::Integer=1, verbose::Bool)

Generate animation from raw simulator output. Adjustable animation `stride`. See simulator parameters and output for more details.
"""
function animate!(xyθ; parameters::NamedTuple, wall_condition::String, collision_correction::Bool, animation_filename=nothing, stride::Integer=1, verbose::Bool)
    iter = 1:stride:parameters.Nt
    # Progress meter if verbose 
    if verbose
        iter = ProgressBar(iter)
        # progress bar description
        set_description(iter, "ANIMATION")
    end
    
    # Animation loop
    anim = @animate for i in iter
        title = "$(parameters.Np) particles, steps $i, "
        if(wall_condition == "elliptical")
            title *= "ellipse a=L/2, b= L/4"
        end
        plot(title=title)

        # Draw particles (circles)
        scatter!(xyθ[1][i][:,1], xyθ[1][i][:,2], aspect_ratio=:equal, lims=(-parameters.L/2, parameters.L/2), markersize=350*parameters.R/parameters.L, marker =:circle,legend=false, title = title)
        
        # Draw the ellipse (for elliptical condition)
        if(wall_condition == "elliptical")
            plot!(parameters.L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π))
        
        # Draw the square (for squared condition)
        elseif(wall_condition == "squared")
            plot!([parameters.L/2], seriestype="vline", color=:red) 
            plot!([-parameters.L/2], seriestype="vline", color=:red)
            plot!([parameters.L/2], seriestype="hline", color=:red) 
            plot!([-parameters.L/2], seriestype="hline", color=:red)
        end

        # Draw particle directions (arrows)
        quiver!(xyθ[1][i][:,1], xyθ[1][i][:,2], quiver=(4*cos.(xyθ[2][i]),4*sin.(xyθ[2][i])), color=:red)
    end

    # Define animation filepath
    if(animation_filename === nothing) # if no name given for animation output
        # Use `instance_marker` to generate filename (see ABP file.jl)
        animation_filename = instance_marker(parameters=parameters, wall_condition=wall_condition, collision_correction=collision_correction) * "_animation.gif"
    end
    # `plotsdir()` gives plots folder path (provided by DrWatson package). We add to the string folder/file names provided in arguments.
    animation_filepath = plotsdir(animation_filename)

    gif(anim, animation_filepath)
end


"""
    plot_trajectory_trace!(xyθ; parameters::NamedTuple, wall_condition::String, plot_filename=nothing, stride::Integer=1, verbose::Bool)

To vizualise particle trajectories in one plot (trace). See simulator parameters and output for more details.
"""
function plot_trajectory_trace!(xyθ; parameters::NamedTuple, wall_condition::String, plot_filename=nothing, stride::Integer=1, verbose::Bool)
    iter = 1:stride:parameters.Nt
    # Progress meter if verbose 
    if verbose
        iter = ProgressBar(iter)
        # progress bar description  
        set_description(iter, "TRAJECTORY TRACE PLOT")
    end

    title = "$(parameters.Np) particles"
    plot(title=title)

    # Same as animate
    if(wall_condition == "elliptical")
        # title *= "ellipse a=L/2, b= L/4"
        plot!(parameters.L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π), title=title) 
    
    # Same as animate
    elseif(wall_condition == "squared")
        plot!([parameters.L/2], seriestype="vline", color=:red)  
        plot!([-parameters.L/2], seriestype="vline", color=:red)
        plot!([parameters.L/2], seriestype="hline", color=:red) 
        plot!([-parameters.L/2], seriestype="hline", color=:red)
    end

    for i in iter
        # Particles (circles)
        scatter!(xyθ[1][i][:,1], xyθ[1][i][:,2], aspect_ratio=:equal, lims=(-parameters.L/2, parameters.L/2), markersize=350*parameters.R/parameters.L, marker =:circle,legend=false, color=:lightblue, alpha=0.2)
        # Directions (arrows)
        quiver!(xyθ[1][i][:,1], xyθ[1][i][:,2], quiver=(4*cos.(xyθ[2][i]),4*sin.(xyθ[2][i])), color=:red, alpha=0.2)
    end

    # Define animation filepath
    if(plot_filename === nothing) # if no name given for plot output
        plot_filename = instance_marker(parameters=parameters, wall_condition=wall_condition) * "_trace_plot.svg"
    end
    # `plotsdir()` gives plots folder path (provided by DrWatson package). We add to the string folder/file names provided in arguments.
    plot_filepath = plotsdir(plot_filename)

    savefig(plot_filepath)
end