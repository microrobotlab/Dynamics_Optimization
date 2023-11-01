using DrWatson
using Plots 
using CSV, DataFrames
using .Threads
using ProgressBars

"""
    animate!(xyθ; parameters::NamedTuple, wall_condition::String, animation_filename=nothing, stride::Integer=1)

Animate simulator output from direct output. Adjustable animation `stride`.
By default, animation filename based on given simulation parameters.
"""
function animate!(xyθ; parameters::NamedTuple, wall_condition::String, animation_filename=nothing, stride::Integer=1, verbose::Bool)
    iter = 1:stride:parameters.Nt
    # progress meter if verbose 
    if verbose
        iter = ProgressBar(iter)
        set_description(iter, "ANIMATION")
    end
    
    anim = @animate for i in iter
        title = "$(parameters.Np) particles, steps $i, "
        if(wall_condition == "elliptical")
            title *= "ellipse a=L/2, b= L/4"
        end
        plot(title=title)

        # for the circles
        scatter!(xyθ[1][i][:,1], xyθ[1][i][:,2], aspect_ratio=:equal, lims=(-parameters.L/2, parameters.L/2), markersize=350*parameters.R/parameters.L, marker =:circle,legend=false, title = title)
        
        # frame for ellipse
        if(wall_condition == "elliptical")
            display(plot!(parameters.L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π))) # ellipse
        # frame for square
        elseif(wall_condition ∈ ["periodic", "squared"])
            display(plot!([parameters.L/2], seriestype="vline", color=:red))  #square
            display(plot!([-parameters.L/2], seriestype="vline", color=:red))
            display(plot!([parameters.L/2], seriestype="hline", color=:red))  #square
            display(plot!([-parameters.L/2], seriestype="hline", color=:red))
        end
        
        # for the arrows
        quiver!(xyθ[1][i][:,1], xyθ[1][i][:,2], quiver=(4*cos.(xyθ[2][i]),4*sin.(xyθ[2][i])), color=:red)
    end
    #marker_z=graph_wall[2][i,1], color=:rainbow, for 

    # define animation filepath
    if(animation_filename === nothing) # if no name given for animation output
        animation_filename = instance_marker(parameters=parameters, wall_condition=wall_condition) * "_animation.gif"
    end
    animation_filepath = plotsdir(animation_filename)

    gif(anim, animation_filepath)
end


function plot_trajectory_trace!(xyθ; parameters::NamedTuple, wall_condition::String, plot_filename=nothing, stride::Integer=1, verbose::Bool)
    iter = 1:stride:parameters.Nt
    # progress meter if verbose 
    if verbose
        iter = ProgressBar(iter)
        set_description(iter, "TRAJECTORY TRACE PLOT")
    end

    title = "$(parameters.Np) particles"
    plot(title=title)
    
    # frame for ellipse
    if(wall_condition == "elliptical")
        # title *= "ellipse a=L/2, b= L/4"
        plot!(parameters.L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π), title=title) # ellipse
    # frame for square
    elseif(wall_condition ∈ ["periodic", "squared"])
        plot!([parameters.L/2], seriestype="vline", color=:red)  #square
        plot!([-parameters.L/2], seriestype="vline", color=:red)
        plot!([parameters.L/2], seriestype="hline", color=:red)  #square
        plot!([-parameters.L/2], seriestype="hline", color=:red)
    end

    for i in iter
        # for the circles
        scatter!(xyθ[1][i][:,1], xyθ[1][i][:,2], aspect_ratio=:equal, lims=(-parameters.L/2, parameters.L/2), markersize=350*parameters.R/parameters.L, marker =:circle,legend=false, color=:lightblue, alpha=0.2)
        # for the arrows
        quiver!(xyθ[1][i][:,1], xyθ[1][i][:,2], quiver=(4*cos.(xyθ[2][i]),4*sin.(xyθ[2][i])), color=:red, alpha=0.2)
    end
    #marker_z=graph_wall[2][i,1], color=:rainbow, for 

    # define animation filepath
    if(plot_filename === nothing) # if no name given for plot output
        plot_filename = instance_marker(parameters=parameters, wall_condition=wall_condition) * "_trace_plot.svg"
    end
    plot_filepath = plotsdir(plot_filename)

    savefig(plot_filepath)
end



# function animate!(data_filepath::String; parameters::NamedTuple, wall_condition::String, animation_filename=nothing, stride::Integer=1)
#     df = CSV.read(data_filepath, DataFrame)
#     gdf = groupby(df, :Time, sort=true)
#     # number of particles given by array size for one time steps
#     Np = nrow(gdf[1])
#     time = unique(df[!,:Time])
#     # as mentionned in docstring, the animation stride is combined with a potential stride on the datafile
#     time_steps = time[1:stride:end]

#     anim = @animate for i in time_steps
#         title = "$(Np) particles, steps $i, "
#         if(wall_condition == "elliptical")
#             title *= "ellipse a=L/2, b= L/4"
#         end
#         # scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(Np) particles, steps $i, ")
#         scatter(gdf[i][!,:xpos], gdf[i][!,:ypos], aspect_ratio=:equal, lims=(-parameters.L/2, parameters.L/2), markersize=350*parameters.R/parameters.L, marker =:circle,legend=false, title = title)
        
#         # frame for ellipse
#         if(wall_condition == "elliptical")
#             plot!(L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π)) # ellipse
#         # frame for square
#         elseif(wall_condition ∈ ["periodic, squared"])
#             plot!([parameters.L/2], seriestype="vline", color=:red)  #square
#             plot!([-parameters.L/2], seriestype="vline", color=:red)
#             plot!([parameters.L/2], seriestype="hline", color=:red)  #square
#             plot!([-parameters.L/2], seriestype="hline", color=:red)
#         end
        
#         # quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
#         quiver!(gdf[i][!,:xpos], gdf[i][!,:ypos],quiver=(4*cos.(gdf[i][!,:orientation]),4*sin.(gdf[i][!,:orientation])), color=:red)
#     end
#     #marker_z=graph_wall[2][i,1], color=:rainbow, for 
    
#     # f1= pathf*".gif"
#     gif(anim, animation_filepath)
# end    

    
#     #------------------------------------------------------------------------------for ellipse-------------------------------------------------------------------------
#     anim = @animate for i = 1:100:Nt
#         scatter(graph_wall[1][i][:,1], graph_wall[1][i][:,2], aspect_ratio=:equal, lims=(-L/2, L/2),markersize=350R/L,marker =:circle,legend=false, title = "$(inside_Np) particles, steps $i, ellipse a=L/2, b= L/4")
#         #plot!(L/2*cos.(-π:0.01:π), L/4*sin.(-π:0.01:π)) # ellipse
#         plot!([L/2], seriestype="vline")  #sqaure
#         plot!([-L/2], seriestype="vline")
#         quiver!(graph_wall[1][i][:,1],graph_wall[1][i][:,2],quiver=(4*cos.(graph_wall[2][i,1]),4*sin.(graph_wall[2][i,1])), color=:red)
    
#     #marker_z=graph_wall[2][i,1], color=:rainbow, for 
    
#     f1= pathf*".gif"
#     gif(anim, f1)
# end

# exp_file = "2023-09-20_21:28:31_periodic_L=100.0_Np=850_Nt=5000_R=1.5_v=10.0"
# data_filepath = datadir("sims", exp_file, exp_file*"_run_1.csv")
# # corresponding parameters, wall condition
# parameters = (L=100.0, Np=100, Nt=10000, R=1.5, v=10.0)
# # and wall condition
# wall_condition = "periodic"
# animation_filepath = plotsdir("last_anim_test.gif")
# animate(data_filepath, animation_filepath, stride=100, parameters=parameters, wall_condition=wall_condition)