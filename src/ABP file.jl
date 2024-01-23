# FILE WRITING
# PURPOSE: Store the output of main code ABP output.jl
# METHOD: Dataframes are used 
# INPUT: array of particles positions and orientation from simulation
# OUTPUT: saved excel file (f1)

using CSV, DataFrames, DrWatson
using Dates


"""
    file_store(graph_wall, Nt, pathf, stride=1)

To store simulator output. Stride can be modified so that particle position and orientation are not stored at each time step.
"""
function file_store(graph_wall, Nt, pathf, stride=1)
    f1= pathf*".csv"  # destination file name

    # particle index, timestep, position and orientation
    pnumber=[]
    time=[]
    x=[]
    y=[]
    θ=[]

    for i = 1:stride:Nt
        for j = 1:length(graph_wall[1][i][:,1])
            push!(pnumber,j)
            push!(time, i)
            push!(x, graph_wall[1][i][j,1])
            push!(y, graph_wall[1][i][j,2])
            push!(θ, graph_wall[2][i,1][j])
        end
    end

    touch(f1)

    efg = open(f1, "w")
    
    # Creating DataFrame
    data = DataFrame(
    N= pnumber,
    Time= time,
    xpos= x,
    ypos= y,
    orientation=θ) 
    
    CSV.write(f1, data)
  
    # println("I am in ABP file")
    return nothing
end


"""
    instance_marker(;parameters=nothing, wall_condition=nothing, collision_correction=nothing)

Utility function to help generating discriminant filenames based on time and experiment parameters.
"""
function instance_marker(;parameters=nothing, wall_condition=nothing, collision_correction=nothing)
    # At least time signature to differentiate
    name = Dates.format(now(),"YYYY-mm-dd_HH:MM:SS:sss") 
    # Add wall_condition if provided (see "ABP main_parallel.jl")
    if(wall_condition !== nothing) name *= "_" * wall_condition end 
    # Add collistion correction indication if provided (see "ABP main_parallel.jl")
    if(collision_correction !== nothing) name *= "_" * "collision-correction=$collision_correction" end 
    # Add experiment parameters if provided
    # `savename`, provided by DrWatson package, to generate description string 
    # based on parameter values (added to `name` considered as prefix)
    if(parameters !== nothing) name = savename(name, parameters; digits=3) end  
    return name
end