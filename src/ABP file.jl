# FILE WRITING
# PURPOSE:Store the output of main code ABP output.jl
# METHOD: Dataframes are used 
# INPUT: array of postions and orientation
# OUTPUT: saved excel file (f1)

using CSV, DataFrames
using Dates


function file_store_inc(graph_wall,Nt,pathf, stride=1)
    f1= pathf*"_inc.csv"               # destination file name
    
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

    # touch(f1)

    efg = open(f1, "w")
    close(efg)
    # creating DataFrame
    data = DataFrame(
    N= pnumber,
    Time= time,
    xpos= x,
    ypos= y,
    orientation=θ) 

    # for row in CSV.RowWriter(data)
    #     print(efg, row)
    # end
    # close(efg)

    # # incremental file writing to avoid memory explosion
    # MAX_LINES = 100_000
    # for (i,sub_data) in enumerate(Iterators.partition(data, MAX_LINES))
    #     if i < 2
    #         CSV.write(f1, Tables.columntable(sub_data), append=true, writeheader=(i==1 ? true : false))       
    #     end 
    # end

    # empty!(data)
    
    # println("I am in ABP file")
    return nothing
end


function file_store(graph_wall,Nt,pathf, stride=1)
    f1= pathf*".csv"               # destination file name

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
    
    #creating DataFrame
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
    instance_marker(parameters::NamedTuple, wall_condition::String)

Utility function to help generating discriminant filenames based on time and experiment parameters.
"""
function instance_marker(;parameters=nothing, wall_condition=nothing, collision_correction=nothing)
    # time signature to differentiate at least
    name = Dates.format(now(),"YYYY-mm-dd_HH:MM:SS:sss") 
    # add wall_condition if provided
    if(wall_condition !== nothing) name *= "_" * wall_condition end 
    # add collistion correction indication if provided
    if(collision_correction !== nothing) name *= "_" * "collision-correction=$collision_correction" end 
    # add experiment parameters if provided
    if(parameters !== nothing) name = savename(name, parameters; digits=3) end  
    return name
end