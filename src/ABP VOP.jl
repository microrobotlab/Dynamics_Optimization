using  Plots, LaTeXStrings, Statistics, CSV, DataFrames,CategoricalArrays
 
"""
    polarization_factor(filename::String; averaged::Bool=false)

Compute polarization factor defined as |(1/N).∑(vi/|vi|)|, where N is the number of particles and vi the velocity of particle i. Take simulation output stored in `filename` as input.

`averaged` determines whether the result is averaged over time or given for each timestep. 
"""
function polarization_factor(filename::String; averaged::Bool=false)
    df = CSV.read(filename, DataFrame)
    time = df[!,:Time]
    x = df[!,:xpos]
    y = df[!,:ypos]
    df[!,:N] = categorical(df[!,:N],compress=true) # it sorts out time step data 
    # Group dataframe by values in categorical column
    gdf = groupby(df,:N,sort=true) 

    # Derive velocity from position
    vel_x = [diff(g[!,:xpos]) for g in gdf]
    vel_y = [diff(g[!,:ypos]) for g in gdf]
    # Velocity norms
    vel = [sqrt.(vx.^2 .+ vy.^2) for (vx,vy) in zip(vel_x, vel_y)]
    
    # Polarization factor defined as : |(1/M).∑(vi/|vi|)| for each time step
    # where is the N number of particles, vi the speed of particle i, |.| the norm
    # mean over particles -> (1/M).∑(vi_x/|vi|)
    mean_vel_x = mean([vx./v for (vx,v) in zip(vel_x, vel)]) 
    # mean over particles -> (1/M).∑(vi_y/|vi|)
    mean_vel_y = mean([vy./v for (vy,v) in zip(vel_y, vel)])
    # norm for the final polarization factor 
    # -> √[(1/M).∑(vi_x/|vi|)² + (1/M).∑(vi_y/|vi|)²] = |(1/M).∑(vi/|vi|)|
    pf = sqrt.(mean_vel_x.^2 .+ mean_vel_y.^2)
    
    # Whether we want the p.f. for each time step or averaged over time
    averaged ? (return mean(pf)) : (return pf)
end