# this funcion will calculate the velocity and velocity polarization of particles , i.e., velocity order parameter
# Date: 29-05-2023
# Method: Absoulte velocity calculation at each instant
# Input: File from the output.jl, which has positions and orientation at each instant
# Output: Velocity of each particle at every instant and average velocity of the ensemble

using  Plots, LaTeXStrings, Statistics, CSV, DataFrames,CategoricalArrays
 
function polarization_factor(filename::String; averaged::Bool=false)
    df= CSV.read(filename, DataFrame)
    #steps= df[!, :N]
    time= df[!,:Time]
    x= df[!,:xpos]
    y= df[!,:ypos]
    df[!,:N] = categorical(df[!,:N],compress=true) # it sorts out time step data 
    ## Group dataframe by values in categorical column
    gdf = groupby(df,:N,sort=true) # only 1000 data groups because I have omitted 100 time steps means 1 s
    # vel = [mean(sqrt.((diff(g[!,:xpos])).^2+(diff(g[!,:ypos]).^2))) for g in gdf]  # dr vector and velocity magnitude at each time second
    # vel_x = [mean(diff(g[!, :xpos])) for g in gdf] 
    # vel_y= [mean(diff(g[!,:ypos])) for g in gdf]

    # vp= mean((vel_x./vel).+ (vel_y./vel))         # polarization vector

    # polarization factor |(1/N).∑(vi/|vi|)| for each time step
    # N number of particle, vi speed of particle i, |.| norm
    vel_x = [diff(g[!,:xpos]) for g in gdf]
    vel_y = [diff(g[!,:ypos]) for g in gdf]
    vel = [sqrt.(vx.^2 .+ vy.^2) for (vx,vy) in zip(vel_x, vel_y)]
    # mean over particles -> (1/N).∑(vi_x/|vi|)
    mean_vel_x = mean([vx./v for (vx,v) in zip(vel_x, vel)]) 
    # mean over particles -> (1/N).∑(vi_y/|vi|)
    mean_vel_y = mean([vy./v for (vy,v) in zip(vel_y, vel)])
    # norm for the final polarization factor 
    # -> √[(1/N).∑(vi_x/|vi|)² + (1/N).∑(vi_y/|vi|)²] = |(1/N).∑(vi/|vi|)|
    pf = sqrt.(mean_vel_x.^2 .+ mean_vel_y.^2)
    
    # whether we want the p.f. for each time step or averaged over time
    averaged ? (return mean(pf)) : (return pf)
end


#plot(vel[100],vel[1000])