# To run optimization with different configurations

using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")

include("objective_functions.jl")
include("utils.jl")
include(projectdir("src","ABP output.jl"))

objective_function = mean_pf
nb_runs = 1
# initializer = 
# optimization_function = 

run((Nt=100,Np=packing_fraction_to_Np(0.2, 1.5, 100.),L=100.,R=1.5,v=10.), nb_runs; wall_condition="periodic", stride=100)