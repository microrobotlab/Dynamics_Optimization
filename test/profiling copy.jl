using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")
using Test, ProfileView, BenchmarkTools

include(projectdir("test", "ABP main_original.jl"))

#---------------------------------------------------------------------------------------------------------------------
# TESTS RUNNER WITH ARGUMENTS 

# simulation parameters
Nt = 10000; Np = 200; L = 100.; R = 1.5; v = 10.


# # FLAMGRAPH
ProfileView.@profview multiparticleE_wall(Np, L, R, v, Nt);
