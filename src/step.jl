# DEFINE DIFFERENT STEP FUNCTIONS TO BE USED IN THE SIMULATOR TO VARY PARTICLE DYNAMICS.
# USEFUL FOR TESTING PARTICLES DYNAMICS INFERENCE METHODS WITH DIFFERENT KNOWN DYNAMICS (see inference folder).  


# Original Active Brownian Particles model
function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    if size(position(abpe),2) == 2
        δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)]
        δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)
    else
        println("No step method available")
    end
    return (δp, δθ)
end

# # Simplistic model with straight particle trajectories
# function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
#     if size(position(abpe),2) == 2
#         Np = abpe.Np; x = abpe.x; y = abpe.y; θ = abpe.θ
#         x_dot = 5.0*ones(Np) 
#         y_dot = -2.0*ones(Np)
#         θ_dot = 3.0*ones(Np)

#         δp = δt*[x_dot y_dot]
#         δθ = δt*θ_dot
#     else
#         println("No step method available")
#     end
#     return (δp, δθ)
# end


# # Lorenz attractor model
# function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
#     if size(position(abpe),2) == 2
#         Np = abpe.Np; x = abpe.x; y = abpe.y; θ = abpe.θ
#         σ = 10; ρ = 28; β = 8/3
#         x_dot = σ*y .- σ*x
#         y_dot = -x.*θ .+ ρ*x .- y
#         θ_dot = x.*y - β*θ

#         δp = δt*[x_dot y_dot]
#         δθ = δt*θ_dot
#     else
#         println("No step method available")
#     end
#     return (δp, δθ)
# end

# # Model with deterministic rotational and translational diffusion to imitate random components ?
# function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
# end