# Simplistic model with straight particle trajectories
function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}

    if size(position(abpe),2) == 2
        # δp = abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)]
        Np = abpe.Np; x = abpe.x; y = abpe.y; θ = abpe.θ
        v0 = 5.0
        x_dot = v0*ones(Np)
        y_dot = v0*ones(Np)
        θ_dot = zeros(Np)

        δp = δt*[x_dot y_dot]
        δθ = δt*θ_dot
    else
        println("No step method available")
    end
    #if nt == 1 
        #println("lo step vero di questo giro è: $δp") 
    #end
    return (δp, δθ)
end

# # Lorenz model, to check results correspondance
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

# function deterministic_rotational_translational_diffusion_step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    
#     if size(position(abpe),2) == 2
#         δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)]
#         δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)
#     else
#         println("No step method available")
#     end
#     #if nt == 1 
#         #println("lo step vero di questo giro è: $δp") 
#     #end
#     return (δp, δθ)
# endplot(X[1,:], label="x")