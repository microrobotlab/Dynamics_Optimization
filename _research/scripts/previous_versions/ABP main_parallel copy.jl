using CalculusWithJulia, ForwardDiff
using Random

Random.seed!(3)

# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
mutable struct ABPE2 <: ABPsEnsemble
    Np::Int64                      # number of particles®®
    L::Float64                      # size of observation space (μm)
	R::Float64  # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Float64 	# velocity (μm/s)                               --> Vector{Float64}(undef,Np)
	DT::Float64 # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64 # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}    # x position (μm)
	y::Vector{Float64}    # y position (μm)
	θ::Vector{Float64}    # orientation (rad)
end


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initialize ABP ensemble (CURRENTLY ONLY 2D)

#------------------------------------------------------------For ellipse ---------------------------------------------------------------------------------------------------------
function initABPE_ellipse(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    DT, DR = diffusion_coeff(1e-6R)

    # ONLY 2D!
    k=0.5
    xyθ1 = (rand(Np,3).-k).*repeat([L L 2π],Np) # 3 dim matrix with x, y and θ 
    r = (xyθ1[:,1]).*(xyθ1[:,1]) + (xyθ1[:,2]).*(xyθ1[:,2])
  
    rₚ = sqrt.(r)   
    α =atan.(xyθ1[:,2], xyθ1[:,1]) 
    a= L/2
        b= L/4
    #rₑ = (a*b)./(sqrt.(((a*sin.(xyθ1[:,3])).^2) .+ (b*cos.((xyθ1[:,3]))).^2))  # r value for boundary
    rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))  # r value for boundary
    #rₑ = b/sqrt.(1 .-((e*cos.(rθ)).^2))
    id = (rₚ .< (rₑ))
    xyθ = [xyθ1[id,1] xyθ1[id,2] xyθ1[id,3]]

    Np1= size(xyθ,1)    # number of particles inside the boundary while Np is total number of particles
    #xyθ = (rand(Np,3).-0.0).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2],R) #xyθ[:,1:2] gives x and y positions of intitial particles
    abpe = ABPE2( Np1, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

    # println("INIT ELLIPSE")
    return abpe, (dists, superpose, uptriang)
end

#------------------------------------------------------------For square ---------------------------------------------------------------------------------------------------------
function initABPE_square(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3)
    # translational diffusion coefficient [m^2/s] & rotational diffusion coefficient [rad^2/s] - R [m]
    # Intial condition will be choosen as per the geometry under study
    DT, DR = diffusion_coeff(1e-6R)

    # ONLY 2D!
    k=0.5
    xyθ = (rand(Np,3).-k).*repeat([L L 2π],Np) # 3 dim matrix with x, y and θ 
   

    Np= size(xyθ,1)    # number of particles inside the sqaure
    #xyθ = (rand(Np,3).-0.0).*repeat([L L 2π],Np)
    xyθ[:,1:2], dists_superpose_list, uptriang = hardsphere(xyθ[:,1:2],R) #xyθ[:,1:2] gives x and y positions of intitial particles
    abpe = ABPE2( Np, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])

    # println("INIT SQUARE")
    return abpe, (dists_superpose_list, uptriang)
end


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
##Calculate diffusion coefficient

function diffusion_coeff(R::Float64, T::Float64=300.0, η::Float64=1e-3)
    # Boltzmann constant [J/K]
    kB = 1.38e-23
    # friction coefficient [Ns/m]
    γ = 6*pi*R*η
    # translational diffusion coefficient [m^2/s]
    DT = kB*T/γ
    # rotational diffusion coefficient [rad^2/s]
    DR = 6*DT/(8*R^2)
    return DT, DR
end;


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to simulate multiple spherical particles

function multiparticleE(;Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1) # Nt is number of time steps
    ABPE[1], matrices = initABPE_square( Np, L, R, v ) # including initial hardsphere correction
    
    simulate!(ABPE, matrices, Nt, δt)

    return position.(ABPE), orientation.(ABPE)
end


function simulate!(ABPE, matrices, Nt, δt)
    # PΘ = [ (position(abpe), orientation(abpe)) ]
    # pθ = PΘ[1]
    
    for nt in 1:Nt
        # if(nt%1000 == 0) println("$nt / $Nt") end
        ABPE[nt+1] = update(ABPE[nt],matrices,δt) #updating information at every step
        # println("Step $nt")
    end
    return nothing
end

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to update particles for the next step

function update(abpe::ABPE, matrices::Tuple{Vector{Tuple{Matrix{Float64}, BitMatrix}}, BitMatrix}, δt::Float64) where {ABPE <: ABPsEnsemble}
    pθ = ( position(abpe), orientation(abpe) ) .+ step(abpe,δt)

    periodic_BC_array!(pθ[1],abpe.L, abpe.R)
    #circular_wall_condition!(pθ[1],L::Float64, R, step_mem::Array{Float64,2})
    hardsphere!(pθ[1], matrices[1], matrices[2], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end


function step(abpe::ABPE, δt::Float64) where {ABPE <: ABPsEnsemble}
    
    if size(position(abpe),2) == 2
        δp = sqrt.(2*δt*abpe.DT)*randn(abpe.Np,2) .+ abpe.v*δt*[cos.(abpe.θ) sin.(abpe.θ)]
        δθ = sqrt(2*abpe.DR*δt)*randn(abpe.Np)
    else
        println("No step method available")
    end
    #if nt == 1 
        #println("lo step vero di questo giro è: $δp") 
    #end
    return (δp, δθ)
end


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for the hard sphere corrections

function indices_per_cell(xy, N, M)
    part = @. ceil(Int, (xy / [L H]) * [N M])
    indices_partition = fill(Array{Int}([]), (N,M))
    for i=1:N
        for j=1:M
            indices_partition[i,j] = findall(r->r==[i,j], eachrow(part))
        end
    end
    return indices_partition
end

function local_hardsphere_correction!(xy::Array{Float64,2}, dists_superpose_list::Vector{Tuple{Matrix{Float64}, BitMatrix}}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
    dists_superpose_list[1][1] .= pairwise(Euclidean(),xy,xy,dims=1)
    dists_superpose_list[1][2] .= (dists_superpose_list[1][1] .< 2R*(1-tol)) .* uptriang
    superpositions = sum(dists_superpose_list[1][2])
    if(superpositions > 0)
        Np = size(dists_superpose_list[1][2],1) # number of particles in the local cell 
        for np1 in 1:Np
            if any(dists_superpose_list[1][2][np1,:]) #if at least one value of the np1 row is true 
                np2 = findfirst(dists_superpose_list[1][2][np1,:])
                Δp = (xy[np1,:] - xy[np2,:]) .* ( ( (1+tol)*2R / dists_superpose_list[1][1][np1,np2] - 1 ) / 2 )
                xy[np1,:] += Δp
                xy[np2,:] -= Δp
                dists_superpose_list[1][1] = pairwise(Euclidean(), xy[np2:np2,:], xy[np2+1:Np,:], dims=1 )  # distances for the row wise pair operation
                dists_superpose_list[1][2] = (dists_superpose_list[1][1][np2,np2+1:Np] .< 2R*(1-tol))   
                # dists[np1+1:np2,np2] = pairwise(Euclidean(), xy[np2:np2,:], xy[np1+1:np2,:], dims=1 )  # distances for the row wise pair operation
                # superpose[np1+1:np2,np2] = (dists[np1+1:np2,np2] .< 2R*(1-tol))   
            end
        end
    end

    # Np = size(superpose,1) ##gives me the number of particles locally (in the considered cell)
    # for i in 1:Np
    #     if any(superpose[i,:]) #if at least one value of the i row is true 
    #         j = findfirst(superpose[i,:])
    #         Δp = (xy[indices[i],:] - xy[indices[j],:]) .* ( ( (1+tol)*2R / dists[i,j] - 1 ) / 2 )
    #         # the lock prevents concurrent access (read/write) to xy from different threads
    #         xy[indices_with_shared[i],:] += Δp
    #         xy[indices_with_shared[j],:] -= Δp

    #         # update horizontal/vertical part of distances with particle corresponding to index j in the matrix
    #         dists[j,j+1:Np] = pairwise(Euclidean(), xy[indices[j]:indices[j],:], xy[indices[j]+1:Np,:], dims=1 )  # distances for the row wise pair operation
    #         dists[i+1:j,j] = pairwise(Euclidean(), xy[j:j,:], xy[i+1:j,:], dims=1 )  # distances for the row wise pair operation
    #         # update superpose similarly
    #         superpose[j,j+1:Np] = (dists[j,j+1:Np] .< 2R*(1-tol)) 
    #         superpose[i+1:j,j] = (dists[i+1:j,j] .< 2R*(1-tol)) 
    #     end
    # end

    return superpositions
end

function hardsphere!(xy::Array{Float64,2}, dists_superpose_list::Vector{Tuple{Matrix{Float64}, BitMatrix}}, uptriang::BitArray{2}, R::Float64; tol::Float64=1e-3)
    superpositions = 1
    counter = 0
    # partition the particles w.r.t. the cells
    # indices_partition = indices_per_cell(xy, N, M)
    # @time begin
    while superpositions > 0
        # @show(superpositions)
        superpositions = local_hardsphere_correction!(xy, dists_superpose_list, uptriang, R; tol=tol)
        counter += 1
        (counter >= 100) && println("$superpositions superpositions remaining after 100 cycles"); break
        # @show(counter)
    end
    return nothing
end

function hardsphere(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3) # called in initABPE
    Np = size(xy,1)
    dists_superpose_list = [(zeros(Np,Np), falses(Np,Np)) for _=1:N*M]
    uptriang = falses(Np,Np)
    for i = 1:Np-1
        uptriang[i,i+1:Np] .= true
    end
    hardsphere!(xy, dists_superpose_list, uptriang, R; tol=tol)
    return xy, dists_superpose_list, uptriang
end
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function periodic_BC_array!(xy::Array{Float64,2},L::Float64, R)   #when a particle crosses an edge it reappears on the opposite side
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> L/2 + R #I create vector idx in which I have 1 where the absolute value of the x coordinate of the particles is outside the observation area
	if any(idx)
		xy[idx,1] .-= sign.(xy[idx,1]).*L   #where I have uni in idx I make the particle reappear on the opposite side of x with respect to 0
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> L/2 + R
	if any(idy)
		xy[idy,2] .-= sign.(xy[idy,2]).*L
	end
    # println("I am in periodic")
	return nothing
end


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Functions for updating reflective boundary AND WALL UPDATE

function multiparticleE_wall(;Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3, wall_condition::String)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1)
    if(wall_condition=="elliptical")
        ABPE[1], matrices = initABPE_ellipse( Np, L, R, v ) # including initial hardsphere correction
    elseif(wall_condition in ["periodic", "squared"])
        ABPE[1], matrices = initABPE_square( Np, L, R, v ) 
    else
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    simulate_wall!(ABPE, matrices, Nt, δt, wall_condition)
    # println("I am in multiwall update")
    return position.(ABPE), orientation.(ABPE)
end

function simulate_wall!(ABPE, matrices, Nt, δt, wall_condition)
    # PΘ = [ (position(abpe), orientation(abpe)) ]
    # pθ = PΘ[1]
    
    for nt in 1:Nt
        ABPE[nt+1] = update_wall(ABPE[nt],matrices,δt, wall_condition)
        #println("Step $nt")
      
    end
    return nothing
end

function update_wall(abpe::ABPE, matrices::Tuple{Matrix{Float64}, BitMatrix, BitMatrix}, δt::Float64, wall_condition::String) where {ABPE <: ABPsEnsemble}
    memory_step = step(abpe,δt)
  
    pθ = ( position(abpe), orientation(abpe) ) .+ memory_step
    
    # wall condition choice
    if(wall_condition == "squared")
        wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])    
    elseif(wall_condition == "elliptical")
        elliptical_wall_condition!(pθ[2],pθ[1],abpe.L, abpe.R, memory_step[1])
    else
        throw(ArgumentError("please provide a correct argument for wall condition"))
    #elliptical_wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])
    end
    hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.R)
    # @btime hardsphere!($p[:,1:2], $matrices[1], $matrices[2], $matrices[3], $params.R)
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end

function wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for square reflective boundary
  
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> (L/2 - R)
	if any(idx)
		xy[idx,1] .-= 2*sign.(xy[idx,1]).*(abs.(xy[idx,1]) .- (L/2 - R)) 
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> (L/2 - R)
	if any(idy)
        xy[idy,2] .-= 2*sign.(xy[idy,2]).*(abs.(xy[idy,2]) .- (L/2 - R))
	end
    # println("I am in square wall")
	return nothing
end

function circular_wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
  # here the condition is calculated w.r.t to r value of the particle and have no edges here
  # this is first method used 
     
    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])

    rr = sqrt.(r)

    rθ= atan.(xy[:,2], xy[:,1])   # angle which the particle make with the origin 
     
    id = rr.> (L/2 - R)        # checking the condition for r vector id = 1 when true, 0 when false
    
    Δr = zeros(length(r),1)
    #println("$id")
             
        Δr[id] = rr[id].- (L/2 - R)
        rr.-= 2*(Δr)
        #println(size(rr)) # gives size of an array
        xy[id,1] = rr[id].*(cos.(rθ[id]))
        xy[id,2] = rr[id].*(sin.(rθ[id]))


	
	#println("I am in circular wall")
	return nothing
end

function circular_wall_condition1g!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
    # here the condition is calculated w.r.t to the normal at the intersection point of the radial distance with the wall
    # this is second method used 
       
      r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
      rr = sqrt.(r)
  
      rθ= atan.(xy[:,2], xy[:,1])   # angle which the particle make with the origin 
       
      id = rr.> (L/2 - R)        # checking the condition for r vector id = 1 when true, 0 when false
      
      correction = zeros(length(r),1)
      projection = zeros(length(r),1)
      inside = zeros(length(r),1)
      hat_normal = zeros(length(r),1)
      x = [[p[1],p[2]] for p in eachrow(xy)]
      
      function grad(x::Array{Float64},R::Float64)
       f(x) = x[1]^2 + x[2]^2 - (L*L)/4
       df=ForwardDiff.gradient(f, [x[1],x[2]])
       return df
    end
    normal = grad.(x,R)
     
    hat_normal = normal./norm.(normal)  # unit normal vector
   
    correction = (rr.-(L/2 - 2*R)).*[[cos.(θ),sin.(θ)] for θ in rθ]
    
    projection = dot.(correction,hat_normal)
                                  
    
    inside = (0.5*L.-1*projection).*[[cos.(θ),sin.(θ)] for θ in rθ]         # the vector which is inside the boundry now
      #println("$id")

      xpos= [p[1] for p in inside]

      ypos= [p[2] for p in inside]
          
          xy[id,1] = xpos[id,1]
          xy[id,2] = ypos[id,1]
  
        # println("I am in circular wall 1g")
    
      return nothing
  end

function elliptical_wall_condition!(orientation::Array{Float64,1},xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
    # here the condition is calculated w.r.t to the normal at the intersection point of the radial distance with the wall
    # this is second method used 
    # orientation of particle will also change   
        a= L/2
        b= L/4

        a1= a-R
        b1= b-R
         
        #e = sqrt(1-(b/a)^2)
      r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
  
      rₚ = sqrt.(r)                # particle r co ordinate
  
      rθ= atan.(xy[:,2], xy[:,1])   # angle which the particle make with the origin 

      rₒ = ((a*b)./(sqrt.(((a*sin.(rθ)).^2) .+ (b*cos.(rθ)).^2)))

      rₑ = rₒ .- R
   
      #rₑ = b/sqrt.(1 .-((e*cos.(rθ)).^2))
      id = (rₚ .> (rₑ))
    
      x = [[p[1],p[2]] for p in eachrow(xy)]
      
      function grad(x::Array{Float64},a::Float64,b::Float64)
        
       f(x) = (x[1]^2)*b^2 + (x[2]^2)*a^2 - (a*b)^2
       df=ForwardDiff.gradient(f, [x[1],x[2]])
       return df
       end
    normal = grad.(x,a1,b1)     # gradient calculated from each particle position onto the ellipse boundary 
     
    hat_normal = normal./norm.(normal)  # unit normal vector

    
    correction_x = (rₚ .- rₑ) .*[(cos.(θ)) for θ in rθ]

    correction_y = (rₚ .- rₑ) .*[(sin.(θ)) for θ in rθ]
    
    c= [correction_x correction_y]

    correction = [[p[1],p[2]] for p in eachrow(c)]
  
    projection = dot.(correction,hat_normal)


    cᵥ = 2*(projection.+0*R).* hat_normal #vector of vector
    # to access this vector I am breaking it in x and y in the following lines

    cᵥx =  [p[1] for p in cᵥ]

    cᵥy =  [p[2] for p in cᵥ]
############################# calculations for orientation vector###############################
###########################3 conversion of vector into matrix was job of this part###########################
        xy[id,1] .-= cᵥx[id]
        xy[id,2] .-= cᵥy[id]

    # println("I am in elliptical wall")
      return nothing
  end
