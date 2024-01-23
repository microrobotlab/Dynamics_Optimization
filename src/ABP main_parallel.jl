using CalculusWithJulia, ForwardDiff
using Distances
using Plots 
using ProgressBars
using Random
using .Threads


# Define an "ABPsEnsemble" Type
abstract type ABPsEnsemble end

# Define a specific type for 2D ABPsEnsembles (CURRENTLY ASSUMING ALL PARTICLES ARE EQUAL)
struct ABPE2 <: ABPsEnsemble
    Np::Int64           # number of particles®®
    L::Float64          # size of observation space (μm)
	R::Float64          # Radius (μm)                                   --> Vector{Float64}(undef,Np)
	v::Float64 	        # velocity (μm/s)                               --> Vector{Float64}(undef,Np)
	DT::Float64         # translational diffusion coefficient (μm^2/s)  --> Vector{Float64}(undef,Np)
	DR::Float64         # rotational diffusion coefficient (rad^2/s)    --> Vector{Float64}(undef,Np)
	x::Vector{Float64}  # x position (μm)
	y::Vector{Float64}  # y position (μm)
	θ::Vector{Float64}  # orientation (rad)
end

# The step function that describes the particles dynamics is provided in an external file `step.jl` for modularity
include("step.jl")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initialize ABP ensemble (CURRENTLY ONLY 2D)

# For ellipse
function initABPE_ellipse(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3, collision_correction::Bool, N::Integer, M::Integer)
    # Translational diffusion coefficient [m^2/s] and rotational diffusion coefficient [rad^2/s] - R [m]
    # Initial condition will be chosen as per the geometry under study
    DT, DR = diffusion_coeff(1e-6R)

    # ONLY 2D!
    k=0.5
    xyθ1 = (rand(Np,3).-k).*repeat([L L 2π],Np) # 3 dim matrix with x, y and θ 
    r = (xyθ1[:,1]).*(xyθ1[:,1]) + (xyθ1[:,2]).*(xyθ1[:,2])
  
    rₚ = sqrt.(r)   
    α =atan.(xyθ1[:,2], xyθ1[:,1]) 
    # ellipse parameters
    a= L/2
    b= L/4
    #rₑ = (a*b)./(sqrt.(((a*sin.(xyθ1[:,3])).^2) .+ (b*cos.((xyθ1[:,3]))).^2))  # r value for boundary
    rₑ = (a-R)*(b-R)./(sqrt.((((a-R)*sin.(α)).^2) .+ ((b-R)*cos.((α))).^2))  # r value for boundary
    #rₑ = b/sqrt.(1 .-((e*cos.(rθ)).^2))

    # To keep particles inside the ellipse only
    id = (rₚ .< (rₑ))
    xyθ = [xyθ1[id,1] xyθ1[id,2] xyθ1[id,3]]

    Np1= size(xyθ,1)    # number of particles inside the boundary while Np is total number of particles
    #xyθ = (rand(Np,3).-0.0).*repeat([L L 2π],Np)

     # Apply hardsphere correction if required
     if(collision_correction)
        # Call hardsphere correction on initial state (modify xyθ) 
        xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2], L, R; N=N, M=M) #xyθ[:,1:2] gives x and y positions of initial particles
    else
        # Otherwise the matrices used for hardsphere correction will remain `nothing`
        dists, superpose, uptriang = nothing, nothing, nothing
    end

    abpe = ABPE2( Np1, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])
    return abpe, (dists, superpose, uptriang)
end

# For square
function initABPE_square(Np::Int64, L::Float64, R::Float64, v::Float64; T::Float64=300.0, η::Float64=1e-3, collision_correction::Bool, N::Integer, M::Integer)
    # Translational diffusion coefficient [m^2/s] and rotational diffusion coefficient [rad^2/s] - R [m]
    # Initial condition will be chosen as per the geometry under study
    DT, DR = diffusion_coeff(1e-6R)

    # ONLY 2D!
    k=0.5
    xyθ = (rand(Np,3).-k).*repeat([L L 2π],Np) # 3 dim matrix with x, y and θ 
    Np= size(xyθ,1)    # number of particles inside the square

    # Apply hardsphere correction if required
    if(collision_correction)
        # Call hardsphere correction on initial state (modify xyθ) 
        xyθ[:,1:2], dists, superpose, uptriang = hardsphere(xyθ[:,1:2], L, R; N=N, M=M) #xyθ[:,1:2] gives x and y positions of intitial particles
    else
        # Otherwise the matrices used for hardsphere correction will remain `nothing`
        dists, superpose, uptriang = nothing, nothing, nothing
    end

    abpe = ABPE2( Np, L, R, v, 1e12DT, DR, xyθ[:,1], xyθ[:,2], xyθ[:,3])
    return abpe, (dists, superpose, uptriang)
end


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Compute diffusion coefficient

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


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to simulate multiple spherical particles in interaction

function multiparticleE(;Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3, wall_condition::String, collision_correction::Bool, N::Integer, M::Integer, verbose::Bool)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1) # Nt is the number of time steps
    ABPE[1], matrices = initABPE_square(Np, L, R, v; collision_correction=collision_correction, N=N, M=M) # including initial hardsphere correction
    
    simulate!(ABPE, matrices, Nt, δt, wall_condition; collision_correction=collision_correction, N=N, M=M, verbose=verbose)

    return position.(ABPE), orientation.(ABPE)
end

function simulate!(ABPE, matrices, Nt, δt, wall_condition; collision_correction, N, M, verbose)
    iter = 1:Nt
    # Progress meter if verbose 
    if verbose
        iter = ProgressBar(iter)
        # for progress bar description
        set_description(iter, "SIMULATION")
    end

    for nt in iter
        ABPE[nt+1] = update(ABPE[nt],matrices,δt, wall_condition; collision_correction=collision_correction, N=N, M=M) # update information at every step
    end
    return nothing
end

position(abpe::ABPE2) = [ abpe.x abpe.y ]
orientation(abpe::ABPE2) = abpe.θ


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions to update particles for the next step

function update(abpe::ABPE, matrices, δt::Float64, wall_condition::String; collision_correction::Bool, N::Integer, M::Integer) where {ABPE <: ABPsEnsemble}
    # Step
    pθ = ( position(abpe), orientation(abpe) ) .+ step(abpe,δt)

    # Boundary conditions
    if(wall_condition == "open") # no modification in "open" mode (no walls)
    elseif(wall_condition == "periodic")
        periodic_BC_array!(pθ[1], abpe.L, abpe.R)
    else
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    # Particle collisions correction if required
    if(collision_correction)
        hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.L, abpe.R; N=N, M=M)
    end

    # Updated system of particles
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )

    return new_abpe
end


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions for particle collisions correction. The space where particles evolve will 
# be divided in cells and the calculation will be performed in parallel on these cells.

function borders(L) 
    # Give borders for the total area.
    # L is the size of the area and it is centered around 0.
    # (FIXED in the current version but could depend upon wall condition) 
    
    # Add a small margin for safety
    ϵ = 0.5
    x_min = -L/2 - ϵ
    x_max = L/2 + ϵ
    y_min = -L/2 - ϵ
    y_max = L/2 + ϵ
    return x_min, x_max, y_min, y_max
end

# Partition particles (xy gives their positions) in order to perform parallel computation
function indices_per_cell(xy, R, tol, N, M, x_min, x_max, y_min, y_max)
    indices_partition = [Array{Integer}([]) for _=1:N*M]

    # Cell dimensions (without overlapping)
    cell_height = (y_max - y_min) / N
    cell_width = (x_max - x_min) / M

    # Distribute particles into the different cells
    for i in axes(xy, 1)
        for n=1:N
            for m=1:M
                # We overlap regions with border 2R*(1-tol), otherwise collisions at the border of the cells will never be detected 
                if((xy[i,1] > x_min + (m-1)*cell_width - 2R*(1-tol)) && (xy[i,1] < x_min + m*cell_width + 2R*(1-tol)) && (xy[i,2] > y_min + (n-1)*cell_height - 2R*(1-tol)) && (xy[i,2] < y_min + n*cell_height + 2R*(1-tol)))
                    push!(indices_partition[(n-1)*M + m], i)
                end 
            end
        end
    end

    return indices_partition
end

# Distance / superpositions matrices update (see `hardsphere` function)
function update_superpositions(xy::Array{Float64,2}, indices, dists::Matrix{Float64}, superpose::BitMatrix, uptriang::BitArray{2}, R::Float64; tol::Float64)
    dists[indices, indices] .= pairwise(Euclidean(),xy[indices, :],xy[indices, :],dims=1)
    superpose[indices, indices] .= (dists[indices, indices] .< 2R*(1-tol)) .* uptriang[indices, indices]
end

# Hardsphere correction for one cell; `indices` provides the indices of the particles in the cell in which we want to apply collisions correction 
function local_hardsphere_correction!(xy::Array{Float64,2}, indices, dists::Matrix{Float64}, superpose::BitMatrix, R::Float64; tol::Float64)
    # Number of superpositions in considered cell
    superpositions = sum(superpose[indices, indices])

    if(superpositions > 0)
        for np1 in indices
            # If at least one value of the np1 row is true, i.e. there is a superposition with the corresponding particle
            if any(superpose[np1,indices]) 
                # Take the first pair of superposed particles (np1, np2)
                # /!\ `findfirst(superpose_list[cell_index][np1,indices])` gives the np2 position in the matrix
                # restricted to [indices, indices] so we have to take the element at this relative position in indices, hence:
                np2 = indices[findfirst(superpose[np1,indices])]

                # superposition correction
                Δp = (xy[np1,:] - xy[np2,:]) .* ( ( (1+tol)*2R / dists[np1,np2] - 1 ) / 2 )
                xy[np1,:] += Δp
                xy[np2,:] -= Δp
                
                # update distances and superpositions matrices (see `hardsphere` function)
                dists[np2,indices[indices .> np2]] = pairwise(Euclidean(), xy[np2:np2,:], xy[indices[indices .> np2],:], dims=1 )  # distances for the row wise pair operation
                superpose[np2,indices[indices .> np2]] = dists[np2,indices[indices .> np2]] .< 2R*(1-tol)
            end
        end
    end
    return superpositions
end

# Main function for hardsphere correction
function hardsphere!(
    xy::Array{Float64,2},
    dists::Matrix{Float64}, superpose::BitMatrix, uptriang::BitArray{2}, 
    L::Float64, 
    R::Float64;
    tol::Float64=1e-3, 
    N::Integer, M::Integer,
    )

    # Borders of the area where the particles evolve
    x_min, x_max, y_min, y_max = borders(L)

    # We introduced overlapping regions (see `indices_per_cell` function) of size 2R(1-tol)
    # If those regions are bigger than half a cell, they start overlapping themselves
    # Vertically, (y_max - y_min) gives the height of the area where particles evolve and N is 
    # the number of vertical divisions. ((y_max - y_min) / N) / 2 thus gives half the height of one cell.
    # Finally, we want: 2R(1-tol) <= ((y_max - y_min) / N) / 2, or: N >= (y_max - y_min) / (4R * (1-tol))
    # and horizontally: M >= (x_max - x_min) / (4R * (1-tol)) where M is the number of horizontal divisions.
    if(N >= (y_max - y_min) / (4R * (1-tol)))
        throw("N must be smaller, possible overlapping between parallel-computed regions")
    elseif(M >= (x_max - x_min) / (4R * (1-tol)))
        throw("M must be smaller, possible overlapping between parallel-computed regions")
    end

    # Partition the particles w.r.t. the cells
    indices_partition = indices_per_cell(xy, R, tol, N, M, x_min, x_max, y_min, y_max)

    # Due to overlaps, neighboring cells will have particles in common and will therefore perform simultaneous accesses 
    # during collision corrections. To avoid it, cells are separated in 4 groups : considering cells division as a checkerboard 
    # (here indices_per_cell is just a 1 dimension list), they correspond to even/even, even/odd, odd/even, odd/odd indices of the cells. 
    # For proper choice of N and M (see above), cells grouped this way by parity  won't have particles in common because they are spaced 
    # by exactly one cell horizontally and vertically. As mentioned below, collision corrections are applied in a "semi-parallel" way: 
    # they are performed in parallel within the 4 groups, but sequentially between them.
    # Here we compute the corresponding indices (of the cells not the particles); indices_per_cell has a linear indexing:

    even_even_cells_indices = []
    even_odd_cells_indices = []
    odd_even_cells_indices = []
    odd_odd_cells_indices = []
    for i in eachindex(indices_partition)
        # to switch from linear to row / colum indexing and then check parity of vertical / horizontal index
        if (iseven(((i-1) ÷ M) + 1) && iseven(((i-1) % M) + 1)); push!(even_even_cells_indices,i);
        elseif (iseven(((i-1) ÷ M) + 1) && isodd(((i-1) % M) + 1)); push!(even_odd_cells_indices,i);
        elseif (isodd(((i-1) ÷ M) + 1) && iseven(((i-1) % M) + 1)); push!(odd_even_cells_indices,i);
        else push!(odd_odd_cells_indices,i); end
    end

    # Keep track of the number of superpositions for each cell of the partition (the threads will access it separately)
    superposition_partition = zeros(Int, length(indices_partition))
    # initialized as 1 to pass the while loop condition at first iteration
    superpositions = 1
    # set a limit to avoid too much hardsphere correction iterations
    counter = 0
  
    while superpositions > 0
        # Reset superpositions count
        superposition_partition = zero(superposition_partition)
        # THREADING REGION: calculation will be carried out in parallel within the four groups, 
        # but separately between them (sequentially). See explanations above. 
        for parity_cell_indices in [even_even_cells_indices, even_odd_cells_indices, odd_even_cells_indices, odd_odd_cells_indices]
            @threads for cell_index in parity_cell_indices
                # update distances and superposition matrices
                update_superpositions(xy, indices_partition[cell_index], dists, superpose, uptriang, R; tol=tol)
                # local hardsphere correction
                superposition_partition[cell_index] += local_hardsphere_correction!(xy, indices_partition[cell_index], dists, superpose, R; tol=tol)
            end
        end

        # Sum superpositions for each cell   
        superpositions = sum(superposition_partition)  
        counter += 1
        # To avoid spending too much time in hardsphere correction
        if counter >= 100
            println("$superpositions superpositions remaining after 100 cycles")
            break
        end
    end
        
    # QUALITY TEST TO SEE THE ACTUAL AMOUNT OF SUPERPOSITIONS AFTER CORRECTION
    # (COMPUTED ON THE OVERALL SYSTEM, NOT PER CELL) 
    # dists_total = pairwise(Euclidean(),xy,xy,dims=1)
    # superpose_total = (dists_total .< 2R*(1-tol)) .* uptriang
    # superpositions_total = sum(superpose_total)
    # println("SUPERPOSITIONS REMAINING (REAL TOTAL) : $superpositions_total")

    return nothing
end
 
# (Called in initABPE) initialize distance and superposition matrices that will store distances and superpositions between particles
# during simulation (they are continuously updated). Also run `hardsphere!` once.
function hardsphere(
    xy::Array{Float64,2},
    L::Float64,
    R::Float64; 
    tol::Float64=1e-3, 
    N::Integer, M::Integer,
    ) 

    Np = size(xy,1)
    dists = zeros(Np,Np) 
    superpose = falses(Np,Np)
    uptriang = falses(Np,Np)
    for i = 1:Np-1
        uptriang[i,i+1:Np] .= true
    end

    hardsphere!(xy, dists, superpose, uptriang, L, R; tol=tol, N=N, M=M)
    return xy, dists, superpose, uptriang
end

# To simulate periodic boundaries
function periodic_BC_array!(xy::Array{Float64,2},L::Float64, R)   # when a particle crosses an edge it reappears on the opposite side
	# Boundary conditions: horizontal edge
	idx = abs.(xy[:,1]) .> L/2 + R # I create vector idx in which I have 1 where the absolute value of the x coordinate of the particles is outside the observation area
	if any(idx)
		xy[idx,1] .-= sign.(xy[idx,1]).*L   # where I have uni in idx I make the particle reappear on the opposite side of x with respect to 0
	end
	# Boundary conditions: vertical edge
	idy = abs.(xy[:,2]) .> L/2 + R
	if any(idy)
		xy[idy,2] .-= sign.(xy[idy,2]).*L
	end
	return nothing
end


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Functions for the version of the simulator that implements different boundary conditions such as squared, elliptical.

function multiparticleE_wall(;Np::Integer, L::Float64, R::Float64, v::Float64, Nt::Int64=2, δt::Float64=1e-3, wall_condition::String, collision_correction::Bool, N::Integer, M::Integer, verbose::Bool)
    (Nt isa Int64) ? Nt : Nt=convert(Int64,Nt)
    
    ABPE = Vector{ABPE2}(undef,Nt+1)
    if(wall_condition=="elliptical")
        ABPE[1], matrices = initABPE_ellipse( Np, L, R, v; collision_correction=collision_correction, N=N, M=M) # including initial hardsphere correction
    elseif(wall_condition in ["periodic", "squared"])
        ABPE[1], matrices = initABPE_square( Np, L, R, v; collision_correction=collision_correction, N=N, M=M) 
    else
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    simulate_wall!(ABPE, matrices, Nt, δt, wall_condition; collision_correction=collision_correction, N=N, M=M, verbose=verbose)
    return position.(ABPE), orientation.(ABPE)
end

function simulate_wall!(ABPE, matrices, Nt, δt, wall_condition; collision_correction, N, M, verbose)
    iter = 1:Nt
    # Progress meter if verbose 
    if verbose
        iter = ProgressBar(iter)
        set_description(iter, "SIMULATION")
    end

    for nt in iter
        ABPE[nt+1] = update_wall(ABPE[nt],matrices,δt, wall_condition; collision_correction=collision_correction, N=N, M=M)
    end
    return nothing
end

function update_wall(abpe::ABPE, matrices, δt::Float64, wall_condition::String; collision_correction::Bool, N::Integer, M::Integer) where {ABPE <: ABPsEnsemble}
    # Step
    memory_step = step(abpe,δt)
    pθ = ( position(abpe), orientation(abpe) ) .+ memory_step
    
    # Boundary condition
    if(wall_condition == "squared")
        wall_condition!(pθ[1],abpe.L, abpe.R, memory_step[1])    
    elseif(wall_condition == "elliptical")
        elliptical_wall_condition!(pθ[2],pθ[1],abpe.L, abpe.R, memory_step[1])
    else
        throw(ArgumentError("please provide a correct argument for wall condition"))
    end

    # Hardsphere correction if required
    if(collision_correction)
        hardsphere!(pθ[1], matrices[1], matrices[2], matrices[3], abpe.L, abpe.R; N=N, M=M)
    end

    # Updated system of particles
    new_abpe = ABPE2( abpe.Np, abpe.L, abpe.R, abpe.v, abpe.DT, abpe.DR, pθ[1][:,1], pθ[1][:,2], pθ[2] )
    return new_abpe
end

# To simulate square reflective boundary
function wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2})
	# Horizontal edge
	idx = abs.(xy[:,1]) .> (L/2 - R)
	if any(idx)
		xy[idx,1] .-= 2*sign.(xy[idx,1]).*(abs.(xy[idx,1]) .- (L/2 - R)) 
	end
	# Vertical edge
	idy = abs.(xy[:,2]) .> (L/2 - R)
	if any(idy)
        xy[idy,2] .-= 2*sign.(xy[idy,2]).*(abs.(xy[idy,2]) .- (L/2 - R))
	end
	return nothing
end

# To simulate circular reflective boundary (generalized by `elliptical_wall_condition!`)
function circular_wall_condition!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2})
    # the condition is calculated w.r.t to r value of the particle and have no edges here, 
    # this is first method used    
    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
    rr = sqrt.(r)
    rθ= atan.(xy[:,2], xy[:,1])     # angle which the particles make with the origin 
    id = rr.> (L/2 - R)             # checking the condition for r vector id = 1 when true, 0 when false
    Δr = zeros(length(r),1)
    #println("$id")
                
    Δr[id] = rr[id].- (L/2 - R)
    rr.-= 2*(Δr)
    #println(size(rr)) # gives size of an array
    xy[id,1] = rr[id].*(cos.(rθ[id]))
    xy[id,2] = rr[id].*(sin.(rθ[id]))

    return nothing
end

function circular_wall_condition1g!(xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
    # Here the condition is calculated w.r.t to the normal at the intersection point of the radial distance with the wall,
    # this is second method used 
    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
    rr = sqrt.(r)
    rθ= atan.(xy[:,2], xy[:,1])     # angle which the particles make with the origin 
    id = rr.> (L/2 - R)             # checking the condition for r vector id = 1 when true, 0 when false

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

    xpos= [p[1] for p in inside]
    ypos= [p[2] for p in inside]

    xy[id,1] = xpos[id,1]
    xy[id,2] = ypos[id,1]
      
    return nothing
end

function elliptical_wall_condition!(orientation::Array{Float64,1},xy::Array{Float64,2},L::Float64, R, step_mem::Array{Float64,2}) # this condition is for cicular reflective boundary
    # Here the condition is calculated w.r.t to the normal at the intersection point of the radial distance with the wall,
    # this is second method used 
    # Orientation of particle will also change   
    a= L/2
    b= L/4

    a1= a-R
    b1= b-R

    #e = sqrt(1-(b/a)^2)
    r = (xy[:,1]).*(xy[:,1]) + (xy[:,2]).*(xy[:,2])
    rₚ = sqrt.(r)                   # particle r coordinate
    rθ= atan.(xy[:,2], xy[:,1])     # angle which the particles make with the origin 
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

    normal = grad.(x,a1,b1)             # gradient calculated from each particle position onto the ellipse boundary 
    hat_normal = normal./norm.(normal)  # unit normal vector

    correction_x = (rₚ .- rₑ) .*[(cos.(θ)) for θ in rθ]
    correction_y = (rₚ .- rₑ) .*[(sin.(θ)) for θ in rθ]

    c= [correction_x correction_y]
    correction = [[p[1],p[2]] for p in eachrow(c)]

    projection = dot.(correction,hat_normal)

    cᵥ = 2*(projection.+0*R).* hat_normal # vector of vector
    # to access this vector I am breaking it in x and y in the following lines
    cᵥx =  [p[1] for p in cᵥ]
    cᵥy =  [p[2] for p in cᵥ]

    ############################# calculations for orientation vector ###############################
    ########################### 3 conversion of vector into matrix was job of this part ############
    xy[id,1] .-= cᵥx[id]
    xy[id,2] .-= cᵥy[id]

    return nothing
end
