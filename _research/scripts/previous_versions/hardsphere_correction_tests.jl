# CODE FOR PARALLEL COMPUTATION OF COLLISION (SUPERPOSITIONS) CORRECTIONS IN DYNAMIC 
# PARTICLES SIMULATION INTEGRATED TO 'ABP main_hardsphere_test.jl'


using .Threads
using CellListMap


mutable struct Atomic{T}; @atomic A::T; end

# TEST : here we use CellListMap package for fast computation of particles vicinity
# THIS FUNCTION WAS FINALLY NOT USED
function hardsphere_correction_cell!(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3)
    # Flag to indicate if superpositions might remain among particles
    superpositions = true
    # Superposition corrections stop after to much iterations
    counter = 0
    # visited = Atomic(Array{Integer}([]))
    # Keep track of visited particles, see below
    visited = Array{Integer}([])

    while(superpositions && (counter < 100))    
        # Compute list of pairs of particles which are closer to each other than cutoff 
        # (= 2R*(1-tol) -> approximately means that particles are superposed)
        close_particles = neighborlist(xy', 2R*(1-tol))
        # if superpositions remain
        superpositions = !isempty(close_particles)
        # @threads for pair in close_particles
        for pair in close_particles
            # elements of the list have format (index1, index2, d(particle_index1, particle_index2) where d is a distance)
            i,j = pair[1:2]
            # if((i ∉ visited.A) && (j ∉ visited.A))
            # eliminate visited particles
            if((i ∉ visited) && (j ∉ visited))
                # pair[3] gives the distance between the two particles
                # correct superposition
                Δp = (xy[i,:] - xy[j,:]) .* (((1+tol)*2R / pair[3] - 1) / 2) 
                xy[i,:] += Δp
                xy[j,:] -= Δp
                # @atomic append!(visited.A, [i, j])
                append!(visited, i, j)
                # println("visited : $visited")
            end 
        end

        # empty!(visited.A)
        empty!(visited)
        counter += 1
        # println("counter : $counter")
        # println("remaining superpositions : $(length(close_particles))")
    end
    # return xy_corrected
end


# see `indices_per_cell` in 'ABP main_parallel.jl'  
function indices_per_cell(xy, R, tol, N, M)
# function indices_per_cell(xy, R, tol, x_min, x_max, y_min, y_max, N, M)
    # partition xy in order to perform parallel computation
    # the limits of the region to be divided are the maximum values of xy for both dimensions
    # which avoids empty cells (if the region was the entire space)
    indices_partition = [Array{Integer}([]) for _=1:N*M]
    # add little margin
    ϵ = 0.1
    x_max = maximum(xy[:,1]) + ϵ
    y_min = minimum(xy[:,2]) - ϵ
    y_max = maximum(xy[:,2]) + ϵ
    x_min = minimum(xy[:,1]) - ϵ

    # We introduced overlapping regions of size 2R(1-tol) (see below)
    # If those regions are bigger than half a cell, they start overlapping themselves
    # Vertically, (y_max - y_min) gives the height of the area where particles evolve and N is 
    # the number of vertical divisions. ((y_max - y_min) / N) / 2 thus gives half the height of one cell.
    # Finally, we want: 2R(1-tol) <= ((y_max - y_min) / N) / 2, or: N >= (y_max - y_min) / (4R * (1-tol))
    # and horizontally: M >= (x_max - x_min) / (4R * (1-tol)) where M is the number of horizontal divisions.
    # Here we select the floor of these values as upper bounds for the number of horizontal / vertical divisions. 
    N_cap = min(N, floor(Int, floor((y_max - y_min) / (4R*(1-tol)))))
    M_cap = min(M, floor(Int, floor((x_max - x_min) / (4R*(1-tol)))))
    # corresponding cell dimensions
    cell_height = (y_max - y_min) / N_cap
    cell_width = (x_max - x_min) / M_cap

    for i in axes(xy, 1)
        for n=1:N
            for m=1:M
                # We overlap cells with an extension of size 2R*(1-tol), otherwise collisions at the border of the cells will never be detected 
                # if((xy[i,1] > x_min + (m-1)*cell_width) && (xy[i,1] < x_min + m*cell_width) && (xy[i,2] > y_min + (n-1)*cell_height) && (xy[i,2] < y_min + n*cell_height))
                if((xy[i,1] > x_min + (m-1)*cell_width - 2R*(1-tol)) && (xy[i,1] < x_min + m*cell_width + 2R*(1-tol)) && (xy[i,2] > y_min + (n-1)*cell_height - 2R*(1-tol)) && (xy[i,2] < y_min + n*cell_height + 2R*(1-tol)))
                    push!(indices_partition[(n-1)*M + m], i)
                end 
            end
        end
    end

    # lx = [x_min + m*cell_width for m=0:M]
    # ly = [y_min + n*cell_height for n=0:N]

    return indices_partition#, lx, ly
end



    # see `hardsphere!` in 'ABP main_parallel.jl'
    function hardsphere!(
        xy::Array{Float64,2},
        dists::Matrix{Float64}, superpose::BitMatrix, uptriang::BitArray{2}, 
        R::Float64; 
        tol::Float64=1e-3, 
        N::Integer, M::Integer,
        x_min, x_max, y_min, y_max
        )
    
        # Partition particles w.r.t. the cells
        indices_partition =  indices_per_cell(xy, R, tol, N, M)
        # println(indices_partition)
    
        # Keep a trace of the number of superpositions per cell in the partition
        # (the threads will access it separately)
        superposition_partition = zeros(Int, length(indices_partition))
        # initialized as 1 to pass the while loop condition at first iteration
        superpositions = 1
        # set a limit to avoid to much hardsphere correction iterations
        counter = 0
        # @time begin


        # Due to overlaps, neighboring cells will have particles in common and will therefore perform simultaneous accesses 
        # during collision corrections. To avoid it, cells are separated in 4 groups : considering cells division as a checkerboard 
        # (here indices_per_cell is just a 1 dimension list), they correspond to even/even, even/odd, odd/even, odd/odd indices of the cells. 
        # For proper choice of N and M (see `indices_per_cell`), cells grouped this way by parity  won't have particles in common because they are spaced 
        # by exactly one cell horizontally and vertically. As mentioned below, collision corrections are applied in a "semi-parallel" way: 
        # they are performed in parallel within the 4 groups, but sequentially between them.
        # Here we compute the corresponding indices (of the cells not the particles); indices_per_cell has a linear indexing:
        even_even_cells_indices = []
        even_odd_cells_indices = []
        odd_even_cells_indices = [] 
        odd_odd_cells_indices = []
        
        for i in eachindex(indices_partition)
            # to switch from linear to row / colum indexing and then check parity of vertical / horizontal index
            if (iseven(((i-1) ÷ M) + 1) && iseven(((i-1) % M) + 1)) 
                push!(even_even_cells_indices,i)
            elseif (iseven(((i-1) ÷ M) + 1) && isodd(((i-1) % M) + 1)) 
                push!(even_odd_cells_indices,i)
            elseif (isodd(((i-1) ÷ M) + 1) && iseven(((i-1) % M) + 1)) 
                push!(odd_even_cells_indices,i)
            else push!(odd_odd_cells_indices,i) 
            end
        end
    
        # # Particles indices in the different partity groups 
        # ee = collect(Iterators.flatten(indices_partition[even_even_cells_indices]))
        # eo = collect(Iterators.flatten(indices_partition[even_odd_cells_indices]))
        # oe = collect(Iterators.flatten(indices_partition[odd_even_cells_indices]))
        # oo = collect(Iterators.flatten(indices_partition[odd_odd_cells_indices]))

        # println(even_even_cells_indices)
        # println(even_odd_cells_indices)
        # println(odd_even_cells_indices)
        # println(odd_odd_cells_indices)
        # println(indices_partition)

        # # verify that there are no particles in to regions with different parities at the same time 
        # println("intersection ee eo : $(intersect(ee, eo))\n")
        # println("intersection ee oe : $(intersect(ee, oe))\n")
        # println("intersection ee oo : $(intersect(ee, oo))\n")
        # println("intersection eo oe : $(intersect(eo, oe))\n")
        # println("intersection eo oo : $(intersect(eo, oo))\n")
        # println("intersection oe oo : $(intersect(oe, oo))\n")
        
        # to visualize cells division
        # scatter(xy[:,1], xy[:,2], aspect_ratio=:equal, lims=(-100/2, 100/2), markersize=350*R/100., marker =:circle, series_annotations = text.(1:length(xy[:,1]), :bottom), legend=false, title = nothing, show=true)
        # scatter(xy[:,1], xy[:,2], aspect_ratio=:equal, lims=(-100/2, 100/2), markersize=350*R/100., legend=false, title = nothing, show=true)
        # hline!(ly, color=:red)
        # hline!(ly .+ 2R*(1-tol), color=:orange)
        # hline!(ly .- 2R*(1-tol), color=:orange)
        # vline!(lx, color=:red)
        # vline!(lx .- 2R*(1-tol), color=:orange)
        # vline!(lx .+ 2R*(1-tol), color=:orange)
        
        # display(current())
        # sleep(0.3)
    
        while superpositions > 0
            # @show(superpositions)
            # reset superpositions count
            superposition_partition = zero(superposition_partition)
            # for cell_index in eachindex(indices_partition) 
            #     # ~ THREADING for the update of the matrices
            #     update_superpositions(xy, cell_index, indices_partition[cell_index], dists_list, superpose_list, uptriang, R; tol=tol)
            # end

            # THREADING REGION: the parallelization is done within the four regions
            for parity_cell_indices in [even_even_cells_indices, even_odd_cells_indices, odd_even_cells_indices, odd_odd_cells_indices]
                @threads for cell_index in parity_cell_indices
                    # update distances and superposition matrices
                    update_superpositions(xy, indices_partition[cell_index], dists, superpose, uptriang, R; tol=tol)
                    # local hardsphere correction
                    superposition_partition[cell_index] += local_hardsphere_correction!(xy, indices_partition[cell_index], dists, superpose, R; tol=tol)
                end
            end
            # sum superpositions for each cell   
            superpositions = sum(superposition_partition)  
            counter += 1
            if counter >= 100
                println("$superpositions superpositions remaining after 100 cycles")
                break
            end
        end
        # if(Inf ∈ xy[:,1] || Inf ∈ xy[:,2] || -Inf ∈ xy[:,1] || -Inf ∈ xy[:,2])
        #     println("INF!")
        # end
            
        # # QUALITY TEST TO SEE THE ACTUAL AMOUNT OF SUPERPOSITIONS AFTER CORRECTION
        # # (COMPUTED ON THE OVERALL SYSTEM, NOT PER CELL) 
        # dists_total = pairwise(Euclidean(),xy,xy,dims=1)
        # superpose_total = (dists_total .< 2R*(1-tol)) .* uptriang
        # superpositions_total = sum(superpose_total)
        # println("SUPERPOSITIONS REMAINING (REAL TOTAL) : $superpositions_total")
    
        return nothing
    end
    