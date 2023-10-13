using .Threads
using CellListMap


mutable struct Atomic{T}; @atomic A::T; end

function hardsphere_correction_cell!(xy::Array{Float64,2}, R::Float64; tol::Float64=1e-3)
    superpositions = true
    counter = 0
    # visited = Atomic(Array{Integer}([]))
    visited = Array{Integer}([])

    while(superpositions && (counter < 100))        
        close_particles = neighborlist(xy', 2R*(1-tol))
        superpositions = !isempty(close_particles)
        # @threads for pair in close_particles
        for pair in close_particles
            i,j = pair[1:2]
            # if((i ∉ visited.A) && (j ∉ visited.A))
            if((i ∉ visited) && (j ∉ visited))
                # pair[3] gives the distance between the two particles
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


function indices_per_cell(xy, R, tol, N, M)
    # function indices_per_cell(xy, R, tol, x_min, x_max, y_min, y_max, N, M)
        # partition xy in order to perform parallel computation
        # the limits of the region to be divided are the maximum values of xy for both dimensions
        # which avoids empty cells (if the region was the entire space)
        indices_partition = [Array{Integer}([]) for _=1:N*M]
        # little margin
        ϵ = 0.1
        x_max = maximum(xy[:,1]) + ϵ
        y_min = minimum(xy[:,2]) - ϵ
        y_max = maximum(xy[:,2]) + ϵ
        x_min = minimum(xy[:,1]) - ϵ
        N_cap = min(N, floor(Int, floor((y_max - y_min) / (4R*(1-tol)))))
        M_cap = min(M, floor(Int, floor((x_max - x_min) / (4R*(1-tol)))))
        cell_height = (y_max - y_min) / N_cap
        cell_width = (x_max - x_min) / M_cap
    
        for i in axes(xy, 1)
            for n=1:N
                for m=1:M
                    # we overlap the regions with border 2R*(1-tol), otherwise collisions at the border of the cells will never be detected 
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






    function hardsphere!(
        xy::Array{Float64,2},
        dists::Matrix{Float64}, superpose::BitMatrix, uptriang::BitArray{2}, 
        R::Float64; 
        tol::Float64=1e-3, 
        N::Integer, M::Integer,
        x_min, x_max, y_min, y_max
        )
    
        # partition the particles w.r.t. the cells
        indices_partition = indices_per_cell(xy, R, tol, N, M, x_min, x_max, y_min, y_max)
        # println(indices_partition)
    
        # keep a trace of the number of superpositions per cell in the partition
        # (the threads will access it separately)
        superposition_partition = zeros(Int, length(indices_partition))
        # initialized as 1 to pass the while loop condition at first iteration
        superpositions = 1
        # set a limit to avoid to much hardsphere correction iterations
        counter = 0
        # @time begin
    
        even_even_cells_indices = [] 
        even_odd_cells_indices = []
        odd_even_cells_indices = [] 
        odd_odd_cells_indices = []
        
        for i in eachindex(indices_partition)
            if (iseven(((i-1) ÷ M) + 1) && iseven(((i-1) % M) + 1)) 
                push!(even_even_cells_indices,i)
            elseif (iseven(((i-1) ÷ M) + 1) && isodd(((i-1) % M) + 1)) 
                push!(even_odd_cells_indices,i)
            elseif (isodd(((i-1) ÷ M) + 1) && iseven(((i-1) % M) + 1)) 
                push!(odd_even_cells_indices,i)
            else push!(odd_odd_cells_indices,i) 
            end
        end
    
        # ee = collect(Iterators.flatten(indices_partition[even_even_cells_indices]))
        # eo = collect(Iterators.flatten(indices_partition[even_odd_cells_indices]))
        # oe = collect(Iterators.flatten(indices_partition[odd_even_cells_indices]))
        # oo = collect(Iterators.flatten(indices_partition[odd_odd_cells_indices]))
        # if(maximum([count(x->x==i, ee) for i in unique(ee)]) > 1) 
        #     throw("Intersections ee") 
    
        # elseif(maximum([count(x->x==i, eo) for i in unique(eo)]) > 1) 
        #     throw("Intersections eo") 
        
        # elseif(maximum([count(x->x==i, oe) for i in unique(oe)]) > 1) 
        #     throw("Intersections oe") 
        
        # elseif(maximum([count(x->x==i, oo) for i in unique(oo)]) > 1) 
        #     throw("Intersections oo") 
        # end
    
        # println(even_even_cells_indices)
        # println(even_odd_cells_indices)
        # println(odd_even_cells_indices)
        # println(odd_odd_cells_indices)
        # println(indices_partition)
    
        # ee = collect(Iterators.flatten(indices_partition[even_even_cells_indices])); println("ee : $ee")
        # eo = collect(Iterators.flatten(indices_partition[even_odd_cells_indices])); println("eo : $eo")
        # oe = collect(Iterators.flatten(indices_partition[odd_even_cells_indices])); println("oe : $oe")
        # oo = collect(Iterators.flatten(indices_partition[odd_odd_cells_indices])); println("oo : $oo")
        # println("intersection : $(intersect(ee, oe, oe, oo))\n")
        
        # scatter(xy[:,1], xy[:,2], aspect_ratio=:equal, lims=(-100/2, 100/2), markersize=350*R/100., marker =:circle,series_annotations = text.(1:length(xy[:,1]), :bottom), legend=false, title = nothing, show=true)
        # scatter(xy[:,1], xy[:,2], aspect_ratio=:equal, lims=(-100/2, 100/2), markersize=350*R/100., marker =:circle,series_annotations = text.(1:length(xy[:,1]), :bottom), legend=false, title = nothing, show=true)
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
            for parity_cell_indices in [even_even_cells_indices, even_odd_cells_indices, odd_even_cells_indices, odd_odd_cells_indices]
                @threads for cell_index in parity_cell_indices
                    update_superpositions(xy, indices_partition[cell_index], dists, superpose, uptriang, R; tol=tol)
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
            # @show(counter)        
        end
        # if(Inf ∈ xy[:,1] || Inf ∈ xy[:,2] || -Inf ∈ xy[:,1] || -Inf ∈ xy[:,2])
        #     println("INF!")
        # end
            
        # # REAL QUALITY TEST
        # dists_total = pairwise(Euclidean(),xy,xy,dims=1)
        # superpose_total = (dists_total .< 2R*(1-tol)) .* uptriang
        # superpositions_total = sum(superpose_total)
        # println("SUPERPOSITIONS REMAINING (REAL TOTAL) : $superpositions_total")
    
        return nothing
    end
    