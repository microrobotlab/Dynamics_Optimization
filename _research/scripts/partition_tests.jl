# TESTS FOR `indices_per_cell` FUNCTION FOR PARALLEL VERSION OF COLLISIONS CORRECTION
# SEE `indices_per_cell` IN "src/ABP main_parallel.jl"


function indices_per_cell1(xy, N, M, H, L)
    part = @. ceil(Int, (xy / [L H]) * [N M])
    indices_partition = fill(Array{Int}([]), (N,M))
    for i=1:N
        for j=1:M
            indices_partition[i,j] = findall(r->r==[i,j], eachrow(part))
        end
    end
    return indices_partition
end


function indices_per_cell2(xy, N, M)
    # partition xy in order to perform parallel computation
    # the limits of the region to be divided are the maximum values of xy for both dimensions
    # which avoids empty cells (if the region was the entire space)
    indices_partition = [Array{Integer}([]) for _=1:N*M]
    # indices_partition = [[Array{Integer}([]) for _=1:M] for _=1:N]
    # to have a little margin
    ϵ = 0.1
    x_min = minimum(xy[:,1]) - ϵ
    x_max = maximum(xy[:,1]) + ϵ
    y_min = minimum(xy[:,2]) - ϵ
    y_max = maximum(xy[:,2]) + ϵ
    cell_height = (y_max - y_min) / N
    cell_width = (x_max - x_min) / M
    lx = []
    ly = []
    for i=1:size(xy,1)
        for n=1:N
            for m=1:M
                if((xy[i,1] > x_min + (m-1)*cell_width) && (xy[i,1] < x_min + m*cell_width) && (xy[i,2] > y_min + (n-1)*cell_height) && (xy[i,2] < y_min + n*cell_height))
                    push!(indices_partition[(n-1)*M + m], i)
                end 
                push!(lx, x_min + (m-1)*cell_width)
                push!(lx, x_min + m*cell_width)
                push!(ly, y_min + n*cell_height)
                push!(ly, y_min + (n-1)*cell_height)
            end
        end
    end
    return indices_partition, lx, ly
end


function indices_per_cell(xy, R, tol, N, M)
    # partition xy in order to perform parallel computation
    # the limits of the region to be divided are the maximum values of xy for both dimensions
    # which avoids empty cells (if the region was the entire space)
    indices_partition = [Array{Integer}([]) for _=1:N*M]
    # little margin
    # ϵ = 0.1
    # x_max = maximum(xy[:,1]) + ϵ
    # y_min = minimum(xy[:,2]) - ϵ
    # y_max = maximum(xy[:,2]) + ϵ
    # x_min = minimum(xy[:,1]) - ϵ
    # N_cap = min(N, floor(Int, floor((y_max - y_min) / (4R*(1-tol)))))
    # M_cap = min(M, floor(Int, floor((x_max - x_min) / (4R*(1-tol)))))
    # cell_height = (y_max - y_min) / N_cap
    # cell_width = (x_max - x_min) / M_cap

    # size of the cells given by the considered region and their number
    ϵ = 3
    x_min = -50 - ϵ
    x_max = 50 + ϵ
    y_min = -50 - ϵ
    y_max = 50 + ϵ
    cell_height = (y_max - y_min) / N
    cell_width = (x_max - x_min) / M

    for i=1:size(xy,1)
        for n=1:N
            for m=1:M
                # we overlap cells with an extension of size 2R*(1-tol), otherwise collisions at the border of the cells will never be detected 
                # if((xy[i,1] > x_min + (m-1)*cell_width) && (xy[i,1] < x_min + m*cell_width) && (xy[i,2] > y_min + (n-1)*cell_height) && (xy[i,2] < y_min + n*cell_height))
                if((xy[i,1] > x_min + (m-1)*cell_width - 2R*(1-tol)) && (xy[i,1] < x_min + m*cell_width + 2R*(1-tol)) && (xy[i,2] > y_min + (n-1)*cell_height - 2R*(1-tol)) && (xy[i,2] < y_min + n*cell_height + 2R*(1-tol)))
                    push!(indices_partition[(n-1)*M + m], i)
                end 
            end
        end
    end
    lx = [x_min + m*cell_width for m=0:M]
    ly = [y_min + n*cell_height for n=0:N]
    return indices_partition#, lx, ly
end