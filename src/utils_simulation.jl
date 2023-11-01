# Utilility functions for simulation

"""
# Give maximum number of divisions of the space for parallel computing
# with respect to the size of the space `L`, radius of the particles `R`
# and tolerance of the simulator for collisions `tol` (see main_parallel.jl)
"""
function max_space_division(L, R, tol) 
    # (What we explain for N in vertical dimension is the same for M in horizontal one)
    # In parallel computing we introduce a necessary overlapping of the cells of size 2R(1-tol)
    # so to avoid overlapping of the overlapping we must have 2R(1-tol) < (L/N)/2
    # i.e. N < L/4R(1-tol) -> we take (N is an integer) the nearest lower integer floor(L/4R(1-tol))
    # if L/4R(1-tol) is already round then floor(L/4R(1-tol)) == L/4R(1-tol) then N == L/4R(1-tol) (not <)
    # so to keep strict inferiority, we take in this case floor(L/4R(1-tol)) - 1
    tmp = L/4R*(1-tol)
    println(tmp)
    bound = floor(tmp)
    if bound==tmp 
        return Int(bound) - 1
    else
        return Int(bound)
    end
end