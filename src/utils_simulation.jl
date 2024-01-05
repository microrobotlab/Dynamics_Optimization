# Utility functions for simulation

"""
    max_space_division(L, R, tol)

Give maximum number of divisions of the space for parallel computing with respect to the size of the space `L` (height or width), particle radius `R` and collision tolerance `tol` (see main_parallel.jl)
"""
function max_space_division(L, R, tol) 
    # Corresponds to the strictly lower integer closest to the limit size defined in `hardsphere!` (see `ABP main_parallel.jl`)
    # where we found the condition N >= (y_max - y_min) / (4R * (1-tol)) vertically and M >= (x_max - x_min) / (4R * (1-tol)) vertically.
    # L can be height or width
    tmp = L/4R*(1-tol)
    bound = floor(tmp)
    if bound==tmp 
        return Int(bound) - 1
    else
        return Int(bound)
    end
end