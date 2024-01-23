# UNIT TESTS FOR SIMULATOR CODE VALIDATION

using ProgressBars
using LinearAlgebra

# List of simulator parameters for testing 
Np_list = [1,2,5,10,20,50,100,200,300,500,650,850]
L_list = collect(range(20., 1200, 5))
R_list = collect(range(0.2, 10., 10))
v_list = collect(range(1., 20., 5))


# `initABPE_ellipse` function unit test
test_initABPE_ellipse()
    iter = ProgressBar(Iterators.product(Np_list, L_list, R_list, v_list))
    set_description(iter, "UNIT TEST : initABPE_ellipse")
    for (Np,L,R,v) in Iterators.product(Np_list, L_list, R_list, v_list)
        abpe, (dists, superpose, uptriang) = initABPE_ellipse(Np, L, R, v; N, M)
        flag = true
        # test particles pairwise distance matrix dimensions (must be number of particles Np)
        # and have distance 0 on diagonal (distance of particles to themselves on diag)
        # (see ABP main_parallel.jl)
        flag &= (size(dists) == (Np,Np) & diag(dists) == zeros(Np))
        # likewise, test particles superposition matrix dimensions 
        flag &= (size(superpose) == (Np,Np))
    end