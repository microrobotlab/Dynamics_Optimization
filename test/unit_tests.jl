using ProgressBars
using LinearAlgebra


Np_list = [1,2,5,10,20,50,100,200,300,500,650,850]
L_list = collect(range(20., 1200, 5))
R_list = collect(range(0.2, 10., 10))
v_list = collect(range(1., 20., 5))


test_initABPE_ellipse()
    iter = ProgressBar(Iterators.product(Np_list, L_list, R_list, v_list))
    set_description(iter, "UNIT TEST : initABPE_ellipse")
    for (Np,L,R,v) in Iterators.product(Np_list, L_list, R_list, v_list)
        abpe, (dists, superpose, uptriang) = initABPE_ellipse(Np, L, R, v; N, M)
        flag = true
        flag &= (size(dists) == (Np,Np) & diag(dists) == zeros(Np)) # tests dists
        flag &= (size(superpose) == (Np,Np))
    end