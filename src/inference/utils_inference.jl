# UTILITY FUNCTIONS FOR AUTOMATIC INFERENCE OF PARTICLE DYNAMICS 

"""
    tup2matrix_states(sim_output::Tuple)

Format simulator output : from tuple ( [ [x₀ y₀], [x₁ y₁], ..., [xₙ yₙ] ], [ [θ₀], [θ₁], [θₙ] ] ) matrix with time in column and state x y θ in rows (inverse of `matrix2tup_states` function)
"""
function tup2matrix_states(sim_output::Tuple)
    # concatenate x,y with θ for each timestep
    concat_matrices =  hcat.(sim_output[1], sim_output[2])
    # concatenate the rows of the matrices to have flat state vectors [x₁,y₁,θ₁,...,xₙ, yₙ, θₙ] 
    # -> one matrix contains x,y,θ for each particle; we want one vector per timestep
    flat_states_list = [vcat(eachrow(m)...) for m in concat_matrices]
    # finally concatenate flat states horizontally in one matrix
    return hcat(flat_states_list...)
end


"""
    matrix2tup_states(matrix_output, nb_states=3)

Format simulator output : from matrix with time in column and state x y θ in rows to tuple ( [ [x₀ y₀], [x₁ y₁], ..., [xₙ yₙ] ], [ [θ₀], [θ₁], [θₙ] ] ) (inverse of `tup2matrix_states` function)
"""
function matrix2tup_states(matrix_output, nb_states=3) # by default we have x,y,θ for the states
    # same steps as in tup2matrix_states in reverse order
    flat_states_list = eachcol(matrix_output)
    concat_matrices = [Matrix(reshape(v, (3,:))') for v in flat_states_list]
    return [m[:,1:2] for m in concat_matrices], [m[:,3] for m in concat_matrices]
end