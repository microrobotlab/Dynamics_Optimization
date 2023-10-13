s = " * simplex_values: [-17.75725646101408, -17.757256494055728, -17.75725649581862]
* time: 0.0007789134979248047
* simplex: [[-1.4275593475262305, -0.0002642646225773252], [-1.4276091871477263, 3.876771027395157e-5], [-1.4274949110934978, -2.0336777782967104e-5]]
   39    -1.775726e+01     9.505832e-09
* simplex_values: [-17.75725651504421, -17.757256494055728, -17.75725649581862]
* time: 0.0007979869842529297
* simplex: [[-1.4275556983234212, -0.00012752457816591646], [-1.4276091871477263, 3.876771027395157e-5], [-1.4274949110934978, -2.0336777782967104e-5]]"

function extract_trace(s)
    output = []
    s_copy = s
    while(true)
        res = match(r"(?<=simplex: )(.*)(?<=])", s_copy)
        if(typeof(res) == Nothing)
            return output
        end
        push!(output, res.match)
        s_copy = s_copy[res.offset + length(res.match):end]
    end
    return output
end

extract_trace(s)