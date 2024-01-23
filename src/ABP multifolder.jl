using FilePaths


"""
    multipledir(path, ICS)

Create folders (`ICS` number of folder created) at `path` to store simulator outputs.
"""
function multipledir(path, ICS)
    rundirs = []
    for i=1:ICS
        rundir = joinpath(path, "run$i")
        if !isdir(rundir)
            mkdir(rundir)
            push!(rundirs, rundir) # Save the path of the created directory.
        else
            println("Directory already exists: ", rundir)
        end
    end
    if isempty(rundirs)
        error("All directories already exist, time for a cleanup!")
    end
    return rundirs # Return the list of created directories.
end
