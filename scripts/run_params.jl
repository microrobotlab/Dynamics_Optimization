# To call run function (see ABP output) with given parameters directly

using DrWatson
@quickactivate "active-brownian-particles"
println("Currently in $(projectdir()) environment !")

include(srcdir("ABP output.jl"))

#---------------------------------------------------------------------------------------------------------------------
# CALL RUNNER WITH ARGUMENTS 

# simulation parameters
Nt = 3000; Np = 1; L = 100.; R = 1.5; v = 10.

# macro-parameters
wall_condition = "open"; nb_runs = 1

# CSV export parameters: `save_stride` is the stride for timesteps while saving simulator output 
# /!\ SAVE WON'T BE TAKEN INTO ACCOUNT FOR nb_runs > 1
save = false; save_stride = 1 

# animation parameters: `animation_stride` stride for the animation, by default `animation_filename` will have the same marker as 
# simulation file output (see file ), no choice in the case of file export (will be based on the data filename) 
animate = true; animation_filename = nothing; animation_stride = 10

# parallelization parameters: (N, M) number of cells (rows, columns)
N = 16; M = 16;

# simulation progress bar display
verbose = true


if(nb_runs==1 && !save)
    # in this case return the simulator output directly
    run((Nt=Nt,Np=Np,L=L,R=R,v=v); wall_condition=wall_condition, animate=animate, animation_stride=animation_stride, N=N, M=M, verbose=verbose)
else
    # in this case return the folder path containing the saved runs
    run_multiple((Nt=Nt,Np=Np,L=L,R=R,v=v); wall_condition=wall_condition, 
        nb_runs=nb_runs, 
        save_stride=save_stride, 
        animate=animate, animation_stride=animation_stride,
        N=N, M=M,
        verbose=verbose);
end