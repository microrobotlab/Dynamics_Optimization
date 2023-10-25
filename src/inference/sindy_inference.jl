include("model_setup.jl")


# ----- INFERENCE WITH SINDY
full_problem = ContinuousDataDrivenProblem(Xₙ, t)
ideal_problem = DirectDataDrivenProblem(X̂, Ȳ)
nn_problem = DirectDataDrivenProblem(X̂, Ŷ)

# SINDy method
DataDrivenSparse.STLSQ