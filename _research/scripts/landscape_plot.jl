using DrWatson
@quickactivate "active-brownian-particles"
using CSV, DataFrames
using Plots


# just to plot the landscape exploration results
df = CSV.read(datadir("processed", "saved_valuations.csv"), DataFrame)

plot()
# separate between the 3 different values of radii
for df_Ri in groupby(df, :R)
    # plot mean polarization factor w.r.t. the velocity on a same plot, indicating the corresponding radius
    # we sort v and mean_polarization_factor w.r.t. v
    order = sortperm(df_Ri[!, :v])
    plot!(df_Ri[!, :v][order], df_Ri[!, :mean_polarization_factor][order], label="R=$(df_Ri[1, :R])")
end

ylims!(0, 1.0) # pf ∈ [0,1] so plot with this range to stay objective
# here we know it was averaged over 200 runs but if we want cleaner code 
# we should make the value data-dependent in title
NB_RUNS = 200
title!("Polarization factor averaged over time and $NB_RUNS runs w.r.t. velocity (for different radii)", titlefontsize=8)
xlabel!("velocity")
ylabel!("mean polarization factor")