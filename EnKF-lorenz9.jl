using DifferentialEquations
using LinearAlgebra
using Plots
using RecursiveArrayTools

include("models/lorenzSystem.jl")
include("kalman_filters/EnKF.jl")
include("kalman_filters/ParKF.jl")
include("plot_framework.jl")

# Definition of integration parameters

int_start = 0.
int_end = 600.
filter_start = 500.
filter_end = 600.

tol = 1e-6

σ_SDE = 5e-3
En_size = 500

filter_dt = 1e-1

first_obs = 6
last_obs = 6


#Definition of plot parameters

plot_start = Int64(filter_start)
plot_end = Int64(filter_end)

fplot_end = Int64((plot_end-plot_start)/filter_dt)+1


# Definition of observation function



#ℋ( u, params ) = [0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0 0]*u
#ℋ( u, params ) = u
function ℋ( u, params )

    a = zeros(1,9)
    index = Int16(params[4])

    a[index] = 1
    a*u
end

# Definition of initial values and parameters

u₀ = [ 0.01, 0, 0.01, 0, 0, 0, 0, 0, 0.01 ]

σ = 0.5
a = 0.5
r = 15.1


par_lorenz9 = [ 0.5, r, 0.5, first_obs, σ_SDE]
# Integration of Lorenz9-Problem

prob_lorenz9 = SDEProblem( lorenz9!, σ_lorenz9!, u₀, ( int_start, int_end ), par_lorenz9 )

sol_lorenz9 = solve(prob_lorenz9, saveat = filter_start:filter_dt:filter_end)



# Refinement of data for Filtering

time_series = [filter_start:filter_dt:filter_end;]

y = [sol_lorenz9(t)[first_obs:last_obs] for t ∈ time_series]

u₀ = [0.,0.,0.,0.,0.,1.,0.,0.,0.]#rand(Float64,9)


# Filtering the given observation with the model

@time res, var = EnKFIntegrator( EnKF!, u₀, time_series, y, lorenz9, ℋ, par_lorenz9, En_size)


# Plotting the filtered results with the actual observations


plot_tit = "\\sigma = $σ_SDE, r = $r, input: C$first_obs - C$last_obs, dt = $filter_dt, N = $En_size"
var
for i ∈ 1:9
    plot(time_series, res[i,:],
    label = "filtered data",
    fill = (res[i,:],res[i,:].+var[i,i,:], 0.2, :blue),
    color=:blue,
    title=plot_tit,
    yaxis="C$i",
    xlims=(plot_start,plot_end),
    linewidth=3)

    plot!(time_series, res[i,:],
    label = nothing,
    fill = (res[i,:],res[i,:].-var[i,i,:], 0.2, :blue),
    color=:blue,
    linewidth = 3)

    plot!(sol_lorenz9,
    vars=(i),
    tspan=(plot_start,plot_end),
    label="model data",
    linewidth = 3)
    plot_attributes!()
    savefig("../results/Enl9_$(first_obs)$(last_obs)_$(r)_$(plot_start)-$(plot_end)_$(filter_dt)_$(σ_SDE)_$i.pdf")
end
sol_lorenz9.u
VA = VectorOfArray(sol_lorenz9.u)
a = convert(Array,VA)

rmsd = reshape(norm_by_dim(a-res,1),fplot_end)
time_series
plot(time_series,rmsd, xaxis="time", yaxis="RMSD", title=plot_tit, legend = nothing,xlims = (plot_start,plot_end))
plot_attributes!(false)

savefig("../results/RMSE_Enl9_$(first_obs)$(last_obs)_$(r)_$(plot_start)-$(plot_end)_$(filter_dt)_$σ_SDE.pdf")
