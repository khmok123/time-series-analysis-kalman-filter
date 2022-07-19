using DifferentialEquations
using Plots
using LinearAlgebra
using BenchmarkTools
include("plot_framework.jl")

include("kalman_filters/EnKF.jl")

include("models/lorenzSystem.jl")



# the observation function
ℋ(r, params) = r

# initial values for parameters and starting point
u₀ = [1.0,0.0,0.0]

p_lorenz = [10.0,28.0,8.0/3.0, 1.]

prob_sde_lorenz = SDEProblem(lorenzequations!,σ_lorenzequations!,u₀,(0.0,2.0),p_lorenz)

sol = solve(prob_sde_lorenz,dt=1e-2,adaptive=false)
sol.retcode
plot(sol,vars=(1,2,3))

# setting up the lorenz equations as an ODEProblem and calculate them retrieving observational data
#lorenzprob = ODEProblem(lorenzequations, u₀, (0.0,5.0), p_lorenz )
#sol = solve(lorenzprob,dt=1e-2,adaptive=false)

# try to replicate the observational data by using the ensemble kalman filter
res = EnKFIntegrator( EnKF!, u₀, sol.t, sol.u, lorenzequations, ℋ, p_lorenz )

# plotting the results
plot(res[1,:],res[2,:],res[3,:], color="blue", label="kalman-filtered data", xlabel="x", ylabel="y", zlabel="z", title="lorenz-attractor")
plot!(sol,vars=(1,2,3), label="original data")
plot_attributes!()

savefig("lorenz.pdf")
