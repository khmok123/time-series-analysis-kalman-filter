using DifferentialEquations
using Plots
using LinearAlgebra
using LaTeXStrings
using Statistics

include("plot_framework.jl")
include("models/infectionModels.jl")
include("kalman_filters/EnKF.jl")

function ℋ(r,params)

    H = [0. 0. 1. 0.; 0. 0. 0. 1.]

    H*r
end


# parameters for the infection models
S₀ = 70e3
E₀ = 400
I₀ = 500
R₀ = 0

β = 1.5
βₜ(t) = t > 14 ? .3 : .7

σ = 1/10
γ = 1/24
μₗ = .1

par = [βₜ,σ,γ,10.0]
u₀ = [S₀,E₀,I₀,R₀]


SEIRprob = ODEProblem(SEIR, u₀, (0.0,42.0),par)
sol = solve(SEIRprob,saveat=[0:42;],save_end=true, timeseries_errors = true)



scatter(sol)

sol[:][3:4]

y = [sol[i][3:4] for i ∈ 1:length(sol)]

res = EnKFIntegrator( EnKF!, u₀, sol.t, y, SEIR, ℋ, par, 100 )

plot(sol.t,res', label=[L"S_{kf}" L"E_{kf}" L"I_{kf}" L"R_{kf}"],ylabel="number of people" )

scatter!(sol, label=[L"S_{obs}" L"E_{obs}" L"I_{obs}" L"R_{obs}"])

plot_attributes!()

#ylims!(0,1e5)
