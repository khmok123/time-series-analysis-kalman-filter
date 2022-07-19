using DifferentialEquations
using Plots
using LinearAlgebra
using LaTeXStrings
using CSV
using DataFrames
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

par = [βₜ,σ,γ,0,0]
u₀ = [S₀,E₀,I₀,R₀]


confirmed = DataFrame(CSV.File("../data/Hubei-confirmed.csv"))
cured = DataFrame(CSV.File("../data/Hubei-cured.csv"))
dead = DataFrame(CSV.File("../data/Hubei-dead.csv"))


recovered = Float64.(cured.cases.+dead.cases)
infected = Float64.(confirmed.cases .- recovered)

recovered = recovered[end:-1:1,end:-1:1]
infected = infected[end:-1:1,end:-1:1]
dates = confirmed.date[end:-1:1,end:-1:1]
covariance_array = []
observations = []

for k in 1:length(infected)
    V = [Float64(infected[k]), Float64(recovered[k])]
    push!(observations, V)
end


# correlation between the recovered and infected cases
corr = cor(recovered, infected)

# correlation matrix
corM = [1 corr; corr 1]


# acceptance and efficiency for infected
acceptanceᵢ = 1-.2
efficiencyᵢ = .99

σ_acceptanceᵢ = .1
σ_efficiencyᵢ = .005


# acceptance and efficiency for recovered
acceptanceᵣ = 1-.2
efficiencyᵣ = .999
σ_acceptanceᵣ = .1
σ_efficiencyᵣ = .001

# true numbers of infected and recovered
infected ./= (acceptanceᵢ*efficiencyᵢ)
recovered ./= (acceptanceᵣ*efficiencyᵣ)

σᵢ = infected.*sqrt((σ_acceptanceᵢ/acceptanceᵢ)^2+(σ_efficiencyᵢ/efficiencyᵢ)^2)
σᵣ = recovered.*sqrt((σ_acceptanceᵣ/acceptanceᵣ)^2+(σ_efficiencyᵣ/efficiencyᵣ)^2)

sigmas = [[σᵢ[k] σᵣ[k]] for k in 1:length(σᵣ)]

D = [diagm([σᵢ[k]; σᵣ[k]]) for k in 1:length(σᵣ)]

covariance_array = [D[k]*corM*D[k] for k in 1:length(σᵣ)]

covariance_array

observations

time_series = collect(1.:length(observations))

y = observations

res = EnKFIntegrator( EnKF!, u₀, time_series, y, SEIR, ℋ, par, 100, covariance_array)

plot(dates,res', label=[L"S_{kf}" L"E_{kf}" L"I_{kf}" L"R_{kf}"],ylabel="number of people" )

scatter!(dates, infected, label=L"I_{obs}", marker =:d)
scatter!(dates, recovered, label=L"R_{obs}",marker =:d)

plot(time_series, res', label=[L"S_{kf}" L"E_{kf}" L"I_{kf}" L"R_{kf}"], color = ["blue" "red" "green" "fuchsia"],ylabel="number of people", ylims=(0,1e5) )
scatter!(time_series, vec(recovered), yerr = σᵣ, label=L"R_{obs}", color = "fuchsia")
scatter!(time_series, vec(infected),yerr = σᵢ, label=L"I_{obs}", color = "green")
plot_attributes!()
plot!(legend=:topleft)
ylims!(0,1e5)
