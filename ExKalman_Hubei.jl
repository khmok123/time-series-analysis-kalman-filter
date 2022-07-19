using DataFrames
using CSV
include("Kalman_extended_second_version.jl")


# Hubei
confirmed = DataFrame(CSV.File("../data/Hubei-confirmed.csv"))
cured = DataFrame(CSV.File("../data/Hubei-cured.csv"))
dead = DataFrame(CSV.File("../data/Hubei-dead.csv"))

recovered = cured.red.+dead.red
infected = confirmed.red .- recovered


recovered = recovered[end:-1:1,end:-1:1]
infected = infected[end:-1:1,end:-1:1]
confirmed

covariance_array = []
observations = []

for k in 1:length(infected)
    M = [infected[k] sqrt(infected[k])*sqrt(recovered[k]); sqrt(infected[k])*sqrt(recovered[k]) recovered[k]]
    push!(covariance_array, M)
    V = [infected[k], recovered[k]]
    push!(observations, V)
end

covariance_array

observations

β = 1.4
σ = 1/20
γ = 1/10

params = [β,σ,γ]

H = [0 0 1 0; 0 0 0 1]
size(H')
H*[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]*H'

length(observations)
implement_extended_kalman_filter(params, I, SEIR_update_exp, SEIR_jac, observations, H, covariance_array[1] ,I, length(observations), 1, length(observations))
