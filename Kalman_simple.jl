using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics

function forcast(xₖ, Fₖ,Pₖ, Qₖ = I)

    xₖ = Fₖ*xₖ
    Pₖ = Fₖ*Pₖ*Fₖ' + Fₖ*Qₖ*Fₖ'

    return xₖ,Pₖ
end

function update(xₖ,yₖ,Fₖ,Pₖ,Hₖ,Rₖ = I)
   
    K = Pₖ*Hₖ'*(Hₖ*Pₖ*Hₖ'+Rₖ)^(-1)
    xₖ += K*(yₖ-Hₖ*xₖ)
    Pₖ -= K*Hₖ*Pₖ
    
    return xₖ,Pₖ
    
end

function KalmanStep!(x,P,F, y,H,R, Q=0*I) 
    x = F*x #+ 
    P = F*P*F' + Q 

    K = P*H'*(H*P*H'+ R)^(-1)
    x = x + K*(y-H*x)
    P = P - K*H*P
    return x,P
end  

function SEIRS_paper(du, u, p,t)
    β,ν,μₗ,α = p
    S,E,I,R = u
    
    du[1] = -β*S*I
    du[2] = β*S*I-ν*E
    du[3] = ν*E-(1-α)*I
    du[4] = (1-μₗ-α)*I
    
end

# parameters of SEIRS_paper
S₀ = 1000.
E₀ = 10.
I₀ = 1.
R₀ = 2.

ν = .4
β = .5
α = .3
μₗ = .1

par = [β,ν,μₗ,α];
init = [S₀,E₀,I₀,R₀]
tspan = (0.0,10.0)

# Integration of nonlinear dynamics
prob = ODEProblem(SEIRS_paper,init,tspan,par)
sol = solve(prob,reltol=1e-8, abstol=1e-8, dense = false);
plot(sol,tspan=(0.0,10.0),label = ["suceptible" "exposed" "infectious" "recovered"], xlabel = "t [days]", ylabel = "% of population")

ρ₀ = S₀*β
A(ρ₀) = [-ν ρ₀; ν 1-α]

H = [1 0; 0 1]

noise = [randn(Float64,4) for k in 1:length(sol)]
observations = sol.u .+ noise*1

obs = observations[1]

for k ∈ 2:length(observations)
    obs = hcat(obs,observations[k])
end

plot(sol.t, obs[1,:],label = ["suceptible" "exposed" "infectious" "recovered"], xlabel = "t [days]", ylabel = "% of population")
plot!(sol.t, obs[2,:],label = ["exposed" "infectious" "recovered"], xlabel = "t [days]", ylabel = "% of population")
plot!(sol.t, obs[3,:],label = ["infectious" "recovered"], xlabel = "t [days]", ylabel = "% of population")
plot!(sol.t, obs[4,:],label = "recovered", xlabel = "t [days]", ylabel = "% of population")

F = A(ρ₀)
Δt = sol.t[2]-sol.t[1]
x, P = forcast([E₀, I₀],exp(F*Δt),I)

x_opt = x

for (t,i) in zip(sol.t, 2:length(sol.t))
    Δt = sol.t[i]-sol.t[i-1]
    ρ₀ = β*obs[1,i]
    F = A(ρ₀)
    x, P = KalmanStep!(x,P,exp(F*Δt), obs[2:3,i], H, I)
    #x, P = forcast(x, exp(F*Δt), P)
    #x, P = update(x, [obs[3,i]], exp(F*Δt), P, H, I*Δt)
    x_opt = hcat(x_opt,x)
end

#sol.u[2,1:5]

plot(sol.t,x_opt[1,1:length(sol.t)], label = "exposed")
plot!(sol.t, obs[2,:],label = ["exposed (model)" "infectious" "recovered"], xlabel = "t [days]", ylabel = "% of population")
plot!(sol.t,x_opt[2,1:length(sol.t)], label = "infected")
plot!(sol.t, obs[3,:],label = "infected in model", xlabel = "t [days]", ylabel = "% of population")




