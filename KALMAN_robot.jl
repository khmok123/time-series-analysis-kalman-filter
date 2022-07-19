using LinearAlgebra
using DifferentialEquations
using StaticArrays
using Plots
include("plot_framework.jl")
include("kalman_filters/ClKF.jl")

#solution of the model
function robot(u, nothing, t)
    x,v = u
    dx = v
    dv = 0
    return [dx, dv]
end
prob = ODEProblem(robot, [0.0,1.], (0.0, 10.0), Δt)
alg = Vern9()
sol = solve(prob; alg = alg)
plot(sol)

function KalmanStep!(x,P,F, y,H,R, Q=0*I)
    x = F*x #+
    P = F*P*F' + Q

    K = P*H'*(H*P*H'+ R)^(-1)
    x = x + K*(y-H*x)
    P = P - K*H*P

    return [x,P]
end


prob = ODEProblem(robot, [0.0,1.], (0.0, 10.0))
sol = solve(prob)


##generating observed data
t = 0:0.5:10

obs = zeros(size(t,1))
for (i, τ) in enumerate(t)
    obs[i] = sol(τ)[1] + sigma * randn()
end


##Kalmanfitler
F = [1. Δt;0. 1.] #forcing

H = [1. 0.] #observation matrix
R = [(sigma)^(1/2)] #covarianz matrix

u₀ = [0.,1.] #starting state
P₀ = [1. 0.; 0. 1.] #starting covariance




x, P = ClKFIntegrator(ClKF!, u₀, P₀,t,obs,F,H,R,1* I)



##plotting data


#y = [x[:,1],x[:,2],obs]
#plot(t,y,label=["x_Kalman" "v_Kalman" "x_observed" "v_observed"], layout = (2,1))
#plot_attributes!()
#plotting data with errorbars
plot(t,x[:,1], fill = (x[:,1].-P[:,1,1],0.3, :red),color=:red,label="Kalman",
    legend=:topleft,xlabel="time t",xlims=(0,10),ylabel="position x",title="position of the robot")
plot!(t,x[:,1], fill = (x[:,1].+P[:,1,1],0.3, :red),label=false)

scatter!(t,obs, yerror=[R[1,1],R[1,1]],color=:blue, label ="observed")

#plot!(t,obs[:,1], fill = (obs[:,1].-R[1,1],0.1,:blue),color=:blue, label ="observed")
#plot!(t,obs[:,1], fill = (obs[:,1].+R[1,1],0.1,:blue),label=false)
plot_attributes!()

#savefig("classic_Kalman:position_robot.pdf")


#plot_attributes!()
