
using LinearAlgebra
using DifferentialEquations
using Plots
using Statistics

# Assistant functions

# The function that switches the dimension of an array of arrays
function switch_arr_dim(arr_of_arr)
    dim1 = size(arr_of_arr)[1]
    dim2 = size(arr_of_arr[1])[1]
    return [[arr_of_arr[i][j] for i in 1:dim1] for j in 1:dim2]
end

# A function that adds Gaussian noise to a solution from solver, given time step and time span
function add_noise(sol,sigma,dt,T)
    y_arr = [[y+randn()*sigma for y in sol(t)] for t in 0:dt:T]
    return y_arr, switch_arr_dim(y_arr)
end

# A function that accepts input of time (Float) and output the corresponding element
get_element(arr, t) = arr[trunc(Int,t/dt)+1]

# A function that serves to plot scatter points at intervals of n, returns the range and the array
function scatter_every_n(arr, n::Int)
    chosen =[]
    for i in 1:size(arr)[1]
        push!(chosen,arr[i][1:n:end])
    end
    return 0:dt*n:T, chosen
end

# Extended Kalman filter - prediction and updating steps

function predict(xₖ,f,Fₖ,Pₖ,Qₖ,repeat)
    for i in 1:repeat
        xₖ = f(xₖ)
        Pₖ = Fₖ*Pₖ*Fₖ' + Qₖ
    end
    return xₖ,Pₖ
end

function update(xₖ,yₖ,Fₖ,Pₖ,Hₖ,Rₖ)
    K = Pₖ*Hₖ'*(Hₖ*Pₖ*Hₖ'+Rₖ)^(-1)
    xₖ += K*(yₖ-Hₖ*xₖ)
    Pₖ -= K*Hₖ*Pₖ
    

    return xₖ,Pₖ
end

# A function that combines the prediction and updating steps

function Kalman(x,p,P,f,y,H,R,Q,repeat)
    jac = SEIR_jac(x,p)
    x,P = predict(x,f,jac,P,Q,repeat)
    x,P = update(x,y,jac,P,H,R)

    return x, P
end

# This is used for the ODE solver
function SEIR(du,u,p,t)
    β,σ,γ = p
    S,E,I,R = u

    du[1] = -β*S*I
    du[2] = β*S*I-σ*E
    du[3] = σ*E-γ*I
    du[4] = γ*I
end

# This function is to be thrown into the Kalman filter (in prediction step)
function SEIR_update(u)
    β,σ,γ = par
    S,E,I,R = u
    u[1] = u[1] - β*S*I*dt
    u[2] = u[2] + (β*S*I-σ*E)*dt
    u[3] = u[3] + (σ*E-γ*I)*dt
    u[4] = u[4] + γ*I*dt
    return u
end

# This function is an alternative to be thrown into the prediction step, but I haven't figured out how to cope with the jacobian
# and fit it into the prediction step
function SEIR_update_exp(u,jac)
    S,E,I,R = u
    u = exp(jac*dt)
    return u
end

# It is the Jacobian matrix of the SEIR model
function SEIR_jac(u,p)
    β,σ,γ = p
    S,E,I,R = u
    jac = [-β*I 0 -β*S 0 ; β*I -σ β*S 0 ; 0 σ -γ 0 ; 0 0 γ 0]
    return jac
end

# State variables and parameters initialization

S₀ = 0.9
E₀ = 0.1
I₀ = 0.
R₀ = 0.

σ = .4
β = .5
γ = .3

dt = 0.002  # time step for both ODE solver and Kalman filter
T = 50.0 # The algorithm will be run up to time T

par = [β,σ,γ];
init = [S₀,E₀,I₀,R₀]
tspan = (0.0,T)

# ODE solver, plots the solution of the model (no noise)
prob = ODEProblem(SEIR,init,tspan,par)
sol = solve(prob,reltol=1e-8, abstol=1e-8, dense = false);
plot(sol,tspan=(0.0,T),label = ["suceptible" "exposed" "infectious" "recovered"], xlabel = "t [days]", ylabel = "% of population")

sigma = 0.02   # The magnitude of Gaussian noise

# Generation of data with noise based on the model (simulated data)
y_arr, y_arr_noised = add_noise(sol,sigma,dt,T)
plt = plot(legend=false)
for i in 1:4
    plot!(0:dt:T, y_arr_noised[i])
end
display(plt)

# Kalman filter implementation

x = y_arr[1]
P = I
f = SEIR_update
H = I
R = I
Q = 0*I
repeat = 1

x_series = []  # This will save the estimation of x, in time step dt
P_series = []  # This will save the covariance matrix P

# Kalman(x,p,P,f,y,H,R,Q,repeat)

# Running the Kalman filter
for t in 0:dt:T
    x,P = Kalman(x,par,P,SEIR_update,get_element(y_arr,t),H,R,Q,1)
    push!(x_series,x)
    push!(P_series,P)
end

x_series_2 = switch_arr_dim(x_series)  # Just for easier plotting
plt2 = plot(legend=false)
for i in 1:4
    plot!(0:dt:T, x_series_2[i])
end
for i in 1:4
    range, chosen_scatter = scatter_every_n(y_arr_noised, 10)
    scatter!(range, chosen_scatter,markersize=0.01)
end
display(plt2)
