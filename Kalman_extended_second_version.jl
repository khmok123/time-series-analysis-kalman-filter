
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

# A function that accepts an array of covariance matrices and the dimension, and output the list of standard deviation
function get_sd_from_P_series(arr, i)
    length = size(arr)[1]
    return [arr[j][i,i] for j in 1:length].^0.5
end

###########################################################

# Extended Kalman filter - prediction and updating steps

function predict!(xₖ,p,f,jac,Pₖ,Qₖ,repeat)
    Fₖ = jac(xₖ,p)
    for i in 1:repeat
        xₖ = f(xₖ)
        Pₖ = Fₖ*Pₖ*Fₖ' + Qₖ
    end
    return xₖ,Pₖ
end

function update!(xₖ,p,yₖ,jac,Pₖ,Hₖ,Rₖ)
    K = Pₖ*Hₖ'*(Hₖ*Pₖ*Hₖ'+Rₖ)^(-1)
    xₖ += K*(yₖ-Hₖ*xₖ)
    Pₖ -= K*Hₖ*Pₖ

    return xₖ,Pₖ
end

# A function that combines the prediction and updating steps

function Kalman!(x,p,P,f,jac,y,H,R,Q,repeat)
    x,P = predict!(x,p,f,jac,P,Q,repeat)
    x,P = update!(x,p,y,jac,P,H,R)

    return x, P
end

# A function that inputs systems of equations, parameters, initial conditions and timespan to produce analytical solution (sol)
function solve_system(system,init,T,p,label;plotstyle="2D",xlabel="",ylabel="")
    prob = ODEProblem(system,init,(0.0,T),p)
    sol = solve(prob,reltol=1e-8, abstol=1e-8, dense = false);
    if plotstyle=="2D"
        display(plot(sol,tspan=(0.0,T),label=label,xlabel=xlabel,ylabel=ylabel))
    elseif plotstyle=="3D"
        display(plot(sol,vars=(1,2,3),tspan=(0.0,T),legend=false))
    end
    return sol
end

function generate_noised_data(sol,sigma,dt,T,xlabel,ylabel;plotstyle="2D")
    y_arr, y_arr_noised = add_noise(sol,sigma,dt,T)
    if plotstyle == "2D"
        plt = plot(legend=false,xlabel=xlabel,ylabel=ylabel)
        for i in 1:size(sol[1])[1]
            plot!(0:dt:T, y_arr_noised[i])
        end
        display(plt)
    elseif plotstyle == "3D"
        plt = plot(y_arr_noised[1], y_arr_noised[2], y_arr_noised[3], legend=false)
        display(plt)
    end
    return y_arr, y_arr_noised
end

function implement_extended_kalman_filter(p,P,f,jac,y_arr,H,R,Q,repeat,dt,T)
<<<<<<< HEAD
<<<<<<< HEAD
    x = [80e6, 4000.,500.,0.]
=======
    x = y_arr[1]
>>>>>>> 1a603e02bc7657b5867662cba60031aec18fac63
=======
    x = [80e6, 4000.,500.,0.]
>>>>>>> 6386eb8ac0226d5b265a8cd2a9c484764171d3b3
    x_series = []  # This will save the estimation of x, in time step dt
    P_series = []  # This will save the covariance matrix P

    # Kalman(x,p,P,f,jac,y,H,R,Q,repeat)

    # Running the Kalman filter
    for t in 0:dt:T
<<<<<<< HEAD
<<<<<<< HEAD
        x,P = Kalman!(x,p,P,f,jac,get_element(y_arr,t),H,R,Q,repeat)
=======
        x,P = Kalman(x,p,P,f,jac,get_element(y_arr,t),H,R,Q,repeat)
>>>>>>> 1a603e02bc7657b5867662cba60031aec18fac63
=======
        x,P = Kalman!(x,p,P,f,jac,get_element(y_arr,t),H,R,Q,repeat)
>>>>>>> 6386eb8ac0226d5b265a8cd2a9c484764171d3b3
        push!(x_series,x)
        push!(P_series,P)
    end
    return x_series, P_series
end

function plot_kalman(x_series,P_series,dt,T;scatter_interval=20,plotstyle="2D")
    x_series_2 = switch_arr_dim(x_series)  # Just for easier plotting
    plt2 = plot(legend=false)
    if plotstyle == "2D"
        for i in 1:size(sol[1])[1]
            plot!(0:dt:T, x_series_2[i])
        end
        for i in 1:size(sol[1])[1]
            range, chosen_scatter = scatter_every_n(y_arr_noised, scatter_interval)
            scatter!(range, chosen_scatter,markersize=0.01)
        end
        display(plt2)
    elseif plotstyle == "3D"
        plt2 = plot(sol,vars=(1,2,3), legend=false, title="Model and filtered")  # model
        plot!(y_arr_noised[1], y_arr_noised[2], y_arr_noised[3], legend=false)
        display(plt2)
        plt3 = plot(x_series_2[1], x_series_2[2], x_series_2[3], legend=false, title="Noised and filtered")
        plot!(y_arr_noised[1], y_arr_noised[2], y_arr_noised[3], legend=false)
        display(plt3)
    end
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
function SEIR_update_exp(u)
    S,E,I,R = u
    jac = SEIR_jac(u,par)
    u = exp(jac.*dt)*u
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

system = SEIR
f = SEIR_update
jac = SEIR_jac

S₀ = 0.9
E₀ = 0.1
I₀ = 0.
R₀ = 0.

σ = .4
β = .5
γ = .3

dt = 0.002  # time step for both ODE solver and Kalman filter
T = 50.0 # The algorithm will be run up to time T

sigma = 0.02

P = I
f = SEIR_update
H = I
R = I
Q = 0.01*I
repeat = 1
scatter_interval = 20

par = [β,σ,γ];
init = [S₀,E₀,I₀,R₀]
label =  ["suceptible" "exposed" "infectious" "recovered"]
xlabel = "t [days]"
ylabel = "% of population"
sol = solve_system(system,init,T,par,label;plotstyle="2D",xlabel=xlabel,ylabel=ylabel)
y_arr, y_arr_noised = generate_noised_data(sol,sigma,dt,T,xlabel,ylabel;plotstyle="2D")
x_series, P_series = implement_extended_kalman_filter(par,P,f,jac,y_arr,H,R,Q,repeat,dt,T)
plot_kalman(x_series,P_series,dt,T;plotstyle="2D", scatter_interval=20)

# Now experiment with Lorenz system

function lorenz(du,u,p,t)
    du[1] = a*(u[2]-u[1])
    du[2] = u[1]*(b-u[3]) - u[2]
    du[3] = u[1]*u[2] - c*u[3]
end

function lorenz_update(u)
    a,b,c = par_lorenz
    x,y,z = u
    u[1] = u[1] + a*(u[2]-u[1])*dt
    u[2] = u[2] + (u[1]*(b-u[3]) - u[2])*dt
    u[3] = u[3] + (u[1]*u[2] - c*u[3])*dt
    return u
end

function lorenz_jac(u,p)
    a,b,c = p
    x,y,z = u
    jac = [-a a 0 ; b-z -1 -x ; y x -c]
    return jac
end

system = lorenz
f = lorenz_update
jac = lorenz_jac


# State variables and parameters initialization

x₀ = 1.0
y₀ = 0.0
z₀ = 0.0

a = 10.
b = 28.
c = 8/3

dt = 0.02  # time step for both ODE solver and Kalman filter
T = 50.0 # The algorithm will be run up to time T

par_lorenz = [a,b,c];
init = [x₀,y₀,z₀]
tspan = (0.0,T)
label=""

sigma = 0.5

P = I
f = lorenz_update
H = I
R = I
Q = 0*I
repeat = 1
jac = lorenz_jac

# ODE solver, plots the solution of the model (no noise)
sol = solve_system(system,init,T,par_lorenz,label;plotstyle="3D")
y_arr, y_arr_noised = generate_noised_data(sol,sigma,dt,T,xlabel,ylabel;plotstyle="3D")
x_series, P_series = implement_extended_kalman_filter(par_lorenz,P,f,jac,y_arr,H,R,Q,repeat,dt,T)
plot_kalman(x_series,P_series,dt,T;plotstyle="3D")
