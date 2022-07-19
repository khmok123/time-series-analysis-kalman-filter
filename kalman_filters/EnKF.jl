using LinearAlgebra
using DifferentialEquations
using Distributions


#=struct FilterProblem

    u₀
    Q
    M

    t
    T

    N



end=#


"""
    EnKF!( E, y, t, Δt, ℳ, ℋ, params, N, M, R=I; alg = Tsit5() )

calculates the ensemble kalman filter for a given ensemble of points `E`, some measurements `y`, model `ℳ` and observation `ℋ`. `M` and `N` are height and length of the Ensemble Array `E`, `R` is the covariance of the measurement `y` and `alg` defines which algorithm to integrate the model with.

"""
function EnKF!( E, y, t, Δt, ℳ, ℋ, params, N, M, R=I; alg = Tsit5() )


    # Calculates the forecast step
    prob = ODEProblem(ℳ,E,(t,t+Δt),params)
    E = solve(prob,alg,saveat=Δt,save_start=false).u[1]+1e-2*randn(Float64,M,N)

    # Initialization the anomalities (distances between points and mean) for model and observation
    X = E - repeat(sum(E,dims=2),1,N)/N
    Y = ℋ(E,params)-repeat(sum(ℋ(E,params),dims=2),1,N)/N

    # Calculates the kalman gain
    K = X*Y'*(Y*Y'+(N-1)*R)^(-1)

    # Calculates the analysis step
    E = E + K*(y*ones(N)' - ℋ(E,params))

    return E

end


"""
    EnKFIntegrator(EnKF!,u₀,t,y,ℳ,ℋ,params,N=1000,R=I,Q=I; alg = Tsit5() )

Integrates the ensemble Kalman filter for a given start point `u₀` in a given array of times `t` and observational data `y`, with model `ℳ` and  observation operator `ℋ`, for given Parameters `params`. `N` is the number of ensemble points used and `R` the covariance of the measurements, `Q` the covariance of the initial value, `alg` defines which algorithm is to be used by `DifferentialEquations`.
"""
function EnKFIntegrator(EnKF!,u₀,t,y,ℳ,ℋ,params,N=1000,R=I,Q=I; alg = Tsit5() )

    # Initialization of used lenghts and the solution and ensemble arrays
    T = length(t)
    M = length(u₀)
    sol = zeros(M,T)
    covariance = zeros(M,M,T)
    E = rand(MvNormal(u₀, Q),N)

    # Initialization of an ODE covering the Ensemble instead of a single state
    enODE( ϵ, ρ, τ = 0 ) = ensembleODE( ϵ, ρ, ℳ, τ )

    # Saving the start value
    sol[:,1] = ones(1,N)*E'/N

    # Iteration of the solution for a given length T
    for i ∈ 1:T-1

        # Initializing the time step and observation value
        Δt = t[i+1]-t[i]
        obs = y[i+1]

        # The single Ensemble Kalman Step
        E = EnKF!(E,obs,t[i+1],Δt,enODE,ℋ,params,N,M,R; alg = Tsit5())

        # Saving the mean and the covariance of the mean in their respective arrays
        sol[:,i+1] = sum(E,dims=2)/N
        X = E - repeat(sum(E,dims=2),1,N)/N
        covariance[:,:,i] = X*X'/(N-1)

        # Returns the timestep of the integration
        println(i)
    end

    return sol, covariance

end



"""
    ensembleODE( E, params, ℳ, t )

Calculates the differential equation `ℳ` for a given Ensemble of points `E`.
"""
function ensembleODE( E, params, ℳ, t )

    # Initializing value of Ensemble length and Ensemblederivative
    N = size(E,2)
    dE = zeros(size(E))

    # calculates the equations for each point solely
    for i in 1:N
        dE[:,i] = ℳ(E[:,i],params,t)
    end

    return dE

end
