using LinearAlgebra
using DifferentialEquations
using StatsBase
using Distributions

norm_by_dim(A,d) = mapslices(norm, A, dims=d)

function Nₑ(w)
    sq = w.*w
    1.0/sum(sq)
end

function resample(particles, w, meandist, M,N)

    # create probability distribution of particles. Since DiscreteNonParametric does not work with multidimensional arrays, the distribution is created for the first row (first variable of ensemble here)
    distribution = DiscreteNonParametric(particles[1,:], w)
    # Here, the resampled particles will be entered
    resampled_particles = zeros(M,N)
    resampled_weights = zeros(N)

    μ = mean(particles, Weights(w), dims=2)
    σ² = mean((particles.-μ).^2,Weights(w), dims = 2)
    if maximum(σ²) < .5
        σℳ = 0.01+0.0001*meandist
    elseif maximum(σ²) < 1
        σℳ = 0.05+0.0001*meandist
    else
        σℳ = 0.1+0.0001*meandist
    end

    for i ∈ 1:N

        # draw random variable from distribution
        val = rand(distribution)

        # check index of random variable in corresponding array
        ind = findall(x -> x == val, particles[1,:])[1]

        # put particle vector with that index into the resampled_particles
        resampled_particles[:,i] = particles[:,ind].+σℳ*randn(M)
        resampled_weights[i] = w[ind]
    end
    particles = resampled_particles
    w = resampled_weights/sum(resampled_weights)

    return particles, w
end

function strat_resample(particles, weights)

    M = size(particles,1)
    N = size(particles,2)
    resampled_particles = zeros(M,N)
    resampled_weights = zeros(N)
    r = rand(Uniform(),N)/N
    c = w[1]
    k = 1
    for i ∈ 1:N
        U = r[i]+(i-1)/N
        while U > c
            k += 1
            c += w[k]
        end
        resampled_particles[:,i] = particles[:,k]
        resampled_weights[:,i] = weights[k]
    end
    return resampled_particles, resampled_weights
end


"""
    ParKF!(E, w, y, t, Δt, ℳ, ℋ, params, σ, σℳ)

-------------------parameters--------------------\n
`E`     :   Ensemble of particles\n
`w`     :   Array of weights of particles\n
`y`     :   Array of observations\n
`t`     :   Array of times at which observations take place\n
`Δt`    :   Time step between observations\n
`ℳ`    :   Model of dynamics\n
`ℋ`    :   Observation operator, projecting from state space to observation space\n
`params`:   parameters of model\n
`σ`

"""
function ParKF!(E,w, y, t, Δt, ℳ, ℋ, params, M, N)

    prob = ODEProblem(ℳ,E,(t,t+Δt),params)

    # prediction
    E = solve(prob,saveat=t+Δt,save_start=false).u[1]#+10*randn(Float64,M,N)
    x = mean(E, Weights(w), dims=2)

    μ = mean(E, Weights(w), dims=2)
    σ² = mean((E.-μ).^2,Weights(w), dims = 2)

    # update/ Sequential importance Sampling
    distance = norm_by_dim(ℋ(E,params).-y,1)
    μdist = mean(distance, Weights(w))
    if μdist < 0.01
        σ = 0.05
    elseif μdist < 0.05
        σ = 0.2
    else
        σ = 0.4
    end
    w .*= pdf.(Normal(0., σ),distance)[1,:]
    w /= sum(w)

    # Resampling
    if Nₑ(w) < N/4
        E, w = resample(E, w, μdist, M, N)
    end

    return E,w

end

"""
    ParKFIntegrator(ParKF!, u₀, t, y, ℳ, ℋ, params, N=1000, σ = [1 1], σℳ = 100)

-------------------parameters--------------------\n
`ParKF!`  :   Particle Filter\n
`u₀`      :   initial condition\n
`t`       :   array of times where observations take place\n
`y`       :   array of observations\n
`ℳ`      :   Model for dynamics\n
`ℋ`      :   observation operator, projects from signal space to
            observation space\n
`params`  :   paramters of `ℳ`\n
`N`       :   Number of ensemble points\n
`σ`       :   Standard deviation of observations, either `size`(`σ`) =
            `dim`(`y`) or `size`(`σ`) = `dim`(`y`) x `length`(`t`)\n
`σℳ`     :   standard deviation of model noise\n
"""
function ParKFIntegrator(ParKF!,u₀,t,y,ℳ,ℋ,params;N=10000,σℳ = 2)


    # Initialization of used lenghts and the solution and ensemble array
    T = length(t)
    M = length(u₀)
    sol = zeros(M,T)
    E = σℳ*randn(Float64,M,N) .+ u₀
    enODE( ϵ, ρ, τ = 0 ) = ensembleODE( ϵ, ρ, ℳ, τ )
    sol[:,1] = ones(1,N)*E'/N
    variance = zeros(M,T)
    w = ones(N)/N           # initialize weights as microcanonical ensemble
    x = u₀
    # Iteration of the solution for a given length T
    for i ∈ 1:T-1
        Δt = t[i+1]-t[i]
        obs = y[i+1]
        E,w = ParKF!(E,w,obs,t[i+1],Δt,enODE,ℋ,params,M,N)
        μ = mean(E, Weights(w), dims=2)
        σ² = mean((E.-μ).^2,Weights(w), dims = 2)
        sol[:,i+1] = μ
        variance[:,i+1] = σ²
    end

    return sol, sqrt.(variance)

end


function ensembleODE( E, params, ℳ, t )

"""
    ensembleODE( E, params, ℳ, t = 0 )

    Calculates the differential equation `ℳ` for a given Ensemble of points `E`
"""
    N = size(E,2)
    dE = zeros(size(E))

    # calculates the lorenz equations for each point solely
    for i in 1:N
        dE[:,i] = ℳ(E[:,i],params,t)
    end

    return dE

end
