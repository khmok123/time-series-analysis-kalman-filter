function ClKF!(x, P, F, y, H, R, Q)
"""
    KF!(x, P, F, y, H, R, Q)

    calculates the classic kalman filter for a given state`x`, some measurements `y`, linear model `F` and observation `H` with covariance `R`
    in the forcast Step gaussian noise with covariance `Q` is added.
"""


    x = F*x + Q*randn(length(x))
    P = F*P*F' + Q

    K = P*H'*(H*P*H'+ R)^(-1)
    x = x + K*(y.-H*x)
    P = P - K*H*P

    return x,P
end


function ClKFIntegrator(ClKF!, u₀, P₀, t, y, F, H, R, Q=0*I)
"""
    ClKFIntegrator(ClKF!, u₀, P₀, t, y, F, H, R, Q=0*I)

    Integrates the classic Kalman filter `ClKF` for a given start point `u₀` with covariance `P₀`
         in a given array of times `t` and observational data `y`, with linear model `F` and  observation operator `H` with covariance R`.
         Q is the covariance of the random gaussian forecing, for Q=0(default setting) no random force is added.
"""

    u_array = zeros(size(t)[1],2) #array for states
    P_array = zeros(size(t)[1],2,2) # array for covariaces

    u_array[1,:] = u₀
    P_array[1,:,:] = P₀
    u = u₀
    P = P₀

    for i in 1:length(t)
        u,P = ClKF!(u,P,F,y[i],H,R,Q)
        u_array[i,:] = u
        P_array[i,:,:] = P

    end
    return u_array, P_array
end
