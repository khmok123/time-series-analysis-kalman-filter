using LinearAlgebra
using Documenter


"""
    lorenzequations( u, params, t = 0 )

Calculates the differential lorenz equations at a given point.
"""
function lorenzequations( u, params, t = 0 )

    # initializes the variables and parameters
    α, β, γ = params
    x, y, z = u

    #calculates the equations
    ẋ = α*(y-x)
    ẏ = x*(β-z)-y
    ż = x*y - γ*z

    [ẋ,ẏ,ż]

end


"""
    lorenzequations!( du, u, params, t = 0 )

Writes the differential lorenz equations at a given point into `du`.
"""
function lorenzequations!( du, u, params, t = 0 )

    # initializes the variables and parameters
    α, β, γ, σ = params
    x, y, z = u

    #calculates the equations
    du[1] = α*(y-x)
    du[2] = x*(β-z)-y
    du[3] = x*y - γ*z

end


"""
    σ_lorenzequations!( du, u, p, t )

The standard noise-function for the lorenz equations. Writes a given noise `σ` into the variable `du`
"""
function σ_lorenzequations!( du, u, p, t )

    α, β, γ, σ = p

    du[1] = σ
    du[2] = σ
    du[3] = σ

end


"""
    lorenz9( u, p, t )

Returns the expansion of the lorenz system to a system of nine equations as done in Reiterer(1998).
"""
function lorenz9( u, p, t )

    # initializes variables and parameters
    σ, r, a = p
    b1, b2, b3, b4, b5, b6 = lorenz9params(a)
    C1, C2, C3, C4, C5, C6, C7, C8, C9 = u

    # calculates the lorenz1998 equations
    dC1dt = - σ*b1*C1 - C2*C4 + b4*C4^2 + b3*C3*C5 - σ*b2*C7
    dC2dt = - σ*C2 + C1*C4 - C2*C5 + C4*C5 - σ*C9/2
    dC3dt = - σ*b1*C3 + C2*C4 - b4*C2^2 - b3*C1*C5 + σ*b2*C8
    dC4dt = - σ*C4 - C2*C3 - C2*C5 + C4*C5 + σ*C9/2
    dC5dt = - σ*b5*C5 + C2^2/2 - C4^2/2
    dC6dt = - b6*C6 + C2*C9 - C4*C9
    dC7dt = - b1*C7 - r*C1 + 2*C5*C8 - C4*C9
    dC8dt = - b1*C8 + r*C3 - 2*C5*C7 + C2*C9
    dC9dt = - C9 - r*C2 + r*C4 - 2*C2*C6 + 2*C4*C6 + C4*C7 - C2*C8

    [dC1dt,dC2dt,dC3dt,dC4dt,dC5dt,dC6dt,dC7dt,dC8dt,dC9dt]

end



"""
    lorenz9( u, p, t )

Writes the expansion of the lorenz system to a system of nine equations as done in Reiterer(1998) into `du`.
"""
function lorenz9!( du, u, p, t )

    # initializes variables and parameters
    σ, r, a = p
    b1, b2, b3, b4, b5, b6 = lorenz9params(a)
    C1, C2, C3, C4, C5, C6, C7, C8, C9 = u

    # calculates the lorenz1998 equations
    du[1] = - σ*b1*C1 - C2*C4 + b4*C4^2 + b3*C3*C5 - σ*b2*C7
    du[2] = - σ*C2 + C1*C4 - C2*C5 + C4*C5 - σ*C9/2
    du[3] = - σ*b1*C3 + C2*C4 - b4*C2^2 - b3*C1*C5 + σ*b2*C8
    du[4] = - σ*C4 - C2*C3 - C2*C5 + C4*C5 + σ*C9/2
    du[5] = - σ*b5*C5 + C2^2/2 - C4^2/2
    du[6] = - b6*C6 + C2*C9 - C4*C9
    du[7] = - b1*C7 - r*C1 + 2*C5*C8 - C4*C9
    du[8] = - b1*C8 + r*C3 - 2*C5*C7 + C2*C9
    du[9] = - C9 - r*C2 + r*C4 - 2*C2*C6 + 2*C4*C6 + C4*C7 - C2*C8

end

"""
    σ_lorenzequations!( du, u, p, t )

The standard noise-function for the lorenz1998 equations. Writes a given noise `σ` into the variable `du`.
"""
function σ_lorenz9!(du,u,p,t)

    σ = p[end]

    du .= σ

end


"""
    lorenz9params( a )

Calculates the parameters for the expansion to the lorenz equations from `a`.
"""
function lorenz9params( a )

    # calculates the six parameters depending on a
    b1 = 4*(1+a^2)/(1+2*a^2)
    b2 = (1+2*a^2)/(2*(1+a^2))
    b3 = 2*(1-a^2)/(1+a^2)
    b4 = a^2/(1+a^2)
    b5 = 8*a^2/(1+2*a^2)
    b6 = 4/(1+2*a^2)

    [b1,b2,b3,b4,b5,b6]

end
