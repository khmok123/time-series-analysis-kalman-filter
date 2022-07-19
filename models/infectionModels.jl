function SEIRS!(du, u, p, t = 0 )
"""
    SEIRS!( du, u, p, t = 0 )

    writes the differential solution for the SEIRS-model into `du`.
"""

    β,ν,μₗ,α = p
    S,E,I,R = u

    du[1] = -β*S*I
    du[2] = β*S*I-ν*E
    du[3] = ν*E-(1-α)*I
    du[4] = (1-μₗ-α)*I

end


function SEIRS( u, p, t = 0 )
"""
    SEIRS( u, p, t = 0 )

    returns the differential solution for the SEIRS-model.
"""

    β,ν,μₗ,α = p
    S,E,I,R = u

    Ṡ = -β*S*I
    Ė = β*S*I-ν*E
    İ = ν*E-(1-α)*I
    Ṙ = (1-μₗ-α)*I

    [Ṡ, Ė, İ, Ṙ]

end


function SEIR!( du, u, p, t = 0 )
"""
    SEIR!( du, u, p, t = 0 )

    writes the differential soultion for the SEIR-model to `du`.
"""

   S,E,I,R = u
   β,σ,γ = p
   N = S + E + I + R

   du[1] = - β*S*I/N
   du[2] =   β*S*I/N - σ*E
   du[3] =   σ*E     - γ*I
   du[4] =   γ*I

end


function SEIR( u, p, t = 0 )
"""
    SEIR( u, p, t = 0 )

    returns the differential solution for the SEIR-model.
"""

   S,E,I,R = u
   β,σ,γ = p
   N = S + E + I + R

   Ṡ = - S*β(t)*E/N
   Ė =   β(t)*S*E/N   - σ*E
   İ =   σ*E       - γ*I
   Ṙ =   γ*I

   [Ṡ, Ė, İ, Ṙ]

end


function SCIRS!( du, u, p, t = 0 )
"""
    SCIRS!( du, u, p, t = 0 )

    writes the differential soultion for the SCIRS-model to `du`.
"""

    S,C,I,R = u
    β, σ, α, ρ, γ = p
    N = S + C + I + R

    du[1] = - β*S*( C+I )/N       + ρ*R
    du[2] =   β*S*( C+I )/N - γ*C - σ*C
    du[3] = - α*I                 + σ*C
    du[4] =   α*I           + γ*C - ρ*R

end


function SCIRS( u, p, t = 0 )
"""
    SCIRS( du, u, p, t = 0 )

    returns the differential soultion for the SCIRS-model.
"""

    S,C,I,R = u
    β, σ, α, ρ, γ = p
    N = S + C + I + R

    Ṡ = - β*S*( C)/N       + ρ*R
    Ċ =   β*S*( C)/N - γ*C - σ*C
    İ = - α*I              + σ*C
    Ṙ =   α*I        + γ*C - ρ*R

    [Ṡ, Ċ, İ, Ṙ]

end

function σ_SEIR!( du, u, p, t = 0 )

"""
    σ_SEIR!( du, u, p, t = 0 )

    writes the `σ` which is given by the 4th element of params into `du`.
"""

    σ = p[4]
    du .= σ

end
