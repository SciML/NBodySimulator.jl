# This is an attempt to include the required for n-body simulations fields into structures
# Indeed, tools of ConcreteAbstractions.jl seem to be more suitable for the fields inhereting
DiffEqBase.@def position_velocity_mass begin
    r::SVector{3,cType}
    v::SVector{3,cType}
    m::mType
end

abstract type Body
end

struct MassBody{cType<:Real, mType<:Real} <: Body
    @position_velocity_mass
end

struct ChargedParticle{cType<:Real, mType<:Real, qType<:Real} <: Body
    @position_velocity_mass
    q::qType
end

struct MagneticParticle{cType<:Real, mType<:Real, mmType<:Real} <: Body
    @position_velocity_mass
    mm::SVector{3,mmType}
end

struct WaterMolecule <: Body
    O::MassBody
    H1::MassBody
    H2::MassBody
end

function generate_bodies_in_cell_nodes(n::Integer, m::Real, v_dev::Real, L::Real; T=Float64)
    rng = MersenneTwister(n)
    velocities = v_dev * randn(rng, T, (3,n))
    bodies = MassBody[]
    count = 1
    dL = T(L / (ceil(n^(1/3))))
    for x = dL/2:dL:L, y = dL/2:dL:L, z = dL/2:dL:L
        if count > n
            break
        end
        r = SVector(x, y, z)
        v = SVector{3,T}(velocities[:,count])
        body = MassBody{T,T}(r, v, m)
        push!(bodies, body)
        count += 1
    end
    return bodies
end
