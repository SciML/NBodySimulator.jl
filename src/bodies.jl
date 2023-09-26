# This is an attempt to include the required for n-body simulations fields into structures
# Indeed, tools of ConcreteAbstractions.jl seem to be more suitable for the fields inheriting
DiffEqBase.@def position_velocity_mass begin
    r::SVector{3, cType}
    v::SVector{3, cType}
    m::mType
end
"""
Bodies or Particles are the objects that will interact with each other and for which the equations of Newton's 2nd law are solved during the simulation process.
"""
abstract type Body end

struct MassBody{cType <: Real, mType <: Real} <: Body
    @position_velocity_mass
end

struct ChargedParticle{cType <: Real, mType <: Real, qType <: Real} <: Body
    @position_velocity_mass
    q::qType
end

struct MagneticParticle{cType <: Real, mType <: Real, mmType <: Real} <: Body
    @position_velocity_mass
    mm::SVector{3, mmType}
end

struct WaterMolecule <: Body
    O::MassBody
    H1::MassBody
    H2::MassBody
end
"""
Places similar particles in the nodes of a cubic cell with their velocities distributed in accordance with the Maxwell–Boltzmann law
"""
function generate_bodies_in_cell_nodes(n::Int, m::Real, v_dev::Real, L::Real;
                                       rng = MersenneTwister(n))
    velocities = v_dev * randn(rng, Float64, (3, n))
    bodies = MassBody[]

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x in (dL / 2):dL:L, y in (dL / 2):dL:L, z in (dL / 2):dL:L
        if count > n
            break
        end
        r = SVector(x, y, z)
        v = SVector{3}(velocities[:, count])
        body = MassBody(r, v, m)
        push!(bodies, body)
        count += 1
    end
    return bodies
end
