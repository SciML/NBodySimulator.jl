"""
Bodies or Particles are the objects that will interact with each other
and for which the equations of Newton's 2nd law are solved during the simulation process.
"""
abstract type NBodySystem end

abstract type BasicPotentialSystem <: NBodySystem end

struct ChargedParticles{bType <: ChargedParticle, kType <: Real} <: BasicPotentialSystem
    bodies::Vector{bType}
    k::kType
end

struct GravitationalSystem{bType <: MassBody, gType <: Real} <: BasicPotentialSystem
    bodies::Vector{bType}
    G::gType
end

struct CustomAccelerationSystem{bType <: Body} <: NBodySystem
    bodies::Vector{bType}
    acceleration::Any # f(u, v, i, system)
    parameters::Vector{<:Number}
end

struct PotentialNBodySystem{bType <: Body} <: NBodySystem
    bodies::Vector{bType}
    potentials::Dict{Symbol, <:PotentialParameters}
end
"""
Structure that represents systems with a custom set of potentials.
In other words, the user determines the ways in which the particles are allowed to interact.
One can pass the bodies and parameters of interaction potentials into that system.
IF the potential parameters are not set, the particles will move at constant velocities without acceleration during the simulation.
"""
function PotentialNBodySystem(bodies::Vector{<:Body}; potentials::Vector{Symbol} = [])
    parameters = Dict{Symbol, PotentialParameters}()

    if :lennard_jones ∈ potentials
        parameters[:lennard_jones] = LennardJonesParameters()
    end

    if :electrostatic ∈ potentials
        parameters[:electrostatic] = ElectrostaticParameters()
    end

    if :gravitational ∈ potentials
        parameters[:gravitational] = GravitationalParameters()
    end

    if :magnetostatic ∈ potentials
        parameters[:magnetostatic] = MagnetostaticParameters()
    end

    return PotentialNBodySystem(bodies, parameters)
end

function Base.show(stream::IO, s::PotentialNBodySystem)
    println(stream, "Potentials: ")

    ordered_list = [:lennard_jones, :electrostatic, :magnetostatic, :gravitational]
    for potential in ordered_list
        if potential ∈ keys(s.potentials)
            show(stream, s.potentials[potential])
        end
    end
    return
end

function PotentialNBodySystem(system::PotentialNBodySystem)
    return system
end

function PotentialNBodySystem(system::ChargedParticles)
    pp = ElectrostaticParameters(system.k)
    potential = Dict{Symbol, PotentialParameters}(:electrostatic => pp)
    return PotentialNBodySystem(system.bodies, potential)
end

function PotentialNBodySystem(system::GravitationalSystem)
    pp = GravitationalParameters(system.G)
    potential = Dict{Symbol, PotentialParameters}(:gravitational => pp)
    return PotentialNBodySystem(system.bodies, potential)
end

struct WaterSPCFw{bType <: Body, pType <: Real} <: NBodySystem
    bodies::Vector{bType}
    mH::pType
    mO::pType
    qH::pType
    qO::pType
    lj_parameters::LennardJonesParameters{pType}
    e_parameters::ElectrostaticParameters{pType}
    scpfw_parameters::SPCFwParameters{pType}
end
