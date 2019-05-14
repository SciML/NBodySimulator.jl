abstract type NBodySystem
end

abstract type BasicPotentialSystem <: NBodySystem
end

struct ChargedParticles{bType <: ChargedParticle,kType <: Real} <: BasicPotentialSystem
    bodies::Vector{bType}
    k::kType
end

struct GravitationalSystem{bType <: MassBody,gType <: Real} <: BasicPotentialSystem
    bodies::Vector{bType}
    G::gType
end

struct CustomAccelerationSystem{bType <: Body} <: NBodySystem
    bodies::Vector{bType}
    acceleration # f(u, v, i, system)
    parameters::Vector{<:Number}
end

struct PotentialNBodySystem{bType <: Body} <: NBodySystem
    bodies::Vector{bType}
    potentials::Dict{Symbol,<:PotentialParameters}
end

function PotentialNBodySystem(bodies::Vector{<:Body}; potentials::Vector{Symbol}=[])

    paramteres = Dict{Symbol,PotentialParameters}()
    
    if :lennard_jones ∈ potentials
        paramteres[:lennard_jones] = LennardJonesParameters()
    end

    if :electrostatic ∈ potentials 
        paramteres[:electrostatic] = ElectrostaticParameters()
    end

    if :gravitational ∈ potentials   
        paramteres[:gravitational] = GravitationalParameters()
    end

    if :magnetostatic ∈ potentials   
        paramteres[:magnetostatic] = MagnetostaticParameters()
    end

    PotentialNBodySystem(bodies, paramteres)
end

function Base.show(stream::IO, s::PotentialNBodySystem)
    println(stream, "Potentials: ")

    ordered_list = [:lennard_jones, :electrostatic, :magnetostatic,:gravitational]
    for potential in ordered_list
        if potential ∈ keys(s.potentials)
            show(stream, s.potentials[potential])
        end
    end
end

function PotentialNBodySystem(system::PotentialNBodySystem) 
    return system
end

function PotentialNBodySystem(system::ChargedParticles) 
    pp = ElectrostaticParameters(system.k)
    potential = Dict{Symbol,PotentialParameters}(:electrostatic => pp)
    PotentialNBodySystem(system.bodies, potential)
end

function PotentialNBodySystem(system::GravitationalSystem)
    pp = GravitationalParameters(system.G)
    potential =  Dict{Symbol,PotentialParameters}(:gravitational => pp)
    PotentialNBodySystem(system.bodies, potential)
end

struct WaterSPCFw{bType <: Body,pType <: Real} <: NBodySystem
    bodies::Vector{bType}
    mH::pType
    mO::pType
    qH::pType
    qO::pType
    lj_parameters::LennardJonesParameters{pType}
    e_parameters::ElectrostaticParameters{pType}
    scpfw_parameters::SPCFwParameters{pType}
end