"""
Abstract supertype for systems passed to [`NBodySimulation`](@ref).
"""
abstract type NBodySystem end

abstract type BasicPotentialSystem <: NBodySystem end

"""
    ChargedParticles(bodies, k)

System of charged particles interacting through an electrostatic constant `k`.

# Arguments

- `bodies`: Charged particles in the system.
- `k`: Coulomb constant.

# Fields

- `bodies`: Charged particles.
- `k`: Coulomb constant.

# Examples

```julia
using NBodySimulator, StaticArrays

body = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0, 1.0)
system = ChargedParticles([body], 1.0)
```
"""
struct ChargedParticles{bType <: ChargedParticle, kType <: Real} <: BasicPotentialSystem
    bodies::Vector{bType}
    k::kType
end

"""
    GravitationalSystem(bodies, G)

System of massive particles interacting through gravitational constant `G`.

# Arguments

- `bodies`: Massive particles in the system.
- `G`: Gravitational constant.

# Fields

- `bodies`: Massive particles.
- `G`: Gravitational constant.

# Examples

```julia
using NBodySimulator, StaticArrays

body = MassBody(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0)
system = GravitationalSystem([body], 1.0)
```
"""
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
    PotentialNBodySystem(bodies, potentials)
    PotentialNBodySystem(bodies; potentials = Symbol[])

System of bodies governed by a selected set of potential parameters.

# Arguments

- `bodies`: Bodies in the system.
- `potentials`: A `Dict{Symbol, <:PotentialParameters}` or a vector of built-in
  potential names.

# Fields

- `bodies`: Bodies in the system.
- `potentials`: Potential parameters keyed by interaction name.

# Examples

```julia
using NBodySimulator, StaticArrays

bodies = [MassBody(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0)]
system = PotentialNBodySystem(bodies; potentials = [:gravitational])
```
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

"""
    WaterSPCFw(bodies, mH, mO, qH, qO, lj_parameters, e_parameters, scpfw_parameters)

System representation for SPC/Fw water molecules and their interaction parameters.

# Arguments

- `bodies`: Oxygen sites or explicit water molecules.
- `mH`, `mO`: Hydrogen and oxygen masses.
- `qH`, `qO`: Hydrogen and oxygen charges.
- `lj_parameters`, `e_parameters`, `scpfw_parameters`: Interaction parameters.

# Fields

- `bodies`: Water molecule definitions.
- `mH`, `mO`, `qH`, `qO`: Per-site physical constants.
- `lj_parameters`, `e_parameters`, `scpfw_parameters`: Interaction parameters.

# Examples

```julia
using NBodySimulator, StaticArrays

bodies = [MassBody(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 18.0)]
system = WaterSPCFw(
    bodies, 1.0, 16.0, 0.41, -0.82, LennardJonesParameters(), ElectrostaticParameters(),
    SPCFwParameters(0.1012, 1.9764, 44315.0, 317.6)
)
```
"""
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
