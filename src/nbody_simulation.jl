include("./bodies.jl")
include("./boundary_conditions.jl")
include("./thermostats.jl")
include("./basic_potentials.jl")
include("./nbody_system.jl")

const kb_SI = 1.38e-23 # J/K

"""
    NBodySimulation(
        system, tspan, boundary_conditions = InfiniteBox(),
        thermostat = NullThermostat(), kb = 1.38e-23
    )

Configuration for an N-body simulation.

# Arguments

- `system`: An `NBodySystem` describing the bodies and interactions.
- `tspan`: Tuple with the simulation start and end times.
- `boundary_conditions`: Boundary-condition model.
- `thermostat`: Thermostat model.
- `kb`: Boltzmann constant in the selected unit system.

# Fields

- `system`, `tspan`: Physical system and integration interval.
- `boundary_conditions`, `thermostat`: Simulation constraints and temperature model.
- `kb`: Boltzmann constant.
- `external_electric_field`, `external_magnetic_field`, `external_gravitational_field`:
  External-field callables.

# Examples

```julia
using NBodySimulator, StaticArrays

body = MassBody(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0)
simulation = NBodySimulation(GravitationalSystem([body], 1.0), (0.0, 1.0))
```
"""
struct NBodySimulation{
        sType <: NBodySystem, bcType <: BoundaryConditions, tType <: Real,
        thermType <: Thermostat,
    }
    system::sType
    tspan::Tuple{tType, tType}
    boundary_conditions::bcType
    thermostat::thermType
    kb::tType
    external_electric_field::Any
    external_magnetic_field::Any
    external_gravitational_field::Any
end

function NBodySimulation(
        system::BasicPotentialSystem,
        tspan::Tuple{tType, tType},
        boundary_conditions::BoundaryConditions,
        thermostat::Thermostat,
        kb::tType,
        external_electric_field,
        external_magnetic_field,
        external_gravitational_field
    ) where {tType <: Real}
    potential_system = PotentialNBodySystem(system)
    return NBodySimulation(
        potential_system, tspan, boundary_conditions, thermostat, kb,
        external_electric_field, external_magnetic_field,
        external_gravitational_field
    )
end

function NBodySimulation(
        system::NBodySystem, tspan::Tuple{tType, tType},
        boundary_conditions::BoundaryConditions,
        kb::tType
    ) where {tType <: Real}
    return NBodySimulation(
        system, tspan, boundary_conditions, NullThermostat(), kb, x -> 0,
        x -> 0, x -> 0
    )
end

function NBodySimulation(
        system::NBodySystem, tspan::Tuple{tType, tType},
        boundary_conditions::BoundaryConditions, thermostat::Thermostat,
        kb::tType
    ) where {tType <: Real}
    return NBodySimulation(
        system, tspan, boundary_conditions, thermostat, kb, x -> 0, x -> 0,
        x -> 0
    )
end

function NBodySimulation(
        system::NBodySystem, tspan::Tuple{tType, tType},
        boundary_conditions::BoundaryConditions,
        thermostat::Thermostat
    ) where {tType <: Real}
    return NBodySimulation(
        system, tspan, boundary_conditions, thermostat, kb_SI, x -> 0, x -> 0,
        x -> 0
    )
end

function NBodySimulation(
        system::NBodySystem, tspan::Tuple{tType, tType},
        boundary_conditions::BoundaryConditions
    ) where {tType <: Real}
    return NBodySimulation(
        system, tspan, boundary_conditions, NullThermostat(), kb_SI, x -> 0,
        x -> 0, x -> 0
    )
end

function NBodySimulation(
        system::NBodySystem,
        tspan::Tuple{tType, tType}
    ) where {tType <: Real}
    return NBodySimulation(
        system, tspan, InfiniteBox(), NullThermostat(), kb_SI, x -> 0, x -> 0,
        x -> 0
    )
end

function Base.show(stream::IO, s::NBodySimulation)
    print(stream, "Timespan: ")
    show(stream, s.tspan)
    println(stream)
    print(stream, "Boundary conditions: ")
    show(stream, s.boundary_conditions)
    println(stream)
    return show(stream, s.system)
end

include("./nbody_to_ode.jl")
include("./nbody_simulation_result.jl")
