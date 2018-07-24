using StaticArrays, DiffEqBase, OrdinaryDiffEq, RecipesBase, DiffEqCallbacks

include("./bodies.jl")
include("./boundary_conditions.jl")
include("./thermostats.jl")
include("./basic_potentials.jl")
include("./nbody_system.jl")

const kb_SI = 1.38e-23 # J/K

# This structure defines conditions under wich we test our system of n-bodies
# With this wrapping we can make such fields as `boundary_conditions` necessary for every simulation
# while allowing one to describe a particular systme of N interacting particles
struct NBodySimulation{sType <: NBodySystem,bcType <: BoundaryConditions,tType <: Real,thermType <: Thermostat}
    system::sType
    tspan::Tuple{tType,tType}
    boundary_conditions::bcType
    thermostat::thermType
    kb::tType
    external_electric_field
    external_magnetic_field
    external_gravitational_field
end

function NBodySimulation(system::BasicPotentialSystem,
    tspan::Tuple{tType,tType},
    boundary_conditions::BoundaryConditions,
    thermostat::Thermostat,
    kb::tType,
    external_electric_field,
    external_magnetic_field,
    external_gravitational_field) where {tType <: Real}
    
    potential_system = PotentialNBodySystem(system)
    NBodySimulation(potential_system, tspan, boundary_conditions, thermostat, kb, external_electric_field, external_magnetic_field, external_gravitational_field)
end

function NBodySimulation(system::NBodySystem, tspan::Tuple{tType,tType}, boundary_conditions::BoundaryConditions, kb::tType) where {tType <: Real}
    NBodySimulation(system, tspan, boundary_conditions, NullThermostat(), kb, x -> 0, x -> 0, x -> 0)
end

function NBodySimulation(system::NBodySystem, tspan::Tuple{tType,tType}, boundary_conditions::BoundaryConditions, thermostat::Thermostat, kb::tType) where {tType <: Real}
    NBodySimulation(system, tspan, boundary_conditions, thermostat, kb, x -> 0, x -> 0, x -> 0)
end

function NBodySimulation(system::NBodySystem, tspan::Tuple{tType,tType}, boundary_conditions::BoundaryConditions, thermostat::Thermostat) where {tType <: Real}
    NBodySimulation(system, tspan, boundary_conditions, thermostat, kb_SI, x -> 0, x -> 0, x -> 0)
end

function NBodySimulation(system::NBodySystem, tspan::Tuple{tType,tType}, boundary_conditions::BoundaryConditions) where {tType <: Real}
    NBodySimulation(system, tspan, boundary_conditions, NullThermostat(), kb_SI, x -> 0, x -> 0, x -> 0)
end

function NBodySimulation(system::NBodySystem, tspan::Tuple{tType,tType}) where {tType <: Real}
    NBodySimulation(system, tspan, InfiniteBox(), NullThermostat(), kb_SI, x -> 0, x -> 0, x -> 0)
end

function Base.show(stream::IO, s::NBodySimulation)
    print(stream, "Timespan: ")
    show(stream, s.tspan)
    println(stream)
    print(stream, "Boundary conditions: ")
    show(stream, s.boundary_conditions)
    println(stream)
    show(stream, s.system)
end

include("./nbody_to_ode.jl")
include("./nbody_simulation_result.jl")