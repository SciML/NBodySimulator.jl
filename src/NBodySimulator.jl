__precompile__()

module NBodySimulator

using Reexport
@reexport using DiffEqBase, OrdinaryDiffEq, RecursiveArrayTools
using StaticArrays, RecipesBase, FileIO
using Random, Printf, LinearAlgebra

include("nbody_simulation.jl")

export NBodySimulation
export MassBody, ChargedParticle, MagneticParticle
export PotentialParameters, LennardJonesParameters, GravitationalParameters, 
       ElectrostaticParameters, MagnetostaticParameters, SPCFwParameters
export PotentialNBodySystem, ChargedParticles, GravitationalSystem, WaterSPCFw
export PeriodicBoundaryConditions, CubicPeriodicBoundaryConditions, InfiniteBox
export AndersenThermostat, BerendsenThermostat, NoseHooverThermostat, LangevinThermostat
export run_simulation, get_position, get_velocity, get_masses, temperature,
       initial_energy, kinetic_energy, potential_energy, total_energy, rdf, msd,
       generate_bodies_in_cell_nodes, load_water_molecules_from_pdb

end # module
