"""
$(DocStringExtensions.README)
"""
module NBodySimulator

import DocStringExtensions
using CommonSolve: solve
using FileIO: @format_str, File, skipmagic
using LinearAlgebra: cross, dot, norm, normalize
using OrdinaryDiffEqSymplecticRK: VelocityVerlet
using OrdinaryDiffEqTsit5: Tsit5
using PrecompileTools: @compile_workload, @setup_workload
using Printf: @sprintf
using Base: rand, randn
using Random: MersenneTwister
using RecipesBase: @recipe, @series
import RecursiveArrayTools
import SciMLBase
using SciMLBase: AbstractTimeseriesSolution, CallbackSet, DECallback, DiscreteCallback,
    ODEProblem, SDEProblem, SecondOrderODEProblem
using StaticArrays: @SVector
using StaticArraysCore: SVector

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
    generate_bodies_in_cell_nodes, get_accelerating_function, load_water_molecules_from_pdb,
    save_to_pdb

include("precompilation.jl")

end # module
