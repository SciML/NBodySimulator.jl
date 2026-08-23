"""
$(DocStringExtensions.README)
"""
module NBodySimulator

import DocStringExtensions
using CommonSolve: solve
using FileIO: @format_str, File, skipmagic
using LinearAlgebra: cross, dot, norm, normalize
using PrecompileTools: @compile_workload, @setup_workload
using Printf: @sprintf
using Base: rand, randn
using Random: MersenneTwister
using RecipesBase: @recipe, @series
import RecursiveArrayTools
import SciMLBase
using SciMLBase: AbstractTimeseriesSolution, DECallback
using StaticArrays: @SVector
using StaticArraysCore: SVector

# ---------------------------------------------------------------------------
# Reexported interface (see the `export` blocks at the bottom of this file).
#
# `using NBodySimulator` has always been enough to write the documented workflow
#
#     simulation = NBodySimulation(system, tspan)
#     result = run_simulation(simulation, VelocityVerlet(), dt = tau)
#
# because the integrator names and the SciML common interface came in through
# `@reexport using DiffEqBase, OrdinaryDiffEq, OrdinaryDiffEqRKN,
# OrdinaryDiffEqSymplecticRK, RecursiveArrayTools`. That blanket reexport is gone;
# these explicit lists keep exactly the surface that normal documented use needs.
# Everything here stays owned and documented upstream.
#
# `run_simulation` always builds a `SecondOrderODEProblem` (see
# `calculate_simulation`), so the algorithms below are the full symplectic and
# Runge-Kutta-Nystrom families -- the two second-order-ODE solver sets this package
# depends on -- plus `Tsit5`, which is `run_simulation`'s default algorithm.
# ---------------------------------------------------------------------------
using OrdinaryDiffEqSymplecticRK: CalvoSanz4, CandyRoz4, KahanLi6, KahanLi8,
    LeapfrogDriftKickDrift, McAte2, McAte3, McAte4, McAte42, McAte5, McAte8,
    PseudoVerletLeapfrog, Ruth3, SofSpa10, SymplecticEuler, VelocityVerlet,
    VerletLeapfrog, Yoshida6
using OrdinaryDiffEqRKN: DPRKN12, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, ERKN4,
    ERKN5, ERKN7, FineRKN4, FineRKN5, IRKN3, IRKN4, Nystrom4,
    Nystrom4VelocityIndependent, Nystrom5VelocityIndependent, RKN4
using OrdinaryDiffEqTsit5: Tsit5
using SciMLBase: CallbackSet, ContinuousCallback, DiscreteCallback, ODEProblem,
    ODESolution, RODESolution, ReturnCode, SDEProblem, SecondOrderODEProblem,
    VectorContinuousCallback, init, remake, solve!, step!, successful_retcode
using RecursiveArrayTools: ArrayPartition

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

# Reexported solver interface; approved via `reexports_allow` in test/qa/qa.jl.
# Second-order-ODE integrators accepted by `run_simulation`.
export CalvoSanz4, CandyRoz4, KahanLi6, KahanLi8, LeapfrogDriftKickDrift, McAte2,
    McAte3, McAte4, McAte42, McAte5, McAte8, PseudoVerletLeapfrog, Ruth3, SofSpa10,
    SymplecticEuler, VelocityVerlet, VerletLeapfrog, Yoshida6
export DPRKN12, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, ERKN4, ERKN5, ERKN7,
    FineRKN4, FineRKN5, IRKN3, IRKN4, Nystrom4, Nystrom4VelocityIndependent,
    Nystrom5VelocityIndependent, RKN4
export Tsit5
# The SciML common interface used to build a problem from an `NBodySimulation`,
# solve it, and inspect the result.
export CallbackSet, ContinuousCallback, DiscreteCallback, ODEProblem, ODESolution,
    RODESolution, ReturnCode, SDEProblem, SecondOrderODEProblem,
    VectorContinuousCallback, init, remake, solve, solve!, step!, successful_retcode
export ArrayPartition

include("precompilation.jl")

end # module
