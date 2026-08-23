# API

```@docs
NBodySimulator
```

## Simulation

```@docs
NBodySimulation
run_simulation
```

## Bodies

```@docs
NBodySimulator.Body
MassBody
ChargedParticle
MagneticParticle
generate_bodies_in_cell_nodes
```

## Potentials

```@docs
PotentialParameters
NBodySimulator.get_accelerating_function
GravitationalParameters
MagnetostaticParameters
ElectrostaticParameters
LennardJonesParameters
SPCFwParameters
```

## Thermostats

```@docs
NBodySimulator.Thermostat
NBodySimulator.NullThermostat
AndersenThermostat
BerendsenThermostat
NoseHooverThermostat
LangevinThermostat
```

### Boundary Conditions

```@docs
CubicPeriodicBoundaryConditions
PeriodicBoundaryConditions
InfiniteBox
```

## Systems

```@docs
PotentialNBodySystem
ChargedParticles
GravitationalSystem
WaterSPCFw
```

# Analyze

```@docs
NBodySimulator.SimulationResult
get_position
get_velocity
get_masses
temperature
rdf
msd
initial_energy
kinetic_energy
potential_energy
total_energy
```

# Protein Database File

```@docs
load_water_molecules_from_pdb
save_to_pdb
```

# Reexported solver interface

`using NBodySimulator` also brings the solver names the documented workflow needs
into scope, so that

```julia
result = run_simulation(simulation, VelocityVerlet(), dt = τ)
```

works without importing a solver package separately. NBodySimulator only
reexports these names — they are owned and documented upstream, at the links
below.

`run_simulation` always builds a `SecondOrderODEProblem`, so the reexported
algorithms are the second-order-ODE families:

  - Symplectic integrators, owned by
    [OrdinaryDiffEqSymplecticRK](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dynamical_solve/#Symplectic-Integrators):
    `SymplecticEuler`, `VelocityVerlet`, `VerletLeapfrog`,
    `LeapfrogDriftKickDrift`, `PseudoVerletLeapfrog`, `McAte2`, `Ruth3`,
    `McAte3`, `CandyRoz4`, `McAte4`, `CalvoSanz4`, `McAte42`, `McAte5`,
    `Yoshida6`, `KahanLi6`, `McAte8`, `KahanLi8`, `SofSpa10`
  - Runge–Kutta–Nyström integrators, owned by
    [OrdinaryDiffEqRKN](https://docs.sciml.ai/DiffEqDocs/stable/solvers/dynamical_solve/#Runge-Kutta-Nystr%C3%B6m):
    `Nystrom4`, `Nystrom4VelocityIndependent`, `Nystrom5VelocityIndependent`,
    `RKN4`, `IRKN3`, `IRKN4`, `ERKN4`, `ERKN5`, `ERKN7`, `FineRKN4`, `FineRKN5`,
    `DPRKN4`, `DPRKN5`, `DPRKN6`, `DPRKN6FM`, `DPRKN8`, `DPRKN12`
  - The default algorithm, owned by
    [OrdinaryDiffEqTsit5](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/):
    `Tsit5`

Together with the parts of the SciML common interface used to build a problem
from an `NBodySimulation`, solve it, and read the result, owned by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/) and
[CommonSolve](https://docs.sciml.ai/CommonSolve/stable/):

  - Problems: `ODEProblem`, `SecondOrderODEProblem`, `SDEProblem` — each has a
    method taking an `NBodySimulation`
  - Solutions: `ODESolution`, `RODESolution`
  - Solving: `solve`, `solve!`, `init`, `step!`, `remake`
  - Return status: `ReturnCode`, `successful_retcode`
  - Callbacks: `ContinuousCallback`, `DiscreteCallback`,
    `VectorContinuousCallback`, `CallbackSet`

and one name from
[RecursiveArrayTools](https://docs.sciml.ai/RecursiveArrayTools/stable/), the
array type a `SecondOrderODEProblem` solution holds:

  - `ArrayPartition`

Anything else from these packages must be imported from them directly. In
particular, NBodySimulator does not reexport the rest of
[OrdinaryDiffEq](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) —
the stiff, implicit, IMEX and DAE solvers do not apply to a
`SecondOrderODEProblem` — nor the integrator-mutation interface (`add_tstop!`,
`u_modified!`, …), since NBodySimulator implements no integrator of its own. Use
[StochasticDiffEq](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)
directly for the SDE algorithms (e.g. `EM()`) used with the `LangevinThermostat`.
