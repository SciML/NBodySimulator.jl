# API

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
GravitationalSystem
WaterSPCFw
```

# Analyze
```@docs
NBodySimulator.SimulationResult
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
