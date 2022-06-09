# NBodySimulator

[![Build Status](https://github.com/SciML/NBodySimulator.jl/workflows/CI/badge.svg)](https://github.com/SciML/NBodySimulator.jl/actions?query=workflow%3ACI)

Simulating systems of N interacting bodies.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://nbodysimulator.sciml.ai/stable/). Use the
[in-development documentation](https://nbodysimulator.sciml.ai/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

```julia
using NBodySimulator
using StaticArrays
using Plots
body1 = MassBody(SVector(0.0, 1.0, 0.0), SVector( 5.775e-6, 0.0, 0.0), 2.0)
body2 = MassBody(SVector(0.0,-1.0, 0.0), SVector(-5.775e-6, 0.0, 0.0), 2.0)
G = 6.673e-11
system = GravitationalSystem([body1,body2], G)
tspan = (0.0, 1111150.0)
simulation = NBodySimulation(system, tspan)
sim_result = run_simulation(simulation)
animate(sim_result, "path_to_animated_particles.gif")
```

<img src="https://user-images.githubusercontent.com/16945627/39958539-d2cf779c-561d-11e8-96a8-ffc3a595be8b.gif" alt="Here should appear a gif of rotating bodies" width="350"/>


