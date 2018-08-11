# NBodySimulator

[![Build Status](https://travis-ci.org/JuliaDiffEq/NBodySimulator.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/NBodySimulator.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/1ofg9ianvcciq26v?svg=true)](https://ci.appveyor.com/project/Mikhail-Vaganov/nbodysimulator-jl)

This project is under development at the moment. The implementation of potential calculations is fairly experimental and has not been extensively verified yet.

You can test simulation of different systems now but be aware of possible changes in future. 

## Gravitational interaction

## Electrostatic interaction

## Magnetic interaction
An n-body system consisting of `MagneticParticle`s can be used for simulation of interacting magnetic dipoles, thoug such dipoles cannot rotate in space. Such a model can represent single domain particles interacting under the influence of a strong external magnetic field.

In order to create a magnetic particle one specifies its location in space, velocity and the vector of its magnetic moment. The following code shows how we can construct an iron particle:

```julia
iron_dencity = 7800 # kg/m^3
magnetization_saturation = 1.2e6 # A/m

mass =  5e-6 # kg
r = SVector(-0.005,0.0,0.0) # m
v = SVector(0.0,0.0,0.0) # m/s
magnetic_moment = SVector(0.0, 0.0, magnetization_saturation * mass / iron_dencity) # A*m^2

p1 = MagneticParticle(r, v, mass, magnetic_moment)
```

For the second particle we will use a shorter form:

```julia
p2 = MagneticParticle(SVector(0.005, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 5e-6, SVector(0.0,0.0,0.00077))
```

To calculate magnetic interactions properly one should also specify the value for the constant μ<sub>0</sub>/4π or its substitute. Having created parameters for the magnetostatic potential, one now can instantiate a system of particles which should interact magnetically. For that pupose we use `PotentialNBodySystem` and pass particles and potential parameters as arguments.

```julia
parameters = MagnetostaticParameters(μ_4π)
system = PotentialNBodySystem([p1, p2], Dict(:magnetic => parameters))
simulation = NBodySimulation(system, (t1, t2))
sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)
```

## Anylizing the Result of Simulation
Once the simulation is completed, one can analyze its result and obtain some useful characteristics of the system. 

Function `run_simulation` returns a structure containig the initial parameters of simulation and the solution of differential equation required for description of the corresponding system of particles. There are different functions which help to intepret solution of DEs into physical quantities.

One of the main charachterestics of a system during molecular dynamics simulations is its thermodynamic temperature. The value of the temperature at a particular time `t` can be obtained via calling this function:

```julia
T = temperature(result, t) 
```

[Radial distribution functions](https://en.wikipedia.org/wiki/Radial_distribution_function) is another popular and essential characteristic of molecules or similar system of particles. It shows the reciprocal location of particles averaged by the time of simulation.

```julia
(rs, grf) = @time rdf(result)
```

The dependence of `grf` on `rs` shows radial distribution of particles at different distances from an average particle in a system.
Here is the radial distribution function for the classic system of liquid argon:
![rdf for liquid argon][https://user-images.githubusercontent.com/16945627/43990348-843b164c-9d74-11e8-8d9e-daaff142c0b7.png]


```julia
(ts, dr2) = @time msd(result)
```
![rdf for liquid argon][https://user-images.githubusercontent.com/16945627/43990362-9a67c0aa-9d74-11e8-9512-08840294d411.png]
