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

To calculate magnetic interactions properly one should also specify the value for the vacuum permeability μ<sub>0</sub>/10<sup>7</sup>. or its substitute:

```julia
parameters = MagnetostaticParameters(μ_4π)
system = PotentialNBodySystem([p1, p2], Dict(:magnetic => parameters))
simulation = NBodySimulation(system, (t1, t2))
sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)
```