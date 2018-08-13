# NBodySimulator

[![Build Status](https://travis-ci.org/JuliaDiffEq/NBodySimulator.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/NBodySimulator.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/1ofg9ianvcciq26v?svg=true)](https://ci.appveyor.com/project/Mikhail-Vaganov/nbodysimulator-jl)

This project is under development at the moment. The implementation of potential calculations is fairly experimental and has not been extensively verified yet.
You can test simulation of different systems now but be aware of possible changes in future. 

## Add Package

In order to start simulating systems of N interacting bodies, it is necessary to add `NBodySimulator` package to Julia and then begin to use it:

```julia
]add NBodySimulator
using NBodySimulator
```

If you cannot wait to start codding, try to run some scripts from `examples` folder. The number of particles `N` and the final timestep of simulations `t2` are the two parameters which will determine the time of script execution.

## Basic Components
There are three basic components required for any simulation of systems of N-bodies: `bodies`, `system` and `simulation`.

**Bodies** or **Particles** are the objects which will interact with each other and for which the equations of Newton's 2nd law are solved during the simulation process. Three parameters of a body are necessary, they are initial location, initial velocity and mass. `MassBody` structure represents such particles:

```julia
using StaticArrays
r = SVector(.0,.0,.0)
v = SVector(.1,.2,.5)
mass = 1.25
body = MassBody(r,v,mass)
```
For the sake of simulation speed it is advised to use [static arrays](https://github.com/JuliaArrays/StaticArrays.jl).  

A **System** covers bodies and necessary parameters for correct simulation of interaction between particles. For example, to create an entity for a system of gravitationally interacting particles, one needs to use `GravitationalSystem` constructor:

```julia
const G = 6.67e-11 # m^3/kg/s^2
system = GravitationalSystem(bodies, G)
```

**Simulation** is an entity determining parameters of the experiment: time span of simulation, global physical constants, borders of the simulation cell, external magnetic or electric fields, etc. The required arguments for `NBodySImulation` constructor are the system to be tested and the time span of simulation.

```julia
tspan = (.0, 10.0)
simulation = NBodySimulation(system, tspan)
```

There are different types of bodies but they are just containers of particle parameters. The interaction and acceleration of particles are defined by the potentials or force fields.

## Generating bodies
The package exports quite a useful function for placing similar particles in the nodes of a cubic cell with their velocities distributed in accordance with the Maxwell–Boltzmann law:

```julia
N = 100 # number of bodies/particles
m = 1.0 # mass of each of them
v = 10.0 # mean velocity
L = 21.0 # size of the cell side

bodies = generate_bodies_in_cell_nodes(N, m, v, L)
```

Molecules for the SPC/Fw water model can be imported from a PDB file:
```julia
molecules = load_water_molecules_from_pdb("path_to_pdb_file.pdb")
```

## Potentials
The potentials or force field determines the interaction of particles and, therefore, their acceleration.

There are several structures for basic physical interactions:

```julia
g_parameters = GravitationalParameters(G)
m_parameters = MagnetostaticParameters(μ_4π)
el_potential = ElectrostaticParameters(k, cutoff_radius)
jl_parameters = LennardJonesParameters(ϵ, σ, cutoff_radius)
spc_water_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
```

The Lennard-Jones potential is used in molecular dynamics simulations for approximating interactions between neutral atoms or molecules. The [SPC/Fw water model](http://www.sklogwiki.org/SklogWiki/index.php/SPC/Fw_model_of_water) is used in water simulations. The meaning of arguments for `SPCFwParameters` constructor will be clarified further in this documentation.

`PotentialNBodySystem` structure represents systems with a custom set of potentials. In other words, the user determines the ways in which the particles are allowed to interact. One can pass the bodies and parameters of interaction potentials into that system. In case the potential parameters are not set, during the simulation particles will move at constant velocities without acceleration.

```julia
system = PotentialNBodySystem(bodies, Dict(:gravitational => g_parameters, electrostatic: => el_potential))
```

### Custom Potential
There exists an [example](http://docs.juliadiffeq.org/latest/models/physical.html) of simulation of an N-body system at absolutely custom potential. 

Here is shown how to create custom acceleration functions using tools of NBodySimulator.

First of all, it is necessary to create a structure for parameters for the custom potential.

```julia
struct CustomPotentialParameters <: PotentialParameters
    a::AbstractFloat
end
```

Next, the acceleration function for the potential is required. The custom potential defined here creates a force acting on all the particles proportionate to their masses. The first argument of the function determines the potential for which the acceleration should be calculated in this method. 

```julia
import NBodySimulator.get_accelerating_function
function get_accelerating_function(p::CustomPotentialParameters, simulation::NBodySimulation)
    ms = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin custom_accel = SVector(0.0, 0.0, p.a); dv .= custom_accel*ms[i] end 
end
```

After the parameters and acceleration function are created, one can instantiate a system of particles interacting with a set of potentials which includes the just created custom potential:

```julia
parameters = CustomPotentialParameters(-9.8)
system = PotentialNBodySystem(bodies, Dict(:custom_potential_params => parameters))
```

### Gravitational Interaction
Using NBodySimulator it is possible to simulate gravitational interaction of celestial bodies.
In fact, any structure for bodies can be used for simulation of gravitational interaction since all those structures are required to have mass as one of their parameters:

```julia
body1 = MassBody(SVector(0.0, 1.0, 0.0), SVector( 5.775e-6, 0.0, 0.0), 2.0)
body2 = MassBody(SVector(0.0,-1.0, 0.0), SVector(-5.775e-6, 0.0, 0.0), 2.0)
```

Solving gravitational problem one needs to specify the gravitational constant G.
```julia
G = 6.673e-11
```

Now we have enough parameters to create a GravitationalSystem object:

```julia
system = GravitationalSystem([body1,body2], G)
```

Usually we solve an N-body problem for a certain period of time:

```julia
tspan = (0.0, 1111150.0)
```

The created objects determine the simulation we want to run:

```julia
simulation = NBodySimulation(system, tspan)
sim_result = run_simulation(simulation)
```

And, finally, we can animate our solution showing two equal bodies rotating on the same orbit:
```julia
using Plots
animate(sim_result, "path_to_animated_particles.gif")
```

<img src="https://user-images.githubusercontent.com/16945627/39958539-d2cf779c-561d-11e8-96a8-ffc3a595be8b.gif" alt="Here should appear a gif of rotating bodies" width="350"/>

### Electrostatic Interaction
Interaction between charged particles obeys Coulomb's law. The movement of such bodies can be simulated using `ChargedParticle` and `ChargedParticles` structures.

The following example shows how to model two oppositely charged particles. If one body is more massive that another, it will be possible to observe rotation of the light body around the heavy one without adjusting their positions in space. The constructor for `ChargedParticles` system requires bodies and Coulomb's constant `k` to be passed as arguments.

```julia
r = 100.0 # m
q1 = 1e-3 # C
q2 = -1e-3 # C
m1 = 100.0 # kg
m2 = 0.1 # kg
v2 = sqrt(abs(k * q1 * q2 / m2 / r)) # m/s - using the centrifugal acceleration
t = 2 * pi * r / v2 # s  - for one rotation
p1 = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0, 0.0), m1, q1)
p2 = ChargedParticle(SVector(100.0, 0.0, 0.0), SVector(0.0, v2, 0.0), m2, q2)
system = ChargedParticles([p1, p2], k)
simulation = NBodySimulation(system, (0.0, t))
sim_result = run_simulation(simulation)
```

### Magnetic Interaction
An N-body system consisting of `MagneticParticle`s can be used for simulation of interacting magnetic dipoles, though such dipoles cannot rotate in space. Such a model can represent single domain particles interacting under the influence of a strong external magnetic field.

In order to create a magnetic particle, one specifies its location in space, velocity and the vector of its magnetic moment. The following code shows how we can construct an iron particle:

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

To calculate magnetic interactions properly one should also specify the value for the constant μ<sub>0</sub>/4π or its substitute. Having created parameters for the magnetostatic potential, one can now instantiate a system of particles which should interact magnetically. For that purpose we use `PotentialNBodySystem` and pass particles and potential parameters as arguments.

```julia
parameters = MagnetostaticParameters(μ_4π)
system = PotentialNBodySystem([p1, p2], Dict(:magnetic => parameters))
simulation = NBodySimulation(system, (t1, t2))
sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)
```

## Molecular Dynamics (MD)
NBodySimulator allows one to conduct molecular dynamic simulations for the Lennard-Jones liquids, SPC/Fw model of water and other molecular systems thanks to implementations of basic interaction potentials between atoms and molecules: 

- Lennard-Jones
- electrostatic and magnetostatic
- harmonic bonds
- harmonic valence angle generated by pairs of bonds 

The comprehensive examples of liquid argon and water simulations can be found in `examples` folder.
Here only the basic principles of the molecular dynamics simulations using NBodySimulator are presented using liquid argon as a classical MD system for beginners.

First of all, one needs to define parameters of the simulation:

```julia
T = 120.0 # °K
T0 = 90.0 # °K
kb = 8.3144598e-3 # kJ/(K*mol)
ϵ = T * kb
σ = 0.34 # nm
ρ = 1374/1.6747# Da/nm^3
m = 39.95# Da
N = 216
L = (m*N/ρ)^(1/3)#10.229σ
R = 0.5*L   
v_dev = sqrt(kb * T / m)
bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

τ = 0.5e-3 # ps or 1e-12 s
t1 = 0.0
t2 = 2000τ
```

Liquid argon consists of neutral molecules so the Lennard-Jones potential runs their interaction:

```julia
parameters = LennardJonesParameters(ϵ, σ, R)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters)); 
```

Then, a thermostat and boundary conditions should be selected and instantiated:

```julia
thermostat = NoseHooverThermostat(T0, 200τ)
pbc = CubicPeriodicBoundaryConditions(L)
simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
result = run_simulation(simulation, VelocityVerlet(), dt=τ)
```
It is recommended to use `CubicPeriodicBoundaryConditions` since cubic boxes are among the most popular boundary conditions in MD. There are different variants of the `NBodySimulation` constructor for MD:

```julia
simulation = NBodySimulation(lj_system, (t1, t2));
simulation = NBodySimulation(lj_system, (t1, t2), pbc);
simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat);
simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
```

The default boundary conditions are `InfiniteBox` without any limits, default thermostat is `NullThermostat` which does no thermostating and default Boltzmann constant `kb` equals its value in SI, i.e. 1.38e-23 J/K.

## Water Simulations
In NBodySImulator the [SPC/Fw water model](http://www.sklogwiki.org/SklogWiki/index.php/SPC/Fw_model_of_water) is implemented. For using this model, one has to specify parameters of the Lennard-Jones potential between the oxygen atoms of water molecules, parameters of the electrostatic potential for the corresponding interactions between atoms of different molecules and parameters for harmonic potentials representing bonds between atoms and the valence angle made from bonds between hydrogen atoms and the oxygen one.

```julia
bodies = generate_bodies_in_cell_nodes(N, mH2O, v, L)
jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
e_parameters = ElectrostaticParameters(k, Rel)
spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
```

For each water molecule here, `rOH` is the equilibrium distance between a hydrogen atom and the oxygen atom, `∠HOH` denotes the equilibrium angle made of those two bonds, `k_bond` and `k_angle` are the elastic coefficients for the corresponding harmonic potentials.

Further, one pass the water system into `NBodySimulation` constructor as a usual system of N-bodies.

```julia
simulation = NBodySimulation(water, (t1, t2), pbc, kb);
```

## Thermostats
Usually during simulation of a system is required to be at a particular temperature. NBodySimulator contains several thermostats for that purpose. Here the thermostating of liquid argon is presented, for thermostating of water one can refer to [this post](https://mikhail-vaganov.github.io/gsoc-2018-blog/2018/08/06/thermostating.html)

### [Andersen Thermostat](http://www.sklogwiki.org/SklogWiki/index.php/Andersen_thermostat)
```julia
τ = 0.5e-3 # timestep of integration and simulation
T0 = 90
ν = 0.05/τ
thermostat = AndersenThermostat(90, ν)
```
![andersen thermostating](https://user-images.githubusercontent.com/16945627/44002487-cd99653a-9e5c-11e8-8481-78945a930d94.png)

### [Berendsen Thermostat](http://www2.mpip-mainz.mpg.de/~andrienk/journal_club/thermostats.pdf)
```julia
τB = 2000τ
thermostat = BerendsenThermostat(90, τB)
```
![berendsen thermostating](https://user-images.githubusercontent.com/16945627/44002495-f07e164a-9e5c-11e8-8db7-16c09a7631cd.png)

### [Nosé–Hoover Thermostat](http://www.sklogwiki.org/SklogWiki/index.php/Nos%C3%A9-Hoover_thermostat)
```julia
τNH = 200τ
thermostat = NoseHooverThermostat(T0, 200τ)
```
![nose-hoover thermostating](https://user-images.githubusercontent.com/16945627/44002501-ffc1aea0-9e5c-11e8-857b-9e3c83197336.png)

### Langevin Thermostat
```julia
γ = 10.0
thermostat = LangevinThermostat(90, γ)
```
![langevin thermostating](https://user-images.githubusercontent.com/16945627/44002505-0683c6b0-9e5d-11e8-8647-5b15b98eb0fa.png)

## Analyzing the Result of Simulation
Once the simulation is completed, one can analyze the result and obtain some useful characteristics of the system. 

Function `run_simulation` returns a structure containing the initial parameters of simulation and the solution of differential equation (DE) required for description of the corresponding system of particles. There are different functions which help to interpret solution of DEs into physical quantities.

One of the main characteristics of a system during molecular dynamics simulations is its thermodynamic temperature. The value of the temperature at a particular time `t` can be obtained via calling this function:

```julia
T = temperature(result, t) 
```

### [Radial distribution functions](https://en.wikipedia.org/wiki/Radial_distribution_function) 
The RDF is another popular and essential characteristic of molecules or similar systems of particles. It shows the reciprocal location of particles averaged by the time of simulation.

```julia
(rs, grf) = rdf(result)
```

The dependence of `grf` on `rs` shows radial distribution of particles at different distances from an average particle in a system.
Here the radial distribution function for the classic system of liquid argon is presented:
![rdf for liquid argon](https://user-images.githubusercontent.com/16945627/43990348-843b164c-9d74-11e8-8d9e-daaff142c0b7.png)


### Mean Squared Displacement (MSD)
The MSD characteristic can be used to estimate the shift of particles from their initial positions.
```julia
(ts, dr2) = msd(result)
```
For a standard liquid argon system the displacement grows with time:
![rdf for liquid argon](https://user-images.githubusercontent.com/16945627/43990362-9a67c0aa-9d74-11e8-9512-08840294d411.png)

### Energy Functions

Energy is a highly important physical characteristic of a system. The module provides four functions to obtain it, though the `total_energy` function just sums potential and kinetic energy:

```julia
e_init = initial_energy(simualtion)
e_kin = kinetic_energy(result, t)
e_pot = potential_energy(result, t)
e_tot = total_energy(result, t)
```

## Plotting Images
Using tools of NBodySimulator one can export results of simulation into a [Protein Database File](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)). [VMD](http://www.ks.uiuc.edu/Research/vmd/) is a well-known tool for visualizing molecular dynamics, which can read data from PDB files. 

```julia
save_to_pdb(result, "path_to_a_new_pdb_file.pdb" )
```

In future it will be possible to export results via FileIO interface and its `save` function.

Using Plots.jl one can draw positions of particles at any time of simulation or create an animation of moving particles, molecules of water:

```julia
using Plots
plot(result)
animate(result, "path_to_file.gif")
```

Makie.jl also has a recipe for plotting results of N-body simulations. The [example](http://makie.juliaplots.org/stable/examples-meshscatter.html#Type-recipe-for-molecule-simulation-1) is presented in the documentation.