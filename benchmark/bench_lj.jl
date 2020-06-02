using NBodySimulator: obtain_data_for_lennard_jones_interaction,
    pairwise_lennard_jones_acceleration!, gather_bodies_initial_coordinates
using StaticArrays

const T = 120.0 # °K
const T0 = 90.0
const kb = 1.38e-23 # J/K
const ϵ = T * kb
const σ = 3.4e-10 # m
const ρ = 1374 # kg/m^3
const m = 39.95 * 1.6747 * 1e-27 # kg
const N = 216#floor(Int, ρ * L^3 / m)
const L = (m*N/ρ)^(1/3)#10.229σ
const R = 0.5*L
const v_dev = sqrt(kb * T / m)
const τ = 0.5e-15 # σ/v
const t1 = 0τ
const t2 = 30000τ

bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

lj_parameters = LennardJonesParameters(ϵ, σ, R)
pbc = CubicPeriodicBoundaryConditions(L)
#thermostat = AndersenThermostat(0.02, T, kb)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => lj_parameters));
simulation = NBodySimulation(lj_system, (t1, t2), pbc, kb)

ms, indxs = obtain_data_for_lennard_jones_interaction(lj_system)
u0, v0, n = gather_bodies_initial_coordinates(simulation)
dv = zero(v0)
i = 1

b = @benchmarkable pairwise_lennard_jones_acceleration!(dv, $u0, $i, $indxs, $ms,
    $lj_parameters, $simulation.boundary_conditions) setup=(dv=zero($v0)) evals=1

SUITE["lennard_jones"] = b
