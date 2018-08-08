using NBodySimulator

const T = 120.0 # °K
const kb = 1.38e-23 # J/K
const ϵ = T * kb
const σ = 3.4e-10 # m
const ρ = 1374 # kg/m^3
const m = 39.95 * 1.6747 * 1e-27 # kg
const N = 216#floor(Int, ρ * L^3 / m)
const L = (m*N/ρ)^(1/3)#10.229σ
const R = 2.25σ   
const v_dev = sqrt(kb * T / m)
const τ = 0.5e-15 # σ/v, fs
const t1 = 0.0
const t2 = 3000τ

const _L = L / σ
const _σ = 1.0
const _ϵ = 1.0
const _m = 1.0
const _v = v_dev / sqrt(ϵ / m)
const _R = R / σ
const _τ = τ * sqrt(ϵ / m) / σ
const _t1 = t1 * sqrt(ϵ / m) / σ
const _t2 = t2 * sqrt(ϵ / m) / σ
bodies = generate_bodies_in_cell_nodes(N, _m, _v, _L)
parameters = LennardJonesParameters(_ϵ, _σ, _R)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
simulation = NBodySimulation(lj_system, (_t1, _t2), CubicPeriodicBoundaryConditions(_L), _ϵ/T);
#result = run_simulation(simulation, Tsit5())
result = @time run_simulation(simulation, VelocityVerlet(), dt=_τ)

#(rs, grf) = rdf(result)
#(ts, dr2) = msd(result)

#using Plots
#plot(rs, grf, xlim=[0, 0.4999_L], label=["Radial distribution function"],ylabel="g(r)", xlabel="r")
#plot(rs/_σ, grf, xlim=[0, 0.4999_L/_σ], label=["Radial distribution function"],ylabel="g(r)", xlabel="r/sigma")


#time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
#Nactual = length(bodies)
#timesteps = round(length(result.solution.t))

#@time save_to_pdb(result, "D:/liquid argon simulation $Nactual molecules and $timesteps steps $time_now.pdb" )