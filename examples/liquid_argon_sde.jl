using NBodySimulator
using StochasticDiffEq

T = 120.0 # °K
T0 = 90
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
t2 = 10000τ

parameters = LennardJonesParameters(ϵ, σ, R)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
#thermostat = AndersenThermostat(T0, 1000)
#thermostat = BerendsenThermostat(T0, 20τ)
thermostat = LangevinThermostat(T0, 10.0)
pbc = CubicPeriodicBoundaryConditions(L)
simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
#result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)
result = @time run_simulation_sde(simulation, EM(),  dt=τ)

#=
time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))


using Plots
t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)
pl=plot(t, temper, ylim=[0,500], xlabel="t, ps", ylabel = "T, °K", label="Temperature, °K", linewidth=2)
plot!(pl, t, T0*ones(length(t)), label = "$T0 °K", linewidth=2)
plot!(pl, title="Langevin thermostat at gamma=$(thermostat.γ)  1/ps")


using JLD
save("d:/$(typeof(thermostat)) thermostat for water $T0 _$time_now.jld", "t",t, "temper", temper, "T0", T0, "thermostat", thermostat)
=#