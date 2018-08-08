using NBodySimulator
using StochasticDiffEq

const T = 120.0 # °K
const T0 = 90.0 # °K
const kb = 8.3144598e-3 # kJ/(K*mol)
const ϵ = T * kb
const σ = 0.34 # nm
const ρ = 1374/1.6747# Da/nm^3
const m = 39.95# Da
const N = 8
const L = (m*N/ρ)^(1/3)#10.229σ
const R = 0.5*L   
const v_dev = sqrt(kb * T / m)
const bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

const τ = 0.5e-4 # ps or 1e-12 s
const t1 = 0.0
const t2 = 3τ

const parameters = LennardJonesParameters(ϵ, σ, R)
const lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
#const thermostat = AndersenThermostat(90, 0.01/τ)
#thermostat = BerendsenThermostat(90, 2000τ)
const thermostat = LangevinThermostat(70, 10.0)
#thermostat = NoseHooverThermostat(T0, 200τ)
const pbc = CubicPeriodicBoundaryConditions(L)
const simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
#result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)
#result = @time run_simulation_sde(simulation, ISSEM(symplectic=true,theta=0.5))
result = @time run_simulation(simulation, EM(), dt=τ)

#(rs, grf) = rdf(result)
#(ts, dr2) = msd(result)

t = t1:τ:result.solution.t[end-1]
temper = @time temperature.(result, t)


#using Plots
#pl=plot(t, temper, ylim=[0,200], xlabel="t, ps", ylabel = "T, °K", label="Temperature, °K", linewidth=2)
#plot!(pl, t, T0*ones(length(t)), label = "90 °K", linewidth=2)
#plot!(pl, title="Andersen thermostat at dt*v=$(thermostat.ν*τ) ")
#plot!(pl, title="Berendsen thermostat at tau=$(thermostat.τ) ps")
#plot!(pl, title="Nose-Hoover thermostat at tau=$(thermostat.τ) ps")


time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

#@time save_to_pdb(result, "D:/liquid argon simulation $Nactual molecules and $timesteps steps $time_now.pdb" )

#using JLD
#save("d:/nosehoover thermostat for water liquid argon $T0 and $(thermostat.τ) _$time_now.jld", "t",t, "temper", temper, "T0", T0, "τ", thermostat.τ)


#=
using MAT

time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

file = matopen("d:/liquid argon $(length(bodies)) omm units at temperature $T0 andgamma $(thermostat.γ)  $timesteps timesteps $T0 _$time_now.mat", "w")
write(file, "t", collect(t))
write(file, "temper", temper)
write(file, "T0", T0)
if thermostat isa LangevinThermostat
    write(file, "gamma", thermostat.γ)
end
#write(file, "rs", rs)
#write(file, "ts", ts)
#write(file, "grf", grf)
#write(file, "dr2", dr2)
close(file)
=#