using NBodySimulator

N = 216
τ = 0.5e-3
T = 300
timesteps = 100

const kb = 8.3144598e-3     # Boltzmann constant, kJ/(K*mol)
const ϵOO = 0.1554253*4.184 # kJ
const σOO = 0.3165492       # nm
const ρ = 997/1.6747        # Da/nm^3
const mO = 15.999           # Da
const mH = 1.00794          # Da
const mH2O = mO+2*mH
const L = (mH2O*N/ρ)^(1/3)
const R = 0.9               # ~3*σOO
const Rel = 0.49*L
const v_dev = sqrt(kb * T /mH2O)
const k_bond = 1059.162*4.184*1e2 # kJ/(mol*nm^2)
const k_angle = 75.90*4.184 # kJ/(mol*rad^2)
const rOH = 0.1012 # nm
const ∠HOH = 113.24*pi/180 # rad
const qH = 0.41
const qO = -0.82
const k = 138.935458 #
#bodies = load_water_molecules_from_pdb("C:/Users/Michael/Desktop/GSoC/pdbs/output_4.pdb")
bodies = generate_bodies_in_cell_nodes(N, mH2O, v_dev, L, T=Float32)
jl_parameters = LennardJonesParameters(Float32(ϵOO), Float32(σOO), Float32(R))
e_parameters = ElectrostaticParameters(Float32(k), Float32(Rel))
spc_paramters = SPCFwParameters(Float32(rOH), Float32(∠HOH), Float32(k_bond), Float32(k_angle))
pbc = CubicPeriodicBoundaryConditions(Float32(L))
water = WaterSPCFw(bodies, Float32(mH), Float32(mO), Float32(qH), Float32(qO),  jl_parameters, e_parameters, spc_paramters)
#thermostat = BerendsenThermostat(T0, 200τ)
#thermostat = NoseHooverThermostat(T0, 200τ)
#thermostat = AndersenThermostat(90, 0.01/τ)
#thermostat = LangevinThermostat(T0, 75)
#simulation = NBodySimulation(water, (t1, t2), pbc, thermostat, kb)
simulation = NBodySimulation(water, (Float32(0τ), Float32(timesteps*τ)), pbc, Float32(kb))
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)
#result = @time run_simulation(simulation, EM(),  dt=τ)
(rs, grf) = @time rdf(result)
(ts, dr2) = @time msd(result)

#@time save_to_pdb(result, "D:/water simulation $Nactual molecules and $timesteps steps $time_now.pdb" )
#using Plots
#import GR
#plot(rs, grf, xlim=[0, 0.4999L], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, nm")
#@time animate(result, "D:/$Nactual H20 particles with $timesteps timesteps $time_now.gif")
#=
t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)
pl=plot(t, temper, ylim=[0,500], xlabel="t, ps", ylabel = "T, °K", label="Temperature, °K", linewidth=2)
plot!(pl, t, T0*ones(length(t)), label = "$T0 °K", linewidth=2, legend = :bottomright)
#plot!(pl, title="Berendsen thermostat at tau=$(thermostat.τ) ps")
plot!(pl, title="Andersen thermostat at dt*v=$(thermostat.ν*τ) ")
#plot!(pl, title="Nose-Hoover thermostat at tau=$(thermostat.τ) ps")
=#
#=
t = t1:τ:result.solution.t[end-1]
temper = @time temperature.(result, t)
time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))
using MAT
file = matopen("d:/water $(length(bodies)) omm units at temperature $T0 and gamma $(thermostat.γ) and $timesteps timesteps with tau $τ _$time_now.mat", "w")
write(file, "t", collect(t))
write(file, "temper", temper)
if thermostat isa LangevinThermostat
    write(file, "gamma", thermostat.γ)
end
#write(file, "rs", rs)
#write(file, "ts", ts)
#write(file, "grf", grf)
#write(file, "dr2", dr2)
close(file)
@time NBodySimulator.save_to_pdb(result, "d:/water from PDB $timesteps timesteps $time_now.pdb" )
=#
