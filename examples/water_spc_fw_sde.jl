using NBodySimulator
using StochasticDiffEq

const T = 370 # °K
const T0 = 275 # °K
const kb = 8.3144598e-3 # kJ/(K*mol)
const ϵOO = 0.1554253*4.184 # kJ 
const σOO = 0.3165492 # nm
const ρ = 997/1.6747# Da/nm^3
const mO = 15.999 # Da
const mH = 1.00794 # Da
const mH2O = mO+2*mH
const N = 216#floor(Int, ρ * L^3 / m)
const L = (mH2O*N/ρ)^(1/3)
const R = 0.9 # ~3*σOO  
const Rel = 0.49*L
const v_dev = sqrt(kb * T /mH2O)# sqrt(3.0/7)*sqrt(kb * T /mH2O)
const τ = 0.5e-3 # ps
const t1 = 0τ
const t2 = 100τ # ps
const k_bond = 1059.162*4.184*1e2 # kJ/(mol*nm^2)
const k_angle = 75.90*4.184 # kJ/(mol*rad^2)
const rOH = 0.1012 # nm
const ∠HOH = 113.24*pi/180 # rad
const qH = 0.41
const qO = -0.82
const k = 138.935458 #
bodies = generate_bodies_in_cell_nodes(N, mH2O, v_dev, L)
jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
e_parameters = ElectrostaticParameters(k, Rel)
spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
pbc = CubicPeriodicBoundaryConditions(L)
water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
thermostat = LangevinThermostat(T0, 0)
simulation = NBodySimulation(water, (t1, t2), pbc, thermostat, kb);
result = @time run_simulation(simulation, EM(),  dt=τ)

time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

#(rs, grf) = @time rdf(result)
#(ts, dr2) = @time msd(result)

#using JLD
#save("D:/water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2)
#save("D:/water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2, "e_tot", e_tot, "e_kin", e_kin, "e_pot", e_pot)
#save("D:/!_good_water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2, "τ", τ, "t1", t1, "t2", t2, "N", N, "L", L)
#@time save_to_pdb(result, "D:/water simulation $Nactual molecules and $timesteps steps $time_now.pdb" )

#using Plots
#import GR
#plot(rs, grf, xlim=[0, 0.4999L], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, nm")

#@time animate(result, "D:/$Nactual H20 particles with $timesteps timesteps $time_now.gif")


#=
t = t1:τ:result.solution.t[end-1]
temper = zeros(length(t))
for i in 1:length(t)
    temper[i] = temperature(result, t[i])
end

using MAT
file = matopen("d:/$(typeof(thermostat)) thermostat for water and $timesteps timesteps $T0 _$time_now.mat", "w")
write(file, "t",collect(t))
write(file, "temper", temper)
write(file, "T0", T0)
write(file, "gamma", thermostat.γ)
close(file)
=#

#using JLD
#save("d:/$(typeof(thermostat)) thermostat for water $T0 _$time_now.jld", "t",t, "temper", temper, "T0", T0, "thermostat", thermostat)

#using Plots
#pl=plot(t, temper, ylim=[0,500], xlabel="t, ps", ylabel = "T, °K", label="Temperature, °K", linewidth=2)
#plot!(pl, t, T0*ones(length(t)), label = "$T0 °K", linewidth=2, legend = :bottomright)
#plot!(pl, title="Berendsen thermostat at tau=$(thermostat.τ) ps")
#plot!(pl, title="Andersen thermostat at dt*v=$(thermostat.ν*τ) ")
#plot!(pl, title="Nose-Hoover thermostat at tau=$(thermostat.τ) ps")

#=
t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)

time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

using MAT
file = matopen("d:/water sde $timesteps timesteps $T0 _$time_now.mat", "w")
write(file, "t", collect(t))
write(file, "temper", temper)
write(file, "T0", T0)
close(file)
=#