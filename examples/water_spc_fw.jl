using NBodySimulator

function generate_bodies_in_line(n::Int, m::Real, v_dev::Real, L::Real)
    dL = L / (ceil(n^(1 / 3)))
    n_line = floor(Int, L / dL)
    rng = MersenneTwister(n);
    velocities = v_dev * randn(rng, Float64, (3, n_line))
    bodies = MassBody[]
    x = y = L / 2
    for i = [3] #1:n_line-1
        r = SVector(x, y, i * dL)
        v = SVector{3}(velocities[:,i])
        body = MassBody(r, v, m)
        push!(bodies, body)  
    end
    return bodies
end

const qe = 1.6e-19
const Na = 6.022e23
const T = 298.16 # °K
const T0 = 275 # °K
const kb = 1.38e-23 # J/K
const ϵOO = 0.1554253*4184/Na
const σOO = 3.165492e-10 # m
const ρ = 216 # kg/m^3
const mO = 15.999 * 1.6747 * 1e-27 # kg
const mH = 1.00794 * 1.6747 * 1e-27#
const mH2O = mO+2*mH
const N = 216#floor(Int, ρ * L^3 / m)
const L = (mH2O*N/ρ)^(1/3)#10.229σ
const R = 9e-10 # 3*σOO  
const Rel = 0.49*L
const v_dev = sqrt(kb * T /mH2O)
const τ = 0.5e-15 # σ/v
const t1 = 0τ
const t2 = 6000τ
const k_bond = 1059.162*4184*1e20/Na # J/m^2
const k_angle = 75.90*4184/Na # J/rad^2
const rOH = 1.012e-10 # m
const ∠HOH = 113.24*pi/180 # rad
const qH = 0.41*qe
const qO = -0.82*qe
const k = 9e9 #
#bodies = generate_bodies_in_cell_nodes(N, mH2O, v_dev, L)
bodies = load_water_molecules_from_pdb("C:/Users/Michael/Desktop/GSoC/pdbs/output_4.pdb")
jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
e_parameters = ElectrostaticParameters(k, Rel)
spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
pbc = CubicPeriodicBoundaryConditions(L)
water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
#thermostat = BerendsenThermostat(T0, 200τ)
simulation = NBodySimulation(water, (t1, t2), pbc, kb);
#result = run_simulation(simulation, Tsit5())
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)


#using JLD
#save("D:/water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2)
#save("D:/water $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2, "e_tot", e_tot, "e_kin", e_kin, "e_pot", e_pot)

#using Plots
#import GR
#plot(rs/1e-10, grf, xlim=[0, 0.4999L/1e-10], label=["Radial distribution function"],ylabel="g(r)", xlabel="r, Å")

#@time animate(result, "D:/$Nactual H20 particles with $timesteps timesteps $time_now.gif")

#=
(rs, grf) = rdf(result)
(ts, dr2) = msd(result)

t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)

time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

using MAT
file = matopen("d:/water real units without thermostat and $timesteps timesteps_$time_now.mat", "w")
write(file, "t", collect(t))
write(file, "temper", temper)
write(file, "rs", rs)
write(file, "ts", ts)
write(file, "grf", grf)
write(file, "dr2", dr2)
close(file)

@time NBodySimulator.save_to_pdb(result, "d:/water in real units from PDB $timesteps timesteps $time_now.pdb" )
=#
