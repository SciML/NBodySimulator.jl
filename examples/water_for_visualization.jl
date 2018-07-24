using NBodySimulator, Dates

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
const v_dev = sqrt(kb * T /mH2O)
const τ = 0.5e-3 # ps
const t1 = 0τ
const t2 = 5τ # ps
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
simulation = NBodySimulation(water, (t1, t2), pbc, kb);
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)

#=
time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

ts = t1:τ:result.solution.t[end-1]

for t in ts
    coordinates = get_position(result, t)
end
=#