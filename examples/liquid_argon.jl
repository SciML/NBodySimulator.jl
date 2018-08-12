using NBodySimulator, StaticArrays

function generate_bodies_in_line(n::Int, m::Real, v_dev::Real, L::Real)
    dL = L / (ceil(n^(1 / 3)))
    n_line = floor(Int, L / dL)
    rng = MersenneTwister(n);
    velocities = v_dev * randn(rng, Float64, (3, n_line))
    bodies = MassBody[]
    x = y = L / 2
    for i = 1:n_line      
        r = SVector(x, y, i * dL)
        v = SVector{3}(velocities[:,i])
        body = MassBody(r, v, m)
        push!(bodies, body)  
    end
    return bodies
end

function generate_random_directions(n::Int)
    theta = acos.(1 - 2 * rand(n));
    phi = 2 * pi * rand(n);
    directions = [@SVector [sin(theta[i]) .* cos(phi[i]), sin(theta[i]) .* sin(phi[i]), cos(theta[i])] for i = 1:n]
end

units = :real
units = :reduced

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
#bodies = generate_bodies_randomly(N, m, v_dev, L)
bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)
#bodies = generate_bodies_in_line(N, m, v_dev, L)
jl_parameters = LennardJonesParameters(ϵ, σ, R)
pbc = CubicPeriodicBoundaryConditions(L)
#thermostat = AndersenThermostat(0.02, T, kb)
lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => jl_parameters));
simulation = NBodySimulation(lj_system, (t1, t2), pbc, kb);
#result = run_simulation(simulation, Tsit5())
result = @time run_simulation(simulation, VelocityVerlet(), dt=τ)

#=
using Plots
import GR
(rs, grf) = rdf(result)
(ts, dr2) = msd(result)
plot(rs/σ, grf, xlim=[0, 0.4999L/σ], label=["Radial distribution function"],ylabel="g(r)", xlabel="r/sigma")

using JLD
time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))
#save("D:/water $Nactual molecules $timesteps steps.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2)
save("D:/liquid argon $Nactual molecules $timesteps steps $time_now.jld", "rs", rs, "grf", grf, "ts", ts, "dr2", dr2, "e_tot", e_tot, "e_kin", e_kin, "e_pot", e_pot)
=#
#=
using Plots
import GR
time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))
@time animate(result, "D:/$Nactual liquid argon particles with $timesteps timesteps $time_now.gif")

#plot(simulation)
=#

#=
(rs, grf) = rdf(result)
(ts, dr2) = msd(result)

t = t1:τ:result.solution.t[end-1]
temper = temperature.(result, t)

time_now = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
Nactual = length(bodies)
timesteps = round(length(result.solution.t))

using MAT
file = matopen("d:/liquid argon SI units and $timesteps timesteps_$time_now.mat", "w")
write(file, "t", collect(t))
write(file, "temper", temper)
write(file, "rs", rs)
write(file, "ts", ts)
write(file, "grf", grf)
write(file, "dr2", dr2)
close(file)
=#