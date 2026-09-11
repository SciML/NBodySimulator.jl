using NBodySimulator, BenchmarkTools
using StaticArrays

const SUITE = BenchmarkGroup()

# Two charged particles in orbit
r = 100.0
q1, q2 = 1.0e-3, -1.0e-3
m1, m2 = 100.0, 0.1
k = 9.0e9
v2 = sqrt(abs(k * q1 * q2 / m2 / r))
t_orbit = 2π * r / v2
p1 = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0, 0.0), m1, q1)
p2 = ChargedParticle(SVector(r, 0.0, 0.0), SVector(0.0, v2, 0.0), m2, q2)
charged_sys = ChargedParticles([p1, p2], k)
sim = NBodySimulation(charged_sys, (0.0, 0.1 * t_orbit))

# Small gravitational N-body
function make_grav_sys(n)
    bodies = [
        MassBody(
            SVector(randn(), randn(), randn()),
            SVector(0.01 .* randn(3)...), 1.0
        ) for _ in 1:n
    ]
    return GravitationalSystem(bodies, 6.674e-11 * 1.0e6)
end
grav_sys = make_grav_sys(30)
grav_sim = NBodySimulation(grav_sys, (0.0, 5.0))

# =============================================================================
# run_simulation
# =============================================================================

SUITE["simulate"] = BenchmarkGroup()

SUITE["simulate"]["charged_orbit"] = @benchmarkable run_simulation(
    $sim, VelocityVerlet(); dt = 0.001 * $t_orbit
)
SUITE["simulate"]["gravitational_30"] = @benchmarkable run_simulation(
    $grav_sim, VelocityVerlet(); dt = 0.01
)

# =============================================================================
# System construction
# =============================================================================

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["nbody_simulation"] = @benchmarkable NBodySimulation(
    $grav_sys, (0.0, 5.0)
)
SUITE["construct"]["grav_system"] = @benchmarkable make_grav_sys(30)
