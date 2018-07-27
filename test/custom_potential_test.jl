struct CustomPotentialParameters <: PotentialParameters
    a::AbstractFloat
end

import NBodySimulator.get_accelerating_function
function get_accelerating_function(p::CustomPotentialParameters, simulation::NBodySimulation)
    (dv, u, v, t, i) -> begin 
        custom_accel = SVector(p.a, 0.0, 0.0); dv .= custom_accel 
    end 
end

@testset "A sysem with custom potential" begin
    m1 = 1.0
    m2 = 1.0

    t1 = 0.0
    t2 = 1.0
    τ = (t2 - t1) / 100

    p1 = MassBody(SVector(1, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1)
    p2 = MassBody(SVector(-1, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2)


    parameters = CustomPotentialParameters(1.5)
    system = PotentialNBodySystem([p1, p2], Dict(:custom_potential_params => parameters))
    simulation = NBodySimulation(system, (t1, t2))
    simResult = run_simulation(simulation, VelocityVerlet(), dt=τ)

    v2 = get_velocity(simResult, t2, 1)
    r2 = get_position(simResult, t2, 1)

    ε = 1e-6
    @test 1.5 ≈ v2[1] atol = ε
    @test 1.75 ≈ r2[1] atol = ε



    potential_system = PotentialNBodySystem([p1, p2]; potentials=[:lennard_jones, :electrostatic, :gravitational, :magnetostatic])
    @test sprint(io -> show(io, potential_system)) == 
"Potentials: \nLennard-Jones:\n\tϵ:1.0\n\tσ:1.0\n\tR:2.5\nElectrostatic:\n\tk:9.0e9\nMagnetostatic:\n\tμ/4π:1.0e-7\nGravitational:\n\tG:6.67408e-11\n"

end