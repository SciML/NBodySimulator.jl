@testset "Magnetostatics Functional Tests" begin
    @testset "Repelling magnetic dipoles" begin
        d1 = 0.01 # m
        m1 = 5e-6 # kg
        m2 = 5e-6 # kg
        ρ = 7800 # kg/m^3
        M = 1.2e6 # A/m
        mm1 = SVector(0.0, 0.0, M * m1 / ρ)
        mm2 = SVector(0.0, 0.0, M * m2 / ρ)
        μ_4π = 1e-7

        t1 = 0.0  # s
        t2 = 1.0 # s
        τ = (t2 - t1) / 100

        p1 = MagneticParticle(SVector(-d1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1, mm1)
        p2 = MagneticParticle(SVector(d1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2, mm2)
        parameters = MagnetostaticParameters(μ_4π)
        system = PotentialNBodySystem([p1, p2], Dict(:magnetic => parameters))
        simulation = NBodySimulation(system, (t1, t2))
        sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        d2 = norm(get_position(sim_result, t2, 2) - get_position(sim_result, t2, 1))
        v_expected = sqrt(μ_4π * ( dot(mm1, mm2) * (1 / d1^3 - 1 / d2^3) ) / m1)
        v_actual = norm(get_velocity(sim_result, t2, 2))

        ε = 0.001
        @test v_expected ≈ v_actual atol = ε
    end

    @testset "Constructing magnetostatic potential parameters entity" begin
        default_potential = MagnetostaticParameters()
        @test 1e-7 == default_potential.μ_4π

        io = IOBuffer()
        potential1 = MagnetostaticParameters()
        potential2 = MagnetostaticParameters(4π)
        @test sprint(io -> show(io, potential1)) == "Magnetostatic:\n\tμ/4π:1.0e-7\n"
        @test sprint(io -> show(io, potential2)) == "Magnetostatic:\n\tμ/4π:12.566370614359172\n"
    end
end