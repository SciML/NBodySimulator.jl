using Test, NBodySimulator, StaticArrays, OrdinaryDiffEq

@testset "Interface Compatibility" begin
    @testset "BigFloat support" begin
        # Create bodies with BigFloat coordinates
        body1 = MassBody(SVector(big"1.0", big"0.0", big"0.0"),
            SVector(big"0.0", big"0.5", big"0.0"), big"1.0")
        body2 = MassBody(SVector(big"-1.0", big"0.0", big"0.0"),
            SVector(big"0.0", big"-0.5", big"0.0"), big"1.0")

        @test eltype(body1.r) == BigFloat
        @test eltype(body1.v) == BigFloat
        @test typeof(body1.m) == BigFloat

        # Create a gravitational system with BigFloat
        system = GravitationalSystem([body1, body2], big"1.0")
        @test typeof(system.G) == BigFloat

        # Create simulation with BigFloat time span and kb
        simulation = NBodySimulation(system, (big"0.0", big"1.0"), InfiniteBox(),
            big"1.38e-23")
        @test typeof(simulation.tspan[1]) == BigFloat

        # Create ODE problem and verify element types are preserved
        prob = SecondOrderODEProblem(simulation)
        @test eltype(prob.u0) == BigFloat

        # Solve and verify solution maintains BigFloat precision
        sol = solve(prob, Tsit5(), dt = big"0.01", adaptive = false, save_everystep = false)
        @test eltype(sol.u[end]) == BigFloat
    end

    @testset "Float32 support" begin
        # Create bodies with Float32 coordinates
        body1 = MassBody(SVector(1.0f0, 0.0f0, 0.0f0),
            SVector(0.0f0, 0.5f0, 0.0f0), 1.0f0)
        body2 = MassBody(SVector(-1.0f0, 0.0f0, 0.0f0),
            SVector(0.0f0, -0.5f0, 0.0f0), 1.0f0)

        @test eltype(body1.r) == Float32
        @test eltype(body1.v) == Float32
        @test typeof(body1.m) == Float32

        # Create a gravitational system with Float32
        system = GravitationalSystem([body1, body2], 1.0f0)
        @test typeof(system.G) == Float32

        # Create simulation with Float32 time span and kb
        simulation = NBodySimulation(system, (0.0f0, 1.0f0), InfiniteBox(), 1.38f-23)

        # Create ODE problem and verify element types are preserved
        prob = SecondOrderODEProblem(simulation)
        @test eltype(prob.u0) == Float32

        # Solve and verify solution maintains Float32 precision
        sol = solve(prob, Tsit5(), dt = 0.01f0, adaptive = false, save_everystep = false)
        @test eltype(sol.u[end]) == Float32
    end

    @testset "Lennard-Jones BigFloat support" begin
        # Test that Lennard-Jones potential also works with BigFloat
        body1 = MassBody(SVector(big"0.0", big"0.0", big"0.0"),
            SVector(big"0.0", big"0.0", big"0.0"), big"1.0")
        body2 = MassBody(SVector(big"2.0", big"0.0", big"0.0"),
            SVector(big"0.0", big"0.0", big"0.0"), big"1.0")

        lj = LennardJonesParameters(big"1.0", big"1.0", big"2.5")
        potentials = Dict{Symbol, PotentialParameters}(:lennard_jones => lj)
        system = PotentialNBodySystem([body1, body2], potentials)

        simulation = NBodySimulation(system, (big"0.0", big"0.1"), InfiniteBox(),
            big"1.38e-23")
        prob = SecondOrderODEProblem(simulation)

        @test eltype(prob.u0) == BigFloat

        sol = solve(prob, Tsit5(), dt = big"0.001", adaptive = false, save_everystep = false)
        @test eltype(sol.u[end]) == BigFloat
    end
end
