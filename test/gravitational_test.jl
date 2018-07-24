using OrdinaryDiffEq 

@testset "Gravitational Functional Test" begin
    G = 1
    @testset "The well-known figure Eight choreography" begin

        m1 = MassBody(SVector(-0.995492, 0.0, 0.0), SVector(-0.347902, -0.53393, 0.0), 1.0)
        m2 = MassBody(SVector(0.995492, 0.0, 0.0), SVector(-0.347902, -0.53393, 0.0), 1.0)
        m3 = MassBody(SVector(0.0, 0.0, 0.0), SVector(0.695804, 1.067860, 0.0), 1.0)
        t1 = 0.0
        t2 = 2pi
        tspan = (t1, t2);
        system = GravitationalSystem([m1, m2, m3], G)
        simulation = NBodySimulation(system, tspan)
        sim_result = run_simulation(simulation)
        solution_simo_3 = sim_result.solution;
        ε = 0.1
        for j = 1:3, i = 1:3
            @test solution_simo_3[1][i,j] ≈ solution_simo_3[end][i,j] atol = ε
        end

        @testset "Analyzing simulation result" begin
            r1 = get_position(sim_result, tspan[1], 1)
            r1_by_index = get_position(sim_result, tspan[1], 1)
            v1 = get_velocity(sim_result, tspan[1], 1)
            masses = get_masses(system)
            for i = 1:3
                @test r1[i] == (-0.995492, 0.0, 0.0)[i]
                @test r1_by_index[i] == (-0.995492, 0.0, 0.0)[i]
                @test v1[i] == (-0.347902, -0.53393, 0.0)[i]
                @test masses[i] == 1.0
            end

            ε = 0.001
            e_kin = 1.218
            @test e_kin ≈ kinetic_energy(sim_result, t1) atol = ε
        end
    

        @testset "Using convertion into SecondOrderODEProblem" begin
            sim_result = run_simulation(simulation, DPRKN6())
            solution_simo_3_2nd = sim_result.solution;
            ε = 0.001
            for i = 1:3, j = 1:3
                @test solution_simo_3_2nd[1][9 + 3(i - 1) + j] ≈ solution_simo_3_2nd[end][9 + 3(i - 1) + j] atol = ε
            end
        end

        @testset "Using symplectic integrators" begin
            sim_result = run_simulation(simulation, VelocityVerlet(), dt=pi / 130)
            solution_simo_3_2nd = sim_result.solution;
            ε = 0.001
            for i = 1:3, j = 1:3
                @test solution_simo_3_2nd[1][9 + 3(i - 1) + j] ≈ solution_simo_3_2nd[end][9 + 3(i - 1) + j] atol = ε
            end

            sim_result = run_simulation(simulation, Yoshida6(), dt=pi / 12)
            solution_simo_3_2nd = sim_result.solution;
            ε = 0.001
            for i = 1:3, j = 1:3
                @test solution_simo_3_2nd[1][9 + 3(i - 1) + j] ≈ solution_simo_3_2nd[end][9 + 3(i - 1) + j] atol = ε
            end
        end
    end

    @testset "5-body choreography" begin
        m1 = MassBody(SVector(1.657666, 0.0, 0.0), SVector(0.0, -0.593786, 0.0), 1.0)
        m2 = MassBody(SVector(0.439775, -0.169717, 0.0), SVector(1.822785, 0.128248, 0.0), 1.0)
        m3 = MassBody(SVector(-1.268608, -0.267651, 0.0), SVector(1.271564, 0.168645, 0.0), 1.0)
        m4 = MassBody(SVector(-1.268608, 0.267651, 0.0), SVector(-1.271564, 0.168645, 0.0), 1.0)
        m5 = MassBody(SVector(0.439775, 0.169717, 0.0), SVector(-1.822785, 0.128248, 0.0), 1.0)
        tspan = (0.0, 2pi);
        system = GravitationalSystem([m1, m2, m3, m4, m5], G)
        simulation = NBodySimulation(system, tspan)
        sim_result = run_simulation(simulation, Tsit5(), abstol=1e-10, reltol=1e-10)
        solution_simo_5 = sim_result.solution;

        ε = 0.01
        for j = 1:5, i = 1:3
            @test solution_simo_5[1][i,j] ≈ solution_simo_5[end][i,j] atol = ε
        end
    end

    @testset "Constructing electorstatic potential parameters entity" begin
        default_potential = GravitationalParameters()
        @test 6.67408e-11 == default_potential.G
    
        io = IOBuffer()
        potential1 = GravitationalParameters()
        potential2 = GravitationalParameters(35.67)
        @test sprint(io -> show(io, potential1)) == "Gravitational:\n\tG:6.67408e-11\n"
        @test sprint(io -> show(io, potential2)) == "Gravitational:\n\tG:35.67\n"
    end
end