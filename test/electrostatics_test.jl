@testset "Electrostatics Functional Tests" begin
    k = 9e9

    @testset "One particle rotates around another" begin
    # small mass with negative charge rotating around more massive object with positive charge

        r = 100.0
        q1 = 1e-3
        q2 = -1e-3
        m1 = 100.0
        m2 = 0.1
        v2 = sqrt(abs(k * q1 * q2 / m2 / r))
        t = 2 * pi * r / v2
        p1 = ChargedParticle(SVector(0.0, 0.0, 0.0), SVector(0.0, 0, 0.0), m1, q1)
        p2 = ChargedParticle(SVector(r, 0.0, 0.0), SVector(0.0, v2, 0.0), m2, q2)
        system = ChargedParticles([p1, p2], k)
        simulation = NBodySimulation(system, (0.0, t))
        sim_result = run_simulation(simulation)


        solution = sim_result.solution;
        ε = 0.1 * r
        for j = 1:2, i = 1:3
            @test solution[1][i,j] ≈ solution[end][i,j] atol = ε
        end

        (qs_act, ms_act, indxs_act, exclude_act) = NBodySimulator.obtain_data_for_electrostatic_interaction(simulation.system)
        @test qs_act[1] == q1 && qs_act[2] == q2
        @test ms_act[1] == m1 && ms_act[2] == m2
        @test length(qs_act) == length(ms_act)
    end

    @testset "Two positive charges repelling from each other" begin
    #   ("               <---⊕-----r1-----⊕--->                ")

        q1 = 1e-3 # C
        q2 = 1e-3 # C
        m1 = 1.0 # kg
        m2 = 1.0 # kg
        r1 = 1 # m

        E0 = k * q1 * q2 / r1 # initial energy
        t1 = 0.0
        t2 = 1.0
        τ = (t2 - t1) / 1000

        p1 = ChargedParticle(SVector(-r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1, q1)
        p2 = ChargedParticle(SVector(r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2, q2)
        system = ChargedParticles([p1, p2], k)
        simulation = NBodySimulation(system, (t1, t2))
        sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        r2 = get_position(sim_result, t2, 2) - get_position(sim_result, t2, 1)
        v_expected = sqrt(k * q1 * q2 / m1 * (1 / norm(r1) - 1 / norm(r2)))
        v_actual = norm(get_velocity(sim_result, t2, 2))

        ε = 0.001 * v_expected
        @test v_expected ≈ v_actual atol = ε
    end

    @testset "Constructing electorstatic potential parameters entity" begin
        default_potential_1 = ElectrostaticParameters()
        @test 9e9 == default_potential_1.k

        default_potential_2 = ElectrostaticParameters(k, 8.0)
        @test 64.0 == default_potential_2.R2

        io = IOBuffer()
        potential1 = ElectrostaticParameters()
        potential2 = ElectrostaticParameters(1.0)
        @test sprint(io -> show(io, potential1)) == "Electrostatic:\n\tk:9.0e9\n"
        @test sprint(io -> show(io, potential2)) == "Electrostatic:\n\tk:1.0\n"
    end

    @testset "Energy conservation test" begin
        n = 8
        bodies = ChargedParticle[]
        L = 1.0
        m = 1.0
        q = 1.0
        count = 1
        dL = L / (ceil(n^(1 / 3)) + 1)
        for x = dL / 2:dL:L, y = dL / 2:dL:L, z = dL / 2:dL:L  
            if count > n
                break
            end
            r = SVector(x, y, z)
            v = SVector(.0, .0, .0)
            body = ChargedParticle(r, v, m, q)
            push!(bodies, body)
            count += 1           
        end

        k = 9e9
        τ = 0.01 * dL / sqrt(2 * k * q * q / (dL * m))
        t1 = 0.0
        t2 = 1000 * τ
        
        potential = ElectrostaticParameters(k, 0.45 * L)
        system = PotentialNBodySystem(bodies, Dict(:electrostatic => potential))
        pbc = CubicPeriodicBoundaryConditions(L)
        simulation = NBodySimulation(system, (t1, t2), pbc)
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        e_tot_1 = total_energy(result, t1)
        e_tot_2 = total_energy(result, t2)

        ε = 0.001
        @test (e_tot_2 - e_tot_1) / e_tot_1 ≈ 0.0 atol = ε 
    end
end