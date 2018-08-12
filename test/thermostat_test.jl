@testset "Testing thermostats on liquid argon " begin 
    T = 120.0 # °K
    T0 = 90 # °K
    kb = 8.3144598e-3 # kJ/(K*mol)
    ϵ = T * kb
    σ = 0.34 # nm
    ρ = 1374/1.6747# Da/nm^3
    m = 39.95# Da
    N = 125
    L = (m * N / ρ)^(1 / 3)#10.229σ
    R = 0.5 * L   
    v_dev = sqrt(kb * T / m)
    @testset "Andersen thermostat" begin 
        bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

        τ = 0.5e-3 # ps or 1e-12 s
        t1 = 0.0
        t2 = 200τ

        parameters = LennardJonesParameters(ϵ, σ, R)
        lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
        thermostat = AndersenThermostat(T0, 0.1/τ)
        simulation = NBodySimulation(lj_system, (t1, t2), PeriodicBoundaryConditions(L), thermostat, kb);
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)
    

        T1 = temperature(result, t1) 
        T2 = temperature(result, t2)
        ε = 0.5
        @test abs(T2 - T0) / T0 ≈ 0.0 atol = ε
    end

    @testset "Berendsen thermostat" begin 
        bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

        τ = 0.5e-3
        t1 = 0.0
        t2 = 200τ

        parameters = LennardJonesParameters(ϵ, σ, R)
        lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
        thermostat = BerendsenThermostat(T0, 10τ)
        pbc = CubicPeriodicBoundaryConditions(L)
        simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)
     
        T2 = temperature(result, t2)
        ε = 0.1
        @test abs(T2 - T0) / T ≈ 0.0 atol = ε
    end

    @testset "Nose-Hoover thermostat" begin 
        bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

        τ = 0.5e-3
        t1 = 0.0
        t2 = 200τ

        parameters = LennardJonesParameters(ϵ, σ, R)
        lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
        thermostat = NoseHooverThermostat(T0, 20τ)
        pbc = CubicPeriodicBoundaryConditions(L)
        simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)
     
        T2 = temperature(result, t2)
        ε = 0.5
        @test abs(T2 - T0) / T0 ≈ 0.0 atol = ε
    end

    @testset "Langevin thermostat" begin 
        bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

        τ = 0.5e-3
        t1 = 0.0
        t2 = 200τ

        parameters = LennardJonesParameters(ϵ, σ, R)
        lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
        thermostat = LangevinThermostat(T0, 10)
        pbc = CubicPeriodicBoundaryConditions(L)
        simulation = NBodySimulation(lj_system, (t1, t2), pbc, thermostat, kb);
        result = run_simulation(simulation, EM(),  dt=τ)
     
        T2 = temperature(result, t2)
        ε = 0.5
        @test abs(T2 - T0) / T0 ≈ 0.0 atol = ε
    end
end