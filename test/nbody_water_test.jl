@testset "Water SPC/Fw test" begin
    qe = 1.6e-19
    Na = 6.022e23
    T = 298.16 # °K
    kb = 1.38e-23 # J/K
    ϵOO = 0.1554253 * 4184 / Na
    σOO = 3.165492e-10 # m
    ρ = 900 # kg/m^3
    mO = 15.999 * 1.6747 * 1e-27 # kg
    mH = 1.00794 * 1.6747 * 1e-27# kg
    mH2O = mO + 2 * mH
    N = 216#floor(Int, ρ * L^3 / m)
    L = (mH2O * N / ρ)^(1 / 3)#10.229σ
    R = 9e-10 # 3*σOO  
    Rel = 0.45 * L
    v_dev = sqrt(kb * T / mH2O)
    τ = 1e-15 # σ/v
    k_bond = 1059.162 * 4184 * 1e20 / Na # J/m^2
    k_angle = 75.90 * 4184 / Na # J/rad^2
    rOH = 1.012e-10 # m
    ∠HOH = 113.24 * pi / 180 # rad
    qH = 0.41 * qe
    qO = -0.84 * qe
    k_el_w = 9e9 #

    @testset "Analyzing simulation result" begin
        N = 216#floor(Int, ρ * L^3 / m)
        L = (mH2O * N / ρ)^(1 / 3)#10.229σ
        t1 = 0τ
        t2 = 10τ 

        r1 = SVector(L / 3, L / 3, 2 * L / 3)
        r2 = SVector(L / 3, 2 * L / 3, L / 3)
        r3 = SVector(2 * L / 3, L / 3, L / 3)
        v1 = SVector(0, 0, -v_dev)
        v2 = SVector(0, -v_dev, 0)
        v3 = SVector(-v_dev, 0, 0)
        p1 = MassBody(r1, v1, mH2O)
        p2 = MassBody(r2, v2, mH2O)
        p3 = MassBody(r3, v3, mH2O)

        bodies = [p1, p2, p3]
        jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
        e_parameters = ElectrostaticParameters(k_el_w, Rel)
        spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
        pbc = CubicPeriodicBoundaryConditions(L)
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
        simulation = NBodySimulation(water, (t1, t2), pbc)

        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        cc = get_position(result, t2)
        ε = 0.01
        for i = 1:3
            indO = 3 * (i - 1) + 1
            @test (rOH - norm(cc[:,indO] - cc[:,indO + 1])) / rOH ≈ 0.0 atol = ε
            @test (rOH - norm(cc[:,indO] - cc[:,indO + 2])) / rOH ≈ 0.0 atol = ε
            ang = acos(dot(cc[:,indO] - cc[:,indO + 1], cc[:,indO] - cc[:,indO + 2]) / (norm(cc[:,indO] - cc[:,indO + 1]) * norm(cc[:,indO] - cc[:,indO + 2])))
            @test (ang - ∠HOH) / ∠HOH ≈ 0.0 atol = ε
        end

        e_tot_1 = total_energy(result, t1)
        e_tot_init = initial_energy(simulation)
        @test e_tot_1 == e_tot_init

        e_tot_2 = total_energy(result, t2)
    #@test (e_tot_1 - e_tot_2)/e_tot_1 ≈ 0.0 atol = ε

        temperature_expected = (dot(v1, v1) + dot(v2, v2) + dot(v3, v3)) * (2 * mH + mO) / (kb * 21)    
        temperature_1 = temperature(result, t1)
        @test temperature_1 == temperature_expected

        count_plotting_data = 0
        for plotting_data ∈ result
            count_plotting_data += 1
            @test plotting_data[2] == result.solution.t[count_plotting_data]
        end
        @test count_plotting_data == length(result.solution.t)

        (ts, mean_square_displacement) = msd(result)
        @test mean_square_displacement[1] < mean_square_displacement[end]
    end

    @testset "Berendsen thermostating" begin
        t1 = 0τ
        t2 = 10τ
        T0 = 275
        r1 = SVector(L / 3, L / 3, 2 * L / 3)
        r2 = SVector(L / 3, 2 * L / 3, L / 3)
        v1 = SVector(0, 0, -v_dev)
        v2 = SVector(0, -v_dev, 0)
        p1 = MassBody(r1, v1, mH2O)
        p2 = MassBody(r2, v2, mH2O)

        bodies = [p1, p2]
        jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
        e_parameters = ElectrostaticParameters(k_el_w, Rel)
        spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
        pbc = CubicPeriodicBoundaryConditions(L)
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
    
        thermostat = BerendsenThermostat(T0, 200τ)
        simulation = NBodySimulation(water, (t1, t2), pbc)
        (ms_ber, kb_ber, n_ber, nc_ber, p_ber) = DiffEqPhysics.obtain_data_for_berendsen_thermostating(simulation)
        @test length(ms_ber) == n_ber
        @test kb_ber == kb
        @test n_ber == 6
        @test nc_ber == 4
    
        thermostat = NoseHooverThermostat(T0, 200τ)
        simulation = NBodySimulation(water, (t1, t2), pbc)
        (ms_nh, kb_nh, n_nh, nc_nh, γind, p_nh) = DiffEqPhysics.obtain_data_for_nosehoover_thermostating(simulation)
        @test length(ms_nh) == n_nh
        @test kb_nh == kb
        @test n_nh == 6
        @test nc_nh == 4
        @test γind == 19
    end

    @testset "Protein Data Bank file compilation" begin
        t1 = 0τ
        t2 = 2τ 

        r1 = SVector(L / 3, L / 3, 2 * L / 3)
        r2 = SVector(L / 3, 2 * L / 3, L / 3)
        v1 = SVector(0, 0, -v_dev)
        v2 = SVector(0, -v_dev, 0)
        p1 = MassBody(r1, v1, mH2O)
        p2 = MassBody(r2, v2, mH2O)

        bodies = [p1, p2]
        jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
        e_parameters = ElectrostaticParameters(k_el_w, Rel)
        spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
        pbc = CubicPeriodicBoundaryConditions(L)
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
        simulation = NBodySimulation(water, (t1, t2), pbc)
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        io = IOBuffer()
        pdb_data = sprint(io -> DiffEqPhysics.write_pdb_data(io, result))
        splitted_data = split(pdb_data, '\n')

        hetatm_count = 0
        timestep_count = 0

        for s in splitted_data
            if length(s) >= 10 && s[1:10] == "REMARK 250"
                timestep_count += 1
            elseif length(s) >= 6 && s[1:6] == "HETATM"
                hetatm_count += 1
            end
        end

        molecule_number = hetatm_count / timestep_count / 3

        @test length(result.solution.t) == timestep_count
        @test molecule_number == length(bodies)
    end
end
