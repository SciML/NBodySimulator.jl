@testset "Water SPC/Fw test" begin
    qp = 1 # charge of a proton
    T = 298.16 # °K
    kb = 8.3144598e-3 # kJ/(K*mol)
    ϵOO = 0.1554253*4.184 # kJ 
    σOO = 0.3165492 # nm
    ρ = 997/1.6747# Da/nm^3
    mO = 15.999 # Da
    mH = 1.00794 # Da
    mH2O = mO+2*mH
    N = 216
    L = (mH2O*N/ρ)^(1/3) # nm
    R = 0.9 # ≈3*σOO  
    Rel = 0.49*L
    v_dev = sqrt(kb * T / mH2O) # nm/ps
    k_bond = 1059.162*4.184*1e2 # kJ/(mol*nm^2)
    k_angle = 75.90*4.184 # kJ/(mol*rad^2)
    rOH = 0.1012 # nm
    ∠HOH = 113.24 * pi / 180 # rad
    qH = 0.41 * qp
    qO = -0.84 * qp
    k_el = 138.935458 #
    τ = 0.5e-3 # ps

    jl_parameters = LennardJonesParameters(ϵOO, σOO, R)
    e_parameters = ElectrostaticParameters(k_el, Rel)
    spc_paramters = SPCFwParameters(rOH, ∠HOH, k_bond, k_angle)
    pbc = CubicPeriodicBoundaryConditions(L)

    @testset "Analyzing simulation result" begin
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
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
        simulation = NBodySimulation(water, (t1, t2), pbc, kb)

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
        @test (e_tot_1 - e_tot_2)/e_tot_1 ≈ 0.0 atol = ε

        temperature_expected = (dot(v1, v1) + dot(v2, v2) + dot(v3, v3)) * (2 * mH + mO) / (kb * 21)    
        temperature_1 = temperature(result, t1)
        @test (temperature_1 - temperature_expected)/temperature_expected ≈ 0.0 atol = ε

        count_plotting_data = 0
        for plotting_data ∈ result
            count_plotting_data += 1
            @test plotting_data[2] == result.solution.t[count_plotting_data]
        end
        @test count_plotting_data == length(result.solution.t)

        (ts, mean_square_displacement) = msd(result)
        @test mean_square_displacement[1] < mean_square_displacement[end]
    end

    @testset "Preparing data for thermostating" begin
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
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
    
        thermostat = BerendsenThermostat(T0, 200τ)
        simulation = NBodySimulation(water, (t1, t2), pbc, kb)
        (ms_ber, kb_ber, n_ber, nc_ber, p_ber) = NBodySimulator.obtain_data_for_berendsen_thermostating(simulation)
        @test length(ms_ber) == n_ber
        @test kb_ber == kb
        @test n_ber == 6
        @test nc_ber == 4
    
        thermostat = NoseHooverThermostat(T0, 200τ)
        simulation = NBodySimulation(water, (t1, t2), pbc, kb)
        (ms_nh, kb_nh, n_nh, nc_nh, γind, p_nh) = NBodySimulator.obtain_data_for_nosehoover_thermostating(simulation)
        @test length(ms_nh) == n_nh
        @test kb_nh == kb
        @test n_nh == 6
        @test nc_nh == 4
        @test γind == 19
    end

    @testset "Water thermostating" begin
        ττ = 0.5e-4 # ps
        t1 = 0ττ
        t2 = 200ττ
        T0 = 275
        bodies = generate_bodies_in_cell_nodes(N, mH2O, v_dev, L)
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters);
        thermostat = LangevinThermostat(T0, 100)
        simulation = NBodySimulation(water, (t1, t2), pbc, thermostat, kb)
        result = run_simulation(simulation, EM(),  dt=τ)
        
        T2 = temperature(result, t2)
        ε = 1.0
        @test abs(T2 - T0) / T0 ≈ 0.0 atol = ε
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
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
        simulation = NBodySimulation(water, (t1, t2), pbc, kb)
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        io = IOBuffer()
        pdb_data = sprint(io -> NBodySimulator.write_pdb_data(io, result))
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

    @testset "Load water molecules from a PDB file" begin

        pdb_data = IOBuffer("MODEL     1
REMARK 250 time=0.0000 picoseconds
HETATM    1  O   HOH     1      18.642  18.642   0.000  1.00  0.00
HETATM    2  H1  HOH     1      17.685  18.642   0.000  1.00  0.00
HETATM    3  H2  HOH     1      18.882  17.715   0.000  1.00  0.00            
HETATM    4  O   HOH     2       0.000   0.000   3.107  1.00  0.00            
HETATM    5  H1  HOH     2       0.957   0.000   3.107  1.00  0.00            
HETATM    6  H2  HOH     2      -0.240   0.927   3.107  1.00  0.00            
HETATM    7  O   HOH     3      18.642  18.642   6.214  1.00  0.00            
HETATM    8  H1  HOH     3      17.685  18.642   6.214  1.00  0.00            
HETATM    9  H2  HOH     3      18.882  17.715   6.214  1.00  0.00       
ENDMDL
MODEL     2
REMARK 250 time=0.0005 picoseconds
HETATM    1  O   HOH     1      18.642  18.642  -0.000  1.00  0.00            
HETATM    2  H1  HOH     1      17.685  18.642   0.000  1.00  0.00            
HETATM    3  H2  HOH     1      18.882  17.715   0.000  1.00  0.00            
HETATM    4  O   HOH     2      -0.000  -0.000   3.107  1.00  0.00            
HETATM    5  H1  HOH     2       0.958  -0.000   3.107  1.00  0.00            
HETATM    6  H2  HOH     2      -0.240   0.928   3.107  1.00  0.00            
HETATM    7  O   HOH     3      18.642  18.642   6.214  1.00  0.00            
HETATM    8  H1  HOH     3      17.685  18.642   6.214  1.00  0.00            
HETATM    9  H2  HOH     3      18.881  17.715   6.214  1.00  0.00           
ENDMDL")

        bodies = NBodySimulator.extract_from_pdb(pdb_data)
        water = WaterSPCFw(bodies, mH, mO, qH, qO,  jl_parameters, e_parameters, spc_paramters)
        t1 = 0τ
        t2 = 2τ
        simulation = NBodySimulation(water, (t1, t2), pbc, kb)
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        
        ε = 0.00005

        @test bodies[1].O.r[1] - 1.8642 ≈ 0.0 atol = ε
        @test bodies[1].O.r[2] - 1.8642 ≈ 0.0 atol = ε
        @test bodies[1].O.r[3] - 0.0 ≈ 0.0 atol = ε

        @test bodies[2].H1.r[1] - 0.0957 ≈ 0.0 atol = ε
        @test bodies[2].H1.r[2] - 0.0 ≈ 0.0 atol = ε
        @test bodies[2].H1.r[3] - 0.3107 ≈ 0.0 atol = ε

        @test bodies[3].H2.r[1] - 1.8882 ≈ 0.0 atol = ε
        @test bodies[3].H2.r[2] - 1.7715 ≈ 0.0 atol = ε
        @test bodies[3].H2.r[3] - 0.6214 ≈ 0.0 atol = ε

        @test bodies[3].H2.v[1] - 0.0 ≈ 0.0 atol = ε
        @test bodies[3].H2.v[2] - 0.0 ≈ 0.0 atol = ε
        @test bodies[3].H2.v[3] - 0.0 ≈ 0.0 atol = ε
    end    
end
