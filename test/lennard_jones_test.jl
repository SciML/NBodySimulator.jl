@testset "Lennard-Jones potential" begin
    @testset "Two interacting particles" begin

        m1 = 1.0
        m2 = 1.0
        r1 = 1.3

        t1 = 0.0
        t2 = 1.0
        τ = (t2 - t1) / 1000

        p1 = MassBody(SVector(-r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m1)
        p2 = MassBody(SVector(r1 / 2, 0.0, 0.0), SVector(0.0, 0.0, 0.0), m2)

        σ = 1.0
        ϵ = 1.0
        parameters = LennardJonesParameters(ϵ, σ, Inf)
        system = PotentialNBodySystem([p1, p2], Dict(:lennard_jones => parameters))
        simulation = NBodySimulation(system, (t1, t2))
        sim_result = run_simulation(simulation, VelocityVerlet(), dt=τ)


        r2 = get_position(sim_result, t2, 2) - get_position(sim_result, t2, 1)
        v_expected = sqrt(4ϵ / m1 * ( ((σ / norm(r1))^12 - (σ / norm(r2))^12) - ((σ / norm(r1))^6 - (σ / norm(r2))^6 ) ))
        v_actual = norm(get_velocity(sim_result, t2, 2))

        ε = 0.001 * v_expected
        @test v_expected ≈ v_actual atol = ε

        io = IOBuffer()
    
        @test sprint(io -> show(io, system)) == 
    "Potentials: \nLennard-Jones:\n\tϵ:1.0\n\tσ:1.0\n\tR:Inf\n"

        @test sprint(io -> show(io, simulation)) == 
    "Timespan: (0.0, 1.0)\nBoundary conditions: InfiniteBox{Float64}([-Inf, Inf, -Inf, Inf, -Inf, Inf])\nPotentials: \nLennard-Jones:\n\tϵ:1.0\n\tσ:1.0\n\tR:Inf\n"

        @test sprint(io -> show(io, sim_result)) == 
    "N: 2\nTimespan: (0.0, 1.0)\nBoundary conditions: InfiniteBox{Float64}([-Inf, Inf, -Inf, Inf, -Inf, Inf])\nPotentials: \nLennard-Jones:\n\tϵ:1.0\n\tσ:1.0\n\tR:Inf\nTime steps: 1001\nt: 0.0, 1.0\n"
    end

    @testset "Three particles of liquid argon and their \"temperature\"" begin 
        T = 120.0 # °K
        kb = 8.3144598e-3 # kJ/(K*mol)
        ϵ = T * kb
        σ = 0.34 # nm
        ρ = 1374/1.6747# Da/nm^3
        m = 39.95# Da
        L = 5σ # 10.229σ
        N = 3 # floor(Int, ρ * L^3 / m)
        R = 2.25σ   
        v_dev = sqrt(3 * kb * T / m)
        r1 = SVector(L / 3, L / 3, 2 * L / 3)
        r2 = SVector(L / 3, 2 * L / 3, L / 3)
        r3 = SVector(2 * L / 3, L / 3, L / 3)
        v1 = SVector(0, 0, -v_dev)
        v2 = SVector(0, -v_dev, 0)
        v3 = SVector(-v_dev, 0, 0)
        p1 = MassBody(r1, v1, m)
        p2 = MassBody(r2, v2, m)
        p3 = MassBody(r3, v3, m)

        τ = 0.5e-3 # ps or 1e-12 s
        t1 = 0.0
        t2 = 100τ

        parameters = LennardJonesParameters(ϵ, σ, R)
        lj_system = PotentialNBodySystem([p1, p2, p3], Dict(:lennard_jones => parameters));
        simulation = NBodySimulation(lj_system, (t1, t2), PeriodicBoundaryConditions(L), kb);
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        T1 = temperature(result, t1) 
        ε = 1e-6
        @test T1 ≈ 120.0 atol = ε

        e_kin_1 = m * (dot(v1, v1) + dot(v2, v2) + dot(v3, v3)) / 2
        @test e_kin_1 == kinetic_energy(result, t1)

        e_tot_1 = total_energy(result, t1)
        ε = 0.1 * e_tot_1
        e_tot_2 = total_energy(result, t2)
        @test e_tot_1 ≈ e_tot_2 atol = ε

        simulation = NBodySimulation(lj_system, (t1, t2), CubicPeriodicBoundaryConditions(L));
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)
        e_tot_1 = total_energy(result, t1)
        ε = 0.1 * e_tot_1
        e_tot_2 = total_energy(result, t2)
        @test e_tot_1 ≈ e_tot_2 atol = ε

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

        molecule_number = hetatm_count / timestep_count

        @test length(result.solution.t) == timestep_count
        @test molecule_number == length(lj_system.bodies)

    end

    @testset "Testing RDF and MSD calculation" begin 
        T = 120.0 # °K
        kb = 8.3144598e-3 # kJ/(K*mol)
        ϵ = T * kb
        σ = 0.34 # nm
        ρ = 1374/1.6747# Da/nm^3
        m = 39.95# Da
        L = 5σ # 10.229σ
        N = 125 # floor(Int, ρ * L^3 / m)
        R = 0.5L   
        τ = 0.5e-3 # ps or 1e-12 s
        t1 = 0.0
        t2 = 400τ
        v_dev = sqrt(kb * T / m)
        bodies = generate_bodies_in_cell_nodes(N, m, v_dev, L)

        parameters = LennardJonesParameters(ϵ, σ, R)
        lj_system = PotentialNBodySystem(bodies, Dict(:lennard_jones => parameters));
        simulation = NBodySimulation(lj_system, (t1, t2), CubicPeriodicBoundaryConditions(L), kb);
        result = run_simulation(simulation, VelocityVerlet(), dt=τ)

        (ts, mean_square_displacement) = msd(result)
        @test mean_square_displacement[1] < mean_square_displacement[end]

        (rs, grs) = rdf(result)
        (val, ind) = findmax(grs)
        @test (rs[ind] / σ) ≈ 1.0 atol = 1.0
    end

    @testset "Constructing electorstatic potential parameters entity" begin
        default_potential = LennardJonesParameters()
        @test 1.0 == default_potential.ϵ
        @test 1.0 == default_potential.σ
        @test 2.5 == default_potential.R
        @test 1.0 == default_potential.σ2
        @test 6.25 == default_potential.R2
    
        io = IOBuffer()
        potential1 = LennardJonesParameters()
        potential2 = LennardJonesParameters(2, 5, 10)
        @test sprint(io -> show(io, potential1)) == "Lennard-Jones:\n\tϵ:1.0\n\tσ:1.0\n\tR:2.5\n"
        @test sprint(io -> show(io, potential2)) == "Lennard-Jones:\n\tϵ:2\n\tσ:5\n\tR:10\n"
    end

    @testset "Periodic boundar conditions" begin
        pbc = PeriodicBoundaryConditions(15e-6)
        boundary = [0.0, 15e-6, 0.0, 15e-6, 0.0, 15e-6]
        ind = 1
        for b in pbc
            b = boundary[ind]
            ind += 1
        end
        for i = 1:6
            @test boundary[i] == pbc[i]
        end
    end
end