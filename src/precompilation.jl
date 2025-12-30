using PrecompileTools

@setup_workload begin
    using StaticArrays

    @compile_workload begin
        # Precompile gravitational system workflow
        bodies_grav = [
            MassBody(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 1.0),
            MassBody(SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), 1.0),
        ]
        system_grav = GravitationalSystem(bodies_grav, 1.0)
        simulation_grav = NBodySimulation(system_grav, (0.0, 0.1))
        result_grav = run_simulation(simulation_grav)

        # Access result data
        get_position(result_grav, 0.05)
        get_velocity(result_grav, 0.05)
        kinetic_energy(result_grav, 0.05)
        total_energy(result_grav, 0.05)
        temperature(result_grav, 0.05)

        # Precompile Lennard-Jones system workflow with periodic boundaries
        bodies_lj = generate_bodies_in_cell_nodes(4, 1.0, 0.1, 2.0)
        lj_params = LennardJonesParameters(1.0, 0.5, 1.5)
        pbc = CubicPeriodicBoundaryConditions(2.0)
        system_lj = PotentialNBodySystem(bodies_lj, Dict(:lennard_jones => lj_params))
        simulation_lj = NBodySimulation(system_lj, (0.0, 0.1), pbc)

        # Run with default solver (Tsit5)
        result_lj = run_simulation(simulation_lj)

        # Run with VelocityVerlet (common for molecular dynamics)
        result_lj_vv = run_simulation(simulation_lj, VelocityVerlet(), dt = 0.01)

        # Access LJ result data
        get_position(result_lj, 0.05)
        get_velocity(result_lj, 0.05)
        kinetic_energy(result_lj, 0.05)
        potential_energy(result_lj, 0.05)
        total_energy(result_lj, 0.05)
    end
end
