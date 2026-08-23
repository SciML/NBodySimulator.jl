using SciMLTesting, NBodySimulator, JET, Test

# The solver interface NBodySimulator deliberately reexports so that `using
# NBodySimulator` on its own is enough to write the documented workflow
# (`run_simulation(simulation, VelocityVerlet(), dt = τ)`), build a problem from a
# simulation, solve it, and inspect the result. Owned and documented upstream; kept in
# sync with the reexport `export` blocks in src/NBodySimulator.jl.
const REEXPORTS = (
    # Second-order-ODE integrators accepted by `run_simulation` (symplectic).
    :CalvoSanz4, :CandyRoz4, :KahanLi6, :KahanLi8, :LeapfrogDriftKickDrift, :McAte2,
    :McAte3, :McAte4, :McAte42, :McAte5, :McAte8, :PseudoVerletLeapfrog, :Ruth3,
    :SofSpa10, :SymplecticEuler, :VelocityVerlet, :VerletLeapfrog, :Yoshida6,
    # Second-order-ODE integrators accepted by `run_simulation` (Runge-Kutta-Nyström).
    :DPRKN12, :DPRKN4, :DPRKN5, :DPRKN6, :DPRKN6FM, :DPRKN8, :ERKN4, :ERKN5, :ERKN7,
    :FineRKN4, :FineRKN5, :IRKN3, :IRKN4, :Nystrom4, :Nystrom4VelocityIndependent,
    :Nystrom5VelocityIndependent, :RKN4,
    # `run_simulation`'s default algorithm.
    :Tsit5,
    # The SciML common interface.
    :ArrayPartition, :CallbackSet, :ContinuousCallback, :DiscreteCallback, :ODEProblem,
    :ODESolution, :RODESolution, :ReturnCode, :SDEProblem, :SecondOrderODEProblem,
    :VectorContinuousCallback, :init, :remake, :solve, :solve!, :step!,
    :successful_retcode,
)

run_qa(NBodySimulator; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using NBodySimulator`, so
    # the allow-list cannot drift into approving names the package no longer provides.
    # `isdefined(@__MODULE__, ...)` tests the property directly: this file's `using
    # NBodySimulator` is what has to bring the name into scope.
    @testset "$name" for name in REEXPORTS
        @test name in names(NBodySimulator)
        @test isdefined(@__MODULE__, name)
    end
end
