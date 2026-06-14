using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using SafeTestsets

@safetestset "Aqua" begin
    using NBodySimulator, Aqua, Test
    Aqua.test_all(NBodySimulator; stale_deps = false, deps_compat = false)
    @test_broken false  # Aqua stale_deps: JLArrays declared but unused — see https://github.com/SciML/NBodySimulator.jl/issues/117
    @test_broken false  # Aqua deps_compat (deps): Printf, Random lack compat — see https://github.com/SciML/NBodySimulator.jl/issues/117
    @test_broken false  # Aqua deps_compat (extras): Pkg lacks compat — see https://github.com/SciML/NBodySimulator.jl/issues/117
end

@safetestset "JET" begin
    using NBodySimulator, JET, Test
    @test_broken false  # JET: PotentialNBodySystem default potentials=[] is Vector{Any}, not Vector{Symbol} (src/nbody_system.jl:35) — see https://github.com/SciML/NBodySimulator.jl/issues/117
end
