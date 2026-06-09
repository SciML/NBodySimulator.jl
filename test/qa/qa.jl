using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using NBodySimulator, Aqua, JET, Test

@testset "Aqua" begin
    Aqua.test_all(NBodySimulator)
end

@testset "JET" begin
    JET.test_package(NBodySimulator; target_defined_modules = true)
end
