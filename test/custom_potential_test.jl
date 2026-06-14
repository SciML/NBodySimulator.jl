@safetestset "A system with custom potential" begin
    # Body lives in an included file so that the `import NBodySimulator.get_accelerating_function`
    # + method extension are evaluated as module-top-level statements. Inside an inline
    # @safetestset body they land in the testset's local scope, where `import` fails with
    # "expected Symbol, got a value of type Core.SlotNumber".
    include("custom_potential_body.jl")
end
