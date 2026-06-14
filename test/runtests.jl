using Test, SafeTestsets

# `using NBodySimulator` must stay at Main top-level: several tests assert on
# `show`/`string` of types defined in NBodySimulator (e.g. InfiniteBox). Julia
# qualifies a type's module in `show` based on visibility in `Base.active_module()`
# (== Main), not the module the test body runs in. Without the bare names in Main,
# the printed strings become "NBodySimulator.InfiniteBox{...}" and the assertions fail.
using NBodySimulator

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Core"
    include("lennard_jones_test.jl")
    include("electrostatics_test.jl")
    include("gravitational_test.jl")
    include("custom_potential_test.jl")
    include("magnetostaic_test.jl")
    include("thermostat_test.jl")
    include("water_test.jl")
    include("interface_test.jl")
end

if GROUP == "QA"
    include("qa/qa.jl")
end
