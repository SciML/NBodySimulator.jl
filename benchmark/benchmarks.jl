using BenchmarkTools, NBodySimulator

const SUITE = BenchmarkGroup()

include("bench_lj.jl")
