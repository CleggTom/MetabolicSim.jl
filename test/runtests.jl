using SafeTestsets

@safetestset "parameter_tests" begin include("parameter_tests.jl") end
@safetestset "simulation_tests" begin include("simulation_tests.jl") end
