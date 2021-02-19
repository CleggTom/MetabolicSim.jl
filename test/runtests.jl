using SafeTestsets

@safetestset "parameter tests:" begin include("parameter_tests.jl") end
@safetestset "simulation tests:" begin include("simulation_tests.jl") end
@safetestset "assembly tests:" begin include("assembly_tests.jl") end
