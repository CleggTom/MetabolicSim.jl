using MetabolicSim
using Test

#generate some parameters
N,M = 2,2
Rm,ρ = [1.0, 1.0],[0.4, 0.2]
u = [1.0 1.0 ; 1.0 0.0]
l = [0.1 0.0 ; 0.2 0.2]
l_sum = sum(l, dims=2)[:]
#make parameter object
p_1 = MetabolicSim.Parameters(N,M,u,Rm,l,l_sum,ρ)
p_2 = MetabolicSim.make_parameters(N,M,u,Rm,l,ρ)

#function to test if parameter structs are equal
function isequal(a::MetabolicSim.Parameters, b::MetabolicSim.Parameters)
    itt = fieldnames(MetabolicSim.Parameters)
    for i = itt
        if getfield(a, i) != getfield(b, i)
            return false
        end
    end
    return true
end

@testset begin
    #test Parameters and make_parameters give the same objects
    @test isequal(p_1,p_2)

    #test make_parameters assertion error for leakage
    l = ones(M,M)
   @test_throws AssertionError MetabolicSim.make_parameters(N,M,u,Rm,l,ρ)
end