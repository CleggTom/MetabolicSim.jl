using MetabolicSim
using Test
using DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks

#generate parameters
N,M = 2,2
Rm,ρ = [1.0, 1.0],[0.1, 0.1]
u = [1.0 0.2 ; 0.2 1.0]
l = [0.5 0.0 ; 0.0 0.2]
p = make_parameters(N,M,u,Rm,l,ρ)

#generate integrator
prob = ODEProblem(MetabolicSim.dx!, ones(N+M), (0,1.0), p)
integrator = DiffEqBase.init(prob, Tsit5())

#test default functions
@testset begin
    @test typeof(MetabolicSim.biomass_func_default(integrator)) == Float64
    
    new_p = MetabolicSim.new_parameters_default(integrator)
    @test size(new_p[1]) == (N+1, M)
    @test typeof(new_p[2]) == Float64
end


f_biomass = MetabolicSim.biomass_func_default
f_new_p = MetabolicSim.new_parameters_default
#test callback function assertion errors
@testset begin
    f_wrong_args(a) = Float64(1.0)
    f_wrong_result(integrator) = [1.0]

    @test_throws AssertionError MetabolicSim.assert_assembly_callbacks(f_wrong_args, f_new_p)
    @test_throws AssertionError MetabolicSim.assert_assembly_callbacks(f_biomass, f_wrong_args)
    @test_throws AssertionError MetabolicSim.assert_assembly_callbacks(f_wrong_result, f_new_p)
    @test_throws AssertionError MetabolicSim.assert_assembly_callbacks(f_biomass, f_wrong_result)
end

#test consumer_equilibirum assembly function
@testset begin
    #test function
    @test_nowarn add_at_equilibrium(f_biomass,f_new_p, tmax = 100.0, abstol = 1e-8, reltol = 1e-6)

    #test abs_warning
    abstol = 1.0
    abs_warn = "abstol is set greater than `1e-5` (`$abstol`). Are you sure you dont want tighter tolerances?"
    @test_logs (:warn, abs_warn) add_at_equilibrium(f_biomass,f_new_p, tmax = 100.0, abstol = abstol, reltol = 1e-6)

    #test rel_warning
    reltol = 1.0
    rel_warn = "reltol is set greater than `1e-5` (`$reltol`). Are you sure you dont want tighter tolerances?"
    @test_logs (:warn, rel_warn) add_at_equilibrium(f_biomass,f_new_p, tmax = 100.0, abstol = 1e-8, reltol = reltol)

    #test simulalation
    tmax = 500.0
    cb = add_at_equilibrium(f_biomass,f_new_p, tmax = tmax, abstol = 1e-6, reltol = 1e-6)
    
    @test cb[1] == :equilibrium
    @test isa(cb[2],SciMLBase.DECallback)

    sol = simulate(ones(N+M), p, cb; t_start=0.0, t_end= 1000.0)

    @test_nowarn sol = simulate(ones(N+M), p, cb; t_start=0.0, t_end= 1000.0)
    @test length(sol.u[end]) == 5
    #test assembly after tmax
    tmax_ind = findfirst(sol.t .> tmax)
    @test all(length.(sol.u[tmax_ind:end]) .== length(sol.u[tmax_ind]))
end

#test add_at_t assembly function
@testset begin
    t_stops = [1.0, 20.0, 30.0]
    @test_nowarn add_at_t(f_biomass,f_new_p,t_stops)

    #test t_tstops are the same
    cb = add_at_t(f_biomass,f_new_p,t_stops)
    cb[1] == :at_time_stop
    isa(cb[2], SciMLBase.DECallback)
    cb[3] == t_stops

    sol = simulate(ones(N+M), p, cb; t_start=0.0, t_end= 100.0)

    #test species are added at correct times
    @test length.(sol.(t_stops)) == (M + N .+ collect(0:2)) #at t_stop
    @test length.(sol.(t_stops .+ 0.1)) == (M + N .+ collect(1:3)) #after t_stops

    #test sol parameters are updated
    p_ = sol.prob.p
    @test p_.N == N+3
    @test p_.M == M
    @test size(p_.u) == (N+3, M)
    @test length(p_.Rm) == N+3

end

