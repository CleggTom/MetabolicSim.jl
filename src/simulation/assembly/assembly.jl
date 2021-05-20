#fucntion to get argument names
function method_argnames(m::Method)
    argnames = ccall(:jl_uncompress_argnames, Vector{Symbol}, (Any,), m.slot_syms)
    isempty(argnames) && return argnames
    return argnames[1:m.nargs]
end

function assert_assembly_callbacks(biomass_func::Function, new_parameters::Function)
    #assert functions have correct args and return types
    @assert Base.return_types(biomass_func)[1] == Float64 "`biomass_func` must return biomass of new consumer as a Float64)"
    @assert method_argnames(collect(methods(biomass_func))[1]) == [Symbol("#self#"), :integrator] "`biomass_func` must take only `integrator` as argument"

    @assert Base.return_types(new_parameters)[1] == Tuple{Float64,Any} "`new_parameters` must return biomass of new uptake matrix and maintenance cost as a `Tuple{Any,Float64}`"
    @assert method_argnames(collect(methods(new_parameters))[1]) == [Symbol("#self#"), :integrator] "`new_parameters` must take only `integrator` as argument"
end


"""
    add_consumer!(integrator, biomass_func::Function, new_parameters::Function)

Function to add consumer to system within a callback. Takes two functions `biomass_func` and `new_parameters` which give the biomass of the new consumer and the updated parameters respectively. Both functions should take `integrator` as the argument and return either a single `Float64` for biomass or a new `u` array and `Rm` value for new parameters. 
"""
function add_consumer!(integrator, biomass_func::Function, new_parameters::Function)


    #find extinct species
    extant = findall(integrator.u[integrator.p.M+1 : end] .> 1e-5)
    # println(extant)
    # println(integrator.u)

    #move biomasses
    integrator.u[(integrator.p.M+1) : (integrator.p.M+length(extant))] .= integrator.u[extant .+ integrator.p.M]
    # #resize cache to add new species ....
    resize!(integrator, length(extant) + integrator.p.M + 1)
    integrator.u[end] = biomass_func(integrator)
    # println(integrator.u)
    #alter parameters
    new_p = new_parameters(integrator)
    #N - get from extant species
    integrator.p.N = length(extant) + 1

    #Rm - resize to add new and then remove
    integrator.p.Rm[1:length(extant)] .= integrator.p.Rm[extant]
    resize!(integrator.p.Rm, integrator.p.N)
    integrator.p.Rm[end] = new_p[1]

    # #U
    new_u = Array{Float64,2}(undef, integrator.p.N, integrator.p.M)
    new_u[1:(end-1), :] .= integrator.p.u[extant,:]
    new_u[end,:] .= new_p[2]
    integrator.p.u = new_u
end

"""
     consumer_equilibrium(integrator, tmax::Float64, abstol::Float64, reltol::Float64)

function to detect equilbirium based on an absolute tolerance `abstol` and relative tolerance `reltol`. Will not detect equilibrium past a given time `tmax` 
"""
function consumer_equilibrium(integrator, tmax::Float64, abstol::Float64, reltol::Float64)
    #test time - dont return true when over Tmax
    # if integrator.t > tmax
    #    return(false)     
    # end

    #get derivative - from DiffEqCallbacks TerminateSteadyState
    if DiffEqBase.isinplace(integrator.sol.prob)
        testval = first(DiffEqBase.get_tmp_cache(integrator))
        DiffEqBase.get_du!(testval, integrator)
        if typeof(integrator.sol.prob) <: DiffEqBase.DiscreteProblem
            @. testval =  testval - integrator.u
        end
    else
        testval = get_du(integrator)
        if typeof(integrator.sol.prob) <: DiffEqBase.DiscreteProblem
            testval =  testval - integrator.u
        end
    end
    #set resource du to 0

    testval[1:integrator.p.M] .= 0.0
    # println(testval)

    any(abs(d) > abstol && abs(d) > reltol*abs(u) for (d,abstol, reltol, u) = zip(testval, Iterators.cycle(abstol), Iterators.cycle(reltol), integrator.u)) && (return false)

    return(true)

end

"""

"""
function add_at_equilibrium(biomass_func::Function, new_parameters::Function; tmax::Float64 = Inf, abstol::Float64 = 1e-8, reltol::Float64 = 1e-6)
    # #warn if tolerances are high
    # if abstol > 1e-2
    #     @warn "abstol is set greater than `1e-2` (`$abstol`). Are you sure you dont want tighter tolerances?"
    # elseif reltol > 1e-2
    #     @warn "reltol is set greater than `1e-2` (`$reltol`). Are you sure you dont want tighter tolerances?"
    # end

    #assert functions have correct args and return types
    assert_assembly_callbacks(biomass_func, new_parameters)

    #define condtion and affect!
    condition = (u, t, integrator) -> consumer_equilibrium(integrator, tmax, abstol, reltol)
    affect! = (integrator) -> add_consumer!(integrator, biomass_func, new_parameters)
    return(:equilibrium , DiffEqCallbacks.DiscreteCallback(condition, affect!; save_positions = (true, true)))
end

function add_at_t(biomass_func::Function, new_parameters::Function, times::Vector{Float64})
    #assert functions have correct args and return types
    assert_assembly_callbacks(biomass_func, new_parameters)

    affect! = (integrator) -> add_consumer!(integrator, biomass_func, new_parameters)
    return(:at_time_stop, DiffEqCallbacks.PresetTimeCallback(times, affect!), times)
end