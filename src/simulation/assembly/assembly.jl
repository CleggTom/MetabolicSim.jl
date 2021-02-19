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

    @assert Base.return_types(new_parameters)[1] == Tuple{Any,Float64} "`new_parameters` must return biomass of new uptake matrix and maintenance cost as a `Tuple{Any,Float64}`"
    @assert method_argnames(collect(methods(new_parameters))[1]) == [Symbol("#self#"), :integrator] "`new_parameters` must take only `integrator` as argument"
end

#default callback functions
function biomass_func_default(integrator)
    return(rand())
end

function new_parameters_default(integrator)
    p = integrator.p
    θ = rand()
    return( [p.u ; θ * rand(p.M)'], (1-θ)^0.5 )
end

"""
    add_consumer!(integrator, biomass_func::Function, new_parameters::Function)

Function to add consumer to system within a callback. Takes two functions `biomass_func` and `new_parameters` which give the biomass of the new consumer and the updated parameters respectively. Both functions should take `integrator` as the argument and return either a single `Float64` for biomass or a new `u` array and `Rm` value for new parameters. 
"""
function add_consumer!(integrator, biomass_func::Function, new_parameters::Function)
    #get arrays
    u = integrator.u 
    p = integrator.p
    
    #resize the actual integrator
    resize!(integrator,length(u)+1)
    #move masses around
    @views u[(p.N + 2):(p.N + p.M + 1)] .= u[(p.N + 1):(p.N + p.M)]
    @views u[p.N + 1] = biomass_func(integrator)

    #update params object
    p.N += 1

    new_sp = new_parameters(integrator)
    p.u = new_sp[1]
    append!(p.Rm, new_sp[2])
end


function consumer_equilibrium(integrator, tmax::Float64, abstol::Float64, reltol::Float64)
    #get dt
    testval = first(DiffEqBase.get_tmp_cache(integrator))
    DiffEqBase.get_du!(testval, integrator)
    #set resource du to 0
    testval[(integrator.p.N+1) : end] .= 0.0
    #test
    if integrator.t < tmax
        if any(abs.(testval) .> abstol) & any(abs.(testval) .> reltol * integrator.u) 
            return(false)
        else
            return(true)
        end
    else
        return(false)
    end
end

function add_at_equilibrium(biomass_func::Function, new_parameters::Function; tmax::Float64 = Inf, abstol::Float64 = 1e-8, reltol::Float64 = 1e-6)
    #warn if tolerances are high
    if abstol > 1e-5
        @warn "abstol is set greater than `1e-5` (`$abstol`). Are you sure you dont want tighter tolerances?"
    elseif reltol > 1e-5
        @warn "reltol is set greater than `1e-5` (`$reltol`). Are you sure you dont want tighter tolerances?"
    end

    #assert functions have correct args and return types
    assert_assembly_callbacks(biomass_func, new_parameters)

    #define condtion and affect!
    condition = (u, t, integrator) -> consumer_equilibrium(integrator, tmax, abstol, reltol)
    affect! = (integrator) -> add_consumer!(integrator, biomass_func, new_parameters)
    return(:equilibrium , DiffEqCallbacks.DiscreteCallback(condition, affect!; save_positions = (true, false)))
end

function add_at_t(biomass_func::Function, new_parameters::Function, times::Vector{Float64})
    #assert functions have correct args and return types
    assert_assembly_callbacks(biomass_func, new_parameters)

    affect! = (integrator) -> add_consumer!(integrator, biomass_func, new_parameters)
    return(:at_time_stop, DiffEqCallbacks.PresetTimeCallback(times, affect!), times)
end