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

function add_at_equilibrium(biomass_func::Function, new_parameters::Function, tmax::Float64 = Inf; abstol::Float64 = 1e-8, reltol::Float64 = 1e-6)
    condition = (u, t, integrator) -> consumer_equilibrium(integrator, tmax, abstol, reltol)
    affect! = (integrator) -> add_consumer!(integrator, biomass_func, new_parameters)
    return(:equilibrium , DiffEqCallbacks.DiscreteCallback(condition, affect!; save_positions = (true, false)))
end

function add_at_t(biomass_func::Function, new_parameters::Function, times::Vector{Float64})
    affect! = (integrator) -> add_consumer!(integrator, biomass_func, new_parameters)
    return(:at_time_stop, DiffEqCallbacks.PresetTimeCallback(times, affect!), times)
end