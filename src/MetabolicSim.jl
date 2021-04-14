module MetabolicSim

    import DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks

    include("parameters/parameters.jl")
    include("simulation/assembly//assembly.jl")
    include("simulation/dxdt.jl")
    include("simulation/simulate.jl")
    include("analysis/analysis.jl")
    
    export make_parameters
    export add_at_equilibrium, add_at_t
    export simulate
    export get_timeseries, get_Nsp
end