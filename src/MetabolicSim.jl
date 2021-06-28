module MetabolicSim

    import DiffEqBase, OrdinaryDiffEq, DiffEqCallbacks

    include("parameters/parameters.jl")
    include("simulation/assembly//assembly.jl")
    include("simulation/dxdt.jl")
    include("simulation/simulate.jl")
    include("analysis/analysis.jl")
    
    #parameter functions
    export make_parameters

    #callback functions 
    # export add_at_equilibrium, add_at_t
    
    #simulate functions
    export dx!
    
    #analysis
    export get_timeseries, get_Nsp
    export get_resource_a, get_overlap_competition, get_overlap_facilitation 
    export get_R,get_ω
end