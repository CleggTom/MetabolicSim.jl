"""
    Parameters(N::Int64, M::Int64, u::Array{Float64,2}, Rm::Vector{Float64}, l::Array{Float64,2}, l_sum::Vector{Float64}, ρ::Vector{Float64})

Type containing the parameters for a simuation.

# Arguments

- `N::Int64` : Number of consumers
- `M::Int64`: Number of Resourcces
- `u::Array{Float64,2}`: Array of uptake rates of resource i by consumer j
- `Rm::Array{Float64,1}`: Vector of maintennce respiration costs for consumers
- `l::Array{Float64,2}`: Array of the proportion of uptake of resource i released as resource j
- `l_sum::Array{Float64,1}`: Array of the total proportion of uptake of resource j released
- `ρ::Array{Float64,1}`: Vector of inflow rates for resources
- `ω::Array{Float64,1}`: Vector of outflow rates for resources
"""
mutable struct Parameters
    N::Int64
    M::Int64
    u::Array{Float64,2}
    Rm::Array{Float64,1}
    l::Array{Float64,2}
    l_sum::Array{Float64,1}
    ρ::Array{Float64,1}
    ω::Array{Float64,1}
end

"""
    make_parameters(N::Int64, M::Int64, u::Array{Float64,2}, Rm::Vector{Float64}, l::Array{Float64,2}, ρ::Vector{Float64}, ω::Array{Float64,1})

Helper function used internally. Takes values for parameters and returns a `Parameters`
object. Also does checks internally to make sure the values are correct.
"""
function make_parameters(N::Int64, M::Int64,
        u::Array{Float64,2}, Rm::Vector{Float64},
        l::Array{Float64,2}, ρ::Vector{Float64},ω::Array{Float64,1})
    
    #assert parameters are the right size
    @assert all([N,M] .> 0) "Number of consumers and resources must both be greater than 0"
    @assert size(u) == (N,M) "Uptake array is the wrong size"
    @assert length(Rm) == N "Rm vector is the wrong length"
    @assert size(l) == (M,M) "Leakage array is the wrong size"
    @assert length(ρ) == M "Inflow vector is the wrong length"

    #assert that the leakage proportions add up to 1
    @assert all(sum(l,dims = 2) .<= 1) "One or more l values are > 1"

    return( Parameters(copy(N), copy(M), copy(u), copy(Rm), copy(l), sum(copy(l), dims=2)[:], copy(ρ), copy(ω)))
end


"""

"""
function make_parameters(N::Int64, M::Int64,
        u::Array{Float64,2}, Rm::Vector{Float64},
        l::Array{Float64,2}, ρ::Vector{Float64})
    
    #assert parameters are the right size
    @assert all([N,M] .> 0) "Number of consumers and resources must both be greater than 0"
    @assert size(u) == (N,M) "Uptake array is the wrong size"
    @assert length(Rm) == N "Rm vector is the wrong length"
    @assert size(l) == (M,M) "Leakage array is the wrong size"
    @assert length(ρ) == M "Inflow vector is the wrong length"

    #assert that the leakage proportions add up to 1
    @assert all(sum(l,dims = 2) .<= 1) "One or more l values are > 1"

    return( Parameters(copy(N), copy(M), copy(u), copy(Rm), copy(l), sum(copy(l), dims=2)[:], copy(ρ)))
end