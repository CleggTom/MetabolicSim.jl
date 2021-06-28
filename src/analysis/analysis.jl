function get_timeseries(sol)
    sol_copy = deepcopy(sol)

    for i = eachindex(sol_copy.u)
        while length(sol_copy.u[i]) < length(sol_copy.u[end])
            ind = length(sol_copy.u[i]) - sol_copy.prob.p.M + 1
            insert!(sol_copy.u[i], ind, 0.0)
        end
    end

    return(hcat(sol_copy.u...)')
end

function get_Nsp(sol, abstol::Float64 = 1e-10)
    series = get_timeseries(sol)
    return(mapslices(x -> sum(x .> abstol), series[: , 1:sol.prob.p.N], dims = 2)) 
end

#respiration callback
function save_func(u, t, integrator)
    integrator.p.Rm .* integrator.u[(integrator.p.M+1 : end)]
end

###
#Interactions
###

#get resource explicit competition 
function get_a_ij_comp(y,u,l_sum,i,j,k)
    y[k] * u[i,k] * (u[j,k] * (1 - l_sum[k]))
end

function get_a_ij_facil(y,u,l,l_sum,i,j,k,M)
    sum_var = 0.0
    for l_ = 1:M
        sum_var += y[l_] * u[i,l_] * l[l_,k]
    end
    return(sum_var * (u[j,k] * (1 - l_sum[k])))
end

function get_resource_a(sol, t, type::Symbol)
    #assert the interaction type is valid
    @assert type ∈ [:competition, :facilitation, :both] "type of interaction must be competition, facilitation or both"
    #get resource enviroment and parameter object
    p = sol.prob.p
    y = sol(t)[1:p.M]

    #loop over all species pairs and resources
    a = zeros(p.N,p.N)
    for i = 1:p.N, j = 1:p.N
        for k = 1:p.M
            if type == :competition
                a[i,j] -= get_a_ij_comp(y,p.u,p.l_sum,i,j,k)
            elseif type == :facilitation
                a[i,j] += get_a_ij_facil(y,p.u,p.l,p.l_sum,i,j,k,p.M)
            elseif type == :both
                a[i,j] += get_a_ij_facil(y,p.u,p.l,p.l_sum,i,j,k,p.M) - get_a_ij_comp(y,p.u,p.l_sum,i,j,k)
            end
            
        end
    end
    return(a)
end

#naive interaction measures
#competition - degree to which the uptake of i overlaps with the uptake of j
function get_overlap_competition(sol)
    p = sol.prob.p
    a = zeros(p.N,p.N)
    for i = 1:p.N, j = 1:p.N
        t = 0.0
        b = 0.0
        for k = 1:p.M
            t += min(p.u[i,k],p.u[j,k])
            b += max(p.u[i,k],p.u[j,k])
        end
        a[i,j] = t / b
    end
    return(a)
end

#facilitation - degree to which the leaked uptake of i overlaps with the uptake of j
function get_overlap_facilitation(sol)
    p = sol.prob.p
    a = zeros(p.N,p.N)
    for i = 1:p.N, j = 1:p.N
        t = 0.0
        b = 0.0
        for k = 1:p.M
            t += min(sum(p.u[i,:] .* p.l[:,k]), p.u[j,k])
            b += max(sum(p.u[i,:] .* p.l[:,k]), p.u[j,k])
        end
        a[i,j] = t / b
    end
    return(a)
end

#get respiration
function get_R(sol,t)
    sum(sol(t)[(sol.prob.p.M+1):end] .* sol.prob.p.Rm)
end

#get outflow
function get_ω(sol,t)
    sum(sol(t)[1:sol.prob.p.M] .* sol.prob.p.ω)
end