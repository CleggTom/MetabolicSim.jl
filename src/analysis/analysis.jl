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
