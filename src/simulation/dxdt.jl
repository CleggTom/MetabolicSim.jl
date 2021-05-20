
"""
    dx!(dx,x,p::Parameters,t)

Function defining change in amount of consumers and resources over time. 

Dynamics occur as described by the set of equations:

``\\frac{dC_i}{dt} = C_i \\left(∑\\right)``

"""
function dx!(dx,x,p::Parameters,t)
  #loop over consumers i = 1:N
  for i = 1:p.N
    dx[p.M + i] = -x[p.M + i] * p.Rm[i]
    #loop over metabolites
    for j = 1:p.M
      dx[p.M + i] +=  x[p.M + i] * x[j] * p.u[i,j] * (1 - p.l_sum[j]) # u_ij * (1-l_j)
    end
  end

 #loop over metabolites j = 1:M
  for j = 1:p.M
    dx[j] = p.ρ[j] - x[j] * p.ω[j]#inflow
    #loop over species
    for i = 1:p.N
      dx[j] += -p.u[i,j] * x[p.M+ i] * x[j] #minus uptake
      #loop over other metabolites
      for k = 1:p.M
        dx[j] += p.u[i,k] * x[p.M + i] * x[k] * p.l[k,j] # +leakage
      end
    end
  end

  # dx[(x .+ dx) .< 0] .= -x[(x .+ dx) .< 0]
end
