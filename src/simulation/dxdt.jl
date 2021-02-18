
"""
    dx!(dx,x,p::Parameters,t)

Function defining change in amount of consumers and resources over time. 

Dynamics occur as described by the set of equations:

``\\frac{dC_i}{dt} = C_i \\left(∑\\right)``

"""
function dx!(dx,x,p::Parameters,t)
  for i = 1:p.N
    dx[i] = -x[i] * p.Rm[i] # - Rm
    #loop over metabolites
    for j = 1:p.M
      dx[i] +=  x[i] * x[p.N + j] * p.u[i,j] * (1-p.l_sum[j]) # u_ij * (1-l_j)
    end
  end

 #loop over metabolites j = 1 : M
  for j = 1:p.M
    dx[p.N + j] = p.ρ[j] # inflow
    #loop over species
    for i = 1:p.N
      dx[p.N + j] += -p.u[i,j] * x[i] * x[p.N + j] # -uptake
      #loop over other metabolites
      for k = 1:p.M
        dx[p.N + j] += p.u[i,k] * x[i] * x[p.N + k] * p.l[k,j] # +leakage
      end
    end
  end
end
