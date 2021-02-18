using MetabolicSim
using Test

#function to create 1x1 matricies
mat(sc) = reshape([sc],1,1)

#set up basic simulation: 1 consumer and 1 resources
M,N = 1,1
u,Rm = mat(1.0), [1.0]
l,ρ = mat(0.2), [0.2]

#Simulate
p = make_parameters(N,M,u,Rm,l,ρ)
@test_nowarn sol = simulate(ones(N+M),p,t_end = 100)

#test vs analytical prediction
C_eq = ρ[1] / Rm[1]
R_eq = -Rm[1] / (u[1] * (l[1] - 1))
@test all(isapprox.(sol[end],[C_eq,R_eq], rtol = 1e-2))
