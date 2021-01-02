# SampledCG

The 'Penalty.m' code usage is as follows:

[v,xbest,sol_MC,constr_viol] = Penalty(InputType, InputData, e1, max_time, tolerance)

InputType = ‘R’ (generates random instance) or ‘G’ (takes GSet as input)

InputData = integer value if Type = ‘R’ or ‘path/filename’ if Type = ‘G’

e1 - value of epsilon (approximation parameter); typically set to 0.1

max_time = maximum run time (in seconds)

tolerance = tolerance level set for computing eigenvector using eigs (set to 0.1)

Output:

v = mapped vector with v(1) objective function value of MaxCut SDP  and v(2) to v(n+1) is constraint violation

xbest_sol = binary vector representing cut

sol_MC = cut value generated for MaxCut

constr_viol = average constraint violation
