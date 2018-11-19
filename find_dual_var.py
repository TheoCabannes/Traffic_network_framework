import numpy as np
import cvxpy as cvx

#Given link flows, find a feasible solution to TAP-C

l = 10 #Number of links
r = 20 #number of routes
d = np.zeros(l,r) # link-route incidence matrix
f = np.ones(l,1)

h = cvx.Variable((r,1))
obj = cvx.Minimize(cvx.sum_squares(f-d*h))
con = [h>=0]
prob = cvx.Problem(obj,con)
result = prob.solve()
h_opt = h.value
