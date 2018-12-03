import numpy as np
import cvxpy as cvx

def lower_bound(f, h, t, Delta, route2od, u_old, u_new, nk):
    """ 
    function that compute a lower bound with updated capacity given the original solution of a TAP-C
    input:
    f(vector): primal solution of link flow
    h(vector): primal solution of route flow
    t(ventor): link travel time of f
    Delta(matrix)[r,a]: 1 if route r use link a; 0 otherwise. 
    route2od(vector): mapping the route index to the index of o-d pairs
    u_old(vector): original capacity of links
    u_new(vector): new capacity of links   
    nk(int): number of OD pairs
    """
    nr, na = Delta.shape # number of non-zero routes and number of links 
    c = Delta.dot(t)
    
    lambda_link = cvx.Variable(na)
    pi = cvx.Variable(nk)
    
    obj = cvx.Maximize((u_old - u_new) * lambda_link)
    
    tolerance = 1e-3 # since Frank-Wolfe alg might not give the best solution, we set the tolerance 
    
    constr_1 = [lambda_link[i]==0 for i in range(na) if f[i] < u_old[i] - tolerance]
    constr_2 = [lambda_link[i]>=0 for i in range(na)]
    constr_3 = [c[r]+Delta[r]*lambda_link == pi[route2od[r]] for r in range(nr) if h[r] > tolerance]
    constr_4 = [c[r]+Delta[r]*lambda_link >= pi[route2od[r]] for r in range(nr)]
    
    constr = constr_1 + constr_2 + constr_3 + constr_4
    
    problem = cvx.Problem(obj, constr)
    bound = problem.solve(solver=cvx.CPLEX, verbose=True, feastol=tolerance, reltol=tolerance, abstol=tolerance)
    
    return bound, lambda_link.value, pi.value

