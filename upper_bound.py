# Given link flows and capacity constraints, find a feasible set of link flows

# Inputs:
# f(vector): given (possibly infeasible) link flows 
# u_old(vector): link capacity constraints to be met
# delta(matrix): link*route incidence matrix
# route2od(vector): index of od pair for a given route
# od2route: gives indices of all routes for a given od pair
# K(vector): indices of OD pairs

# Output: 
# f_feas(vector): feasible link flows for new capacity constraints


import numpy as np
l = 10 #number of links in network
K = 15 #number of od pairs in network
r = 20 #number of routes in network
max_r = 5 #maximum number of routes for any od pair
l_nb = np.zeros((l,1)) #link number for each link in nw
f = np.ones((l,1)) #Given flows
u = np.zeros((l,1)) #link capacities
d = np.zeros((l,r))
delta = np.zeros((l,r))
od2route = np.zeros((K,max_r))
f_feas = f #vector of feasible flows satisfying capacity constraints
flag = [f>u]
for k in range(K):
    print('OD pair: %d'%k)
    r_k = (od2route[k,:]).astype(int) #vector of indices of all routes bw od pair k
    not_done = 1
    while not_done:
        flag = [f_feas>u] 
        sat_f = np.nonzero(flag)[1] #indices of saturated flows
        d = (delta>0) #binary link route incidence matrix
        #update the flag matrix of saturated flows
        for r in r_k:
            for i in sat_f:
                if d[i,r]>0:
                    d[i,r] = d[i,r]+1 #matrix indicating routes with saturated links
                    #0: link not in route, 1: unsaturated link, 2: saturated link
        sat = np.where(np.amax(d[:,r_k],axis=0)==2)[0] #saturated routes
        unsat = np.where(np.amax(d[:,r_k],axis=0)==1)[0] #unsaturated routes
        if np.size(sat)==0:#if no saturated route is present for a given od pair
            not_done = 0 
        else:#presence of at least one saturated route
            if np.size(unsat)==0:#if no unsaturated route is present for a given od pair
                print('Cannot find feasible solution for OD pair %d'%k)
            else:
                us = unsat[0] #unsaturated route to which flow is to be moved
                for r in sat: #route has saturated link/s
                    pos_sat_f = np.where(d[:,r]==2)[0] #position of saturated link/s
                    t = np.where(d[:,us]==1)[0] #position of unsaturated links in unsaturated route
                    f_feas[pos_sat_f] = u[pos_sat_f]
                    for p in pos_sat_f:
                        i = 0
                        t = t[i]
                        f_feas[t] = f_feas[t]+f_feas[p]-u[p]
                        i = i+1
                
    