import numpy as np
f = np.matrix([3.5,1,2.5,1,3.5])
f = f.reshape((5,1))
u = np.matrix([10.0,10.0,2.0,10.0,10.0])
u = u.reshape((5,1))
delta = np.matrix([[1,1,0],[0,0,1],[0,1,0],[1,0,0],[0,1,1]])
od2route = np.array([0,1,2]).reshape((1,3))
#def find_feasible(f,u,delta,od2route):

#import numpy as np
l = np.shape(f)[0] #number of links in network
K = np.shape(od2route)[0] #number of od pairs in network
r = np.shape(delta)[1] #number of routes in network
d = np.zeros((l,r))
f_feas = f #vector of feasible flows satisfying capacity constraints
for k in range(K):
    print('OD pair: %d'%(k+1))
    r_k = (od2route[k,:]).astype(int) #vector of indices of all routes bw od pair k
    flag1 = [f_feas>u] 
    sat_f = np.nonzero(flag1)[1] #indices of saturated flows
    d = (delta>0) #binary link route incidence matrix
    d = d*1
    #update the incidence matrix to find saturated flows
    
    for R in r_k:
        for i in sat_f:
            if d[i,R]>0:
                d[i,R] = d[i,R]+1 #matrix indicating routes with saturated links
                #0: link not in route, 1: unsaturated link, 2: saturated link
    sat = np.where(np.max(d[:,r_k],axis=0)==2)[1] #saturated routes
    unsat = np.where(np.max(d[:,r_k],axis=0)==1)[1] #unsaturated routes
    if np.size(sat)==0:#if no saturated route is present for a given od pair
        print('No saturated route for OD pair: %d'%(k+1))
        
    else:#presence of at least one saturated route
        f_r = 1.0
        for r in sat:#For each saturated route
            if (np.max(d[:,r],axis=0)==2)&(f_r>0):
                l_r = np.where(d[:,r]>0)[0] #links in saturated route r
                f_r = max(f_feas[l_r]-u[l_r]) #flow to be rerouted from saturated route r
                flag = 1
                while flag: # as long as the route r remains saturated
                    unsat = np.where(np.max(d[:,r_k],axis=0)==1)[1] #unsaturated routes
                    cap = np.zeros(np.size(unsat))
                    c = 0
                    for us in unsat: #compute extra capacity for unsaturated routes
                        l_u = np.where(d[:,us]>0)[0] #links in unsaturated route u
                        cap[c] = min(u[l_u]-f_feas[l_u]) #extra capacity of unsaturated route u
                        c = c+1
                    i_unsat_sort = np.flip(np.argsort(cap),0) #unsaturated routes sorted in decreasing order of capacities
                    cap_sort = cap[i_unsat_sort] 
                    if f_r<=cap_sort[0]: #all flow to be rerouted can be put in one unsat route
                        f_feas[np.where(d[:,r]>0)[0]] = f_feas[np.where(d[:,r]>0)[0]]-f_r #remove f_r from sat route r
                        f_feas[np.where(d[:,unsat[i_unsat_sort[0]]]>0)[0]] = f_feas[np.where(d[:,unsat[i_unsat_sort[0]]]>0)[0]]+f_r #move f_r to route with highest capacity
                        gl = np.where(d[:,r]==2)[0]
                        d[np.where(d[:,r]>0)[0],r] = 1
                        d[gl,np.where(d[gl,:]>0)[1]] = 1
                        flag = 0
                    else: #split flow to be rerouted
                        f_feas[np.where(d[:,r]>0)[0]] = f_feas[np.where(d[:,r]>0)[0]]-cap_sort[0] #remove cap_sort[0] from sat route r
                        f_feas[np.where(d[:,unsat[i_unsat_sort[0]]]>0)[0]] = f_feas[np.where(d[:,unsat[i_unsat_sort[0]]]>0)[0]]+cap_sort[0] #move f_r to route with highest capacity
                        f_r = f_r-cap_sort[0]
                        flag = 1
    #return f_feas