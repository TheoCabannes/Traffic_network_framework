
# coding: utf-8

# #### Remarks
# First, we tried to solve the STA with CVX.
# This does not work because the values of the travel time are to big. 
# 
# Secondly, we solve the STA using links flow. 
# This does not help us in our project because we need the path flow to compute the dual of the TAP-C.
# 
# Thirdly, we solve the STA using paths flow. 
# The issue here is that finding every path in a network in a NP-complete problem. 
# But if we can compute every path (we only need to do it one time), then it is really fast to compute the shortest path. We do not need to do a Dikjstra's algorithm to compute the shortest path, only a big matrix multiplication and a array sorting are enough. This might be more difficult on a laptop. But it might be faster on a HPC.
# 
# Finally, we solve the STA using links flow and we compute the paths during the all or nothing step of the Frank Wolf algorithm.

I210 = 'data/I210'
Chic = 'data/Chicago'
Anah = 'data/Anaheim'
Siou = 'data/SiouxFalls'
Brae = 'data/braess'

network_name = Siou
debug = True


# First we load the network
graph, demand = load_network(network_name)


# ## 1. We load the graph and the demand
# Both graph and demand are in csv file, we load them

import numpy as np
import scipy.sparse


# ### We clean the data

"""
We need to demand to be an array like:
[
[o1, d1, demand from o1 to d1],
...,
[on, dn, demand from on to dn]
]

We need the graph to be an array like:
[
[id link 1, node origin link 1, node destination link 1, a0, a1, a2, a3, a4],
...,
[id link n, node origin link n, node destination link n, a0, a1, a2, a3, a4]
]
where the node are indexed from 0 to nb_nodes-1,
    the links are indexed from 0 to nb_links-1,
    and a0, ..., a4 are the coeficient to calculate the travel time of one link as a function
    of the flow on the link: t = a0 + a1 * f + a2 * f**2 + a3 * f**3 + a4 * f**4

One can add some checks here to be sure that demand and graph respect this format!

Some demand and graph data can be found on:
    https://github.com/bstabler/TransportationNetworks
    Here one need to perprocess the data before using this code
"""

def cleaning_input(graph, demand):
    # in the case where there is only one o-d, then demand is interpret as a single row and not as a matrix (2d array)
    try:
        demand.shape[1]
    except:
        demand = np.array([demand])
    nb_ods = int(demand.shape[0])

    # in the case where the index of the od pairs does not begin by 0, we rename the od pairs
    first_index_od = min(np.min(graph[:,1]), np.min(graph[:,2]))
    graph[:,1] = graph[:,1]-first_index_od
    graph[:,2] = graph[:,2]-first_index_od
    demand[:,0] = demand[:,0] - first_index_od
    demand[:,1] = demand[:,1] - first_index_od
    # WE SHOULD ALSO CHECK THAT EVERY DEMAND IS WELL DEFINED: >0
    return graph, demand

def load_network(network):
    graph = np.loadtxt(network + '_net.csv', delimiter=',', skiprows=1)
    demand = np.loadtxt(network + '_od.csv', delimiter=',', skiprows=1)
    graph, demand = cleaning_input(graph, demand)
    if debug:
        print("The network " + network + " has been charged. The network characteristics are:")
        nb_links = int(np.max(graph[:,0])+1)
        nb_nodes = int(max(np.max(graph[:,1]), np.max(graph[:,2]))+1)
        nb_ods = int(demand.shape[0])
        print("\t nb nodes = " + str(nb_nodes))
        print("\t nb links = " + str(nb_links))
        print("\t nb ods = " + str(nb_ods))
    return graph, demand

debug_local = False
if debug_local:
    graph, demand = load_network(network_name)
    if network_name == Brae:
        demand[0][2] = 10


# Then, we define the function which gives the travel time as a function of the flow


## Here Max float is fixed to not overlap the maximum float number when we compute the potential function
max_float = 1e+10
float_appro = 1e-3
if debug:
    print("The max capacity and to max travel time are " + str(max_float))

"""
We define the travel time of each link as a function of the flow on each link.
This function is given by the network topology (in the graph file) and by the capacity of each link.
One row of the graph is:
[id link 1, node origin link 1, node destination link 1, a0, a1, a2, a3, a4]
Each row correspond to one link.
The travel time of the link is t = a0 + a1 * f + a2 * f**2 + a3 * f**3 + a4 * f**4 if f < c 
    It is t = + inf (max_float) is f >= c
"""
def travel_time(graph, f, c = -1):
    if c == -1:
        c = [max_float for i in range(len(f))]
    # here we need to have the same indexation for graph and f. That why we need to store graph in a dictionnary
    tt_tmp = graph[:,3] + graph[:,4]*f + graph[:,5]*(f**2) + graph[:,6]*(f**3) + graph[:,7]*(f**4)
    if len(c) != len(f):
        print("ERROR: in travel_time(graph, f, c = -1)\n******* the dimension of the capacity vector does not match the number of links *******")
    tt_tmp = [tt_tmp[i] if f[i]<c[i]-float_appro else max_float for i in range(len(f))]
    return tt_tmp


# ## 2. We compute the all or nothing flow allocation

# To use the Dijkstra's algorithm class of scipy we need to define the adjacent matrix of the graph


"""
We store the travel time of each link in a sparce adjacency matrix G.
Given a flow allocation f and a capacities constraints c, we compute 
the travel time of each link and store it in G.

This allow us after to us the Disjtra's algorithm of the class scipy.sparse.csgraph
"""

def update_travel_time_from_tt(graph, tt, nb_nodes):
    G = np.zeros(shape=(nb_nodes,nb_nodes))
    G = scipy.sparse.lil_matrix(G)
    for i in range(graph.shape[0]):
        G[int(graph[i][1]),int(graph[i][2])] = tt[i]
    G = G.tocsr()
    return G

def update_travel_time_from_flow(graph, f, nb_nodes, c = -1):
    return update_travel_time_from_tt(graph,travel_time(graph, f, c), nb_nodes)

debug_local = False
if debug_local:
    graph, _ = load_network(network_name)
    nb_links = int(np.max(graph[:,0])+1)
    nb_nodes = int(max(np.max(graph[:,1]), np.max(graph[:,2]))+1)
    G = update_travel_time_from_flow(graph, np.zeros(nb_links), nb_nodes) # , [1,1,-1,1,1]))
    print(G)

from scipy.sparse.csgraph import dijkstra


"""
Because the classic Frank Wolf algorithm does not take into account the paths,
we need to store the path in memory. We do it by using Object and a dictionnary 
called paths_used.

This idea is that we do not want to compute all the possible paths of the network
before the beginning of the algorithm. So we compute the paths used during the algorithm,
and we store them in memory.

We build at the same time the incidence matrix.

CAREFUL. We use the number of links of the graph in the calcul of the hash of every path.
This variable should be a parameter.
"""

class path:
    __nb_links = -1
    def set_nb_links(nb_l):
        path.__nb_links = nb_l
    def __init__(self,links):
        self.links = links
        self.__flow = 0
        if path.__nb_links == -1:
            raise Exception('The number of links as not be defined inside the path class. We need it for defining the hash of a link. Please use the function path.set_nb_links(nb_links) to define it')
    def set_flow(self, flow):
        self.__flow = flow
    def get_flow(self,):
        return self.__flow
    def __eq__(self, other):
        return np.all(self.links == other.links)
    # I am not sure about the hash table structure in Python. I am used to Java.
    def __hash__(self):
        return hash(np.sum([hash(self.links[i]*(path.__nb_links**i)) for i in range(len(self.links))]))
    def __str__(self):
        return str(self.__flow) + " is on " + str(self.links) 

debug_local = False
if debug_local:
    graph, _ = load_network(network_name)
    nb_links = int(np.max(graph[:,0])+1)
    path.set_nb_links(nb_links)
    paths_used = {}
    print(paths_used)
    p = path([0, 1, 2, 3])
    paths_used[hash(p)] = (p, 0)
    print(paths_used)
    p = path([0, 1, 1, 3])
    paths_used[hash(p)] = (p, 0)
    print(paths_used)
    p = path([0, 0, 2, 3])
    paths_used[hash(p)] = (p, 0)
    print(paths_used)


# Now let's compute the all or nothing allocation

"""
This cell is the main point of the entire solver.
It is the all or nothing allocation.

"""


# graph_dict gives the line of the graph matrix corresponding to the destination d and the origin o
def build_graph_adjacency(graph):
    graph_dict = {}
    for i in range(graph.shape[0]):
        try: 
            graph_dict[int(graph[i][1])]
        except:
            graph_dict[int(graph[i][1])] = {}
        graph_dict[int(graph[i][1])][int(graph[i][2])] = int(graph[i][0])
    return graph_dict

def put_flow_on_short_path(faon, o_tmp, d_tmp, flow_tmp, return_predecessors, graph_dict, paths_used, paths_used_tmp, k, i, full_dijkstra):
    node_tmp = d_tmp
    links = []
    # using the dijkstra, we build the fastest path and we put the flow on it.
    while node_tmp != o_tmp:
        if debug_local:
            print(o_tmp)
            print(node_tmp)
        if not full_dijkstra:
            node_tmp_d = return_predecessors[i][node_tmp]
        else:
            node_tmp_d = return_predecessors[o_tmp][node_tmp]
        # Here we need the graph_dict to recover the link id from the nodes id.
        link_tmp = graph_dict[node_tmp_d][node_tmp]
        # we recover the path from the predecessor of the 
        links.insert(0, link_tmp)
        faon[link_tmp] += flow_tmp
        node_tmp = node_tmp_d

    p = path(links)
    p.set_flow(flow_tmp)


    if debug_local:
        print(p)
    # If we do not already know p, we add it to the delta matrix (here a dict)
    if not hash(p) in paths_used:
        if debug_local:
            print(k)
        # we add p and his index
        paths_used[hash(p)] = (p, k)
        paths_used_tmp[hash(p)] = (p, k)
        k = k+1
    else:
        # we recover the index of p using the paths_used dict
        paths_used_tmp[hash(p)] = (p, paths_used[hash(p)][1])
    return faon, paths_used, paths_used_tmp, k


# computing the all or nothing flow
def all_or_nothing(demand, G, graph_dict, paths_used, k, full_dijkstra = True):
    # k is the next index to give to a new path.
    
    debug_local = False
    # paths_used_tmp is the set of the path in the all or nothing allocation of this step
    paths_used_tmp = {}
    # faon is the flow allocation of the all or nothing allocation
    nb_links = G.nnz
    faon = np.zeros(nb_links)
    
    # using scipy to compute dijkstra
    """
    HERE ONE CAN CHANGE THE ALGORITHM TO BE FASTER. WE MIGHT REWRITE THE DIJKSTRA ALGORITHM.
    indices : array_like or int, optional
        if specified, only compute the paths for the points at the given indices.
    """
    indices = None
    if not full_dijkstra:
        indices = [int(demand[i,0]) for i in range(len(demand))]
        
    dist_matrix, return_predecessors = dijkstra(G, return_predecessors = True, indices = indices)
    if debug_local:
        print(return_predecessors)
        print(indices)
        
    # for every origin destination pairs (i.e. one line of the demand file)
    nb_ods = int(demand.shape[0])
    for i in range(nb_ods):
        # we compute the shortest path and add the demand on it.
        o_tmp = int(demand[i][0])
        d_tmp = int(demand[i][1])
        flow_tmp = demand[i][2]
        faon, paths_used, paths_used_tmp, k = put_flow_on_short_path(faon, o_tmp, d_tmp, flow_tmp, return_predecessors, graph_dict, paths_used, paths_used_tmp, k, i, full_dijkstra)

    return faon, paths_used_tmp, paths_used, k


# We define the line search

debug_local = False
if debug_local:
    # we load the network
    graph, demand = load_network(network_name)
    nb_links = int(np.max(graph[:,0])+1)
    print(nb_links)
    path.set_nb_links(nb_links)
    nb_nodes = int(max(np.max(graph[:,1]), np.max(graph[:,2]))+1)
    G = update_travel_time_from_flow(graph, np.zeros(nb_links), nb_nodes)
    graph_dict = build_graph_adjacency(graph)
    paths_used = {}
    k = 0 
    _, _ , pp, k = all_or_nothing(demand, G, graph_dict, paths_used, k)
    for p in pp.values():
        print(p[0])
    all_or_nothing(demand, G, graph_dict, paths_used, k, False)
    for p in paths_used.values():
        print(p[0])


"""
The function potential is used to compute the line search between 
the all or nothing allocation and the current flow allocation
The function potential returns the objective function corresponding to
the flow allocation f.

The function line search does a 1D line search.
"""

def potential(graph, f, c=-1):
    # this routine is useful for doing a line search
    # computes the potential at flow assignment f
    if c == -1:
        c = [max_float for i in range(len(f))] # this might be to much, to do once we have everything done
    # here we need to have the same indexation for graph and f.
    pot_tmp = graph[:,3]*f + 1/2*graph[:,4]*(f**2) + 1/3*graph[:,5]*(f**3) + 1/4*graph[:,6]*(f**4) + 1/5*graph[:,7]*(f**5)
    pot_tmp = [pot_tmp[i] if f[i]<c[i] else f[i]*max_float for i in range(len(f))]
    return np.sum(pot_tmp)


def line_search(f, res=20):
    debug_local = False
    # on a grid of 2^res points bw 0 and 1, find global minimum
    # of continuous convex function
    # here we do a bisection
    d = 1. / (2**res - 1)
    l, r = 0, 2**res - 1
    while r - l > 1:
        if f(l * d) <= f(l * d + d):
            return l * d
        if f(r * d - d) >= f(r * d):
            return r * d
        # otherwise f(l) > f(l+d) and f(r-d) < f(r)
        m1, m2 = (l + r) / 2, 1 + (l + r) / 2
        if debug_local:
            print(l * d, end=": ")
            print(f(l * d))
            print(r * d, end=": ")
            print(f(r * d))
            print(str(m1 * d) + " = " + str(f(m1 * d)))
            print(str(m2 * d) + " = " + str(f(m2 * d)))
            print()
        if f(m1 * d) < f(m2 * d):
            r = m1
        if f(m1 * d) > f(m2 * d):
            l = m2
        if f(m1 * d) == f(m2 * d):
            return m1 * d
    return l * d


# Now let run the Frank-Wolf's algorithm with a line search to find alpha

def build_network(graph, c):
    nb_links = int(np.max(graph[:,0])+1)
    path.set_nb_links(nb_links)
    nb_nodes = int(max(np.max(graph[:,1]), np.max(graph[:,2]))+1)
    G = update_travel_time_from_flow(graph, np.zeros(nb_links), nb_nodes, c)
    graph_dict = build_graph_adjacency(graph)
    return nb_links, nb_nodes, G, graph_dict
    
def initialization_FW(demand, G, graph_dict, all_paths_used, k, graph, nb_nodes, c):
    debug_local = True
    print(G)
    f, paths_used_for_this_iter, all_paths_used, k = all_or_nothing(demand, G, graph_dict, all_paths_used, k)
    
    print(f)
    G = update_travel_time_from_flow(graph, f, nb_nodes, c)
    path_flow_matrix = np.zeros(k)
    for val in paths_used_for_this_iter.values():
        path_flow_matrix[val[1]] = val[0].get_flow()

    if debug_local:
        print("Test of initialization_FW")
        print(f)
        for p in paths_used_for_this_iter.values():
            print(p[0])
        print(path_flow_matrix)
    return f, paths_used_for_this_iter, all_paths_used, k, G, path_flow_matrix

def iteration_FW(demand, G, graph_dict, all_paths_used, k, graph, f, nb_nodes, c, path_flow_matrix, i):
    # WE COMPUTE THE ALL OR NOTHING ALGORITHM
    faon, paths_used_for_this_iter, all_paths_used, k = all_or_nothing(demand, G, graph_dict, all_paths_used, k)

    # we find the better convex combinaison of f and faon
    s = line_search(lambda a: potential(graph, (1. - a) * f + a * faon, c))
    # TO DO
    # HERE WE SHOULD BE CAREFUL IN THE CASE WHERE WE HAVE THE CAPACITY CONSTRAINTS
    # THE GRADIENT OF THE FUNCTION IS NOT THE SHORTEST PATH
    # WE SHOULD REMOVE THE PATH THAT SATURATED THE LINKS
    # AND THEN COMPUTE THE SOLUTION WITHOUT THE SATURATION, AND WITHOUT THE 
    # CORRESPONDING DEMAND
    f = (1. - s) * f + s * faon
    G = update_travel_time_from_flow(graph, f, nb_nodes, c)

    # we multiply the previous path flow matrix by the coeficient of the line search
    path_flow_matrix = (1-s) * path_flow_matrix
    # we add the new path at the end of the path flow matrix.
    path_flow_matrix = np.append(path_flow_matrix, np.zeros(len(all_paths_used)-path_flow_matrix.shape[0]))
    # we add the all or nothing path flow (mulitply by the coeficient of the line search) to the path flow matrix.
    for val in paths_used_for_this_iter.values():
        # val[1] is the index of the path, val[0] is the path object
        path_flow_matrix[val[1]] += s * val[0].get_flow()

    if debug and i % (nb_iter / 10) == 0:
        print("Iteration: " + str(i))
        print("s: " + str(s))
        print("The paths used at this iteration are:")
        for p in paths_used_for_this_iter.values():
            print(p[0])
    return path_flow_matrix, f, G, all_paths_used, k, s

def output_FW(all_paths_used, nb_links, graph, f, demand):
    # MAYBE I SHOULD SPLIT THE CODE HERE TO MAKE EASIER TEST
    # At the end we do some work to return a proper output
    nb_paths = len(all_paths_used.keys())
    delta = np.zeros(shape=(nb_paths, nb_links)) # nb_links should be a parameter
    route2od = [0 for _ in range(nb_paths)]# np.zeros(shape=nb_paths)
    delta = scipy.sparse.lil_matrix(delta)
    tt_f = np.array(travel_time(graph, f))

    for p in all_paths_used.values():
        # here I can built route2od matrix at the same time
        try:
            links_tmp = p[0].links
            for l in links_tmp:
                delta[p[1],l] = 1
            route2od[p[1]] = (int(graph[int(links_tmp[0])][1]), int(graph[int(links_tmp[-1])][2]))
        except:
            ("")
    delta = delta.tocsr()
    # I should built the route2od matrix here,
    # route2od should be build from the all_paths_used_dict
    return tt_f, delta, route2od

def Frank_Wolf_solver(graph, demand, eps, nb_iter, c=-1):
    ######### FIRST, WE INITIALIZE THE ALGORITHM #########
    # We initialize the number of paths_used to 0, 
    k = 0
    all_paths_used = {}
    
    # we built the network as a matrix from the graph file
    nb_links, nb_nodes, G, graph_dict = build_network(graph, c)
    
    # The initialization step: we put all the demand on the fastest free flow travel time paths.
    f, paths_used_for_this_iter, all_paths_used, k, G, path_flow_matrix = initialization_FW(demand, G, graph_dict, all_paths_used, k, graph, nb_nodes, c)

    ######### THEN, I RUN THE ITERATION OF THE FRANK-WOLF ALGORITHM #########
    for i in range(nb_iter):
        path_flow_matrix, f, G, all_paths_used, k, s = iteration_FW(demand, G, graph_dict, all_paths_used, k, graph, f, nb_nodes, c, path_flow_matrix, i)
        if s < eps:
            break
    if debug:
        print(path_flow_matrix)
       
    ######### FINALLY, I WORK ON THE OUTPUT TO RETURN #########
    tt_f, delta, route2od = output_FW(all_paths_used, nb_links, graph, f, demand)
    return path_flow_matrix, tt_f, delta, route2od

eps=1e-8
nb_iter = 1000
graph, demand = load_network(network_name)
if network_name == Brae:
    demand[0][2] = 10
        
path_flow_matrix, tt_f, delta, route2od = Frank_Wolf_solver(graph, demand, eps, nb_iter) #, [11,11,2,11,11])

if debug:
    nb_paths = len(path_flow_matrix)
    print(nb_paths)
    print(tt_f)
    print(delta)
    print(delta @ tt_f)
    print(route2od)
    print(delta.shape)


def check_wardrop(j, demand, route2od, delta, path_flow_matrix, tt_f):
    f = delta.T @ path_flow_matrix
    tt_p = delta @ tt_f
    od = (int(demand[j][0]), int(demand[j][1]))
    tab = []
    for i in range(len(route2od)):
        if route2od[i] == od:
            tab.append(i)
    print(od)
    print(tab)
    print(tt_p[tab])
# check_wardrop(231, demand, route2od, delta, path_flow_matrix, tt_f)

