import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt

def read_metadata(file):
    """
    Parameter: 
        file: a string which represent the path to the data file tntp.
    Return:
        The function read the first metadata lines (in the <> tag) and return the numbers of lines read.
    """
    inputfile = open(file)
    line = inputfile.readline()
    line_skip = 0
    # parsing metadata
    while "<" in line:
        line_skip += 1
        # print(line, end="")
        line = inputfile.readline()
    """   
    # skip all the blank line
    while len(line)<2:
        line_skip += 1
        line = inputfile.readline()
    """ 
    return line_skip
    # print(line)
    # return line, line_skip

def read_data(file, line_skip):
    """
    Parameter: 
        file: a string which represent the path to the data file tntp
        line_skip: the number of lines of metadata to skip (should be the output of read_metadata(file))
    Return:
        The function read the tntp file using pandas.
        It drops the columns ~ and ;
        Then it return the legend and the data in a np.array
    """
    table = pd.read_csv(file, sep='\t', header=line_skip)
    
    # remove columns with ~ or ;
    to_remove = set({})
    for i in range(table.columns.shape[0]):
        if '~' in table.columns[i] or ';' in table.columns[i]:
            # print(i)
            to_remove.add(i)
    table = table.drop(table.columns[list(to_remove)], axis=1)
    
    # HERE we should do some data type change (init node, term node and power should be 32 bit integer)
    legend = np.array(table.columns)
    table_tmp = np.array(table)
    del table
    # get the legend and the data with np
    return legend, table_tmp

def init_travel_time_function(network):
    """
    Parameter: 
        network: a string which represent the path to the data file tntp
    Return:
        The table that can be used to compute the travel time of every links.
        The index of the free flow travel time, of the B, of the capacity and of the power columns.
        
        The travel time of a link (one row of the table) is t(f) = t0 * (1 + B*(f/capacity)**power)
    """
    file_type = "net"
    suffix = "_" + file_type + ".tntp"
    line_skip = read_metadata(network + suffix)
    legend_net, table_net = read_data(network + suffix, line_skip)
    
    j = 0
    for legend_name in legend_net:
        # print(legend_name)
        if "B" in legend_name:
            index_B = j
        if "Power" in legend_name:
            index_power = j
        if "Free Flow Time" in legend_name:
            index_fft = j
        if "Capacity" in legend_name:
            index_capacity = j
        if ("Init" in legend_name) or ("Tail" in legend_name):
            index_init = j
        if ("Term" in legend_name) or ("Head" in legend_name):
            index_term = j
        j += 1
    graph_wrapped = table_net, index_fft, index_B, index_capacity, index_power, index_init, index_term
    return graph_wrapped

def flow_cost_solution(network):
    """
    Read and return the flow and the cost of the STA solution using the flow.tntp file
    """
    file_type = "flow"
    suffix = "_" + file_type + ".tntp"
    line_skip = read_metadata(network + suffix)
    legend_flow, table_flow = read_data(network + suffix, line_skip)
    j = 0
    for legend_name in legend_flow:
        if "Volume" in legend_name:
            index_flow = j
        if "Cost" in legend_name:
            index_cost = j
        j += 1

    flow = table_flow[:,index_flow]
    cost_solution = table_flow[:,index_cost]
    return flow, cost_solution

def read_demand(network):
    """
    Parameter: 
        file: a file name that encode the demand under a TNTP format (network_trips.tntp).
    Return:
        A dictionary of keys based sparse matrix encoding the demand.
    """
    file_type = "trips"
    file = network + "_" + file_type + ".tntp"
    inputfile = open(file)
    line = inputfile.readline()
    # parsing metadata
    while "<" in line:
        print(line, end="")
        line = inputfile.readline()

    # skip all the blank line
    while len(line)<2:
        line = inputfile.readline()

    # read the Origin
    assert 'Origin' in line
    data = {}
    # print(line, end="")
    origin = None
    try:
        origin = int(line.split('\t')[1])
    except:
        origin = int(line.split(' ')[1])
    data[origin] = {}
    for line in inputfile:
        # print(line, end="")
        if 'Origin ' in line:
            origin = None
            try:
                origin = int(line.split('\t')[1])
            except:
                origin = int(line.split(' ')[1])
            data[origin] = {}
        # else read the destination
        else:
            dest_array = line.split(";")
            for dest in dest_array[:-1]:
                dest_tmp = dest.split(":")
                data[origin][int(dest_tmp[0])] = float(dest_tmp[1])

    d = pd.DataFrame(data)
    del data

    demand = sparse.csr_matrix(np.array(d))
    del d
    # print(demand)
    return demand.todok()

def get_graph(graph_wrapped):
    """
    Parameter:
        graph_wrapped: should be = (table_net, index_fft, index_B, index_capacity, index_power, index_init, index_term)
    Return:
        A dictionary encoding the graph: graph[origin][destination] = link_id
    """
    table_net, index_fft, index_B, index_capacity, index_power, index_init, index_term = graph_wrapped
    
    graph = {}
    for i in range(table_net.shape[0]):
        link_index = i
        node_init = int(table_net[i][index_init])
        node_term = int(table_net[i][index_term])
        
        if node_init not in graph:
            graph[node_init] = {}
        # if node_term not in graph[node_init]:
        #     graph[node_init][node_term] = {}
        # we assume that there is only one possible directed link between two nodes
        graph[node_init][node_term] = link_index
    return graph

def neighbours(node, g):
    """
    Parameter:
        node: node for which to find neigbours for
        adj: adjacency matrix that describes network
    Return:
        List of nodes that are neighbours to node
    """
    return [k-1 for k in list(g[node+1].keys())]

def add_flow_dijkstra(src, target, faon, flow, g, tt):
    """
    Parameter:
        adj: adjacency matrix that describes network
        src: starting node for OD pair
        target: end node for OD pair
        faon: current all or nothing flow allocation that the algorithm is building
        flow: demand of the current of OD pair
        g: graph dictionary, g[node_init] = {node_term_1: link_from_init_to_term_1, node_term_2: link_from_init_to_term_2}
    Return:
        Returns all or nothing flow allocation, where flow is allocated by shortest path for OD pair (as determined by Dijkstra's algorithm)
    """
    q = [i for i in range(len(g))] # queue of nodes
    dist = [np.inf for i in range(len(g))] # node distances
    prev = [None for i in range(len(g))] # predecessor in shortest path
    
    dist[src] = 0 # set source node's distance to itself as zero
    
    while len(q) != 0:
        u = min(q, key=lambda n:dist[n])
        q.remove(u)
        
        if u == target: # when target is reach, use predecessors to allocate all or nothing flow to links in shortest path
            while u != src:
                prev_node = prev[u]
                curr_link = g[prev_node+1][u+1]
                faon[curr_link] += flow
                u = prev_node
            return faon
        
        for v in neighbours(u, g):
            alt = dist[u] + tt[g[u+1][v+1]]
            if alt < dist[v]:
                dist[v] = alt
                prev[v] = u

# all or nothing function
def all_or_nothing(d, tt, g):
    """
    Parameter:
        d: a demand as sparse matrix, d[origin][dest] = demand between origin and dest
        tt: travel time on every links, tt[l] = travel time of the link l
        g: graph dictionary, g[node_init] = {node_term_1: link_from_init_to_term_1, node_term_2: link_from_init_to_term_2}
    Return:
        The all or nothing flow allocation corresponding to the demand d and the travel time tt.
    """
    nb_links = tt.shape
    faon = np.zeros(nb_links)

    for o_tmp, d_tmp in d.keys():
        #sp = find_sp(o_tmp, d_tmp)
        #faon = add_flow(faon, sp, d[o_tmp, d_tmp])
#         faon = add_flow_dijkstra(o_tmp, d_tmp, faon, d[o_tmp, d_tmp], g, tt) # TO DO: Pass g and tt instead of adj
        faon = add_flow_dijkstra(o_tmp, d_tmp, faon, d[o_tmp, d_tmp], g, tt)
        # flow_tmp = d[o_tmp, d_tmp]
        # faon = put_on_shortest_path(faon, o_tmp, d_tmp, flow_tmp, tt, g, G) # to change when we write our own 
    return faon

#def find_sp(g):

def potential(f, graph_wrapped):
    """
    Parameter:
        flow: a vector that represents the flow of every links
        graph_wrapped: should be = (table_net, index_fft, index_B, index_capacity, index_power, index_init, index_term)
    Return:
        The potential function tha t we want to minimize.

        The potential function is (one row of the table) is sum(int(t_e(s), s in [0, f_e]), e) = sum(t0 * (f + B*(f/capacity)**(power+1)/power)).    
        This function is useful for doing a line search, it computes the potential at flow assignment f.
    """
    table_net, index_fft, index_B, index_capacity, index_power, index_init, index_term = graph_wrapped
    pot_tmp = table_net[:,index_fft] * f * (1 + (table_net[:,index_B] / table_net[:,index_power])*(f/table_net[:,index_capacity])**table_net[:,index_power])
    return np.sum(pot_tmp)

def line_search(f, res=20):
    """
    Parameter:
        f: the function that we want to minimize between 0 and 1
        res: the resolution of the grid ([0, 1] is divided in 2^res points).
    Return:
        On a grid of 2^res points bw 0 and 1, find global minimum for a continuous convex function using bisection 
    """
    d = 1. / (2**res - 1) # d: Size of interval between adjacent points
    l, r = 0, 2**res - 1 # l: first point, r: last point
    while r - l > 1:
        if f(l * d) <= f(l * d + d): # Interval l: [l*d, l*d + d]
            return l * d
        if f(r * d - d) >= f(r * d):
            return r * d
        # otherwise f(l) > f(l+d) and f(r-d) < f(r)
        m1, m2 = (l + r) / 2, 1 + (l + r) / 2
        if f(m1 * d) < f(m2 * d):
            r = m1
        if f(m1 * d) > f(m2 * d):
            l = m2
        if f(m1 * d) == f(m2 * d):
            return m1 * d
    return l * d