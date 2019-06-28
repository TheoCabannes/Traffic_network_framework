#test
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
    print(line)
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

def neighbours(node, adj):
    """
    Parameter:
        node: node for which to find neigbours for
        adj: adjacency matrix that describes network
    Return:
        List of nodes that are neighbours to node
    """
    return np.nonzero(adj[node])[0].tolist()

def add_flow_dijkstra(adj, src, target, faon, flow, g):
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
    q = [i for i in range(len(adj))] # queue of nodes
    dist = [np.inf for i in range(len(adj))] # node distances
    prev = [None for i in range(len(adj))] # predecessor in shortest path
    
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
        
        for v in neighbours(u, adj):
            alt = dist[u] + adj[u][v]
            if alt < dist[v]:
                dist[v] = alt
                prev[v] = u