#test
import numpy as np
import pandas as pd

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
        if "Init" in legend_name:
            index_init = j
        if "Term" in legend_name:
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
    origin = int(line.split('\t')[1])
    data[origin] = {}
    for line in inputfile:
        # print(line, end="")
        if 'Origin ' in line:
            origin = int(line.split('\t')[1])
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