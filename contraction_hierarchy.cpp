#include <cstdlib>
#include <ctime>
#include <memory>
#include <iostream>
#include <limits>
#include <queue>
#include <list>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <set>
#include "graph.h"

#define OMP 0
#define DEBUG 1

#if OMP
    #include <omp.h>
#endif

using namespace std;

int int_max = numeric_limits<int>::max();

// WE ASSUME THAT THE GRAPH IS NOT DIRECTED
unordered_map<node_type, int>* generate_ch_nodes(Graph& g) {
    // we have all the vertices in a set
    set<node_type> vertices = g.vertices();
    int n = vertices.size();

    // contraction hierarchy order
    unordered_map<node_type, int>* heights = new unordered_map<node_type, int>();
    
    int i = 0;
    // for every node ...
    for(node_type it : vertices){
        // ... if they only have one edge we contract them. 
        if(g.successors(it)->size() <= 1){
            (*heights)[it] = i;
        }
        i = i+1;
    }

    int counter = 2 * n;

    // do the contraction until break
    while (true) {
        int most_edges_node = -1;
        int most_edges_num = int_max;

        for (node_type it : vertices) {
            if (!(*heights)[it] && most_edges_num > g.successors(it)->size()) {
                most_edges_node = it;
                most_edges_num = g.successors(it)->size();
            }
        }

        // if there is no more nodes to contract break
        if (most_edges_node == -1)
            break;

        (*heights)[most_edges_node] = counter++;

        #if DEBUG
            cout << "Contracting node " << most_edges_node << " with hierarchy " << counter;
            cout << " and edges " << g.successors(most_edges_node)->size() << endl;
        #endif

        // TO DO? having a thread safe vector structure here (cuckoohash_map from libcuckoo)
        vector<node_type> *neighbors = g.successors(most_edges_node);


        // collapse(2) means that 2 for loops imbricated used openMP 
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < most_edges_num; i++) {
            for (int j = 0; j < most_edges_num; j++) {
                node_type begin_val = neighbors->at(i);
                node_type end_val = neighbors->at(j);

                double begin_len = g.getWeight(most_edges_node, begin_val);
                double end_len = g.getWeight(most_edges_node, end_val);
                if (begin_val >= end_val || (*heights)[begin_val] != 0 || (*heights)[end_val] != 0)
                    continue;

                double new_distance = begin_len + end_len;

                #if DEBUG
                    cout << "\t Examining " << begin_val << " " << end_val << endl;
                #endif
                
                double cur_length = g.getWeight(begin_val, end_val);
                if (cur_length == -1){
                    g.addEdgeWeighted(begin_val, end_val, new_distance);
                    g.addEdgeWeighted(end_val, begin_val, new_distance);
                } else if (cur_length > begin_len + end_len) {
                    g.setWeight(begin_val, end_val, new_distance);
                    g.setWeight(end_val, begin_val, new_distance);
                }
            }
        }
    }

    return heights;
}

