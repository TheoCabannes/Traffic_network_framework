#include <unordered_map>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <vector> 
#include "graph.h"
#include <set>

using namespace std;

#define DEBUG_LOCAL 1

unordered_map<node_type, vector<node_type>> adjacency_list;
unordered_map<node_type, unordered_map<node_type, dist_type>> sparse_adjacency_matrix;

set<node_type> vertices_set;

vector<node_type>* Graph::successors(node_type i){
	return &(adjacency_list[i]);
};

void Graph::addEdge(node_type i, node_type j){
	if(sparse_adjacency_matrix[i][j] != 0){
		#if DEBUG_LOCAL
			cout << "There is already a link between " << i << " and " << j << " with length=" << sparse_adjacency_matrix[i][j] << endl;
		#endif
		return;
	}
	
	vertices_set.insert(i);
	vertices_set.insert(j);

	adjacency_list[i].push_back(j);
};

void Graph::setWeight(node_type i, node_type j, dist_type dist){
	// Maybe we can check that the edge i->j exists
	sparse_adjacency_matrix[i][j] = dist;
};

dist_type Graph::getWeight(node_type i, node_type j){
	try{
		if(sparse_adjacency_matrix[i][j]==0){
			throw -1;
		}
		else{
			return sparse_adjacency_matrix[i][j];
		}
	}
	catch(int e){
		#if DEBUG_LOCAL
			cout << "There is no links between " << i << " and " << j << endl;
		#endif
		return e;
	}
}

set<node_type> Graph::vertices(){
    return vertices_set;
}
