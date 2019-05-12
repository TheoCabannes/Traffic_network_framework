#include <unordered_map>
#include <tuple>
#include <iostream>
#include <stdlib.h>
#include <vector> 
#include "graph.h"
#include <set>
#include "omp.h"

using namespace std;

#define DEBUG_LOCAL 1

// We can add a boolean variable to never do two times the same calculus: memoization

unordered_map<node_type, vector<node_type>> adjacency_list;
unordered_map<node_type, unordered_map<node_type, dist_type>> sparse_adjacency_matrix;

unordered_map<node_type, omp_lock_t> writelock;

set<node_type> vertices_set;

vector<node_type>* Graph::successors(node_type i){
	return &(adjacency_list[i]);
};

// Thread safe addEdge
void Graph::addEdge(node_type i, node_type j){
	if(sparse_adjacency_matrix[i][j]!=0 & sparse_adjacency_matrix[i][j]!=-1){
		#if DEBUG_LOCAL
			cout << "There is already a link between " << i << " and " << j << " with length=" << sparse_adjacency_matrix[i][j] << endl;
		#endif
		return;
	}

	if(!vertices_set.count(i))
		omp_init_lock(& writelock[i]);
	if(!vertices_set.count(j))
		omp_init_lock(& writelock[j]);
	omp_set_lock(&writelock[i]);
	
	vertices_set.insert(i);
	vertices_set.insert(j);

	adjacency_list[i].push_back(j);

	omp_unset_lock(&writelock[i]);
};

// Thread safe setWeight
void Graph::setWeight(node_type i, node_type j, dist_type dist){
	// Maybe we can check that the edge i->j exists
	omp_set_lock(&writelock[i]);
	sparse_adjacency_matrix[i][j] = dist;
	omp_unset_lock(&writelock[i]);
};

// Thread safe get
dist_type Graph::getWeight(node_type i, node_type j){
	try{
		if(sparse_adjacency_matrix.at(i).count(j)==0){
			throw -1;
		}
		else{
			return sparse_adjacency_matrix.at(i)[j];
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
