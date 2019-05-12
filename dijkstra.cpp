#include <unordered_map>
#include <tuple>
#include <iostream>
#include <limits> 
#include <stdlib.h>
#include <queue>
#include <vector> 
#include<algorithm> 
#include "graph.h"
#include "dijkstra.h"

using namespace std;


typedef pair<node_type, dist_type> pair_node_dist;

vector<node_type> recover_path(unordered_map<node_type, node_type> parent, node_type origin, node_type destination){
	vector<node_type> path;
	node_type temp = destination;
	path.push_back(temp);
	while(temp != origin){
		temp = parent[temp];
		path.push_back(temp);
	}
	std::reverse(path.begin(),path.end());
	return path;
}

vector<node_type> Dijkstra::single_source_single_destination_dijkstra(node_type origin, node_type destination){
	// some test
	set<node_type> vert = vertices();
	if((vert.find(origin)==vert.end()) || (vert.find(destination) ==vert.end())){
		cout << "ERROR: Either the origin or the destination are not vertices of the graph" << endl;
		throw -1;
	}

	// The dijkstra algorithm
	int n = vert.size();
	unordered_map<node_type, node_type> parent(n); 

	priority_queue<pair_node_dist, vector<pair_node_dist>, greater<pair_node_dist>> pq;
	unordered_map<node_type, dist_type> dist;

	for(node_type i : vert){
		dist[i] = numeric_limits<dist_type>::max();
		parent[i] = -1;
	}
	pq.push(make_pair(0, origin)); 
	dist[origin] = 0;

	while (!pq.empty()) {
        node_type u = pq.top().second;
        pq.pop();

        vector<node_type>* vec = successors(u);

        // Here we can add some parallelism #pragma omp parallel
        for (auto i = vec->cbegin(); i != vec->cend(); ++i){
            node_type v = *i;
            dist_type weight = getWeight(u, v);

            if (dist[v] > dist[u] + weight) {
                dist[v] = dist[u] + weight;
                parent[v] = u;
                pq.push(make_pair(dist[v], v));
            }
            else if(v == destination){
            	return recover_path(parent, origin, destination);
            }
        }
    }
    return recover_path(parent, origin, destination);
};

unordered_map<node_type, node_type> Dijkstra::single_source_all_nodes_dijkstra(node_type origin){
	// some test
	set<node_type> vert = vertices();
	if(vert.find(origin)==vert.end()){
		cout << "ERROR: Either the origin or the destination are not vertices of the graph" << endl;
		throw -1;
	}

	// The dijkstra algorithm
	int n = vert.size();
	unordered_map<node_type, node_type> parent(n); 

	priority_queue<pair_node_dist, vector<pair_node_dist>, greater<pair_node_dist>> pq;
	unordered_map<node_type, dist_type> dist;

	for(node_type i : vert){
		dist[i] = numeric_limits<dist_type>::max();
		parent[i] = -1;
	}
	pq.push(make_pair(0, origin)); 
	dist[origin] = 0;

	while (!pq.empty()) {
        node_type u = pq.top().second;
        pq.pop();

        vector<node_type>* vec = successors(u);
        // Here we can add some parallelism #pragma omp parallel
        for (auto i = vec->cbegin(); i != vec->cend(); ++i){
            node_type v = *i;
            dist_type weight = getWeight(u, v);

            if (dist[v] > dist[u] + weight) {
                dist[v] = dist[u] + weight;
                parent[v] = u;
                pq.push(make_pair(dist[v], v));
            }
        }
    }
    return parent;
};