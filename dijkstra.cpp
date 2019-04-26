#include <memory>
#include <iostream>
#include <tuple>
#include <limits>
#include <queue>
#include <unordered_map>
#include "parse_csv.h"
#include "graph.h"
#include "dijkstra.h"

typedef std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> graph_type;

std::unique_ptr<dist_mat> dijkstra(int origin, std::unique_ptr<graph>& G) {
    using std::unordered_map;
    using std::vector;
    using std::priority_queue;
    using std::unique_ptr;
    using std::unordered_map;
    using std::cout;
    using std::endl;

    // hard coded -1, fix later
    auto graph = G.get()->graph;
    priority_queue<dist_mat, vector<dist_mat>, std::greater<dist_mat>> pq;

    vector<double> dist(graph.size(), std::numeric_limits<double>::max());

    pq.push(std::make_pair(0.0, origin));
    dist[origin - 1] = 0.0;

    while (!pq.empty()) {
        int u = pq.top().second;
        cout << u << endl;
        pq.pop();

        for (auto neighbour : graph.at(u)) {
            int v = (int) neighbour.first;
            double weight = graph.at(u).at(v);
            if (dist[v - 1] > dist[u - 1] + weight) {
                dist[v - 1] = dist[u - 1] + weight;
                pq.push(std::make_pair(dist[v - 1], v));
            }
        }
    }

    unique_ptr<dist_mat> distances((dist_mat*) malloc(sizeof(dist_mat) * dist.size()));
    cout << "Hi" << endl;
    for (int i = 0; i < dist.size(); i++) {
        *(distances.get() + i) = std::make_pair(dist[i], i + 1); 
    }
    
    for (int i = 0; i < dist.size(); i++) {
        cout << (distances.get() + i)->second;
        cout << " " << endl;
        cout << (distances.get() + i)->first << endl;
    }

    return distances;
}
