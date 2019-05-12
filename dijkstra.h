#include <unordered_map>
#include <tuple>
#include <vector>

#ifndef GRAPH_H
#define GRAPH_H

using namespace std;

class Dijkstra : virtual public Graph{
    public:
        // add initialization? (for CH)

    	unordered_map<node_type, node_type> single_source_all_nodes_dijkstra(node_type origin); //!< Return a map where map[final node] = previous node to go to final node from origin

    	vector<node_type> single_source_single_destination_dijkstra(node_type origin, node_type destination); //!< Return a vector of nodes that forms the path from origin to destination

        unordered_map<node_type, node_type> single_source_multiple_destinations_dijkstra(node_type origin, set<node_type> destination); // TODO

        dist_type path_distance(vector<node_type> path){
            /** Given a sequence of nodes that forms a path, path_distance computes the distance of this path.
            *
            */
            dist_type dist = 0;

            auto i = path.cbegin();
            node_type t = *i;
            for (i = i+1;i != path.cend(); ++i){
                node_type t_1 = *i;
                dist += getWeight(t, t_1);
                t = *i;
            }
            return dist;
        }

    	dist_type distance_dijkstra(node_type origin, node_type destination){
            /** Return the distance computed with dijkstra between the origin and the destination
            *   
            */
    		return path_distance(single_source_single_destination_dijkstra(origin, destination));
    	}
};

#endif