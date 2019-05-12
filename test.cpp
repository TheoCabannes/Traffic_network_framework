#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include "graph.h"
#include "dijkstra.h"
#include "road_network.h"
#include "contraction_hierarchy.h"


using namespace std;

extern unordered_map<node_type, int>* generate_ch_nodes(Graph& g);
class NetworkDikjstra : virtual public Dijkstra, virtual public RoadNetwork {
    /** The Road network class that use Dijkstra as a shortest path algorithm
    *   The inheritances are virtual in order to have a "diamond" inheriance to Graph
    */
};

void test_network_1(){
    // Graph graph;
    Dijkstra graph;
    graph.addEdgeWeighted(1, 2, 1.0);
    graph.addEdgeWeighted(1, 3, 2.0);
    graph.addEdgeWeighted(1, 4, 3.0);
    graph.addEdgeWeighted(2, 5, 4.0);
    graph.addEdgeWeighted(3, 6, 2.0);
    graph.addEdgeWeighted(4, 6, 3.0);
    graph.addEdgeWeighted(4, 7, 3.0);
    graph.addEdgeWeighted(5, 8, 1.0);
    graph.addEdgeWeighted(6, 8, 3.0);
    graph.addEdgeWeighted(7, 8, 4.0);

    /*
    The graph is the following one:
    
       _______2_______5_______
      /  1.0     4.0     1.0  \
    1/________3_______6________\8
     \   2.0     3.0     3.0   /
      \_______4_______7_______/
         3.0     3.0     4.0

    */
    cout << "The graph is the following one:\n\n    _______2_______5______ \n  /  1.0     4.0      1.0 \\ \n1/________3_______6________\\8\n \\   2.0     3.0     3.0   /\n  \\_______4_______7_______/\n     3.0     3.0     4.0\n"; 


    vector<node_type> path = graph.single_source_single_destination_dijkstra(1, 8);
    cout << "The shortest path between 1 and 8 is: ";
    for(auto i : path){
        cout << i << " ";
    }
    cout << endl;

    cout << "The distance is: " << graph.distance_dijkstra(1, 8) << endl;
}

void test_network_2(){
    string filename = "data/ChicagoRegional";
    
    NetworkDikjstra chicago;
    chicago.read_from_tntp(filename);

    clock_t start;
    double duration;

    start = clock();

    chicago.single_source_all_nodes_dijkstra(1);

    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Running: "<< duration << endl;

    vector<node_type> path = chicago.single_source_single_destination_dijkstra(1, 12980);

    for(auto i : path){
        cout << i << " ";
    }
    cout << endl;

    cout << chicago.distance_dijkstra(1, 12980) << endl;
}

void test_network_3(){
    string filename = "data/LA";
    
    NetworkDikjstra losAngeles;
    losAngeles.read_from_csv(filename);

    // return;
    clock_t start;
    double duration;

    start = clock();

    losAngeles.single_source_all_nodes_dijkstra(1);

    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Running: "<< duration << endl;

    cout << losAngeles.vertices().size() << endl;
    vector<node_type> path = losAngeles.single_source_single_destination_dijkstra(1, 12980);

    for(auto i : path){
        cout << i << " ";
    }
    cout << endl;

    cout << losAngeles.distance_dijkstra(1, 12980) << endl;
}

void test_network_4(){
    string filename = "data/ChicagoRegional";
    
    NetworkDikjstra chicago;
    chicago.read_from_tntp(filename);
    set<node_type> all_nodes = chicago.vertices();
    clock_t start;
    double duration;

    cout << all_nodes.size() << endl;
    start = clock();

    for(auto i : all_nodes){
        chicago.single_source_all_nodes_dijkstra(i);

        // duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
        // cout << "Running: "<< duration << " for node " << i << endl;
        if(i>100)
            break;
    }
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Running: "<< duration << endl;

    vector<node_type> path = chicago.single_source_single_destination_dijkstra(1, 12980);

    for(auto i : path){
        cout << i << " ";
    }
    cout << endl;

    cout << chicago.distance_dijkstra(1, 12980) << endl;
}

void test_network_5(){
    string filename = "data/LA";
    
    NetworkDikjstra losAngeles;
    losAngeles.read_from_csv(filename);

    set<node_type> origins = losAngeles.getOrigins();

    int total_od_pairs = 0;
    double total_flow = 0;
    cout << "Total flow demand: " << total_flow << endl;


    for (node_type o : origins){
        vector<node_type> dest = losAngeles.getDestination(o);
        total_od_pairs += dest.size();
        for(node_type odest : dest){
            total_flow += losAngeles.getDemand(o, odest);
        }
    }
    cout << "Total number of origins: " << origins.size() << endl;
    cout << "Total number of origin-destination pairs: " << total_od_pairs << endl;
    cout << "Total flow demand: " << total_flow << endl;
    cout << "Demand between nodes 1 and 635: " << losAngeles.getDemand(1, 635) << endl;
    // return;
    clock_t start;
    double duration;

    start = clock();

    losAngeles.single_source_all_nodes_dijkstra(1);

    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Running: "<< duration << endl;

    cout << losAngeles.vertices().size() << endl;
    vector<node_type> path = losAngeles.single_source_single_destination_dijkstra(1, 12980);

    for(auto i : path){
        cout << i << " ";
    }
    cout << endl;

    cout << losAngeles.distance_dijkstra(1, 12980) << endl;
}

void test_network_6(){
    string filename = "data/ChicagoRegional";
    
    NetworkDikjstra chicago;
    chicago.read_from_tntp(filename);

    ifstream data;
    string line, temp;
    data.open(filename + "_flow.tntp");
    getline(data, temp);
    cout << temp << endl;
    link_type id = 0;
    while(!data.eof()) {
        getline(data, line);
        stringstream buffer(line);
        cout << line << endl;

        getline(buffer, temp, '\t');
        getline(buffer, temp, '\t');
        int tail = strtol(temp.c_str(), NULL, 10);
        getline(buffer, temp, '\t');
        int end = strtol(temp.c_str(), NULL, 10);
        
        getline(buffer, temp, '\t');
        double volume = strtol(temp.c_str(), NULL, 10);
        getline(buffer, temp, '\t');
        double cost = strtol(temp.c_str(), NULL, 10);

        cout << volume << endl;
        chicago.update_flow(id, volume);
        chicago.print_link(id);
        id = id + 1;
        break;
    }
    data.close();

    chicago.update_all_weights();
    cout << "Travel time calculated " << chicago.getWeight(1, 10293) << endl;

    // cout << chicago.distance_dijkstra(1, 12980) << endl;
}

int test_network_7() {
    string filename = "data/ChicagoRegional";
    NetworkDikjstra chicago;
    chicago.read_from_tntp(filename);
    
    cout << "Start..." <<endl;
    clock_t start;
    double duration;

    start = clock();
    unordered_map<node_type, int>* heights = generate_ch_nodes(chicago);
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Preprocessing: "<< duration << endl;
}

int test_network_8() {
    Dijkstra graph;
    graph.addEdgeWeighted(1, 2, 1.0);
    graph.addEdgeWeighted(1, 3, 2.0);
    graph.addEdgeWeighted(1, 4, 3.0);
    graph.addEdgeWeighted(2, 5, 4.0);
    graph.addEdgeWeighted(3, 6, 2.0);
    graph.addEdgeWeighted(4, 6, 3.0);
    graph.addEdgeWeighted(4, 7, 3.0);
    graph.addEdgeWeighted(5, 8, 1.0);
    graph.addEdgeWeighted(6, 8, 3.0);
    graph.addEdgeWeighted(7, 8, 4.0);
    

    clock_t start;
    double duration;

    start = clock();
    unordered_map<node_type, int>* heights = generate_ch_nodes(graph);
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Preprocessing: "<< duration << endl;
}

int main(int argc, char const *argv[])
{
    test_network_7();
    return 0;
}