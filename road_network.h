#include <unordered_map>
#include <vector>
#include <set>
#include <cmath>
#include <set>
#include <iostream>

// #ifndef GRAPH_H
// #define GRAPH_H

using namespace std;
typedef unsigned int link_type; 

class Link{
    public:
        link_type id;
        node_type head, tail;
        double capacity, free_flow_time;

        double travel_time(double flow){
            cout << "Free flow travel time is: ";
            cout << free_flow_time;
            cout << " Flow: ";
            cout << flow;
            cout << "; cap: ";
            cout << capacity;
            cout << " power: ";
            cout << power;
            cout << ", B: ";
            cout << B;
            cout << ", (flow/capacity): ";
            cout << (flow/capacity);
            cout << ", pow((flow/capacity),power): ";
            cout << pow((flow/capacity),power);
            cout << ", (1 + B * pow((flow/capacity),power) ): ";
            cout << (1 + B * pow((flow/capacity),power) );
            cout << ", free_flow_time * (1 + B * pow((flow/capacity),power) ): ";
            cout << free_flow_time * (1 + B * pow((flow/capacity),power) );
            cout << endl;
            return free_flow_time * (1 + B * pow((flow/capacity),power) );
        }

        void define_link(link_type idp, node_type start, node_type end, double capacityp, double lengthp, double freeFlowTTp, double Bp, int powerp){
            id = idp;
            head = start;
            tail = end;
            capacity = capacityp;
            length = lengthp;
            free_flow_time = freeFlowTTp;
            B = Bp;
            power = powerp;
        }

        // to test
        void print_link(){
            cout << id << " " << head << " " << tail << endl;
        }
    private:
        double B;
        int power;
        
        // Not useful now
        double slope, length, speed_limit, toll;
        char* osm_link_type;
};

class RoadNetwork : virtual public Graph {
    /** The Road network class contains the road network
    *   It is an interface to implement in the class road_network.cpp
    */
    private:
        unordered_map<link_type, Link> links;
        unordered_map<link_type, double> link_flows;

    public:
        void read_from_csv(string filename);
        void read_from_tntp(string filename);

        set<node_type> getOrigins();
        vector<node_type> getDestination(node_type origin);
        double getDemand(node_type origin, node_type destination);
        /** INPUTS: Filename in CSV format delimited by tabs (\t) from the github: https://github.com/bstabler/TransportationNetworks
        *   OUTPUTS: Instantiate the graph
        */
        void addLink(Link l){
            links[l.id] = l;
        }

        void update_flow(link_type link_id, double flow){
            cout << "update flow of link " << link_id << " to " << flow << endl;
            link_flows[link_id] = flow;
            cout << "Flow is now " << link_flows[link_id] << endl;
        }

        void update_all_weights(){
            for(auto i: links){
                Link link = i.second;
                setWeight(link.head, link.tail, link.travel_time( (link_flows[i.first]) ) );
                // if(link.head == 1){
                //     cout << "Update weight based on flow " << endl;
                //     cout << "Flow of " << i.first << " is " << link_flows[i.first] << endl;
                //     cout << "The corresponding travel time is then:\n\t " << link.travel_time( (link_flows[i.first]) ) << endl;
                //     link.print_link();
                //     setWeight(link.head, link.tail, link.travel_time( (link_flows[i.first]) ) );   
                // }
            }
        }

        // to test
        void print_link(link_type id){
            links[id].print_link();
        }
};

// #endif


