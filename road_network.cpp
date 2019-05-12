#include <unordered_map>
#include <tuple>
#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>

#include <unordered_map>
#include <tuple>
#include <stdlib.h>
#include <vector>
#include <set> 
#include "graph.h"
#include "road_network.h"

using namespace std;

set<node_type> origins;
unordered_map<node_type, vector<node_type>> demand_list;
unordered_map<node_type, unordered_map<node_type, double>> sparse_demand_matrix;

double RoadNetwork::getDemand(node_type origin, node_type destination){
	return sparse_demand_matrix[origin][destination];
}

set<node_type> RoadNetwork::getOrigins(){
	return origins;
}

vector<node_type> RoadNetwork::getDestination(node_type origin){
	return demand_list[origin];
}

void RoadNetwork::read_from_csv(string filename){
    ifstream data;
    string line, temp;

    // Read the net file
    string file = filename + "_net.csv";
    data.open(file);
    if(data.fail()){
        fprintf(stderr, "The file %s doesn't exist\n", file.c_str());
        throw 3;
    }
    
    getline(data, temp);

    while(!data.eof()) {
        getline(data, line);
        stringstream buffer(line);

        link_type link_id;
        node_type start, end;
        dist_type freeFlowTT;

        getline(buffer, temp, ',');
        link_id = strtol(temp.c_str(), NULL, 10);

        getline(buffer, temp, ',');
        start = strtol(temp.c_str(), NULL, 10);

        getline(buffer, temp, ',');
        end = strtol(temp.c_str(), NULL, 10);
        
        getline(buffer, temp, ',');
        freeFlowTT = strtod(temp.c_str(), NULL);

        addEdge(start, end);
        if(freeFlowTT == 0){
            cout << line << endl;
            freeFlowTT = 0.000001;
        }
        setWeight(start, end, freeFlowTT);
    }

    data.close();

    // Read the demand file
    file = filename + "_od.csv";
    data.open(file);
    if(data.fail()){
        fprintf(stderr, "The file %s doesn't exist\n", file.c_str());
        throw 3;
    }

    
    getline(data, temp);
	cout << temp << endl;
    while(!data.eof()) {
        getline(data, line);
        stringstream buffer(line);

        node_type start, end;
        double dod;

        getline(buffer, temp, ',');
        start = strtol(temp.c_str(), NULL, 10);

        getline(buffer, temp, ',');
        end = strtol(temp.c_str(), NULL, 10);
        
        getline(buffer, temp, ',');
        dod = strtod(temp.c_str(), NULL);

        if(dod==0){
        	cout << line << endl;
        }
        origins.insert(start);
        demand_list[start].push_back(end);
        sparse_demand_matrix[start][end] = dod;
    }
}

void RoadNetwork::read_from_tntp(string filename){
	/** INPUTS:
     *
     * Filename in CSV format delimited by tabs (\t).
     *
     * OUTPUTS:
     *
     * Struct containing the titles of each column
     * and the csv as vector<vector<double>>.
     *
     */

    ifstream data;
    string line, temp;

    data.open(filename + "_net.tntp");
    
    getline(data, temp);
    int n = strtol(temp.c_str(), NULL, 10);

    link_type id = 0;
    while(!data.eof()) {
        getline(data, line);
        stringstream buffer(line);

        node_type start, end;
        dist_type length;
        double capacity, freeFlowTT, B, speed;
        int power;

        getline(buffer, temp, '\t');
        start = strtol(temp.c_str(), NULL, 10);

        getline(buffer, temp, '\t');
        end = strtol(temp.c_str(), NULL, 10);

        getline(buffer, temp, '\t');
        capacity = strtod(temp.c_str(), NULL);

        getline(buffer, temp, '\t');
        length = strtod(temp.c_str(), NULL);

        addEdge(start, end);
        setWeight(start, end, length);

  //       getline(buffer, temp, '\t');
  //       freeFlowTT = strtod(temp.c_str(), NULL);

  //       getline(buffer, temp, '\t');
  //       B = strtod(temp.c_str(), NULL);

  //       getline(buffer, temp, '\t');
  //       power = strtod(temp.c_str(), NULL);

		// getline(buffer, temp, '\t');
  //       speed = strtol(temp.c_str(), NULL, 10);

  //       if(freeFlowTT == 0){
  //           cout << line << endl;
  //           freeFlowTT = length/speed;
  //           cout << freeFlowTT << endl;
  //       }
  //       Link l;
  //       l.define_link(id, start, end, capacity, length, freeFlowTT, B, power);
  //       addLink(l);

        id = id + 1;
    }
    data.close();
}
