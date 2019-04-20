#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <tuple>
#include <stdlib.h>
#include "parse_csv.h"
#include "graph.h"
#include "dijkstra.h"

struct travel_time_function {
    vector<vector<double>> table_net;
    int index_fft;
    int index_B;
    int index_capacity; 
    int index_power;
};

//header_and_csv read_data(string, unique_ptr<graph>&);
travel_time_function init_travel_time_function(string, unique_ptr<graph>&);
vector<double> test_travel_time_function(string, travel_time_function&, unique_ptr<graph>&);
vector<double> travel_time(vector<double>&, travel_time_function&);

int main(int argc, char** argv) { 
    using std::vector;
    using std::string;
    unique_ptr<graph> G;

    string network = "data/SiouxFalls/SiouxFalls";
    /// DEBUG PURPOSES
//    vector<int> keys;
//    keys.reserve(g.get()->graph.size());
//    for (auto k : g.get()->graph)
//        keys.push_back(k.first);
//    
//    for (auto i : keys) {
//        printf("Neighbours of %d:\n", i);
//        for (auto j : g.get()->graph[i]) {
//            cout << get<0>(j);
//            cout << " ";
//        }
//        cout << "" << endl;
//    }
    /// END DEBUG
    auto init_tt_vals_sf = init_travel_time_function(network, G);
    test_travel_time_function(network, init_tt_vals_sf, G);
    for (int i = 1; i < 20; i++) 
    dijkstra(i, G);
    
    cout << "" << endl;

   /* 
    network = "data/Anaheim/Anaheim";
    read_data(network + "_net.tntp", g);
    auto init_tt_vals_an = init_travel_time_function(network, g);
    test_travel_time_function(network, init_tt_vals_an, g);
    */

    return 0;
}

travel_time_function init_travel_time_function(string network, unique_ptr<graph>& G) {

    /* INPUTS:
     *
     * Filename to a network.
     *
     * OUTPUTS:
     *
     * Struct containing CSV and indices of constants in the 
     * travel time formula.
     *
     */

    string file_type = "net";
    string suffix = "_" + file_type + ".tntp";
    auto legend_table_net = read_data(network + suffix, G, file_type);
    travel_time_function ttf;

    int index_B, index_power, index_fft, index_capacity, j = 0;

    for (auto legend_name : legend_table_net.headers) {
        if (legend_name.find("B") != string::npos) {
            index_B = j;
        }

        if (legend_name.find("Power") != string::npos) {
            index_power = j;
        }

        if (legend_name.find("Free Flow Time") != string::npos) {
            index_fft = j;
        }

        if (legend_name.find("Capacity") != string::npos) {
            index_capacity = j;
        }

        ++j;
    }

    ttf.table_net = legend_table_net.csv;
    ttf.index_B = index_B;
    ttf.index_power = index_power;
    ttf.index_fft = index_fft;
    ttf.index_capacity = index_capacity;

    return ttf;

}

vector<double> test_travel_time_function(string network, travel_time_function& ttf, unique_ptr<graph>& G) {

    /* INPUTS:
     *
     * Filename to a network and struct containing travel time formula constants. 
     *
     *
     * OUTPUTS:
     *
     * Vector describing the difference between the travel time of flow_solutions and
     * cost_solutions for each link.
     *
     */

    string file_type = "flow";
    string suffix = "_" + file_type + ".tntp";
    auto legend_table_flow = read_data(network + suffix, G, file_type);;

    int index_flow = -1, index_cost = -1, from = -1, to = -1, j = 0;

    for (auto legend_name : legend_table_flow.headers) {
        if (legend_name.find("From") != string::npos) {
            from = j;
        }

        if (legend_name.find("To") != string::npos) {
            to = j;
        }

        if (legend_name.find("Volume") != string::npos) {
            index_flow = j;
        }

        if (legend_name.find("Cost") != string::npos) {
            index_cost = j;
        }

        ++j;
    }

    // Order matters
    vector<double> flow;
    vector<tuple<unsigned int, unsigned int>> indices_per_flow;
    vector<double> cost_solution;
    vector<double> travel_time_difference;

    for (auto f : legend_table_flow.csv) {
        flow.push_back(f[index_flow]);
        indices_per_flow.push_back(make_tuple(f[from], f[to]));
    }

    for (vector<double> f : legend_table_flow.csv) {
        cost_solution.push_back(f[index_cost]);
    }
    
    vector<double> computed_travel_time = travel_time(flow, ttf);
    for (int i = 0; i < flow.size(); i++) {
        travel_time_difference.push_back(computed_travel_time[i] - cost_solution[i]);
        G.get()->graph.at((int) get<0>(indices_per_flow[i])).at((int) get<1>(indices_per_flow[i])) = computed_travel_time[i];
    }
    /*
    for (auto k : travel_time_difference) {
        printf("%20f\n", k);;
    }
    */

    return travel_time_difference;
}

vector<double> travel_time(vector<double>& flow, travel_time_function& ttf) {

    /* INPUTS:
     *
     * The initial flow allocation, the struct ttf containing the 
     * power m_a, capacity c_a, fft t_a^0, B.
     *
     * OUTPUTS:
     *
     * The travel time calculated using the formula
     * t(f_a) = t_a^0 * (1 + B(f_a/c_a)^m_a)
     * for each link.
     *
     */

    vector<double> aon_fft;
    vector<double> aon_power;
    vector<double> aon_capacity;
    vector<double> aon_B;

    for (auto i : ttf.table_net) {
        aon_fft.push_back(i[ttf.index_fft]);
        aon_power.push_back(i[ttf.index_power]);
        aon_capacity.push_back(i[ttf.index_capacity]);
        aon_B.push_back(i[ttf.index_B]);
    }

    vector<double> aon_retval;

    //#pragma omp parallel for
    for (int i = 0; i < flow.size(); i++) {
        aon_retval.push_back(aon_fft[i] * (1 + aon_B[i] * pow(flow[i] / aon_capacity[i], aon_power[i])));
    }

    cout << "Hello" << endl;
    cout << aon_retval.size() << endl;
    return aon_retval;

}

