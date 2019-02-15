#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <omp.h>
#include "parse_csv.h"

using namespace std;

struct travel_time_function {
    vector<vector<double>> table_net;
    int index_fft;
    int index_B;
    int index_capacity; 
    int index_power;
};

header_and_csv read_data(string);
travel_time_function init_travel_time_function(string);
vector<double> test_travel_time_function(string, travel_time_function&);
vector<double> travel_time(vector<double>&, travel_time_function&);

int main(int argc, char** argv) { 
    string network = "data/SiouxFalls/SiouxFalls";
    read_data(network + "_net.tntp");
    auto init_tt_vals_sf = init_travel_time_function(network);
    test_travel_time_function(network, init_tt_vals_sf);
    
    cout << "" << endl;
    
    network = "data/Anaheim/Anaheim";
    read_data(network + "_net.tntp");
    auto init_tt_vals_an = init_travel_time_function(network);
    test_travel_time_function(network, init_tt_vals_an);
    

    return 0;
}

travel_time_function init_travel_time_function(string network) {

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
    auto legend_table_net = read_data(network + suffix);
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

vector<double> test_travel_time_function(string network, travel_time_function& ttf) {

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
    auto legend_table_flow = read_data(network + suffix);

    int index_flow = -1, index_cost = -1, j = 0;

    for (auto legend_name : legend_table_flow.headers) {
        if (legend_name.find("Volume") != string::npos) {
            index_flow = j;
        }

        if (legend_name.find("Cost") != string::npos) {
            index_cost = j;
        }

        ++j;
    }

    vector<double> flow;
    vector<double> cost_solution;
    vector<double> travel_time_difference;

    for (auto f : legend_table_flow.csv) {
        flow.push_back(f[index_flow]);
    }

    for (vector<double> f : legend_table_flow.csv) {
        cost_solution.push_back(f[index_cost]);
    }
    
    vector<double> computed_travel_time = travel_time(flow, ttf);
    for (int i = 0; i < flow.size(); i++) {
        travel_time_difference.push_back(computed_travel_time[i] - cost_solution[i]);
    }

    for (auto k : travel_time_difference) {
        printf("%20f\n", k);;
    }

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
    return aon_retval;

}

