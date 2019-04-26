#include <vector>
#include <cstring>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <memory>
#include <tuple>
#include "parse_csv.h"
#include "graph.h"

using std::begin;
using std::end;

header_and_csv read_data(std::string file, std::unique_ptr<graph>& G, std::string suffix) {

    /* INPUTS:
     *
     * Filename in CSV format delimited by tabs (\t).
     *
     * OUTPUTS:
     *
     * Struct containing the titles of each column
     * and the csv as vector<vector<double>>.
     *
     */
    using std::ifstream;
    using std::vector;
    using std::unordered_set;
    using std::string;
    using std::unique_ptr;
    using std::unordered_map;
    using std::getline;
    using std::find;
    using std::stringstream;
    using std::exception;

    ifstream data(file);
    string line;
    vector<vector<double>> parsed_csv;
    graph adjacency_lists;
    unordered_set<int> delete_indices;
    string cell;
    do {
        stringstream lineStream(line);
        getline(data, line);
    } while (line.find('<') == 0 
            || line.find(' ') == 0 
            || line.empty()
            || data.fail());

    if (data.fail()) {
        throw "Filename does not exist";    
    }

    vector<string> titles;
    stringstream titleStream(line);

    while (getline(titleStream, cell, '\t')) {
        titles.push_back(cell);
    }

    if (suffix == "net") { 
        int init_node_index = -1, term_node_index = -1;
        auto init_iterator = find_if(titles.begin(), titles.end(), 
                [](const string& str) { 
                    return str.find("Init node") != string::npos; 
                });
        auto term_iterator = find_if(titles.begin(), titles.end(), 
                [](const string& str) { 
                    return str.find("Term node") != string::npos; 
                });

        // HARD CODED. CHANGE THIS
        if(init_iterator != titles.end()) {
            init_node_index = distance(titles.begin(), init_iterator) - 1;
        } else {
            perror("Init node not found");
        }

        if(term_iterator != titles.end()) {
            term_node_index = distance(titles.begin(), term_iterator) - 1;
        } else {
            perror("Term node not found");
        }

        // Perhaps move above previous comment?
        for (int i = 0; i < titles.size(); i++) {
            if ((titles[i].find("~") != string::npos)
                    || (titles[i].find(";") != string::npos) 
                    || (titles[i].find(":") != string::npos)) {
                delete_indices.insert(i);
            }
        }

        for (int i = titles.size() - 1; i >= 0; --i) {
            if (find(delete_indices.begin(), delete_indices.end(), i) 
                    != delete_indices.end()) {
                titles.erase(titles.begin() + i);
            }
        }

        while (getline(data, line)) {
            stringstream lineStream(line);
            string cell; 
            vector<double> parsedRow;
            int k = 0;
            while (getline(lineStream, cell, '\t')) {
                if (find(delete_indices.begin(), delete_indices.end(), k) 
                        == delete_indices.end()) {
                    parsedRow.push_back(stod(cell));
                }
                k++;
            }
            unordered_map<unsigned int, double> end_node_and_tt {{parsedRow[term_node_index], 0.0}};
            adjacency_lists.graph[parsedRow[init_node_index]].insert(end_node_and_tt.begin(), end_node_and_tt.end());
            parsed_csv.push_back(parsedRow);
        }

        G = unique_ptr<graph>(new graph(adjacency_lists));
        
        header_and_csv retval;
        retval.headers = titles;
        retval.csv = parsed_csv;
        return retval;
    } else if (suffix == "flow") {
        int init_node_index = -1, term_node_index = -1;
        auto init_iterator = find_if(titles.begin(), titles.end(), 
                [](const string& str) { 
                    return str.find("From") != string::npos; 
                });
        auto term_iterator = find_if(titles.begin(), titles.end(), 
                [](const string& str) { 
                    return str.find("To") != string::npos; 
                }); 

        // HARD CODED. CHANGE THESE
        if(init_iterator != titles.end()) {
            init_node_index = distance(titles.begin(), init_iterator);
        } else {
            perror("From not found");
        }

        if(term_iterator != titles.end()) {
            term_node_index = distance(titles.begin(), term_iterator);
        } else {
            perror("To node not found");
        }

        for (int i = 0; i < titles.size(); i++) {
            if ((titles[i].find("~") != string::npos)
                    || (titles[i].find(";") != string::npos) 
                    || (titles[i].find(":") != string::npos)) {
                delete_indices.insert(i);
            }
        }

        for (int i = titles.size() - 1; i >= 0; --i) {
            if (find(delete_indices.begin(), delete_indices.end(), i) 
                    != delete_indices.end()) {
                titles.erase(titles.begin() + i);
            }
        }

        while (getline(data, line)) {
            stringstream lineStream(line);
            string cell; 
            vector<double> parsedRow;
            int k = 0;
            while (getline(lineStream, cell, '\t')) {
                if (find(delete_indices.begin(), delete_indices.end(), k) 
                        == delete_indices.end()) {
                    parsedRow.push_back(stod(cell));
                }
                k++;
            }
            unordered_map<unsigned int, double> end_node_and_tt {{parsedRow[term_node_index], 0.0}};
            adjacency_lists.graph[parsedRow[init_node_index]].insert(end_node_and_tt.begin(), end_node_and_tt.end());
            parsed_csv.push_back(parsedRow);
        }

        G = unique_ptr<graph>(new graph(adjacency_lists));
        
        header_and_csv retval;
        retval.headers = titles;
        retval.csv = parsed_csv;
        return retval;
 
    }
}
