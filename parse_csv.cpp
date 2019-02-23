#include <vector>
#include <cstring>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <memory>
#include "parse_csv.h"
#include "graph.h"

using namespace std;

header_and_csv read_data(string file, unique_ptr<graph>& g) {

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
            || line.empty());

    vector<string> titles;
    stringstream titleStream(line);

    while (getline(titleStream, cell, '\t')) {
        titles.push_back(cell);
    }

    int init_node_index = -1, term_node_index = -1;
    auto init_iterator = find_if(titles.begin(), titles.end(), 
            [](const string& str) { 
                return str.find("Init node") != string::npos; 
            });
    auto term_iterator = find_if(titles.begin(), titles.end(), 
            [](const string& str) { 
                return str.find("Term node") != string::npos; 
            });

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
        adjacency_lists.graph[parsedRow[init_node_index]]
            .push_back(parsedRow[term_node_index]);
        parsed_csv.push_back(parsedRow);
    }

    g = unique_ptr<graph>(new graph(adjacency_lists));
    
    header_and_csv retval;
    retval.headers = titles;
    retval.csv = parsed_csv;
    return retval;
}