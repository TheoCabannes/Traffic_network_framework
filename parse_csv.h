#ifndef __PARSE_CSV_H_INCLUDED__
#define __PARSE_CSV_H_INCLUDED__

#include <vector>
#include <cstring>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>

using namespace std;

struct header_and_csv {
    vector<string> headers;
    vector<vector<double>> csv;
};

inline header_and_csv read_data(string file) {

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
        parsed_csv.push_back(parsedRow);
    }
    
    header_and_csv retval;
    retval.headers = titles;
    retval.csv = parsed_csv;
    return retval;
}

#endif
