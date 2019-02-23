#ifndef __PARSE_CSV_H_INCLUDED__
#define __PARSE_CSV_H_INCLUDED__

#include <vector>
#include <cstring>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>
#include <memory>
#include "graph.h"

using namespace std;

struct header_and_csv {
    vector<string> headers;
    vector<vector<double>> csv;
};

inline header_and_csv read_data(string file, unique_ptr<graph>& g);

#endif
