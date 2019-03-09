#ifndef __PARSE_CSV_H_INCLUDED__
#define __PARSE_CSV_H_INCLUDED__

#include <memory>
#include <tuple>
#include "parse_csv.h"
#include "graph.h"

using namespace std;

struct header_and_csv {
    vector<string> headers;
    vector<vector<double>> csv;
};

header_and_csv read_data(string file, unique_ptr<graph>& g, string suffix);

#endif
