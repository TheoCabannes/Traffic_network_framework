#ifndef __PARSE_CSV_H_INCLUDED__
#define __PARSE_CSV_H_INCLUDED__

#include <memory>
#include <tuple>
#include "parse_csv.h"
#include "graph.h"

struct header_and_csv {
    std::vector<std::string> headers;
    std::vector<std::vector<double>> csv;
};

header_and_csv read_data(std::string file, std::unique_ptr<graph>& g, std::string suffix);

#endif
