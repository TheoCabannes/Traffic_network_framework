#ifndef __GRAPH_H_INCLUDED__
#define __GRAPH_H_INCLUDED__

#include <vector>
#include <unordered_map>
#include <tuple>

struct graph {
    std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> graph;
};

#endif
