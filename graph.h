#ifndef __GRAPH_H_INCLUDED__
#define __GRAPH_H_INCLUDED__

#include <vector>
#include <unordered_map>

struct graph {
    std::unordered_map<unsigned int, std::vector<unsigned int>> graph;
};

#endif
