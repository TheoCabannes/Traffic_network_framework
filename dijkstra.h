#ifndef __PARSE_DIJKSTRA_H_INCLUDED__
#define __PARSE_DIJKSTRA_H_INCLUDED__

#include <memory>
#include <tuple>
#include <unordered_map>
#include "parse_csv.h"
#include "graph.h"

typedef std::pair<double, unsigned int> dist_mat;

std::unique_ptr<dist_mat> dijkstra(int origin, std::unique_ptr<graph>& G);

#endif

