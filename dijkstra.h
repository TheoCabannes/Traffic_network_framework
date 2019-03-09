#ifndef __PARSE_DIJKSTRA_H_INCLUDED__
#define __PARSE_DIJKSTRA_H_INCLUDED__

#include <memory>
#include <tuple>
#include <unordered_map>
#include "parse_csv.h"
#include "graph.h"

using namespace std;

typedef pair<double, unsigned int> dist_mat;

dist_mat* dijkstra(int origin, unique_ptr<graph>& G);

#endif

