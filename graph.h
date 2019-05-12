#include <tuple>
#include <vector>
#include <set>

using namespace std;

typedef unsigned int node_type;  //!< Node_type is a macro for the integer type that indexes the nodes
typedef double dist_type; //!< Dist_type is a macro for the real type used to define distances

class Graph{
    /** The graph class contains the graph
    *   It is an interface to implement in the class graph.cpp
    */
    public:
    	virtual void addEdge(node_type i, node_type j); //!< Graph initialization. Adding a edge between the vertices i and j.
    	virtual void setWeight(node_type i, node_type j, dist_type dist); //!< Weight initialization. Setting the weight between the vertices i and j.
    	virtual dist_type getWeight(node_type i, node_type j); //!< Getting the weight between the vertices i and j.
    	vector<node_type>* successors(node_type i); //!< Return all the successors of vertice i.


    	// to test
    	void addEdgeWeighted(node_type i, node_type j, dist_type dist){
    		addEdge(i, j);
    		setWeight(i, j, dist);
    	};

    	set<node_type> vertices();
    	/** This function returns all the vertices that have successors in the networks */
};