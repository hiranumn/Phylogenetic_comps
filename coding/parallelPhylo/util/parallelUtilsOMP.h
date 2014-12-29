//parallelUtilsOMP.h
//Parallelization utility methods; general methods used by multiple algorithms that assist
//with the paralellization of their functionality. Open MP implementations.
#ifndef PARAUTILOMP
#define PARAUTILOMP

#include "../treeBuilders/tree.h"
#include <omp.h>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>

namespace ParallelUtilsOMP
{
    //given a topology and an edge descriptor, returns the differential
    //in subgraph size between the two resulting subgraphs when cutting the 
    //graph along that edge
    int cutDifferential(const Tree& topology, edge_descriptor e);
    
    //Given a topology, picks a good edge to root the topology on.
    edge_descriptor pickRootEdge(const Tree& topology);

    //given a topology and a vertex in that topology,
    //returns the number of verticies reachable from that vertex.
    int reachableSize(const Tree& topology, vertex_descriptor v);

}

#endif
