//parallelUtils.h
//Parallelization utility methods; general methods used by multiple algorithms that assist
//with the paralellization of their functionality
#ifndef PARAUTIL
#define PARAUTIL

#include "../treeBuilders/tree.h"
#include <boost/thread.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>

namespace ParallelUtils
{
    //given a topology and an edge descriptor, assigns the differential
    //in subgraph size between the two resulting subgraphs when cutting the 
    //graph along that edge to result
    void cutDifferential(const Tree& topology, edge_descriptor e, int &result);
    
    //Given a topology, picks a good edge to root the topology on.
    edge_descriptor pickRootEdge(const Tree& topology);

    //given a topology and a vertex in that topology, assigns the
    //number of vertices reachable from that vertex in the topology to result
    void reachableSize(const Tree& topology, vertex_descriptor v, int& result);

}

#endif
