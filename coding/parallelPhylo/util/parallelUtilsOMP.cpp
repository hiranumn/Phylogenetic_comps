//parallelUtils.cpp
//Contains implementations of the methods defined in parallelUtils.h

#include "parallelUtilsOMP.h"

int ParallelUtilsOMP::cutDifferential(const Tree& topology, edge_descriptor e)
{
    Tree myTopology(topology);
    vertex_descriptor v1, v2;
    v1 = source(e, myTopology.g);
    v2 = target(e, myTopology.g);
    remove_edge(v1, v2, myTopology.g);
    int size1, size2;
    
    #pragma omp parallel shared(size1, size2)
    {
	#pragma omp single
	{
	    #pragma omp task shared(size1)
	    size1 = ParallelUtilsOMP::reachableSize(myTopology, v1);
	    #pragma omp task shared(size2)
	    size2 = ParallelUtilsOMP::reachableSize(myTopology, v2);
	    #pragma omp taskwait
	}
    }
    return abs(size1 - size2);
}


edge_descriptor ParallelUtilsOMP::pickRootEdge(const Tree& topology)
{
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(topology.g);
    //start gradient decent at ei
    edge_descriptor curEdge = *ei;
    
    //keep track of values for edges so we don't repeat work
    map<edge_descriptor, int> edgeValues;

    //gradient decent
    while(true)
    {
	vector<edge_descriptor> edgesToSearch;
	edgesToSearch.push_back(curEdge);
	vertex_descriptor v1, v2;
	v1 = source(curEdge, topology.g);
	v2 = target(curEdge, topology.g);
	graph_traits<Tree_graph>::out_edge_iterator oei, oei_end;
	for(tie(oei, oei_end) = out_edges(v1, topology.g); oei != oei_end; ++oei)
	    if(*oei != curEdge) edgesToSearch.push_back(*oei);
	for(tie(oei, oei_end) = out_edges(v2, topology.g); oei != oei_end; ++oei)
	    if(*oei != curEdge) edgesToSearch.push_back(*oei);

	//determine how many searches we are going to be doing
	int numSearches = 0;
	for(int i = 0; i < edgesToSearch.size(); ++i)
	    if(edgeValues[edgesToSearch[i]]==0) ++ numSearches;

	vector<int> myResults(numSearches);
	int curNumSearches = 0;
	
	#pragma omp parallel for
	for(int i = 0; i < edgesToSearch.size(); ++i)
	{
	    if(edgeValues[edgesToSearch[i]]==0)
	    {
		myResults[i] = ParallelUtilsOMP::cutDifferential(topology, edgesToSearch[i]);
	    }
	}
	
	//put the results into the map...
	int curNumResults = 0;
	for(int i = 0; i < edgesToSearch.size(); ++i)
	{
	    if(edgeValues[edgesToSearch[i]]==0)
	    {
		edgeValues[edgesToSearch[i]] = myResults[curNumResults] + 1; //add one so we can still use == 0 as an indicator as to if we've searched
		curNumResults++;
	    }
	}

	//Now see if we can do better than the current edge
	int bestVal = edgeValues[curEdge];
	edge_descriptor bestEdge = curEdge;
	for(int i = 0; i < edgesToSearch.size(); ++i)
	{
	    if(edgeValues[edgesToSearch[i]] < bestVal)
	    {
		bestEdge = edgesToSearch[i];
		bestVal = edgeValues[edgesToSearch[i]];
	    }
	}
	if(bestEdge == curEdge) break;
	curEdge = bestEdge;
    }

    return curEdge;
}

int ParallelUtilsOMP::reachableSize(const Tree& topology, vertex_descriptor v)
{
    vector<vertex_descriptor> reachable;
    boost::breadth_first_search(topology.g, v, visitor(make_bfs_visitor(write_property(
									    identity_property_map(),
									    back_inserter(reachable),
									    on_discover_vertex()))));
    return reachable.size();
}
