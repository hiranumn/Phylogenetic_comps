//maximumParsimonyParallel.cpp
//Given an aligned sequences and options, it creates a tree topology
//that explains the evolutionary history of sequences, using maximumParsimony method.
//parallel implementation using open mp.

#include "maximumParsimonyParallel.h"

Tree MaximumParsimonyParallel::solve(const SequenceData &sd, Options opt)
{
    Tree curTree = getStartingTopology(sd);
    int curLength = getLength(curTree, sd.alignedLength);
    bool done = false;

    while (!done)
    {
	vector<Tree> neighbors = getNeighbors(curTree);
	vector<int> neighborLengths(neighbors.size());
	
	#pragma omp parallel for
	for(int i = 0; i < neighbors.size(); ++i)
	{
	    neighborLengths[i] = getLength(neighbors[i], sd.alignedLength);
	}
	
	vector<int>::iterator minResult = min_element(neighborLengths.begin(), neighborLengths.end());
	int minIndex = distance(neighborLengths.begin(), minResult);

	if (*minResult < curLength)
	{
	    curTree = neighbors[minIndex];
	    curLength = *minResult;
	}
	else done = true;
    }
    return curTree;
}


vector<Tree> MaximumParsimonyParallel::getNeighbors(const Tree & curTree)
{
    vector<Tree> returnVector;

    vector<edge_descriptor> internalEdges;
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(curTree.g); ei != ei_end; ++ei)
    {
	vertex_descriptor vertexOne = source(*ei, curTree.g);
	vertex_descriptor vertexTwo = target(*ei, curTree.g);
	graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
	int neighborCountOne = 0;
	int neighborCountTwo = 0;
	//"Internal edges" connect two nodes with exactly three neighbors each.
	for (tie(ai, ai_end) = adjacent_vertices(vertexOne, curTree.g); ai != ai_end; ++ai)
	    neighborCountOne ++;
	for (tie(ai, ai_end) = adjacent_vertices(vertexTwo, curTree.g); ai != ai_end; ++ai)
	    neighborCountTwo ++;
	if (neighborCountOne == 3 && neighborCountTwo == 3)
	    internalEdges.push_back(*ei);
    }
    //Now that we have vector of "internal edges". For each "internal edge",
    //generate two neighboring trees.
    for (int i=0; i < internalEdges.size(); ++i)
    {
	returnVector.push_back(getSwapNeighbor(curTree, internalEdges[i], 0));
	returnVector.push_back(getSwapNeighbor(curTree, internalEdges[i], 1));
    }
    return returnVector;
}


Tree MaximumParsimonyParallel::getSwapNeighbor(const Tree &curTree,
					       edge_descriptor curEdge,
					       int neighborNum)
{
    Tree returnTree(curTree);
	
    //Construct the first tree by swapping subtrees one and three
    //based on the Coursera video.
    vertex_descriptor vertexOne = source(curEdge, returnTree.g);
    vertex_descriptor vertexTwo = target(curEdge, returnTree.g);
    
    vector<vertex_descriptor> v1Neighbors, v2Neighbors;
    
    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    for(tie(ai, ai_end) = adjacent_vertices(vertexOne, returnTree.g);
	ai != ai_end;
	++ai)
    {
	if(*ai != vertexTwo) v1Neighbors.push_back(*ai);
    }
    
    for(tie(ai, ai_end) = adjacent_vertices(vertexTwo, returnTree.g);
	ai != ai_end;
	++ai)
    {
	if(*ai != vertexOne) v2Neighbors.push_back(*ai);
    }
    
    
    remove_edge(vertexOne, v1Neighbors[neighborNum], returnTree.g);
    remove_edge(vertexTwo, v2Neighbors[neighborNum], returnTree.g);
    
    edge_descriptor newEdgeOne = add_edge(vertexOne, v2Neighbors[neighborNum], returnTree.g).first;
    edge_descriptor newEdgeTwo = add_edge(vertexTwo, v1Neighbors[neighborNum], returnTree.g).first;
    
    put(edge_weight_t(), returnTree.g, newEdgeOne, 1.0);
    put(edge_weight_t(), returnTree.g, newEdgeTwo, 1.0);
    return returnTree;
}


Tree MaximumParsimonyParallel::getStartingTopology(const SequenceData &sd)
{
    Tree returnTree;
    typedef map<string, Species*>::const_iterator MapIterator;
    vector<vertex_descriptor> leafNodes;
    for (MapIterator iter = sd.speciesMap.begin();
	 iter != sd.speciesMap.end();
	 ++iter)
    {
	vertex_descriptor curLeaf = add_vertex(returnTree.g);
	put(Vertex_t(), returnTree.g, curLeaf, iter->second);
	leafNodes.push_back(curLeaf);
    }
    
    vertex_descriptor curJoin = leafNodes[0];
    int numDummy = 0;
    //this loop doesn't apply to the first or last leaf nodes
    for(int i = 1; i < leafNodes.size() - 1; ++i)
    {
	ostringstream convert;
	convert << numDummy;
	string dummyName = string("DUMMY_") + convert.str();
	Species* dummy = new Species(dummyName);
	this->dummies.push_back(dummy);

	vertex_descriptor curInternal = add_vertex(returnTree.g);
	put(Vertex_t(), returnTree.g, curInternal, dummy); 

	edge_descriptor firstEdge = add_edge(curInternal, leafNodes[i], returnTree.g).first;
	edge_descriptor secondEdge = add_edge(curInternal, curJoin, returnTree.g).first;

	put(edge_weight_t(), returnTree.g, firstEdge, 1.0);
	put(edge_weight_t(), returnTree.g, secondEdge, 1.0);

	curJoin = curInternal;
	++ numDummy;
    }

    edge_descriptor lastEdge = add_edge(curJoin, leafNodes[leafNodes.size()-1], returnTree.g).first;
    put(edge_weight_t(), returnTree.g, lastEdge, 1.0);
    return returnTree;
}

int MaximumParsimonyParallel::getLength(const Tree &constTopologyIn, int alignedLength)
{
    //First, we need to root the topology by inserting a node ON an edge.
    //Get first edge and treat it as root

    Tree topologyIn(constTopologyIn);

    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(topologyIn.g);
    vertex_descriptor neighborOne = source(*ei, topologyIn.g);
    vertex_descriptor neighborTwo = target(*ei, topologyIn.g);
    
    remove_edge(neighborOne, neighborTwo, topologyIn.g);
    //Put the "dummy root" into the graph.
    Species* dummyRoot = new Species("DUMMY_ROOT");
    vertex_descriptor nodeId = add_vertex(topologyIn.g);
    put(Vertex_t(), topologyIn.g, nodeId, dummyRoot);
    
    //Add edges from dummy root to the neighbors.
    edge_descriptor firstEdge = add_edge(nodeId, neighborOne, topologyIn.g).first;
    edge_descriptor secondEdge = add_edge(nodeId, neighborTwo, topologyIn.g).first;

    put(edge_weight_t(), topologyIn.g, firstEdge, 1.0);
    put(edge_weight_t(), topologyIn.g, secondEdge, 1.0);

    //Now, we need to make the neighbor dictionary, which maps from
    //vertex_descriptor -> vector<vertex_descriptor> where the key is the parent
    //and the vector is a vector of children.
    map<vertex_descriptor, bool> visited;
    map<vertex_descriptor, vector<vertex_descriptor> > childMap;
    vector<vertex_descriptor> verticesToMap;
    verticesToMap.push_back(nodeId);
    while(verticesToMap.size() != 0)
    {
	vertex_descriptor curVertex = verticesToMap.back();
	verticesToMap.pop_back();
	visited[curVertex] = true;

	graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
	vector <vertex_descriptor> unvisitedNeighbors;
	
	// get the vector of unvisited neighbors, which will either
	// have length of 2 or 0.
	for(tie(ai, ai_end) = adjacent_vertices(curVertex, topologyIn.g); ai != ai_end; ++ai)
	{
	    if(! visited[*ai])
	    {
		unvisitedNeighbors.push_back(*ai);
		verticesToMap.push_back(*ai);
	    }
	}
	childMap[curVertex] = unvisitedNeighbors;
    }


    int totalLength = 0;
    //This can be parallelized.
    #pragma omp parallel for reduction(+:totalLength)
    for(int i = 0; i < alignedLength; ++i)
    {
	vector<char> rootNodeSet;
	int val = getTotalLengthRecursive(topologyIn, nodeId, childMap, rootNodeSet, i);
	totalLength += val;
	
    }
    delete dummyRoot;
    
    return totalLength;
    
}

int MaximumParsimonyParallel::getTotalLengthRecursive(Tree topology,
						      vertex_descriptor rootNode,
						      map<vertex_descriptor, vector<vertex_descriptor> > childMap,
						      vector<char> &myNodeCharSet,
						      int seqIndex)
{   
    //No children = base case...
    if(childMap[rootNode].size() == 0)
    {
	char curChar = topology.speciesmap[rootNode]->alignedSequence.at(seqIndex);
	myNodeCharSet.push_back(curChar);
	return 0;
    }
    
    vertex_descriptor childOne = childMap[rootNode][0];
    vertex_descriptor childTwo = childMap[rootNode][1];
    int childOneLength, childTwoLength;
    vector<char> childOneChars, childTwoChars;

    #pragma omp task shared(childOneLength, childOneChars)
    childOneLength = getTotalLengthRecursive(topology, childOne, childMap,
					     childOneChars, seqIndex);
    #pragma omp task shared(childTwoLength, childTwoChars)
    childTwoLength = getTotalLengthRecursive(topology, childTwo, childMap,
					     childTwoChars, seqIndex);
    #pragma omp taskwait

    sort(childOneChars.begin(), childOneChars.end());
    sort(childTwoChars.begin(), childTwoChars.end());
    vector <char> intersection;
    set_intersection(childOneChars.begin(), childOneChars.end(),
		     childTwoChars.begin(), childTwoChars.end(),
		     std::back_inserter(intersection));
    
    if (intersection.size() == 0)
    {
	for(int i = 0; i < childOneChars.size(); ++i)
	{
	    myNodeCharSet.push_back(childOneChars[i]);
	}
	for(int i = 0; i < childTwoChars.size(); ++i)
	{
	    myNodeCharSet.push_back(childTwoChars[i]);
	}
	return 1 + childOneLength + childTwoLength;
    }
    else
    {
	myNodeCharSet = intersection;
	return childOneLength + childTwoLength;
    }
        
}
