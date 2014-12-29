/*
 * maximumParsimonyProg.cpp
 * Inherits from abstractTreeBuilders.h
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#include "maximumParsimonyProg.h"

Tree MaximumParsimonyProgressive::solve(const SequenceData &sd, Options opt)
{
    int siteNum = sd.speciesMap.at(sd.getFirstSpeciesName())->alignedSequence.size();

    map<string, Species*>::const_iterator iter;
    vector<Species*> species;
    // Place all species into a vector for ease of iteration later
    for(iter = sd.speciesMap.begin(); iter != sd.speciesMap.end(); iter++)
    {
	species.push_back(iter->second);
    }
    
    Tree curTree;
    // Init base tree
    vertex_descriptor v1, v2;
    v1 = add_vertex(curTree.g);
    v2 = add_vertex(curTree.g);
    put(Vertex_t(), curTree.g, v1, species[0]);
    put(Vertex_t(), curTree.g, v2, species[1]);
    edge_descriptor firstEdge = add_edge(v1, v2, curTree.g).first;
    put(edge_weight_t(), curTree.g, firstEdge, 1);


    //build up the tree
    for(int i = 2; i != species.size(); i++)
    {
	// Generate possible topologies of the size + 1, find the best choice
	vector<Tree> toplist = getPossibleTopology(curTree, species[i], i-2);
	int bestNeighborIndex = 0;
	int bestNeighborLength = getLength(toplist[0], siteNum);
	for(int j=1; j != toplist.size(); j++)
	{
	    int compLength = getLength(toplist[j], siteNum);
	    if(compLength < bestNeighborLength)
	    {
		bestNeighborIndex = j;
		bestNeighborLength = compLength;
	    } 
	}
	curTree = toplist[bestNeighborIndex];
    }
    return curTree;
}

vector<Tree> MaximumParsimonyProgressive::getPossibleTopology(const Tree& baseTree,
							      Species* s,
							      int dummyNum)
{
    vector<Tree> ret;
    
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    // For each branch in the given tree.
    for(tie(ei, ei_end) = edges(baseTree.g); ei != ei_end; ++ei)
    {
	// Duplicate the tree
	Tree t(baseTree);
	
	// Get the edge information from the original tree by vertex descriptor
	vertex_descriptor v1 = source(*ei, baseTree.g);
	vertex_descriptor v2 = target(*ei, baseTree.g);
	
	// Put a dummy node on a new edge
	edge_descriptor e = edge(v1,v2,t.g).first;
	
	remove_edge(v1, v2, t.g);
	//remove_edge(*ei, t.g);
	vertex_descriptor dummyNodeId = add_vertex(t.g);

	string dummyName = "DUMMY_"; dummyName += to_string(dummyNum);
	Species* dummy = new Species(dummyName);
	dummy -> dummy = true;
	put(Vertex_t(), t.g, dummyNodeId, dummy);
	dummies.push_back(dummy);

	edge_descriptor newEdge1 = add_edge(dummyNodeId, v1, t.g).first;
	edge_descriptor newEdge2 = add_edge(dummyNodeId, v2, t.g).first;
	
	put(edge_weight_t(), t.g, newEdge1, 1);
	put(edge_weight_t(), t.g, newEdge2, 1);

	//Now attach current species to dummy node with edge weight 5
	vertex_descriptor nodeId = add_vertex(t.g);
	put(Vertex_t(), t.g, nodeId, s);
	
	edge_descriptor newEdgeReal = add_edge(nodeId, dummyNodeId, t.g).first;
	put(edge_weight_t(), t.g, newEdgeReal, 1);

	// Push t into the vector
	ret.push_back(t);
    }
    return ret;
}

int MaximumParsimonyProgressive::getLength(const Tree &constTopologyIn, int alignedLength)
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
    dummyRoot -> dummy = true;
    vertex_descriptor nodeId = add_vertex(topologyIn.g);
    put(Vertex_t(), topologyIn.g, nodeId, dummyRoot);
    
    //Add edges from dummy root to the neighbors.
    edge_descriptor firstEdge = add_edge(nodeId, neighborOne, topologyIn.g).first;
    edge_descriptor secondEdge = add_edge(nodeId, neighborTwo, topologyIn.g).first;

    put(edge_weight_t(), topologyIn.g, firstEdge, 1.0);
    put(edge_weight_t(), topologyIn.g, secondEdge, 1.0);
    
    int totalLength = 0;
    //This can be parallelized.
    #pragma omp parallel for reduction(+:totalLength)
    for(int i = 0; i < alignedLength; ++i)
    {
	map<vertex_descriptor, bool> visited;
	graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
	for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
	{
	    visited[*vi] = false;
	}
	map<vertex_descriptor, vector<char> > nodeCharSets;
	totalLength += getTotalLengthRecursive(topologyIn, nodeId, visited, nodeCharSets, i);
    }

    delete dummyRoot;
    
    return totalLength;
    
}

int MaximumParsimonyProgressive::getTotalLengthRecursive(const Tree &topology,
					      vertex_descriptor rootNode,
					      map<vertex_descriptor, bool> &visited,
					      map<vertex_descriptor, vector<char> > &nodeCharSets,
					      int seqIndex)
{
    visited[rootNode] = true;

    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    vector <vertex_descriptor> unvisited_neighbors;

    // get the vector of unvisited neighbors, which will either
    // have length of 2 or 0.
    for (tie(ai, ai_end) = adjacent_vertices(rootNode, topology.g); ai != ai_end; ++ai)
    {
	if (! visited[*ai])
	    unvisited_neighbors.push_back(*ai);
    }
    
    // base case when looking at a leaf node:
    // add the sequence character to charSet, 
    // and return length 0.
    if (unvisited_neighbors.size() == 0)
    {
	vector <char> charSet;
        string seq = get(Vertex_t(), topology.g, rootNode)->alignedSequence;
	char curChar = seq[seqIndex];
	charSet.push_back(curChar);
	nodeCharSets[rootNode] = charSet;
	return 0;
    }
    vertex_descriptor childOne = unvisited_neighbors[0];
    vertex_descriptor childTwo = unvisited_neighbors[1];
    int childOneLength = getTotalLengthRecursive(topology, childOne, visited,
						 nodeCharSets, seqIndex);
    int childTwoLength = getTotalLengthRecursive(topology, childTwo, visited,
						 nodeCharSets, seqIndex);
    
    vector <char> childOneChars = nodeCharSets[childOne];
    vector <char> childTwoChars = nodeCharSets[childTwo];
    sort(childOneChars.begin(), childOneChars.end());
    sort(childTwoChars.begin(), childTwoChars.end());
    vector <char> intersection;
    set_intersection(childOneChars.begin(), childOneChars.end(),
		     childTwoChars.begin(), childTwoChars.end(),
		     std::back_inserter(intersection));
    
    if (intersection.size() == 0)
    {
	vector <char> charUnion;
	for(int i = 0; i < childOneChars.size(); ++i)
	{
	    charUnion.push_back(childOneChars[i]);
	}
	for(int i = 0; i < childTwoChars.size(); ++i)
	{
	    charUnion.push_back(childTwoChars[i]);
	}
	nodeCharSets[rootNode] = charUnion;
	return 1 + childOneLength + childTwoLength;
    }
    else
    {
	nodeCharSets[rootNode] = intersection;
	return childOneLength + childTwoLength;
    }
        
}
