//maximmLikelihood.cpp
//Given an aligned sequences and options, it creates a tree topology
//that explains the evolutionary history of sequences.

#include "maximumLikelihoodClimb.h"
#include "../util/species.h"            /* MutationMap() */
#include <iostream>
#include <math.h>                       /* exp() */

Tree MaximumLikelihoodClimb::solve(const SequenceData& sd, Options opt)
{
    //set parameters
    this->uRate = opt.uRate;
    this->convergence = opt.convergence;
    this->nuc = opt.nuc;
    this->alphabetSize = nuc.size();
    this->freq = opt.freq;

    int siteNum = sd.speciesMap.at(sd.getFirstSpeciesName())->alignedSequence.size();
    Tree curTree = getStartingTopology(sd);
    adjustTopology(curTree, siteNum);
    double curLikelihood = calculateLikelihood(curTree, sd.alignedLength);
    bool done = false;

    while (!done)
    {
	cout << "loglike = " << curLikelihood << endl;
	vector<Tree> neighbors = getNeighbors(curTree);
	vector<int> neighborLikelihoods(neighbors.size());
	#pragma omp parallel for shared(neighborLikelihoods)
	for (int i = 0; i < neighbors.size(); i++)
	{
	    adjustTopology(neighbors[i], sd.alignedLength);
	    neighborLikelihoods[i] = calculateLikelihood(neighbors[i], sd.alignedLength);
	}
	vector<int>::iterator maxResult = max_element(neighborLikelihoods.begin(),
						      neighborLikelihoods.end());
	int maxIndex = distance(neighborLikelihoods.begin(), maxResult);

	cout << "best neighbor was this: " << *maxResult << endl;

	if (*maxResult > curLikelihood)
	{
	    curTree = neighbors[maxIndex];
	    curLikelihood = *maxResult;
	}
	else done = true;

    }
    return curTree;
}

vector<Tree> MaximumLikelihoodClimb::getNeighbors(const Tree & curTree)
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


Tree MaximumLikelihoodClimb::getSwapNeighbor(const Tree &curTree,
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
    
    put(edge_weight_t(), returnTree.g, newEdgeOne, .5);
    put(edge_weight_t(), returnTree.g, newEdgeTwo, .5);
    return returnTree;
}


Tree MaximumLikelihoodClimb::getStartingTopology(const SequenceData &sd)
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
    put(edge_weight_t(), returnTree.g, lastEdge, .5);
    return returnTree;
}

void MaximumLikelihoodClimb::adjustTopology(Tree& t, int siteNum)
{
    /* Tweaks the topology by adjusting branch lengths one at a time to find improvements */
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    //Try looping over 3 times...

    
    for(tie(ei, ei_end) = edges(t.g); ei != ei_end; ++ei)
    {
	float weight = getOptimalBranchLength(*ei, t, siteNum);
	put(edge_weight_t(), t.g, *ei, weight); 
    }
}


float MaximumLikelihoodClimb::getOptimalBranchLength(edge_descriptor e, const Tree &t, int siteNum)
{
    /* Finds the optimal length for a single branch */
    double weight = get(edge_weight_t(), t.g, e);
    
    double curP = 1 - exp(-weight);
    double oldP = 9999999; // Dramatic increase in likelihood of entering the improvement loop
    
    // Get the source and the target of an edge
    vertex_descriptor v1 = source(e, t.g);
    vertex_descriptor v2 = target(e, t.g);
    
    // Single time calculation of the A & B values which relate to subtree likelihood
    vector<double> vecA;
    vector<double> vecB;
    for(int i=0; i!=siteNum; i++)
    {
	vecA.push_back(A(v1, v2, i, t));
	vecB.push_back(B(v1, v2, i, t));
    }
    
    // While curP does not converge according to our set criteria
    while(abs(curP-oldP) > this->convergence)
    {
	oldP = curP;
	curP = 0;
	for(int i=0; i!=siteNum; i++)
	    curP += (vecB[i]*oldP)/(vecA[i]*(1-oldP)+vecB[i]*oldP);
	curP = (1.0/siteNum)*curP;
    }
    float newBranchLength = -log(1-curP);
    return newBranchLength;
}


double MaximumLikelihoodClimb::A(vertex_descriptor v1, vertex_descriptor v2, int index,
				 const Tree& t)
{
    /* Calculated the likelihood of a subtree given a matching character in a state, "likelihood boost" */
    Tree copy(t);
    remove_edge(v1, v2, copy.g);

    double sum = 0;
    double subprod;
    Species* s1 = get(Vertex_t(), t.g, v1); // Vertex 1 info
    Species* s2 = get(Vertex_t(), t.g, v2); // Vertex 2 info
    
    for(int i=0; i<this->alphabetSize; i++)
    {
	double curTerm = 1;
	curTerm *= this->freq[this->nuc[i]];
	map<vertex_descriptor, bool> visited;
	graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
	for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	    visited[*vi] = false;
	}
	curTerm *= calculateLikelihoodRecursively(copy, v1, this->nuc[i],
						      visited, index);
	
	for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	    visited[*vi] = false;
	}
	curTerm *= calculateLikelihoodRecursively(copy, v2, this->nuc[i],
						      visited, index);
	sum += curTerm;
    }
    return sum;
}


double MaximumLikelihoodClimb::B(vertex_descriptor v1, vertex_descriptor v2, int index,
				 const Tree& t)
{
    //copy the tree...
    Tree copy(t);
    remove_edge(v1, v2, copy.g);

    /* Calculates the likelihood of a subtree given mismatched characters in a state */
    double sum1 = 0;
    double sum2 = 0;

    for(int i = 0; i < this->nuc.size(); ++i)
    {
	map<vertex_descriptor, bool> visited;
	graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
	for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	    visited[*vi] = false;
	}
	sum1 += freq[this->nuc[i]]*calculateLikelihoodRecursively(copy, v1, this->nuc[i],
								  visited, index);
	//reset visited
	for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	    visited[*vi] = false;
	}
	sum2 += freq[this->nuc[i]]*calculateLikelihoodRecursively(copy, v2, this->nuc[i],
								  visited, index);
    }
    double ret = sum1*sum2;
    return ret;
}


double MaximumLikelihoodClimb::mutationProbability(char fromState, char toState, float branchLength)
{
    /* Calculates the probability of transitioning from a state to another in a single site, currently only accounts for freq */  
    double prob;
    if (fromState == toState)
    {
	prob = exp(-(uRate*branchLength)) + (1 - exp(-(uRate*branchLength)))*this->freq[toState];
    }
    else
    {
	prob = (1 - exp(-(uRate*branchLength)))*this->freq[toState];
    }
    return prob;
}

double MaximumLikelihoodClimb::calculateLikelihood(const Tree& inputTree, int siteNum)
{
    //First, we need to root the topology by inserting a node ON an edge.
    //Get first edge and treat it as root

    Tree topologyIn(inputTree);

    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(topologyIn.g);
    vertex_descriptor neighborOne = source(*ei, topologyIn.g);
    vertex_descriptor neighborTwo = target(*ei, topologyIn.g);
    
    float oldEdgeWeight = get(edge_weight_t(),
			      topologyIn.g,
			      edge(neighborOne, neighborTwo, topologyIn.g).first);

    float ew1 = oldEdgeWeight/2;
    float ew2 = oldEdgeWeight - ew1;

    remove_edge(neighborOne, neighborTwo, topologyIn.g);
    //Put the "dummy root" into the graph.
    Species* dummyRoot = new Species("DUMMY_ROOT");
    vertex_descriptor nodeId = add_vertex(topologyIn.g);
    put(Vertex_t(), topologyIn.g, nodeId, dummyRoot);
    
    //Add edges from dummy root to the neighbors.
    edge_descriptor firstEdge = add_edge(nodeId, neighborOne, topologyIn.g).first;
    edge_descriptor secondEdge = add_edge(nodeId, neighborTwo, topologyIn.g).first;

    put(edge_weight_t(), topologyIn.g, firstEdge, ew1);
    put(edge_weight_t(), topologyIn.g, secondEdge, ew2);
    
    double logLikelihood = 0;
    //This can be parallelized.
    for(int i = 0; i < siteNum; ++i)
    {
	double perSiteL = 0;
	for(int j = 0; j < this->alphabetSize; j++)
	{
	    map<vertex_descriptor, bool> visited;
	    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
	    for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
	    {
		visited[*vi] = false;
	    }
	    perSiteL += this->freq[this->nuc[j]] * calculateLikelihoodRecursively(topologyIn,
										  nodeId,
										  this->nuc[j],
										  visited,
										  i);
	}
    logLikelihood += log(perSiteL);
    }

    delete dummyRoot;
    
    return logLikelihood;
}

double MaximumLikelihoodClimb::calculateLikelihoodRecursively(const Tree& t,
							      vertex_descriptor rootNode,
							      char state,
							      map<vertex_descriptor, bool> &visited,
							      int seqIndex)
{
    visited[rootNode] = true;
    
    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    vector <vertex_descriptor> unvisitedNeighbors;

    // get the vector of unvisited neighbors, which will either
    // have length of 2 or 0.
    for (tie(ai, ai_end) = adjacent_vertices(rootNode, t.g); ai != ai_end; ++ai)
    {
	if (! visited[*ai])
	    unvisitedNeighbors.push_back(*ai);
    }
    
    //Leaf case
    if(unvisitedNeighbors.size() == 0)
    {
	// Get sequence
	string seq = get(Vertex_t(), t.g, rootNode)->alignedSequence;
	char nucleotide = seq[seqIndex];
	if(nucleotide == state)
	{
	    //map clean up
	    visited[rootNode] = false;
	    return 1;
	}
	else
	{
	    //map clean up
	    visited[rootNode] = false;
	    return 0;
	}
    }
    else
    {
	vertex_descriptor childOne = unvisitedNeighbors[0];
	vertex_descriptor childTwo = unvisitedNeighbors[1];
	double childP1 = 0;
	double childP2 = 0;
	
	// Check all possible nucleotides that can appear in child nodes and how that will affect the likelihood; 
	for(int j = 0; j<this->alphabetSize; j++)
	{
	    childP1 += 
		mutationProbability(state, this->nuc[j], 
				    get(edge_weight_t(), t.g, edge(rootNode,childOne,t.g).first))
		* calculateLikelihoodRecursively(t, childOne, this->nuc[j], visited, seqIndex);
	    childP2 += 
		mutationProbability(state, this->nuc[j],
				    get(edge_weight_t(), t.g, edge(rootNode,childOne,t.g).first))
		* calculateLikelihoodRecursively(t, childTwo, this->nuc[j], visited, seqIndex);
	}
	//reset the map
	visited[rootNode] = false;
	return childP1*childP2;
    }
}
