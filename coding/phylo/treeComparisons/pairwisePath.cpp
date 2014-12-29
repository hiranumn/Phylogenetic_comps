/* 
 * pairwisePath.cpp
 * Implementations of the functions defined in pairwisePath.h
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/22/2014
 */

#include "pairwisePath.h"


vector<double> PairwisePathLength::getDistances(const Tree& comparison,
						const vector<Tree>& trees)
{
    //Maps from tree index -> species name -> vertex, IE to get the vertex_descriptor
    //for "X" in the third tree, you would do vertexIDMap[3]["X"]
    map<int, map<string, vertex_descriptor> > vertexIDMap;
    //ID map for the comparison tree
    map<string, vertex_descriptor> comparisonIDMap; 
    //get a vector of all non-dummy species names
    vector<string> realSpecies;
    
    //get a vector of all verticies in graph that are actually species
    graph_traits<Tree_graph>::vertex_iterator v1, v2;
    
    //first get them for the comparison tree (TREE 0)
    for(tie(v1, v2) = vertices(comparison.g); v1 != v2; ++v1)
	if(!get(Vertex_t(), comparison.g, *v1)->dummy)
	    comparisonIDMap[get(Vertex_t(), comparison.g, *v1)->name] = *v1;
    
    //now get them for the rest...
    for(int i = 0; i < trees.size(); ++i)
	for(tie(v1, v2) = vertices(trees[i].g); v1 != v2; ++v1)
	    if(!get(Vertex_t(), trees[i].g, *v1)->dummy)
	    {
		vertexIDMap[i][get(Vertex_t(), trees[i].g, *v1)->name] = *v1;
		if(i == 0) realSpecies.push_back(get(Vertex_t(), trees[i].g, *v1)->name);
	    }
    
    //now that we have all species vertex, generate the pairwise ordering
    vector<pair<string, string> > pairOrdering;

    for(int i = 0; i < realSpecies.size(); ++i)
    {
	for(int j = i+1; j < realSpecies.size(); ++j)
	{
	    pair<string, string> curPair(realSpecies[i], realSpecies[j]);
	    pairOrdering.push_back(curPair);
	}
    }
    

    vector<vector<double> > treePathDistanceVectors(trees.size(), vector<double>(pairOrdering.size(), 0));
    for(int i = 0; i < trees.size(); ++i)
    {
	treePathDistanceVectors[i] = getPairwisePathlengths(trees[i], pairOrdering, vertexIDMap[i]);
    }

    vector<double> comparisonVector = getPairwisePathlengths(comparison, pairOrdering, comparisonIDMap);

    vector<double> returnVector;

    for(int i = 0; i < trees.size(); ++i)
    {
	vector<double> curCompare = treePathDistanceVectors[i];
	double curSum = 0;
	for(int j = 0; j < comparisonVector.size(); ++j)
	    curSum += pow(comparisonVector[j] - curCompare[j], 2.0);
	double result = pow(curSum, .5);
	returnVector.push_back(result);
	break;
    }

    return returnVector;
}

vector<double> PairwisePathLength::getPairwisePathlengths(const Tree& myTree,
							  vector<pair<string, string> > stringPairs,
							  map<string, vertex_descriptor> vertexIDMap)
{
    vector<double> returnVector;
    map<vertex_descriptor, map<vertex_descriptor, float> > distances;
    floyd_warshall_all_pairs_shortest_paths(myTree.g, distances);

    for(int i = 0; i < stringPairs.size(); ++i)
    {
	pair<string, string> curPair = stringPairs[i];
	pair<vertex_descriptor, vertex_descriptor> curVertexPair(vertexIDMap[curPair.first], vertexIDMap[curPair.second]);
	returnVector.push_back(distances[curVertexPair.first][curVertexPair.second]);
    }

    double maxEle = *max_element(returnVector.begin(), returnVector.end());
    for(int i = 0; i < returnVector.size(); ++i)
	returnVector[i] /= maxEle;


    return returnVector;
}
