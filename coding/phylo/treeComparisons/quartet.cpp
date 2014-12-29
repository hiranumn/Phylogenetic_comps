/* 
 * quartet.cpp
 * Implementations of the functions defined in quartet.h
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/22/2014
 */

#include "quartet.h"

vector<double> QuartetDistance::getDistances(const Tree& comparison,
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
	if(!get(Vertex_t(), comparison.g, *v1)->dummy && get(Vertex_t(), comparison.g, *v1)->name.substr(0,5).compare(to_string("DUMMY")) != 0)
	    comparisonIDMap[get(Vertex_t(), comparison.g, *v1)->name] = *v1;

    //now get them for the rest...
    for(int i = 0; i < trees.size(); ++i)
	for(tie(v1, v2) = vertices(trees[i].g); v1 != v2; ++v1)
	    if(!get(Vertex_t(), trees[i].g, *v1)->dummy && get(Vertex_t(), trees[i].g, *v1)->name.substr(0,5).compare(to_string("DUMMY")) != 0)
	    {
		vertexIDMap[i][get(Vertex_t(), trees[i].g, *v1)->name] = *v1;
		if(i == 0) realSpecies.push_back(get(Vertex_t(), trees[i].g, *v1)->name);
	    }
    vector<double> returnVec(trees.size());

    for(int i = 0; i < trees.size(); ++i)
	returnVec[i] = getPairwiseQuartetDist(comparison, comparisonIDMap,
					      trees[i], vertexIDMap[i], realSpecies);
    return returnVec;


}


int QuartetDistance::getPairwiseQuartetDist(const Tree &T1, map<string, vertex_descriptor> T1map,
					    const Tree &T2, map<string, vertex_descriptor> T2map,
					    vector<string> realSpecies)
{
    //each vector of strings has 3 elements corresponding to the three
    //species we are currently considering.
    vector<vector<string> > triplets;

    for(int i = 0; i < realSpecies.size(); ++i)
	for(int j = i+1; j < realSpecies.size(); ++j)
	    for(int k = j+1; k < realSpecies.size(); ++k)
	    {
		vector<string> curTriplet;
		curTriplet.push_back(realSpecies[i]);
		curTriplet.push_back(realSpecies[j]);
		curTriplet.push_back(realSpecies[k]);
		triplets.push_back(curTriplet);
	    }

    int count = 0;

    for(int i = 0; i < triplets.size(); ++i)
    {
	vector<string> curTriplet = triplets[i];
	//First, get the centers for this triplet
	vertex_descriptor tree1Center = getCenter(T1, T1map, curTriplet);
	vertex_descriptor tree2Center = getCenter(T2, T2map, curTriplet);

	//Get the subtree sets. These maps map from one species in curtriplet
	//to a set of species in that subtree with respect to the center
	map<string, set<string> > T1subtreeSets;
	map<string, set<string> > T2subtreeSets;

	graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
	vector <vertex_descriptor> unvisited_neighbors;
	
	for (tie(ai, ai_end) = adjacent_vertices(tree1Center, T1.g); ai != ai_end; ++ai)
	{
	    set<string> curSpecies = getSpeciesInSubtree(T1, T1map, tree1Center, *ai);
	    //this could be made a lot more efficient.
	    if(find(curSpecies.begin(), curSpecies.end(), curTriplet[0]) != curSpecies.end())
		T1subtreeSets[curTriplet[0]] = curSpecies;
	    else if(find(curSpecies.begin(), curSpecies.end(), curTriplet[1]) != curSpecies.end())
		T1subtreeSets[curTriplet[1]] = curSpecies;
	    else if(find(curSpecies.begin(), curSpecies.end(), curTriplet[2]) != curSpecies.end())
		T1subtreeSets[curTriplet[2]] = curSpecies;
	    else cout << "no matching species from the triplet were found in the set!!!" << endl;
	}
	//T2 now
	for (tie(ai, ai_end) = adjacent_vertices(tree2Center, T2.g); ai != ai_end; ++ai)
	{
	    set<string> curSpecies = getSpeciesInSubtree(T2, T2map, tree2Center, *ai);
	    	    
	    //this could be made a lot more efficient.
	    if(find(curSpecies.begin(), curSpecies.end(), curTriplet[0]) != curSpecies.end())
		T2subtreeSets[curTriplet[0]] = curSpecies;
	    else if(find(curSpecies.begin(), curSpecies.end(), curTriplet[1]) != curSpecies.end())
		T2subtreeSets[curTriplet[1]] = curSpecies;
	    else if(find(curSpecies.begin(), curSpecies.end(), curTriplet[2]) != curSpecies.end())
		T2subtreeSets[curTriplet[2]] = curSpecies;
	    else cout << "no matching species from the triplet were found in the set!!!" << endl;
	}

	//Now we just need to compute set intersections...
	typedef map<string, set<string> >::iterator MyIter;
	for(MyIter x = T1subtreeSets.begin();
	    x != T1subtreeSets.end();
	    ++x)
	{
	    set<string> intersect;
	    set<string> s1 = T1subtreeSets[x->first];
	    set<string> s2 = T2subtreeSets[x->first];

	    set<string>::iterator myIter;

	    set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),
			     inserter(intersect,intersect.begin()));
	    
	    count += intersect.size();
	    count -= 1;
	}
    }

    int n = realSpecies.size();
    double binomialCoef = boost::math::binomial_coefficient<double>(n,4);
    return binomialCoef - count/4.0;
}


set<string> QuartetDistance::getSpeciesInSubtree(const Tree& tree,
						 map<string, vertex_descriptor> speciesMap,
						 vertex_descriptor center,
						 vertex_descriptor centerNeighbor)
{

    Tree tempTree(tree);
    vector<vertex_descriptor> reachable;
    remove_edge(center, centerNeighbor, tempTree.g);
    breadth_first_search(tempTree.g, centerNeighbor, 
			 visitor(make_bfs_visitor(write_property(
						      identity_property_map(),
						      back_inserter(reachable),
						      on_discover_vertex()))));

    typedef map<string, vertex_descriptor>::iterator MyIter;

    set<string> returnSet;

    for(MyIter x = speciesMap.begin();
	x != speciesMap.end();
	++x)
    {
	if(find(reachable.begin(), reachable.end(), x->second) != reachable.end())
	{
	    returnSet.insert(x->first);
	}
    }
    return returnSet;
}



vertex_descriptor QuartetDistance::getCenter(const Tree& tree,
					     const map<string, vertex_descriptor>& vertexIDMap,
					     const vector<string>& species)
{
    
    //maps from pathID -> path
    map<int, vector<vertex_descriptor> > pathMap;
    typedef map<int, vector<vertex_descriptor> >::iterator PathMapIter;

    pathMap[0] = getPath(vertexIDMap.at(species[0]),
			 vertexIDMap.at(species[1]),
			 tree);

    pathMap[1] = getPath(vertexIDMap.at(species[0]),
			 vertexIDMap.at(species[2]),
			 tree);

    pathMap[2] = getPath(vertexIDMap.at(species[1]),
			 vertexIDMap.at(species[2]),
			 tree);

    //find the common element of the three

    //maps from vertex -> numPaths of that vertex
    map<vertex_descriptor, int> pathCount;
    
    for(PathMapIter x = pathMap.begin();
	x != pathMap.end();
	++x)
    {
	for(int i = 0; i < x->second.size(); ++i)
	{
	    if (pathCount[x->second[i]] == 2) return x->second[i];
	    pathCount[x->second[i]] ++;
	}
    }
    cout << "No vertex was in three paths, this shouldn't happen" << endl;
    return -1;
}

vector<vertex_descriptor> QuartetDistance::getPath(vertex_descriptor start,
						   vertex_descriptor goal,
						   const Tree &tree)
{
    vector<vertex_descriptor> predMap(num_vertices(tree.g));
    vector<vertex_descriptor> path;

    dijkstra_shortest_paths(tree.g, start, predecessor_map(&predMap[0]));

    vertex_descriptor current = goal;
    
    while(current!=start)
    {
	path.push_back(current);
	current=predMap[current];
    }
    path.push_back(start);

    reverse(path.begin(), path.end());

    return path;
} 

