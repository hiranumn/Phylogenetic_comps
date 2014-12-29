/* 
 * quartet.h
 * Contains the QuartetDistance class. Quartet distance is a topological
 * dissimilarity metric. A quartet is a phylogenetic tree with only
 * 4 species, and this class contains functions that compute the percentage
 * of quartet reductions that trees share. This class inherits from
 * the one defined in abstractTreeComparison.h
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/18/2014
 */

#ifndef QUARTET
#define QUARTET

#include "../treeBuilders/tree.h"
#include "abstractTreeComparison.h"

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "boost/graph/adjacency_list.hpp"
#include <boost/graph/breadth_first_search.hpp>

#include <boost/math/special_functions/binomial.hpp>

#include <vector>
#include <algorithm>
#include <set>

using namespace std;
using namespace boost;

class QuartetDistance: public AbstractTreeComparison
{
public:
    /*
     * Given a "root" tree and a vector of trees to compare against, return
     * a vector of doubles corresponding to the ordered pairwise quartet distances between
     * the "root" tree and the given vector of trees.
     * args: a "root" tree to compare against and a vector of trees to compare against the
     * root tree
     * returns: a vector of distances from each of the comparison trees to the root tree.
     * The length of the return vector is equal to the length of the input tree vector.
     */ 
    vector<double> getDistances(const Tree& comparison,
				const vector<Tree>& trees);

private:

    /*
     * Given a tree, a map from {species names -> vertex IDs} for that tree,
     * and a vector of three species,
     * returns the vertex descriptor in the tree that is the "center"
     * of those species (ie the path between any two contains
     * that central node)
     * args: a tree to find the center of, a map from species names to vertex
     * ids in that tree, and a vector of three species in that tree
     * returns: the unique vertex_descriptor contained in the pahts from species[0]
     * to species[1], from species[1] to species[2], and from species[0] to species[2].
     */
    vertex_descriptor getCenter(const Tree& tree,
				const map<string, vertex_descriptor>& vertexIDMap,
				const vector<string>& species);

    
    /*
     * Given a source and a target vertex, return the path between the two
     * in tree T.
     * args: a starting and ending vertex_descriptor in the tree
     * returns: the path from start to goal, including start and goal.
     */
    vector<vertex_descriptor> getPath(vertex_descriptor start,
				      vertex_descriptor goal,
				      const Tree &tree);

    
    /*
     * Given Tree1 and its species -> vertex map, Tree2 and its species -> vertex map,
     * and a vector of the names of the actual species in the tree, return the
     * pairwise quartet distance between the two trees.
     * args: T1, its corresponding {species name -> vertex_descriptor} map,
     * T2, its its corresponding {species name -> vertex_descriptor} map, and
     * a vector of strings representing the species to consider, usually the
     * set of real species in the tree, rather than the internal dummy nodes.
     * returns: the quartet distance between T1 and T2.
     */
    int getPairwiseQuartetDist(const Tree &T1, map<string, vertex_descriptor> T1map,
			       const Tree &T2, map<string, vertex_descriptor> T2map,
			       vector<string> realSpecies);



    /*
     * Given a tree, a species -> vertex map of that tree, a center point, and
     * a neighbor of that centerpoint (to define a subtree, you need a node and
     * which neighbor to consider), returns a vector of strings corresponding
     * to the speices in that particular subtree.
     * args: a tree, its correspoding {species name -> vertex_descriptor} map,
     * a vertex and one of its neighbors that define a subtree (ie the subree
     * is all nodes reachable from centerNeighbor)
     * return: a set of species in a subtree
     */
     set<string> getSpeciesInSubtree(const Tree& tree,
				    map<string, vertex_descriptor> speciesMap,
				    vertex_descriptor center,
				    vertex_descriptor centerNeighbor);

};

#endif
