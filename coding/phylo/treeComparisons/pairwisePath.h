/* 
 * pairwisePath.h
 * Contains the PairwisePathLength class. Pairwise pathlength is a distance
 * metric that first computes all possible unique pairs of species. Next,
 * it computes the distance between those two species in a graph, for each
 * pair. Finally, it calculates the euclidian distance between the vector
 * of pairwise differences between two trees. That euclidean norm is a measure
 * of tree dissimilarity.
 * This class inherits from the one defined in abstractTreeComparison.h
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/19/2014
 */

#ifndef PAIRWISE
#define PAIRWISE

#include "../treeBuilders/tree.h"
#include "abstractTreeComparison.h"

#include "boost/graph/floyd_warshall_shortest.hpp"
#include "boost/graph/adjacency_list.hpp"
#include <boost/graph/exterior_property.hpp>

#include <vector>
#include <math.h>

using namespace std;
using namespace boost;

class PairwisePathLength: public AbstractTreeComparison
{
public:
    /*
     * Given a "root" tree and a vector of trees to compare against, return
     * a vector of doubles corresponding to the ordered pairwise pairwise path distances between
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
     * Given a tree, a vector of string pairs to calculate the distances between,
     * and a {species name -> vertex_descriptor> map for the given tree, returns
     * a vector of doubles representing the distances between all species pairs.
     * If stringPairs = < (species1, species2), (species1, species3), (species2, species3) >
     * for instance, this function returns 
     * < d(species1, species2), d(species1, species3), d(species2, species3) >
     * where d(x,y) is the distance between two species in the tree, accounting for
     * branch lengths
     * args: const reference to a tree, a vector of pairs of species names you would
     * like to calculate the distance between, and a map of {species names -> vertex_descriptors>
     * for the given tree
     * returns: a vector of length equal to stringPairs.size() representing the pairwise
     * distances between the input species pairs
     */
    vector<double> getPairwisePathlengths(const Tree& myTree,
					  vector<pair<string, string> > stringPairs,
					  map<string, vertex_descriptor> vertexIDMap);
};

#endif
