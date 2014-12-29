/*
 * neighborJoining.h
 * Inherits from abstractTreeBuilders.h
 * Given alinged sequences and options, this class performs neighborjoining reconstruction of a phylogenetic tree.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/14/2014
 */

#ifndef NEIGHBOR
#define NEIGHBOR

#include "abstractTreeBuilder.h" 
#include "../util/distanceMatrix.h"
#include <string>
#include <iostream>

class NeighborJoining: public AbstractTreeBuilder
{
private:
    /*
     * Given a constant reference to a distance matrix object, returns
     * {tip name -> "average" dist to all other tips} dictionary
     * args: const reference to the Distance Matrix
     * return: map of {tip name -> "average" distance to all other tips}
     */
    map<string, double> computeAverageTipDists(const DistanceMatrix &distMatrix);
    /* 
     * Given a distanceMatrix and the "average" dist dictionary, returns
     * the pair that minimizes Dij-ui-uj
     * args: const reference to the Distance Matrix, const reference to the map of tip name -> "average" distance
     * return: pair of tip names that minimizes Dij-ui-uj
     */
    pair<string, string> findClosestTips(const DistanceMatrix &distMatrix,
					 const map<string, double> &aveDists);

    /*
     * Given a constant reference to the distance matrix, a pair of strings
     * representing the closest pair of "tips" in that distance matrix, and
     * the map storing the "average distance" from each original node to 
     * all other nodes, returns a pair of doubles of the form (dist from new
     * node to child1, dist from new node to child2) where the input pair is
     * (child1, child2)
     * args: const reference to the Distance Matrix, pair of strings representing the closest pair of "tips" in the matrix, and the const reference to the map of tip name -> "average" distance
     * return: pair of doubles of the form (dist from new node to child1
     * dist from new node to child2) where the input pair is
     * (child1, child2)
     */
    pair<double, double> getDistancesToNewNode(const DistanceMatrix &dm,
					       pair<string, string> minimumDistPair, 
					       const map<string, double> &aveDists);
    
public:
    /*
     * Given a constant reference to the distance matrix and options,
     * return the optimal tree using neighbor joining
     * args: const reference to the Distance Matrix, the Options object, opt.
     * return: the optimal tree created using neighbor joining
     */
    Tree solve(const SequenceData &sd, Options opt);
    /*
     * Alternate solve method: if you want to use your own, precomputed
     * distance matrix instead of the default one from SequenceData, you can
     * use this method.
     */
    Tree solve(const SequenceData &sd, DistanceMatrix distMatrix, Options opt);
}; 

#endif
