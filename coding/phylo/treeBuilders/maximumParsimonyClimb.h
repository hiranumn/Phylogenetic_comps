/*
 * maximumParsimonyClimb.h
 * Inherits from abstractTreeBuilders.h
 * Given an aligned sequences and options, it creates a tree topology that explains the evolutionary history of sequences, progressively, using maximumParsimony method.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#ifndef MAXPARSCLIMB
#define MAXPARSCLIMB

#ifdef _OPENMP
#include <omp.h>
#endif

#include "abstractTreeBuilder.h" 
#include "tree.h"
#include <iostream>
#include <algorithm>

class MaximumParsimonyClimb: public AbstractTreeBuilder
{
public:
     /*
     * Given a constant reference to the distance matrix and options,
     * return the optimal tree using maximum parsimony.
     * args: const reference to the Distance Matrix, the Options object, opt.
     * return: the optimal tree created using maximum parsimony.
     */
    Tree solve(const SequenceData &sd, Options opt);
private:
    /*
     * args: a tree's topology
     * return: the length of that tree
     */
    int getLength(const Tree &topologyIn, int alignedLength);

    /*
     * args: a const reference to a rooted topology, a vertex descriptor of the root node, a map of vertex descriptors to booleans initially all false for all nodes in the beginning, a map of vertex descriptors to sets of characters those vertices could be in accordance with Fitch's algorithm, and the index in the sequence alignment we are considering
     * return: the total length of the tree.
     */
    int getTotalLengthRecursive(const Tree &topology,
				vertex_descriptor rootNode,
				map<vertex_descriptor, bool> &visited,
				map<vertex_descriptor, vector<char> > &nodeCharSets,
				int seqIndex);
    
    /* 
     * Given a SequenceData object, make a new tree that will act as a
     * "starting point" for the hill climber.
     * args: const reference to a SequenceData sd
     * return: a starting topology of the tree.
     */
    Tree getStartingTopology(const SequenceData &sd);

    /*
     * Given a constant reference to a tree, returns a vector of "neighbors"
     * of that tree, in accordance with the coursera video's nearest
     * neighbor interchange function.
     * args: const reference to a tree
     * return: a vector of trees
     */
    vector<Tree> getNeighbors(const Tree & curTree);

    /* 
     * Given a tree, an internal edge, and either 0 or 1 
     * (the "neighbor number" which represents the swapping pattern 
     * to generate a specific neighbor), returns a neighbor of this 
     * tree swapped on the curEdge with respect to the given neighbor number.
     * args: const reference to a tree
     * return: a neighbor tree
     */
    Tree getSwapNeighbor(const Tree &curTree,
			 edge_descriptor curEdge,
			 int neighborNum);
};

#endif
