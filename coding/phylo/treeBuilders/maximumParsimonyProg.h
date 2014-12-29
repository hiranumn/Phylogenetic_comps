/*
 * maximumParsimonyProg.h
 * Inherits from abstractTreeBuilders.h
 * Given an aligned sequences and options, it creates a tree topology that explains the evolutionary history of sequences, progressively, using maximumParsimony method.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#ifndef MAXPARSPROG
#define MAXPARSPROG

#ifdef _OPENMP
#include <omp.h>
#endif

#include "abstractTreeBuilder.h" 
#include "tree.h"
#include <iostream>
#include <algorithm>

class MaximumParsimonyProgressive: public AbstractTreeBuilder
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
    /* Generate a vector containing all possible new trees of size 
     * n+1 from the original structure.
     * args: a const reference to a rooted topology, a pointer to a Species, an integer representing the dummy node
     * return: a vector of Trees
     */
    vector<Tree> getPossibleTopology(const Tree& baseTree,
				     Species* s,
				     int dummyNum);
};

#endif
