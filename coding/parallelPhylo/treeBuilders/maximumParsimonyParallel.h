//maximumParsimonyParallel.h
//Inherits from abstractTreeBuilders.h
//Given alinged sequences and options, this class performs neighborjoining reconstruction of a phylogenetic tree.

#ifndef MAXPARSPARA
#define MAXPARSPARA

#include "abstractTreeBuilder.h" 
#include "tree.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include <omp.h>

class MaximumParsimonyParallel: public AbstractTreeBuilder
{
public:
    Tree solve(const SequenceData &sd, Options opt);
private:
//Given a tree's topology, returns the length of that tree
    int getLength(const Tree &topologyIn, int alignedLength);
    
    //Given a rooted topology, a vertex descriptor
    //of the root node, a map that maps from parents -> children, a reference to
    //a vector<char> to be filled (represents THIS node's charset) and an index
    //in the alignment we're considering, returns the length of the tree
    //from the rootNode down.
    int getTotalLengthRecursive(Tree topology,
				vertex_descriptor rootNode,
				map<vertex_descriptor, vector<vertex_descriptor> > childMap,
				vector<char> &myNodeCharSet,
				int seqIndex);
    
    //Given a SequenceData object, make a new tree that will act as a
    //"starting point" for the hill climber.
    Tree getStartingTopology(const SequenceData &sd);

    //Given a constant reference to a tree, returns a vector of "neighbors" of that
    //tree, in accordance with the coursera video's nearest neighbor interchange
    //function.
    vector<Tree> getNeighbors(const Tree & curTree);

    //Given a tree, an internal edge, and either 0 or 1 (the "neighbor number" which
    //represents the swapping pattern to generate a specific neighbor) returns
    //a neighbor of this tree swapped on the curEdge with respect to the given neighbor
    //number.
    Tree getSwapNeighbor(const Tree &curTree,
			 edge_descriptor curEdge,
			 int neighborNum);
};

#endif
