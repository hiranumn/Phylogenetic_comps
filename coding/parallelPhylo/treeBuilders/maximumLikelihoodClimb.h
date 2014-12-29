//maximumLikelihood.h
//inherits from abstractTreeBuilder.h
//Given an aligned sequences and options, it creates a tree topology that explains the evolutionary history of sequences. 

#ifndef MAXLIKECLIMB
#define MAXLIKECLIMB

#include "abstractTreeBuilder.h"
#include "../util/mutationMap.h" 

using namespace std;
using namespace boost;



class MaximumLikelihoodClimb: public AbstractTreeBuilder
{
public:
    double uRate;
    double convergence;
    int alphabetSize;
    vector<char> nuc;
    map<char, double> freq;
    Tree solve(const SequenceData& sd, Options opt);
    
    //Given a constant reference to a tree, returns a vector of "neighbors" of that
    //tree, in accordance with the coursera video's nearest neighbor interchange
    //function.
    vector<Tree> getNeighbors(const Tree & curTree);

    //Given a SequenceData object, make a new tree that will act as a
    //"starting point" for the hill climber.
    Tree getStartingTopology(const SequenceData &sd);
    
    //Given a tree, an internal edge, and either 0 or 1 (the "neighbor number" which
    //represents the swapping pattern to generate a specific neighbor) returns
    //a neighbor of this tree swapped on the curEdge with respect to the given neighbor
    //number.
    Tree getSwapNeighbor(const Tree &curTree,
			 edge_descriptor curEdge,
			 int neighborNum);


    //Calculates the likelifood of a given tree, starting from a given root.
    //The index number indicates which position of alinged sequence we are looking at.
    //Lets implement this first.
    double calculateLikelihood(const Tree& t, int siteNum); //Pass in a tree object as a copy. 
    
    double calculateLikelihoodRecursively(const Tree& t,
					  vertex_descriptor rootNode,
					  char state,
					  map<vertex_descriptor, bool> &visited,
					  int seqIndex);
    
    double mutationProbability(char fromState, char toState, float branchLength); 
    
    
    //Finding a best tree topology
    //implemented this next.
    //vector<Tree> getPossibleTopology(const Tree& baseTree, Species* s);
    void adjustTopology(Tree& t, int siteNum);
    float getOptimalBranchLength(edge_descriptor e, const Tree &t, int siteNum);
    double A(vertex_descriptor v1, vertex_descriptor v2, int index, const Tree& t);
    double B(vertex_descriptor v1, vertex_descriptor v2, int index, const Tree& t);
    
    //For testing purpose
    Tree makeTestTree();
}; 

#endif
