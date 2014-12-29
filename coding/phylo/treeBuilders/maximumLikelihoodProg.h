/*
 * maximumLikelihoodProg.h
 * Given an aligned sequences and options, it creates a tree topology that explains the evolutionary history of sequences
 * by using maximumLiokelihood method and progressive tree search approach.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */
#ifndef MAXLIKEPROG
#define MAXLIKEPROG

#ifdef _OPENMP
#include <omp.h>
#endif

#include "abstractTreeBuilder.h"
#include "../util/mutationMap.h" 
#include "../util/species.h"
#include <iostream>
#include <math.h>


using namespace std;
using namespace boost;


//Inherits from abstractTreeBuilders.h
class MaximumLikelihoodProgressive: public AbstractTreeBuilder
{
public:
    double uRate;
    double convergence;
    int alphabetSize;
    vector<char> nuc;
    map<char, double> freq;

    /*
     * args: a const reference to alinged sequences, and an option object which stores some parameters. 
     * return: an optimal tree created by a maximum likelihood method with progressive tree topology search. 
     */
    Tree solve(const SequenceData& sd, Options opt);
    
    /* Generate a vector containing all possible new trees of size 
     * n+1 from the original structure.
     * args: a const reference to a rooted topology, a pointer to a Species, an integer representing the dummy node
     * return: a vector of Trees
     */
    vector<Tree> getPossibleTopology(const Tree& baseTree, Species* s, int dummyNum);

    /* 
     * Calculates the likelihood value of a given tree. It is a wrapper function for 
     * the calculateLikelhoodRecursively function. This function also preprocesses the tree 
     * for the recursive function by adding a root node to a random location in the tree.  
     * args: const reference to a tree, and the length of DNA sequences
     * return: likelihood value of the tree
     */
    double calculateLikelihood(const Tree& t, int siteNum); //Pass in a tree object as a copy. 
    
    /* 
     * Calculates the likelihood value of a given tree by recursion.
     * args: a const reference t to a tree, vertex_descriptor of a current node, 
     * a character for a currently assumed state, a map that tells if a node is already visited,
     * and integer for the length of the DNA sequences
     * return: likelihood value of the tree
     */
    double calculateLikelihoodRecursively(const Tree& t,
					  vertex_descriptor rootNode,
					  char state,
					  map<vertex_descriptor, bool> &visited,
					  int seqIndex);
    
    /*
     * args: states of parent and child nodes, and branch length between them. 
     * return: probability of the transition happening
     */
    double mutationProbability(char fromState, char toState, float branchLength);    
    
    /* adjusts topology of a tree by using getOptimalBranchLength function on all branches in the tree.
     * args: a reference to a tree and the length of DNA sequence. 
     * return: does not return anything, but modifies the input tree. 
     */   
    void adjustTopology(Tree& t, int siteNum);

    /* adjusts the length of a given branch by using the A and B values of the branch.
     * A value is the sum of the products of the likelihood values of two subtrees 
     * over all possible states, multiplied by the frequency of each state.
     * A = sum( freq[state] * Likelihood[subtree1, state] * Likelihood[subtree2, state]) over states {A, T, C, G, -}
     * B value is the product of the sum of the likelihood values of two subtrees over all possible states,
     * multiplied by the frequency of each state
     * B = [sum( freq[state] * Likelihood[subtree1, state]) over states {A, T, C, G, -}]
     *     * [sum( freq[state] * Likelihood[subtree2, state]) over states {A, T, C, G, -}]
     *
     * args: an edge to optimize, a const reference to the tree that the edge belongs to, the length of the DNA sequences
     * return: optimal length for the input branch
     */
    float getOptimalBranchLength(edge_descriptor e, const Tree &t, int siteNum);

    /* This function calculates the A value of a given branch. 
     * A value is the sum of the products of the likelihood values of two subtrees 
     * over all possible states, multiplied by the frequency of each state.
     * A = sum( freq[state] * Likelihood[subtree1, state] * Likelihood[subtree2, state]) over states {A, T, C, G, -}
     *
     * args: two nodes that belong an edge whose A value is being calculated, index for the DNA site to look at
     *      and a const reference to the tree. 
     * return: A value 
     */
    double A(vertex_descriptor v1, vertex_descriptor v2, int index, const Tree& t);

    /* This function calculates the B value of a given branch. 
     * B value is the product of the sum of the likelihood values of two subtrees over all possible states,
     * multiplied by the frequency of each state
     * B = [sum( freq[state] * Likelihood[subtree1, state]) over states {A, T, C, G, -}]
     *     * [sum( freq[state] * Likelihood[subtree2, state]) over states {A, T, C, G, -}]
     *
     * args: two nodes that belong an edge whose B value is being calculated, index for the DNA site to look at
     *      and a const reference to the tree. 
     * return: B value
     */
    double B(vertex_descriptor v1, vertex_descriptor v2, int index, const Tree& t);
};

#endif
