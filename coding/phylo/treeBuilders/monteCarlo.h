/* 
 * monteCarlo.h
 * Inherits from abstractTreeBuilders.h
 * Given alinged sequences and options, this class uses Monte Carlo Markov Chain reconstruction of a phylogenetic tree.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#ifndef MCMC
#define MCMC

#ifdef _OPENMP
#include <omp.h>
#endif

#include "abstractTreeBuilder.h"
#include "maximumLikelihoodClimb.h"
#include "tree.h"
#include <iostream>
#include <algorithm>

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "boost/graph/adjacency_list.hpp"
#include <boost/graph/breadth_first_search.hpp>

#include <boost/math/special_functions/binomial.hpp>

class MonteCarlo: public AbstractTreeBuilder
{
public:
    /*
     * Given a constant reference to the distance matrix and options,
     * return the final tree using a Monte-Carlo Markov Chain simulation
     * args: const reference to the Distance Matrix, the Options object, opt.
     * return: the final tree after the simulation.
     */
    Tree solve(const SequenceData &sd, Options opt);
    
private:
    double phi;
    vertex_descriptor root;
    double HastingsRatio;
    double uRate;
    int alphabetSize;
    int iterations;
    vector<char> nuc;
    map<char, double> freq;
        
    /* 
     * Given a SequenceData object, make a new tree that will act as a
     * starting point" for the hill climber.
     * args: const reference to the Distance Matrix, vertex_descriptor representing the root
     * return: the starting tree topology
     */
          
    Tree getStartingTopology(const SequenceData &sd, vertex_descriptor* root);
    
    /*
     * Generates a random number between given min and max doubles
     * args: a double min value, a double max falue
     * return: a random double between the max and min values
     */
    double randPerturb(double fMin, double fMax);

    /* 
     * Given a const reference to a rooted topology, return a vector
     * of all internal edges.
     * args: a const reference to a tree
     * return: a vector of edge_descriptors representing the tree's internal edges
     */
    vector<edge_descriptor> getInternalEdges(const Tree &curTree);
    
    /* 
     * Given a const reference to a rooted topology, return a vector
     * of all internal vertices.
     * args: a const reference to  a tree
     * return: a vector of vertex_descriptors representing the tree's internal nodes
     */
    vector<vertex_descriptor> getInternalNodes(const Tree &curTree);

    /*
     * Given a source and a target vertex, return the path between the two
     * in tree T.
     * args: a vertex_descriptor representing the source node, a vertex_descriptor representing the target node, and a const reference to a tree.
     * return: a vector of vertex_descriptors representing the path between the source and target nodes in the tree.
     */
    vector<vertex_descriptor> getPath(vertex_descriptor start,
				      vertex_descriptor goal,
				      const Tree &tree);

    /* 
     * Given a const reference to a rooted topology and the parameter
     * phi, run Global with a Molecular Clock method to adjust the
     * current tree topology.
     * args: a const reference to a rooted tree topology
     * return: an adjusted tree
     */
    Tree globalWithMol(const Tree &curTree);

    /* 
     * Given a const reference to a rooted topology phi, run Local 
     * with a Molecular Clock method to adjust the current tree topology.
     * args: a const reference to a tree
     * returns: an adjusted tree
     */
    Tree localWithMol(const Tree &curTree);

    /* 
     * Given two vertex descriptors start, goal, return the total path
     * length (weights of branches) between the two nodes.
     * args: a vertex_descriptor representing the start node, a vertex_descriptor representing the target node, a const reference to a tree
     * return: a double representing the length of the path.
     */
    double getPathLength(vertex_descriptor start,
			 vertex_descriptor goal,
			 const Tree &tree);

    /*
     * Given the current and proposed tree topology
     * and the siteNum, return the probability of acceptance for the 
     * proposede new topology.
     * args: a const reference to the current tree, a const reference to the proposed tree, an integer representing the number of sites used
     * return: a double representing the probability of accepting the tree
     */
    double metropolisHastings(const Tree &curTree, const Tree &newTree, int sitenum);

    /* 
     * Calculates the likelihood of a given tree, starting from a given root.
     * The index number indicates which position of aligned sequence we 
     * are looking at.
     * args: a const reference to a tree, an integer representing the current position of the aligned sequence
     * return: a double representing the likelihood of the tree
     */
    double calculateLikelihood(const Tree& t, int siteNum); 

    
    /*
     * Recursively calculates the likelihood of the given tree.
     * args: a const reference to a tree, a vertex_descriptor representing the root, a char representing the state of the node, a map of vertex_descriptor to bools representing its visitedness, and an integer representing the sequence index
     * return: a double representing the likelihood of the tree
     */
    double calculateLikelihoodRecursively(const Tree& t,
					  vertex_descriptor rootNode,
					  char state,
					  map<vertex_descriptor, bool> &visited,
					  int seqIndex);
    
    /*
     * Calculates the mutation probability from the fromState to the toState,
     * taking into account the branch length between the two nodes.
     * args: a char representing the source state, a char representing the two state, a float representing the branch length.
     * return: a double representing the mutation probability
     */
    double mutationProbability(char fromState, char toState, float branchLength); 
    
    /*
     * Returns a prior probability given a tree.
     * for our algorithm, the prior has a uniform distribution
     * args: a const reference to a tree
     * return: a double representing the prior probability of that tree.
     */
    double prior(const Tree &tree);
    
};

    struct vertexAndPathLength{
      vertex_descriptor vert;
      double path;
    };
/*
 * Compares the path of a and b and returns true if a is smaller than b.
 * args: two const references to struct vertexAndPathLength, a and b.
 * return: a bool indicating whether a's pathlength is smaller than b's.
 */
bool cmp(const vertexAndPathLength &a, const vertexAndPathLength &b);
#endif
