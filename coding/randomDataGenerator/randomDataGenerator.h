/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* randomDataGenerator.h is the header file for the randomDataGenerator class, which contains
* an implementation of the random data generation functions
*/

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "../phylo/treeBuilders/tree.h"
#include "../phylo/util/mutationMap.h"
#include "../phylo/util/options.h"
#include <iostream>
#include <fstream>

using namespace boost;
using namespace std;

typedef boost::mt19937 RNGType;
typedef graph_traits<Tree_graph>::vertex_descriptor vd;
typedef graph_traits<Tree_graph>::edge_descriptor ed;

class RandomDataGenerator{
 public:
  //Constructor and associated values
  RandomDataGenerator(Options opt);
  int numS;
  int maxBranchLength;
  int randSpeciesDNALength; 
  string nameBase;
  int scount;
  double uRate;
  double insRate;
  double delRate;
  MutationMap mMap;

  RNGType rng;  //Boost random number generator

  /*
    Creates the topology of a tree
  */
  Tree* makeRandomTree(); 
  /*
    Fills the tree off a starting "seed" species and input parameters like global mutation rate
  */
  void assignSpecies(Tree* t);

  /*
    Workhorse function that mutates species from one sequence to a new one given input like branch lengths etc
  */
  void recursivelyMutate(Tree* t, vertex_descriptor cur,
                        Species* source,
                        map<vertex_descriptor, bool> visited,
                        vector<Species*> *leave);
};
