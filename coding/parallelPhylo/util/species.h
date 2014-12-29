//species.h
//This file contains the non-abstract class species, which holds the data
//relating to species: name, sequence (unmodified), and sequence (aligned).

#ifndef SPECIES
#define SPECIES

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "mutationMap.h"
#include <string>

using namespace std;
using namespace boost;

typedef mt19937 RNGType;

class Species
{
 public:
  //two constructors used in different situations: first one takes in just
  //the raw sequence; second one takes in aligned sequence as well
  Species();
  Species(string nameIn);
  Species(string nameIn, string rawSequenceIn);
  Species(string nameIn, string rawSequenceIn, string alignedSequenceIn);
  string mutate(float weight, double uRate, MutationMap mMap); //weight refers to branching length of the tree for mutation purposes
  string name;
  string rawSequence; //I expect only rawSequence can be mutated on.
  string alignedSequence;
  bool visited;

  //Boost random number generator
  RNGType rng;
};
#endif
