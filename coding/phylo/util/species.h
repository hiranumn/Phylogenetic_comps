/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* species.h
* This file contains the non-abstract class species, which holds the data
* relating to species: name, sequence (unmodified), and sequence (aligned).
*/

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
  /*
    two constructors used in different situations: first one takes in just
    the raw sequence; second one takes in aligned sequence as well
  */
  Species();
  Species(string nameIn);
  Species(string nameIn, string rawSequenceIn);
  Species(string nameIn, string rawSequenceIn, string alignedSequenceIn);

  /* 
    mutation function, converts rawSequence to a new sequence based on inputs
    weight refers to branching length of the tree for mutation purposes
    returns the new sequence
  */
  string mutate(float weight,
                double uRate,
                double insRate,
                double delRate,
                MutationMap mMap); 
  string name;
  string rawSequence;
  string alignedSequence;
  string aminoSequence;
  map<string,double> kmerFreqs;

  // temporary / fabricated species flag (always an internal node on the tree)
  bool dummy;

  // Holds the probability of mutating to the named character
  double probA;
  double probT;
  double probG;
  double probC;

  //Boost random number generator
  RNGType rng;
};
#endif
