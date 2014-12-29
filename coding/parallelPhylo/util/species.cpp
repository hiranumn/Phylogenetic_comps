//species.cpp
//This file contains the implementations of functions defined in
//the Species class (see species.h)
#include "species.h"
#include <iostream>


//Mutation random functionality
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using namespace std;

Species::Species()
{
  /* no input species constructor */
  this->name = "";
  this->rawSequence = "";
  this->alignedSequence = "";
  this->visited = false;
}

Species::Species(string nameIn)
{
  /* named input species constructor */
  this->name = nameIn;
  this->rawSequence = "";
  this->alignedSequence = "";
  this->visited = false;
}

Species::Species(string nameIn, string rawSequenceIn)
{
  /* named, and raw sequence input species constructor */
  this->name = nameIn;
  this->rawSequence = rawSequenceIn;
  //THIS IS JUST FOR NOW THOUGH. CHANGE THIS NEXT LINE AFTER WE MAKE
  //SEQUENCE ALIGNERS.
  this->alignedSequence = rawSequenceIn;
  this->visited = false;
}

Species::Species(string nameIn, string rawSequenceIn, string alignedSequenceIn)
{
  /* named, raw sequence, and alligned sequence input species constructor */
  this->name = nameIn;
  this->rawSequence = rawSequenceIn;
  this->alignedSequence = alignedSequenceIn;
  this->visited = false;
}


// NO INSERTION OR DELETION HANDLING
string Species::mutate(float branchLength, double uRate, MutationMap mMap)
{
  /* Function for mutating a species across a given branch of a set length
     Takes as input the branch, and uses a set global mutation rate and character-character conversion rate */
  string mSequence = ""; //mutated sequence

  rng = RNGType(time(0)); //boost's random value generator, used to calculate values of high precision to allow for miniscule chances
  uniform_int<> intgen( 1, 2147483647 ); //2147483647 = 2^31 - 1 (max number for the boost rnd gen)
  variate_generator< RNGType, uniform_int<> > randint(rng, intgen);  

  // Holds the probability of mutating to the named character
  double probA;
  double probT;
  double probG;
  double probC;


  //Site independant mutation, directs a generated random value to an action
  //Generate probabilities using the branch length, uRate, and target char conversion rate
  //Probabilities are distinct for matching conversions and differing conversions
  for (int i = 0; i < this->rawSequence.length(); ++i){ 
    char currentChar = this->rawSequence[i];
    if (currentChar == 'A'){
      probA = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['A']; //(1-u*dt)*1 + u*dt*pi_j
    }
    else {
      probA = (uRate*branchLength)*mMap.mutationRates[currentChar]['A']; //(1-u*dt)*0 + u*dt*pi_j;
    }
    if (currentChar == 'T'){
      probT = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['T']; //(1-u*dt)*1 + u*dt*pi_j
    }
    else {
      probT = (uRate*branchLength)*mMap.mutationRates[currentChar]['T']; //(1-u*dt)*0 + u*dt*pi_j
    }
    if (currentChar == 'G'){
      probG = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['G']; //(1-u*dt)*1 + u*dt*pi_j
    }
    else {
      probG = (uRate*branchLength)*mMap.mutationRates[currentChar]['G']; //(1-u*dt)*0 + u*dt*pi_j
    }
    if (currentChar == 'C'){
      probC = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['C']; //(1-u*dt)*1 + u*dt*pi_j
    }
    else {
      probC = (uRate*branchLength)*mMap.mutationRates[currentChar]['C']; //(1-u*dt)*0 + u*dt*pi_j
    }

    // Using the probabilities and a random generated value, pick a character to add to the mutated sequence
    float rnd = (float) randint()/2147483647;
    bool changed = false;
    //A conversion
    if (rnd < probA && !changed){
      mSequence = mSequence + 'A';
      changed = true;
    }
    else{
      rnd = rnd - probA;
    }
    //T conversion
    if (rnd < probT && !changed){
      mSequence = mSequence + 'T';
      changed = true;
    }
    else{
      rnd = rnd - probT;
    }
    //G conversion
    if (rnd < probG && !changed){
      mSequence = mSequence + 'G';
      changed = true;
    }
    else{
      rnd = rnd - probG;
    }
    //C conversion
    if (rnd < probC && !changed){
      mSequence = mSequence + 'C';
    }

  }

  return mSequence;
}
