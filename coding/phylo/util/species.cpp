/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* sequenceData.cpp
* This file implements the non-abstract class species, which holds the data
* relating to species: name, sequence (unmodified), and sequence (aligned).
*/

#include "species.h"
#include <iostream>

//Mutation random functionality
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using namespace std;

/*
  no input species constructor 
*/
Species::Species()
{
  this->name = "";
  this->rawSequence = "";
  this->alignedSequence = "";
  this->aminoSequence = "";
  this->dummy = false;
}

/* 
  named input species constructor
*/
Species::Species(string nameIn)
{
  this->name = nameIn;
  this->rawSequence = "";
  this->alignedSequence = "";
  this->aminoSequence = "";
  this->dummy = false;
}

/*
  named, and raw sequence input species constructor 
*/
Species::Species(string nameIn, string rawSequenceIn)
{
  this->name = nameIn;
  this->rawSequence = rawSequenceIn;
  this->alignedSequence = rawSequenceIn;
  this->aminoSequence = "";
  this->dummy = false;
}

/*
  named, raw sequence, and alligned sequence input species constructor
*/
Species::Species(string nameIn, string rawSequenceIn, string alignedSequenceIn)
{
  this->name = nameIn;
  this->rawSequence = rawSequenceIn;
  this->alignedSequence = alignedSequenceIn;
  this->aminoSequence = "";
  this->dummy = false;
}


/* Function for mutating a species across a given branch of a set length
   Takes as input the branch, and uses a set global mutation rate and character-character conversion rate
   Assumes site independant mutation.

   Caveat :  is currently not flexible to alphabets outside 'ATGC-' 
              additionally the mutation should always retain a higher likelihood of no change on a character 
*/
string Species::mutate(float branchLength,
                      double uRate,
                      double insRate,
                      double delRate, 
                      MutationMap mMap)
{

  string mSequence = ""; //mutated sequence

  if (uRate * branchLength > 1){
    cout << "Please alter the mutation rate or the branch lengths, values are too high." << endl;
    return "ERRORSTRING";
  }

  //the following equations are awfully sensitive to input
  //the uRate and branchlength cannot multiply to exceed 1!
  for (int i = 0; i < this->rawSequence.length(); ++i)
  { 
    char currentChar = this->rawSequence[i];

    //Set values to the mismatch probabilities
    // u*dt*pi_j;
    probA = (uRate*branchLength)*mMap.mutationRates[currentChar]['A']; 
    probT = (uRate*branchLength)*mMap.mutationRates[currentChar]['T'];
    probG = (uRate*branchLength)*mMap.mutationRates[currentChar]['G']; 
    probC = (uRate*branchLength)*mMap.mutationRates[currentChar]['C'];
    //Alter the matching value to a higher mutation probability
    // (1-u*dt) + u*dt*pi_j
    switch (currentChar)
    {
      case 'A':
            probA = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['A'];
            break;
      case 'T':
            probT = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['T'];
            break;
      case 'G':
            probG = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['G'];
            break;
      case 'C':
            probC = (1-(uRate*branchLength))+(uRate*branchLength)*mMap.mutationRates[currentChar]['C'];
            break;
    }

   // Deletion event occurance
    float delRnd = (float) rand()/RAND_MAX;
    if (delRnd < delRate*branchLength){
      mSequence = mSequence.substr(0, mSequence.size()-1); //deletes the last item
    }
    // Insertion event occurence
    float insRnd = (float) rand()/RAND_MAX;
    if (insRnd < insRate*branchLength){ 
      int i = rand()%4;
      switch (i)
      {
        case 0:
          mSequence += 'A'; //duplicates the last item
          break;
        case 1:
          mSequence += 'T'; //duplicates the last item
          break;
        case 2:
          mSequence += 'G'; //duplicates the last item
          break;
        case 3:
          mSequence += 'C'; //duplicates the last item
          break;
      }
    }

    // Using the probabilities and a random generated value, pick a character to add to the mutated sequence
    float rnd = (float) rand()/RAND_MAX;
    rnd -= probA;
    if (rnd < 0){
      mSequence  += 'A';
      continue;
    }
    rnd -= probT;
    if (rnd < 0){
      mSequence  += 'T';
      continue;
    }
    rnd -= probG;
    if (rnd < 0){
      mSequence  += 'G';
      continue;
    }
    rnd -= probC;
    if (rnd < 0){
      mSequence  += 'C';
      continue;
    }


  }

  return mSequence;
}
