//distanceMatrix.h
//This file contains the non-abstract function definitions
//for the distance matrix class

#ifndef DISTANCEMATRIX
#define DISTANCEMATRIX

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "sequenceData.h"

using namespace std;

class DistanceMatrix
{
public:
    //Constructor which takes in a precomputed matrix data structure
    DistanceMatrix(map<string, map<string, double> > matrixIn);

    //Given a sequenceData object and a boolean representing
    //whether or not to use the raw data (for the sequence
    //alignment algorihtms, for instance) fills in matrix
    //with appropriate values.
    DistanceMatrix(const SequenceData &data);

    //Distance function for sequences of same length (aligned)
    static double alignedDistance(string seq1, string seq2);
    
    //Distance function for sequences with varying length (pre-aligned)
    //this function is sometimes used in MSA algorithms
    static double rawDistance(string seq1, string seq2);

    //the actual matrix
    //this maps from {speciesA -> {speciesB -> distance(A,B) = distance(B,A)}}
    map<string, map<string, double> > matrix;

    //returns all unique pairs of species in the distanceMatrix
    vector<pair<string, string> > getUniquePairs() const;

    //given a pair of strings, merge the nodes that the represent under
    //the new name. Distance measures is based on the NJ coursera video.
    void mergePair(pair<string, string> mergePair, string newName);

    //prints this distance matrix
    void printMatrix() const;
};

#endif
