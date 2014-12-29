//sequenceData.h
//This file contains the non-abstract class SequenceData, which holds the data
//relating to the input sequences.

#ifndef SEQDATA
#define SEQDATA

#include <map>
#include <string>
#include <vector>
#include "species.h"

using namespace std;

class SequenceData
{
public:
    map<string, Species*> speciesMap;  // maps from {species name -> species object} for all species
    
    bool aligned;
    int alignedLength;  // length of aligned sequences, set this when the alignment done.
    
    SequenceData(string filename);
    
    // prints out the data in the species map
    void printData();
    // prints out the aligned sequences
    void printAlignment();
    
    // returns a vector of all the names of the species in this SequenceData object.
    vector<string> getAllSpeciesNames() const;
    
    // gets the name of the first species
    string getFirstSpeciesName() const;
    
    //Given a filename formatted in the standard sequence format,
    //fills this SequenceData object with the information from
    //that file.
    void readFromFile(string filename);
    ~SequenceData();
};

#endif
