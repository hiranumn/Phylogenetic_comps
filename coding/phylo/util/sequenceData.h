/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* sequenceData.h
* widely used non-abstract class storing species and the sequences associated
*/

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
    /*
        maps from {species name -> species object} for all species
    */
    map<string, Species*> speciesMap; 
    
    //aligned flag and attribute
    bool aligned;
    int alignedLength;             
    
    SequenceData(string filename);
    SequenceData();

    /*
        prints out all the data in the species map
    */
    void printData();

    /*
        prints out the aligned sequences
    */
    void printAlignment();
    
    /*
        returns a vector of all the names of the species in this SequenceData object.
    */
    vector<string> getAllSpeciesNames() const;
    
    /*
        gets the name of the first species
    */
    string getFirstSpeciesName() const;

    /*
        given an output filename, saves the aligned sequences in the standard format
    */
    void saveAlignment(string output);

    /*
        Given a filename formatted in the standard sequence format
        with pre-aligned genomes, fills this sequence object with the
        information from that file.
    */
    void readAlignedFromFile(string filename);
    
    /*
        Given a filename formatted in the standard sequence format,
        fills this SequenceData object with the information from
        that file.
    */
    void readFromFile(string filename);

    /*
        Outputs the value which represents the sequence data's total distance
        higher values are generally better
    */
    double totalDistance();

    ~SequenceData();
};

#endif
