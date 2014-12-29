/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* muscle.h is the header file for the muscle class, which contains
* an implementation of the MUSCLE multiple sequence alignment algorithm
*/

#ifndef MUSCLE
#define MUSCLE

#include "abstractAligner.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include "../util/sequenceData.h"
#include "../util/aminoAlphabet.h"
#include "../util/aminoSub.h"
#include "../util/options.h"
#include "../util/distanceMatrix.h"
#include "../treeBuilders/neighborJoining.h"
#include "../treeBuilders/tree.h"


class Muscle : public AbstractAligner
{
public:
    void align(SequenceData &data, const Options &options);


private:
    /*
        Given some data, finds the SP distance of all sequences as a value to be optimized
    */
    double SP(SequenceData &data);

    AminoAlphabet translate;

    /*
        Given UNALIGNED sequences, compute the distance using the frequencies of observed "words" in the sequence converted to amino acids
    */
    map<string, map<string, double> > kmer(SequenceData &data,
					   const Options &options);


    /*
        Given ALIGNED sequences, compute the distance off a measurement of matching/mistaching base proportions
    */
    map<string, map<string, double> > kimura(SequenceData &data,
					     const Options &options);


    /*
        Given 2 vectors of points to species representing the two alignment
        profiles, modify the aligned sequences in the species to align
        these two profiles.
    */
    void alignProfiles(vector<Species*> profile1,
		       vector<Species*> profile2,
		       const Options &options);
		       
    /*
        Given a constant reference to a Tree and a constant reference to an Options
        object, sets the aligned sequences within the species pointers correctly.
    */
    void alignWithTree(const Tree &constTopologyIn, const Options &options);

    /*
        Recursive function utilized to explore the tree and align from the base up.
    */
    void alignWithTreeRecursive(const Tree &topology,
				vertex_descriptor rootNode,
				map<vertex_descriptor, bool> &visited,
				map<vertex_descriptor, vector<Species*> > &nodeProfiles,
				const Options &options);

    /*
        Given 2 vectors of pointers to species representing the two alignment
        profiles, calculate the SP score
    */
    double logExpecScore(vector<Species*> profile1, 
          vector<Species*> profile2, 
          int profile1Col, int profile2Col,
          const Options &options);
};


#endif
