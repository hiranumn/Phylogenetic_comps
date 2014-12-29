/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* clustalW.h is the header file for the clustalW class, which contains
* an implementation of the clustalW multiple sequence alignment algorithm
*/

#ifndef CLUSTALW
#define CLUSTALW

#include "abstractAligner.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include "../util/sequenceData.h"
#include "../util/options.h"
#include "../util/distanceMatrix.h"
#include "../treeBuilders/neighborJoining.h"
#include "../treeBuilders/tree.h"

class ClustalW : public AbstractAligner
{
public:
    /*
        Main entry point for the clustal W algorithm
    */
    void align(SequenceData &data, const Options &options);

private:
    /*
        Given 2 strings s1 and s2, returns a pair of strings representing their optimal alignment
    */
    pair<string, string> pairwiseAlign(string s1,
				       string s2,
				       const Options &options);

    /*
        Given 2 vectors of pointers to species representing the two alignment
        profiles, calculate the SP score
    */
    double calcSP(vector<Species*> profile1, 
		  vector<Species*> profile2, 
		  int profile1Col, int profile2Col,
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
    void alignWithTree(const Tree &constTopologyIn,
                const Options &options);

    /*
        Recursive function utilized to explore the tree and align from the base up.
    */
    void alignWithTreeRecursive(const Tree &topology,
				vertex_descriptor rootNode,
				map<vertex_descriptor, bool> &visited,
				map<vertex_descriptor, vector<Species*> > &nodeProfiles,
				const Options &options);    
};


#endif
