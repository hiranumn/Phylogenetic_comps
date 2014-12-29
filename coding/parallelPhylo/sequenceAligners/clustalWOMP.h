//clustalW.h is the header file for the clustalW class, which contains
//an implementation of the clustalW multiple sequence alignment algorithm
#ifndef CLUSTALWOMP
#define CLUSTALWOMP

#include "abstractAligner.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include "../util/sequenceData.h"
#include "../util/options.h"
#include "../util/distanceMatrix.h"
#include "../treeBuilders/neighborJoining.h"
#include "../treeBuilders/tree.h"
#include "../util/parallelUtilsOMP.h"

#include <omp.h>

class ClustalWOMP : public AbstractAligner
{
public:
    void align(SequenceData &data, const Options &options);

private:
    //Given 2 strings s1 and s2, returns a pair of strings representing
    //their optimal alignment
    pair<string, string> pairwiseAlign(string s1,
				       string s2,
				       const Options &options);


    //Given two raw sequences and the options, returns the distance between the
    //two.
    double parallelDistanceCalc(string seq1,
				string seq2,
				Options options);


    //Given 2 vectors of pointers to species representing the two alignment
    //profiles, calculate the SP score, based on the red book
    double calcSP(vector<Species*> profile1, 
		  vector<Species*> profile2, 
		  int profile1Col, int profile2Col,
		  const Options &options);

    //Given 2 vectors of points to species representing the two alignment
    //profiles, modify the aligned sequences in the species to align
    //these two profiles.
    void alignProfiles(vector<Species*> profile1,
		       vector<Species*> profile2,
		       const Options &options);
		       
    //Given a constant reference to a Tree and a constant reference to an Options
    //object, sets the aligned sequences within the species pointers correctly.
    void alignWithTree(const Tree &constTopologyIn, const Options &options);

    //Given a reference to a tree, a root node, a reference to a profile
    //to be filled by the function (a list
    //of species this particular call has aligned), a map of node->children, and the
    //options, threadsafe aligns the profiles of rootNode's two children nodes.
    void alignWithTreeRecursive(Tree topology,
				vertex_descriptor rootNode,
				vector<Species*> &myProfile,
				map<vertex_descriptor, vector<vertex_descriptor> > childMap,
				Options options);
};


#endif
