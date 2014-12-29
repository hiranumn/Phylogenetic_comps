//clustalW.h is the header file for the clustalW class, which contains
//an implementation of the clustalW multiple sequence alignment algorithm
#ifndef CLUSTALWPARA
#define CLUSTALWPARA

#include "abstractAligner.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/thread.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include "../util/sequenceData.h"
#include "../util/options.h"
#include "../util/distanceMatrix.h"
#include "../treeBuilders/neighborJoining.h"
#include "../treeBuilders/tree.h"
#include "../util/parallelUtils.h"


class ClustalWParallel : public AbstractAligner
{
public:
    void align(SequenceData &data, const Options &options);

private:
    //Given 2 strings s1 and s2, returns a pair of strings representing
    //their optimal alignment
    pair<string, string> pairwiseAlign(string s1,
				       string s2,
				       const Options &options);


    //Our parallel call -- each of the (numSeq choose 2) threads calls this method
    //exactly once. Given the index of the first sequence, the index of the second
    //sequence, a place to put the result, and constant references to the normal data
    //we need, puts the correct distance in result.
    void parallelDistanceCalc(int index1,
			      int index2,
			      double &result,
			      const SequenceData &data,
			      const vector<string> &names,
			      const Options &options);


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
    void alignWithTreeRecursive(const Tree& topology,
				const vertex_descriptor& rootNode,
				vector<Species*> &myProfile,
				const map<vertex_descriptor, vector<vertex_descriptor> >& childMap,
				const Options &options);
};


#endif
