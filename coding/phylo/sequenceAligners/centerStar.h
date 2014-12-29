/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* centerStar.h is the header file for the muscle class, which contains
* an implementation of the centerStar multiple sequence alignment algorithm
*/

#ifndef CSTAR
#define CSTAR

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


class CenterStar : public AbstractAligner
{
public:
    void align(SequenceData &data, const Options &options);


private:
    /*
        Given 2 strings s1 and s2, returns a pair of strings representing their optimal alignment
    */
    pair<string, string> pairwiseAlign(string s1,
                                string s2,
                                bool fitToStar,
                                const Options &options);
};


#endif