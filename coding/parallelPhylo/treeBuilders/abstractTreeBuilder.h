//abstractTreeBuilder.h
//This file contains the abstract class that the tree reconstruction algorithms
//will subclass.

#ifndef ABSTREEBUILD
#define ABSTREEBUILD

#include "../util/sequenceData.h"
#include "../util/options.h"
#include "../util/species.h"
#include "tree.h"

class AbstractTreeBuilder
{
public:
    //The vector of temporary species this particular solver has created.
    //Each of these "dummy" species needs to be freed, probably when the
    //destructor is called.
    vector<Species*> dummies;

    //Given an aligned SequenceData object and an Options object,
    //constructs a pair<Tree, vector<
    virtual Tree solve(const SequenceData& sd, Options opt) = 0;

    virtual ~AbstractTreeBuilder()
    {
	for(int i = 0; i < dummies.size(); ++i)
	    delete dummies[i];
    };
    
};


#endif
