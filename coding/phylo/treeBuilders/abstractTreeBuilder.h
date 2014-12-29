/* 
 * abstractTreeBuilder.h
 * Contains the abstract tree builder class. Tree builders are used
 * to build phylogenetic trees.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/18/2014
 */

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
    
    /*
     * Given a sequence data object and an options object, return
     * a complete phylogenetic tree based on the data and options given.
     * args: a sequence data object and an options object
     * return: a phylogenetic tree
     */
    virtual Tree solve(const SequenceData& sd, Options opt) = 0;

    /*
     * Destructor, deletes the dummy nodes associated with the returned tree
     * args: none
     * return: none
     */
    virtual ~AbstractTreeBuilder()
    {
	for(int i = 0; i < dummies.size(); ++i)
	    delete dummies[i];
    };
    
};


#endif
