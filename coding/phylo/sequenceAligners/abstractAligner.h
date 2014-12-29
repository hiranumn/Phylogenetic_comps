/*
* abstractAligner.h
* This file contains the abstract class that the sequence alignment algorithm classes
* will subclass.
*/

#ifndef ABSALIGN
#define ABSALIGN

#include "../util/sequenceData.h"
#include "../util/options.h"

class AbstractAligner
{
public:

    /*
      Given a sequence data object, set all alignedSequence
      variables in the data's species objects, and set
      alignedLength and aligned within the sequenceData object.
    */
    virtual void align(SequenceData &data, const Options &options) = 0;
    
    virtual ~AbstractAligner() {};
    
};


#endif
