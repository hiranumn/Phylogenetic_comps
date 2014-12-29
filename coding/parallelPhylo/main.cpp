//main.cpp
#include "main.h"

int main(int argc, char** argv)
{
    Options myOptions;
    if(argc == 2)
    {
	myOptions.readFromFile(argv[1]);
    }

    string filename = "../generatedData/RealData.txt";
    
    SequenceData myData(filename);
    AbstractAligner *aligner = new ClustalWParallel();
    aligner->align(myData, myOptions);
    delete aligner;
    
    myData.printAlignment();

    AbstractTreeBuilder *seqBuilder = new MaximumLikelihoodClimb();
    Tree result = seqBuilder->solve(myData, myOptions);
    result.serializeTree("seq.dot");
    delete seqBuilder;
}
