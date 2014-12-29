//main.cpp
#include "main.h"

int main(int argc, char** argv)
{  
    Options myOptions;

    QuartetDistance curCompare;

    ApeTreeBuilder* at = new ApeTreeBuilder();
    Tree apeTree = at->makeTree();
    string alignmentArray[1] = {"CW"};
    string reconArray[6] = {"MLP", "MPP", "MPH", "MC", "MLP", "MLH"};

    vector<string> alignmentAlgs(alignmentArray, alignmentArray + sizeof(alignmentArray) / sizeof(alignmentArray[0]));
    vector<string> reconstructionAlgs(reconArray, reconArray + sizeof(reconArray) / sizeof(reconArray[0]));
    vector<pair<string, string> > algPairs;

    for(int i = 0; i < alignmentAlgs.size(); ++i)
    {
	for(int j = 0; j < reconstructionAlgs.size(); ++j)
	{
	    pair<string, string> curPair(alignmentAlgs[i], reconstructionAlgs[j]);
	    algPairs.push_back(curPair);
	}
    }

    for(int i = 0; i < algPairs.size(); ++i)
    {
	vector<Tree> comparisonTrees;
	pair<string, string> curPair = algPairs[i];
	for(int j = 1; j < 15; ++j)
	{
	    string curFile = to_string("realExperiment/outputTrees/random5COX1_") + to_string(j) 
		+ "_" + curPair.first + "_" + curPair.second + ".txt";
	    Tree curTree;
	    curTree.readTree(curFile);
	    comparisonTrees.push_back(curTree);
	}
	
	cout << curPair.first << " " << curPair.second << endl;

	vector<double> distances;
	for(int i = 0; i < comparisonTrees.size(); ++i)
	{
	    comparisonTrees[i].serializeTree(to_string(i));
	    distances.push_back(curCompare.getDistance(apeTree, comparisonTrees[i]));
	}
	for(int i = 0; i < distances.size(); ++i) cout << distances[i] << " ";
	cout << endl;
    }
}
