/*
 * neighborJoining.cpp
 * Given an aligned sequences and options, it creates a tree topology
 * that explains the evolutionary history of sequences, using neighborjoining method.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/14/2014
 */

#include "neighborJoining.h"

Tree NeighborJoining::solve(const SequenceData &sd, DistanceMatrix distMatrix, Options opt)
{
    Tree returnTree;
    typedef map<string, Species*>::const_iterator SpeciesMapIterator;
    for(SpeciesMapIterator iter = sd.speciesMap.begin();
	iter != sd.speciesMap.end();
	++iter)
    {
	returnTree.addNode(iter->second);
    }
    
    int numTempSpecies = 0;
    
    //map to keep track of our dummy species.
    map<string, Species*> dummyMap;
    
    while(distMatrix.matrix.size() > 2)
    {
	//Maps from {tip name -> "average" distance to all other tips}
	//these are the "u" values from the video
	map<string, double> aveDists = NeighborJoining::computeAverageTipDists(distMatrix);
	pair<string, string> minimumDistPair = NeighborJoining::findClosestTips(distMatrix, aveDists);

	
	//Create a new species ID for this new "dummy species"
	ostringstream convert;
	convert << numTempSpecies;
	string dummyName = string("DUMMY_") + convert.str();
	++numTempSpecies;

	Species *dummySpecies = new Species(dummyName);
	dummySpecies->dummy = true;
	dummyMap[dummyName] = dummySpecies;
	returnTree.addNode(dummySpecies);
	this->dummies.push_back(dummySpecies);

	pair<double, double> newDistances = getDistancesToNewNode(distMatrix, minimumDistPair, aveDists);

	Species *child1, *child2;
	if(minimumDistPair.first.substr(0,6).compare("DUMMY_") == 0) child1 = dummyMap[minimumDistPair.first];
	else child1 = sd.speciesMap.at(minimumDistPair.first);

	if(minimumDistPair.second.substr(0,6).compare("DUMMY_") == 0) child2 = dummyMap[minimumDistPair.second];
	else child2 = sd.speciesMap.at(minimumDistPair.second);

	returnTree.addEdge(dummySpecies, child1, newDistances.first);
	returnTree.addEdge(dummySpecies, child2, newDistances.second);

	distMatrix.mergePair(minimumDistPair, dummyName);
    }
    

    //Merge the final two nodes. Little messy, but gets job done.
    vector<string> finalKeys;
    for(map<string, map<string, double> >::iterator finalIter = distMatrix.matrix.begin();
	finalIter != distMatrix.matrix.end();
	++finalIter)
    {
	finalKeys.push_back(finalIter->first);
    }

    Species *child1, *child2;
    if(finalKeys[0].substr(0,6).compare("DUMMY_") == 0) child1 = dummyMap[finalKeys[0]];
    else child1 = sd.speciesMap.at(finalKeys[0]);

    if(finalKeys[1].substr(0,6).compare("DUMMY_") == 0) child2 = dummyMap[finalKeys[1]];
    else child2 = sd.speciesMap.at(finalKeys[1]);
    
    returnTree.addEdge(child1, child2, distMatrix.matrix[finalKeys[0]][finalKeys[1]]);
    return returnTree;
}


Tree NeighborJoining::solve(const SequenceData &sd, Options opt)
{
    DistanceMatrix distMatrix(sd);
    return solve(sd, distMatrix, opt);
 
}

pair<double, double> NeighborJoining::getDistancesToNewNode(const DistanceMatrix &dm, 
							    pair<string, string> minimumDistPair,
							    const map<string, double> &aveDists)
{
    double avgDist1 = aveDists.at(minimumDistPair.first);
    double avgDist2 = aveDists.at(minimumDistPair.second);
    double childDist = 0.5 * dm.matrix.at(minimumDistPair.first).at(minimumDistPair.second);
    double childOneDist = childDist + 0.5 * (avgDist1 - avgDist2);
    double childTwoDist = childDist + 0.5 * (avgDist2 - avgDist1);
    pair<double, double> returnPair (childOneDist, childTwoDist);
    return returnPair;
}

pair<string, string> NeighborJoining::findClosestTips(const DistanceMatrix &distMatrix, 
						     const map<string, double> &aveDists)
{
    vector<pair<string, string> > uniquePairs = distMatrix.getUniquePairs();
    pair<string, string> curMinDistPair = uniquePairs[0];
    double curMin = std::numeric_limits<double>::infinity();
    for (int i = 0; i < uniquePairs.size(); ++i)
    {
	double dist;
	dist = distMatrix.matrix.at(uniquePairs[i].first).at(uniquePairs[i].second)
	    - aveDists.at(uniquePairs[i].first) - aveDists.at(uniquePairs[i].second);
	if (dist < curMin)
	{
	    curMin = dist;
	    curMinDistPair = uniquePairs[i];
	}
    }
    
    
    return curMinDistPair;
}

map<string, double> NeighborJoining::computeAverageTipDists(const DistanceMatrix &distMatrix)
{
    map<string, double> aveDists;
    typedef map<string, map<string, double> >::const_iterator MatrixIterator;
    typedef map<string, double>::const_iterator RowIterator;

    int curNumTips = distMatrix.matrix.size();

    for(MatrixIterator matIter = distMatrix.matrix.begin();
	matIter != distMatrix.matrix.end();
	++matIter)
    {
	double runningDistSum = 0;
	for(RowIterator rowIter = matIter->second.begin();
	    rowIter != matIter->second.end();
	    ++rowIter)
	{
	    runningDistSum += distMatrix.matrix.at(matIter->first).at(rowIter->first);
	}
	double averageDist = runningDistSum / (curNumTips-2);
	aveDists[matIter->first] = averageDist;
    }
    return aveDists;
}
