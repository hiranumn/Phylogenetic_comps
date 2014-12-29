/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* distanceMatrix.cpp
* 2d matrix for housing distances between sequence pairs
*/

#include "distanceMatrix.h"

/*
    Creation using another distance matrix
*/
DistanceMatrix::DistanceMatrix(map<string, map<string, double> > matrixIn)
{
    this->matrix = matrixIn;
}

/*
    Creation using a set of sequences, fills in values by the aligned sequences if available
*/
DistanceMatrix::DistanceMatrix(const SequenceData &data)
{
    typedef map<string, Species*>::iterator SpeciesIterator;
    
    //a vector of pairs of species that we need to calculate
    //the distance between
    vector<pair<string, string> > distanceToCalculate;
    
    //a vector of all the species names in our dataset.
    vector<string> species = data.getAllSpeciesNames();
    
    
    for(int i = 0; i < species.size(); ++i)
    {
	for(int j = i+1; j < species.size(); ++j)
	{
	    pair<string, string> newPair;
	    newPair.first = species[i];
	    newPair.second = species[j];
	    distanceToCalculate.push_back(newPair);
	}
    }

    for(int i = 0; i < distanceToCalculate.size(); ++i)
    {
	double curDist;
	Species *specOne, *specTwo;
	specOne = data.speciesMap.at(distanceToCalculate[i].first);
	specTwo = data.speciesMap.at(distanceToCalculate[i].second);
	if(data.aligned)
	{
	    curDist = alignedDistance(specOne->alignedSequence, specTwo->alignedSequence);
	}
	else
	{
	    curDist = rawDistance(specOne->rawSequence, specTwo->rawSequence);
	}
	matrix[distanceToCalculate[i].first][distanceToCalculate[i].second] = curDist;
	matrix[distanceToCalculate[i].second][distanceToCalculate[i].first] = curDist;
    }
}

/*
    Printing utility function
*/
void DistanceMatrix::printMatrix() const
{
    typedef map<string, map<string, double> >::const_iterator OutsideMapIter;
    typedef map<string, double>::const_iterator InsideMapIter;

    for(OutsideMapIter iter1 = this->matrix.begin();
	iter1 != this->matrix.end();
	++iter1)
    {
	cout << iter1->first << " ";
	for(InsideMapIter iter2 = iter1->second.begin();
	    iter2 != iter1->second.end();
	    ++iter2)
	{
	    cout << iter2->second << " ";
	}
	cout << endl;
    }
}

/*
    Returns a vector of all unique species combinations
*/
vector<pair<string, string> > DistanceMatrix::getUniquePairs() const
{
    typedef map<string, map<string, double> >::const_iterator MatrixIterator;
    vector<string> species;
    for(MatrixIterator matIter = this->matrix.begin();
	matIter != this->matrix.end();
	++ matIter)
    {
	species.push_back(matIter->first);
    }
	
    vector<pair<string, string> >speciesPair;
    for(int i = 0; i < species.size(); ++i)
    {
	for(int j = i+1; j < species.size(); ++j)
        {
            pair<string, string> newPair;
            newPair.first = species[i];
            newPair.second = species[j];
            speciesPair.push_back(newPair);
        }
    }
    return speciesPair;
}

/*
    Utility function, measures the distance between two equal length sequences, returns a double
*/
double DistanceMatrix::alignedDistance(string seq1, string seq2)
{
    double dist = 0;
    
    for (int i = 0; i < seq1.size(); ++i)
    {
	if (seq1[i] != seq2[i]) dist += 1;
    }
    
    return dist;
}

/*
    Utility function, measures the distances between two unequal length sequences, returns a double

    Currently derelicht, not completed / in use
*/
double DistanceMatrix::rawDistance(string seq1, string seq2)
{
    if(seq1.size() != seq2.size())
    {
	return 666;
    }
    double dist = 0;
    for (int i = 0; i < seq1.size(); ++i)
    {
	if(seq1[i] != seq2[i]) dist+=1;
    }
    return dist;
}


/*
    Merges two items in the matrix and reduces the size accordingly
*/
void DistanceMatrix::mergePair(pair<string, string> mergePair, string newName)
{
    map<string, double> newDistMap;
    typedef map<string, map <string, double> >::iterator MapIterator;
    for (MapIterator mapIter = this->matrix.begin();
	 mapIter != this->matrix.end();
	 ++mapIter)
    {
	double newDist;
	if (mapIter->first.compare(mergePair.first) != 0 && mapIter->first.compare(mergePair.second) != 0)
	{
	    newDist = (this->matrix[mergePair.first][mapIter->first] +
		       this->matrix[mergePair.second][mapIter->first] -
		       this->matrix[mergePair.first][mergePair.second]) / 2;
	    newDistMap[mapIter->first] = newDist;
	}
    }
    
    this->matrix.erase(mergePair.first);
    this->matrix.erase(mergePair.second);

    for(MapIterator mapIter = this->matrix.begin();
	mapIter != this->matrix.end();
	++mapIter)
    {
	mapIter->second.erase(mergePair.first);
	mapIter->second.erase(mergePair.second);
    }

    //add in the new entries
    for(map<string, double>::iterator myIter = newDistMap.begin();
	myIter != newDistMap.end();
	++myIter)
    {
	this->matrix[newName][myIter->first] = myIter->second;
	this->matrix[myIter->first][newName] = myIter->second;
    }
}
