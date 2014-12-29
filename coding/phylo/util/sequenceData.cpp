/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* sequenceData.cpp
* widely used non-abstract class storing species and the sequences associated
*/

#include "sequenceData.h"
#include <iostream>
#include <fstream>

/* 
    sequence data constructor class, calls helper function to load information in from file
*/
SequenceData::SequenceData(string filename)
{
    this->aligned = false;
    this->readFromFile(filename);
}

/*
    base constructor
*/
SequenceData::SequenceData()
{
    this->aligned = false;
}

/*
    returns a vector of all the names of species stored in the object
*/
vector<string> SequenceData::getAllSpeciesNames() const
{
    vector<string> returnVec;

    typedef map<string, Species*>::const_iterator MyIter;
    for(MyIter iter = this->speciesMap.begin();
	iter != this->speciesMap.end();
	++iter)
    {
	returnVec.push_back(iter->first);
    }
    return returnVec;
}

/*
    returns the first name stored in the object, useful for getting sequence lengths
    or checking if the object is empty
*/
string SequenceData::getFirstSpeciesName() const
{

    typedef map<string, Species*>::const_iterator MyIter;
    for(MyIter iter = this->speciesMap.begin();
	iter != this->speciesMap.end();
	++iter)
    {
	return iter->first;
    }
    return "NONE";
}

/*
    prints the sequence data for all species
*/
void SequenceData::printData()
{
    typedef map<string, Species*>::iterator MyIter;
    for(MyIter iter = this->speciesMap.begin();
	iter != this->speciesMap.end();
	++iter)
    {
	cout << iter->first << endl;
	cout << (iter->second)->rawSequence << endl;
    }
}

/*
    prints the aligned sequences for all species
*/
void SequenceData::printAlignment()
{
    typedef map<string, Species*>::iterator MyIter;
    for(MyIter iter = this->speciesMap.begin();
	iter != this->speciesMap.end();
	++iter)
    {
	cout << (iter->second)->alignedSequence << endl;
    }
}

/*
    saves the sequences to the designated file
*/
void SequenceData::saveAlignment(string output)
{
    using namespace std;
    ofstream out(output.c_str());
    typedef map<string, Species*>::iterator MyIter;
    for(MyIter iter = this->speciesMap.begin();
	iter != this->speciesMap.end();)
    {
	out << (iter->first) << endl;
	out << (iter->second)->alignedSequence;
	++iter;
	if(iter != this->speciesMap.end()) out << endl;
    }
}

/*
    reads in aligned sequences from a designated file
*/
void SequenceData::readAlignedFromFile(string filename)
{
    using namespace std;
    ifstream inputStream(filename.c_str());
    if(inputStream.is_open())
    {
	string speciesName;
	while(getline(inputStream, speciesName))
	{
	    string sequence;
	    getline(inputStream, sequence);
	    Species *newSpecies;
	    this->alignedLength = sequence.size();
	    newSpecies = new Species(speciesName, sequence, sequence);
	    this->speciesMap[speciesName] = newSpecies;
	}
	this->aligned = true;
    }
    else
    {
	cout << filename << " has not been found!!" << endl;
    }
}

/* 
    Load in sequence data from a specified file
    current format is 2 line pairs, started with a species name
    and followed by the sequence of the species
*/
void SequenceData::readFromFile(string filename)
{
    using namespace std;
    ifstream inputStream(filename.c_str());
    if(inputStream.is_open())
    {
	string speciesName;
	while(getline(inputStream, speciesName))
	{
	    string sequence;
	    getline(inputStream, sequence);
	    Species *newSpecies;
	    newSpecies = new Species(speciesName, sequence);
	    this->speciesMap[speciesName] = newSpecies;
	}
    }
    else
    {
	cout << filename << " has not been found!!" << endl;
    }
}

/*
    returns the total distance of all sequences as a double, useful for evaluating alignments
*/
double SequenceData::totalDistance() 
{
    double totalDistance = 0.0;
    vector<string> names = this->getAllSpeciesNames();
    for(int i = 0; i < names.size(); ++i)
    {
        for(int j = i + 1; j < names.size(); ++j)
        {
            int numNonGapPositions = 0;
            int numMatches = 0;
            string seq1 = this->speciesMap[names[i]]->alignedSequence;
            string seq2 = this->speciesMap[names[j]]->alignedSequence;
            for(int k = 0; k < seq1.length(); ++k)
            {
                if(seq1[k] != '-' && seq2[k] != '-') 
                {
                    numNonGapPositions++;
                    if(seq1[k] == seq2[k]) 
                        numMatches++;
                }
            }
            totalDistance += (1 - 1.0*numMatches/numNonGapPositions) * seq1.length();
        }
    }
    return  totalDistance;
}

/*
    Sequence data deconstructor
*/
SequenceData::~SequenceData()
{
    typedef map<string, Species*>::iterator MapIterator;
    for(MapIterator myIter = this->speciesMap.begin();
	myIter != this->speciesMap.end();
	++myIter)
    {
	delete myIter->second;
    }
}
