//sequenceData.cpp
//This file contains implementations of functions defined in
//the SequenceData class (see sequenceData.h)
#include "sequenceData.h"
#include <iostream>
#include <fstream>

SequenceData::SequenceData(string filename)
{
  /* sequence data constructor class, calls helper function to load information in from file */
    this->aligned = false;
    this->readFromFile(filename);
}

vector<string> SequenceData::getAllSpeciesNames() const
{
  /* returns a vector containing all the names of species in the sequence data */
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

void SequenceData::printData()
{
  /* prints the sequence data for all loaded species */
    typedef map<string, Species*>::iterator MyIter;
    for(MyIter iter = this->speciesMap.begin();
	iter != this->speciesMap.end();
	++iter)
    {
	cout << iter->first << endl;
	cout << (iter->second)->rawSequence << endl;
    }
}

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

void SequenceData::readFromFile(string filename)
{
  /* load in sequence data from a specified file
     current format is 2 line pairs, started with a species name
     and followed by the sequence of the species */
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
