//options.h
//This file contains the non-abstract class options, which holds the data
//relating to runtime options: alignment options, reconstruction options,
//etc. The options.h file contains default values for all these variables.

#ifndef OPTIONS
#define OPTIONS

#include <string>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

class Options
{
public:
    Options();
    Options(string optionsFileIn);
    
    void readFromFile(string optionsFileIn);
    
    //General Options
    int verbosity;
    int maxNumThreads;

    //Maximum Likelihood
    double uRate;
    double convergence;
    int alphabetSize;
    vector<char> nuc;
    map<char, double> freq;
    
    //Random Data Generator
    int numS;
    int maxBranchLength;
    int randSpeciesDNALength;
    string nameBase;
    double rndRate;
    double insRate; 
    double delRate; 

    //Maximum parsimony
    
    //Clustal W
    double gapReward; //reward for putting in a gap in pairwise alignment
    double matchReward; //reward for matching in pairwise alignments
    
    double profileGapReward; //reward for inserting a gap when aligning profiles

    vector<char> alphabet; //alphabet of characters
    map<char, map<char, double> > similarities; //character similarities for profile similarities

    //Muscle
    int kmerSize; // size of "words" to count frequencies of in a sequence

private:
    void setValue(string key, string value);
    void setDictionaryValue(string key, string value);
    void setSingularValue(string key, string value);
    void setTwoDDictionaryValue(string key, string value);
    void setVectorValue(string key, string value);

    //given a string, and a character to split on
    //tokenizes txt on ch, returns vector of tokens
    vector<string> split(const string& txt, char ch);
};

#endif
