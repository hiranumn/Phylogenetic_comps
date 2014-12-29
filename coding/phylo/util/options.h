/*
 * options.h
 * This file contains the non-abstract class options, which holds the data relating to runtime options: alignment options, reconstruction options, etc. The options.h file contains default values for all these variables.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/23/2014
 */

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
    
    /*
     * Reads in the options from the options file
     * args: string representing the options file
     * return: nothing
     */
    void readFromFile(string optionsFileIn);
    
    //General Options
    int verbosity;
    int maxNumThreads;
    
    //Monte Carlo
    int iterations;
    double phi;
    double hastingsRatio;
    
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
    double convergenceValue; // convergence for improvements to the data alignment
    double mProfileGapReward;
    
    map<char, map<char, double> > mSimilarities;
    vector<char> mAlphabet;

    //Center Star
    double csGapReward;
    double csMatchReward;

private:
    /*
     * Depending on the key and value taken in, create a dict, a 2d-dict,
     * or a vector of key-value structure.
     * args: string key, string value
     * returns: nothing
     */
    void setValue(string key, string value);

    /*
     * Sets the dictionary key-value pairs
     * args: string key, string value
     * returns: nothing
     */
    void setDictionaryValue(string key, string value);

    /*
     * Sets the singular values
     * args: string key, string value
     * returns: nothing
     */
    void setSingularValue(string key, string value);

    /*
     * Sets the 2-d dictionary key-value pairs
     * args: string key, string value
     * returns: nothing
     */
    void setTwoDDictionaryValue(string key, string value);

    /*
     * Sets the vectors key-value pairs
     * args: string key, string value
     * returns: nothing
     */
    void setVectorValue(string key, string value);

    /*
     * Given a string, and a character to split on
     * tokenizes txt on ch, returns vector of tokens
     * args: const reference to a string, a split character
     * returns: vector of strings representing the tokens
     */
    vector<string> split(const string& txt, char ch);
};

#endif
