//options.cpp
//This file contains implementations of functions defined in
//the Options class (see options.h)
#include "options.h"

Options::Options()
{
    //General
    this->verbosity = 2;
    this->maxNumThreads = 256; //TODO: start using threadpool to manage max num threads

    //Maximum Likelihood
    this->convergence = 0.01;
    this->uRate = 1;//2.5e-08; // constant value set from a paper
    this->nuc.push_back('C');
    this->nuc.push_back('G');
    this->nuc.push_back('A');
    this->nuc.push_back('T');
    this->nuc.push_back('-');
    this->freq['C'] = .2;
    this->freq['G'] = .2;
    this->freq['A'] = .2;
    this->freq['T'] = .2;
    this->freq['-'] = .2;
    
    
    //Random Data Generator
    this->numS = 6; // number of species to create a tree with
    this->maxBranchLength = 10; // maximum branch length
    this->randSpeciesDNALength = 100; // DNA sequence length
    this->nameBase = "Taxa-"; // naming convention for species
    this->rndRate = 2.5e-02; // global mutation rate for random data generation
    this->insRate = 5.0e-03; // insertion event likelihood
    this->delRate = 5.0e-03; // deletion event likelhiood
    
    //Maximum parsimony
    
    //Clustal W
    this->gapReward = -1;
    this->matchReward = 2;
    this->profileGapReward = -.2;
    
    this->similarities['A']['A'] = 1;
    this->similarities['C']['C'] = 1;
    this->similarities['G']['G'] = 1;
    this->similarities['T']['T'] = 1;
    
    this->alphabet.push_back('A');
    this->alphabet.push_back('G');
    this->alphabet.push_back('C');
    this->alphabet.push_back('T');

    //Muscle
    this->kmerSize = 3; // essentially going to lead to 20 choose x possibilities
    
}

Options::Options(string optionsFileIn)
{
    readFromFile(optionsFileIn);
}


void Options::readFromFile(string optionsFileIn)
{
    ifstream myFile(optionsFileIn.c_str());
    if(myFile.is_open())
    {
	string line;
	while(getline(myFile, line))
	{
	    //deal with comments
	    if(line.size() == 0 || line.at(0) == '#')
		continue;
	    istringstream is_line(line);
	    string key;
	    if(getline(is_line, key, '='))
	    {
		string value;
		if(getline(is_line, value)) 
		    setValue(key, value);
	    }
	}
    }
    else
    {
	cout << optionsFileIn << " not found!" << endl;
    }
}

void Options::setValue(string key, string value)
{
    if(value.substr(0, 4).compare("dict") == 0)
	setDictionaryValue(key, value.substr(4));
    else if(value.substr(0, 6).compare("2ddict") == 0)
	setTwoDDictionaryValue(key, value.substr(6));
    else if(value.substr(0, 3).compare("vec") == 0)
	setVectorValue(key, value.substr(3));
    else
	setSingularValue(key, value);
}

void Options::setVectorValue(string key, string value)
{
    //cut off parens
    value = value.substr(1, value.size()-2);
    
    vector<string> tokens;
    tokens = split(value, ' ');

    //We know have a vector of tokens to put in our dictionary!

    //General
    //Maximum Likelihood
    //Random Data Generator
    //Maximum Parsimony
    //Clustal W
    if(key.compare("alphabet") == 0)
    {
	this->alphabet.clear();
	for(int i = 0; i < tokens.size(); ++i)
	{
	    this->alphabet.push_back(tokens[i].at(0));
	}
    }
	
    else
	cout << "invalid option " << value << " for " << key << " specified." << endl;
}



void Options::setTwoDDictionaryValue(string key, string value)
{
    //cut off parens
    value = value.substr(1, value.size()-2);
    
    vector<string> tokens;
    tokens = split(value, ' ');

    vector<vector<string> > entries;
    for(int i = 0; i < tokens.size(); ++i)
    {
	//cut off end parens
	string curToken = tokens[i];
	curToken = curToken.substr(1, curToken.size()-2);
	vector<string> keyVal = split(curToken, ',');
	vector<string> curEntry;
	curEntry.push_back(keyVal[0]);//key1
	curEntry.push_back(keyVal[1]);//key2
	curEntry.push_back(keyVal[2]);//value
	entries.push_back(curEntry);
    }

    //We know have a vector of (key1, key2, val) to put in our dictionary!

    //General
    //Maximum Likelihood
    //Random Data Generator
    //Maximum Parsimony
    //Clustal W
    if(key.compare("similarities") == 0)
    {
	this->similarities.clear();
	for(int i = 0; i < entries.size(); ++i)
	{
	    this->similarities[entries[i][0].at(0)][entries[i][1].at(0)] = atof(entries[i][2].c_str());
	}
    }
	
    else
	cout << "invalid option " << value << " for " << key << " specified." << endl;
}

void Options::setDictionaryValue(string key, string value)
{
    //cut off end parens
    value = value.substr(1, value.size()-2);
    
    vector<string> tokens;
    tokens = split(value, ' ');

    vector<pair<string, string> > entries;
    for(int i = 0; i < tokens.size(); ++i)
    {
	//cut off end parens
	string curToken = tokens[i];
	curToken = curToken.substr(1, curToken.size()-2);
	vector<string> keyVal = split(curToken, ',');
	pair<string, string> curPair(keyVal[0], keyVal[1]);
	entries.push_back(curPair);
    }
    //now we have a vector of (key, val) pairs to put in our dictionary!


    //General
    //Maximum Likelihood
    //Random Data Generator
    //Maximum Parsimony
    //Clustal W

//    else
    cout << "invalid option " << value << " for " << key << " specified." << endl;
}

void Options::setSingularValue(string key, string value)
{    
    //General
    if(key.compare("verbosity") == 0)
	this->verbosity = atoi(value.c_str());
    
    //Maximum Likelihood
    else if(key.compare("convergence") == 0)
	this->convergence = atof(value.c_str());

    //Random Data Generator
    else if(key.compare("numS") == 0)
	this->numS = atoi(value.c_str());
    else if(key.compare("maxBranchLength") == 0)
	this->maxBranchLength = atoi(value.c_str());
    else if(key.compare("randSpeciesDNALength") == 0)
	this->randSpeciesDNALength = atoi(value.c_str());
    else if(key.compare("nameBase") == 0)
	this->nameBase = value;
    else if(key.compare("uRate") == 0)
	this->uRate = atof(value.c_str());
    
    //Maximum Parsimony

    //Clustal W
    else if(key.compare("gapReward") == 0)
	this->gapReward = atof(value.c_str());
    else if(key.compare("matchReward") == 0)
	this->matchReward = atof(value.c_str());
    else if(key.compare("profileGapReward") == 0)
	this->profileGapReward = atof(value.c_str());

    else
	cout << "invalid option " << value << " for " << key << " specified." << endl;

}


vector<string> Options::split(const string &txt, char ch)
{
    vector<string> strs;
    
    unsigned int pos = txt.find(ch);
    unsigned int initialPos = 0;
    
    // Decompose statement
    while(pos < txt.size())
    {
	string curString = txt.substr(initialPos, pos - initialPos + 1);
	curString = curString.substr(0, curString.size()-1);
        strs.push_back(curString);
        initialPos = pos + 1;
        pos = txt.find(ch, initialPos);
    }

    // Add the last one
    strs.push_back(txt.substr(initialPos, min(pos, static_cast<unsigned int>(txt.size())) - initialPos + 1));
    return strs;
}
