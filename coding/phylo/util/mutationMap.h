/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* mutationMap.h
* contains a 2 layer map which identifies evolutionary rates between nucleic characters
* essentially represents our chosen model of evolution
*/

#ifndef MUTMAP
#define MUTMAP

#include <map>

using namespace std;

class MutationMap
{
	public:
  	/*
  		constructor
  	*/
  	MutationMap();

	/*
		the mutation maps
	*/
	map<char,float> adenineRates;
	map<char,float> thymineRates;
	map<char,float> guanineRates;
	map<char,float> cytosineRates;
	map<char,map<char,float> > mutationRates;
};

#endif






