//mutationMap.h
//This file constructs a 2 layer map to hold the conversion rates for each nucleic acid.

#ifndef MUTMAP
#define MUTMAP

#include <map>

using namespace std;

class MutationMap
{
 public:
  //constructor
  MutationMap();

  //mutation maps
  map<char,float> adenineRates;
  map<char,float> thymineRates;
  map<char,float> guanineRates;
  map<char,float> cytosineRates;
  map<char,map<char,float> > mutationRates;

};

#endif






