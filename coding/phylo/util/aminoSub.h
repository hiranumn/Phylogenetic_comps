/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* aminoSub.h
* Implements a 2d matrix housing "scores" of amino acid conversions, used in MUSCLE
*/

#ifndef AMNSUB
#define AMNSUB

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

class aminoSub
{
public:
	/*
  		constructor
	*/
	aminoSub();

  	/*
  		mutation maps
 	*/
 	map<string, map<string, int> > subRates;
};

#endif






