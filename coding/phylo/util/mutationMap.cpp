/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* mutationMap.cpp
* contains a 2 layer map which identifies evolutionary rates between nucleic characters
* essentially represents our chosen model of evolution
*/

#include "mutationMap.h"
#include <iostream>
#include <map>

using namespace std;

/*
  A nested map structure to store nucleotide-nucleotide substitution rates
  Currently using Juke's & Cantor's single variable substitution model (all rates = .25)
*/
MutationMap::MutationMap()
{
  //adenineRates;
  this->adenineRates['A']=.25; //no mutation
  this->adenineRates['T']=.25;
  this->adenineRates['G']=.25;
  this->adenineRates['C']=.25;
  //thymineRates;
  this->thymineRates['A']=.25;
  this->thymineRates['T']=.25; //no mutation
  this->thymineRates['G']=.25;
  this->thymineRates['C']=.25;
  //guanineRates;
  this->guanineRates['A']=.25;
  this->guanineRates['T']=.25;
  this->guanineRates['G']=.25; //no mutation
  this->guanineRates['C']=.25;
  //cytosineRates;
  this->cytosineRates['A']=.25;
  this->cytosineRates['T']=.25;
  this->cytosineRates['G']=.25;
  this->cytosineRates['C']=.25; //no mutation


  //map<char,map<char,float> >
  this->mutationRates['A']=this->adenineRates;
  this->mutationRates['T']=this->thymineRates;
  this->mutationRates['G']=this->guanineRates;
  this->mutationRates['C']=this->cytosineRates;

  // usage example below
  //cout << this->mutationRates['A']['G'] << endl;
  
}
