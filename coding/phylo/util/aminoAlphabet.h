/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* aminoAlphabet.h is the header file for the amino acid conversion class, which contains
* a dictionary translating nucleic characters to amino acids
*/

#ifndef AMNALPH
#define AMNALPH

#include <map>
#include <iostream>

using namespace std;

class AminoAlphabet
{
 public:
  /*
  	constructor
  */
  AminoAlphabet();

  /*
  	conversion map
  */
  map<string,string> amnAlpha;
};

#endif






