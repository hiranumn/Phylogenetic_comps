/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* mutationMap.cpp
* Creates a map translating nucleic characters to amino acids.
*/

#include "aminoAlphabet.h"

using namespace std;


/*
  Maps sequences of nucleotides to the corresponding amino acid.
  all 64 combinations (of 4 letter alphabet, with 3 chosen) are represented so no lookup error should be possible
  only works on the ATGC character set, no '-' characters allowed! or U for like RNA stuff
*/
AminoAlphabet::AminoAlphabet()
{
  //Translations
  //3
  this->amnAlpha["ATT"]="I";
  this->amnAlpha["ATC"]="I";
  this->amnAlpha["ATA"]="I";

  //9
  this->amnAlpha["CTT"]="L";
  this->amnAlpha["CTC"]="L";
  this->amnAlpha["CTA"]="L";
  this->amnAlpha["CTG"]="L";
  this->amnAlpha["TTA"]="L";
  this->amnAlpha["TTG"]="L";

  //13
  this->amnAlpha["GTT"]="V";
  this->amnAlpha["GTC"]="V";
  this->amnAlpha["GTA"]="V";
  this->amnAlpha["GTG"]="V";

  //15
  this->amnAlpha["TTT"]="F";
  this->amnAlpha["TTC"]="F";

  //16
  this->amnAlpha["ATG"]="M"; //STARTING CODON

  //18
  this->amnAlpha["TGT"]="C";
  this->amnAlpha["TGC"]="C";

  //22
  this->amnAlpha["GCT"]="A";
  this->amnAlpha["GCC"]="A";
  this->amnAlpha["GCA"]="A";
  this->amnAlpha["GCG"]="A";

  //26
  this->amnAlpha["GGT"]="G";
  this->amnAlpha["GGC"]="G";
  this->amnAlpha["GGA"]="G";
  this->amnAlpha["GGG"]="G";

  //30
  this->amnAlpha["CCT"]="P";
  this->amnAlpha["CCC"]="P";
  this->amnAlpha["CCA"]="P";
  this->amnAlpha["CCG"]="P";

  //34
  this->amnAlpha["ACT"]="T";
  this->amnAlpha["ACC"]="T";
  this->amnAlpha["ACA"]="T";
  this->amnAlpha["ACG"]="T";

  //40
  this->amnAlpha["TCT"]="S";
  this->amnAlpha["TCC"]="S";
  this->amnAlpha["TCA"]="S";
  this->amnAlpha["TCG"]="S";
  this->amnAlpha["AGT"]="S";
  this->amnAlpha["AGC"]="S";

  //42
  this->amnAlpha["TAT"]="Y";
  this->amnAlpha["TAC"]="Y";

  //43
  this->amnAlpha["TGG"]="W";

  //45
  this->amnAlpha["CAA"]="Q";
  this->amnAlpha["CAG"]="Q";

  //47
  this->amnAlpha["AAT"]="N";
  this->amnAlpha["AAC"]="N";

  //49
  this->amnAlpha["CAT"]="H";
  this->amnAlpha["CAC"]="H";

  //51
  this->amnAlpha["GAA"]="E";
  this->amnAlpha["GAG"]="E";

  //53
  this->amnAlpha["GAT"]="D";
  this->amnAlpha["GAC"]="D";

  //55
  this->amnAlpha["AAA"]="K";
  this->amnAlpha["AAG"]="K";

  //61
  this->amnAlpha["CGT"]="R";
  this->amnAlpha["CGC"]="R";
  this->amnAlpha["CGA"]="R";
  this->amnAlpha["CGG"]="R";
  this->amnAlpha["AGA"]="R";
  this->amnAlpha["AGG"]="R";

  //64
  this->amnAlpha["TAA"]="!";
  this->amnAlpha["TAG"]="!";
  this->amnAlpha["TGA"]="!";
  

  // usage example
  //cout << this->amnAlpha["AAG"] << endl;
  
}
