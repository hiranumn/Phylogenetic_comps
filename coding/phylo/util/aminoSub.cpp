/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* aminoSub.cpp
* reads in a 2d matrix housing "scores" of amino acid conversions, used in MUSCLE
*/

#include "aminoSub.h"
using namespace std;

aminoSub::aminoSub()
  /*
    a 2d subsitution matrix read in from a file of a strict header line and header character format
  */
{
  string filename = "util/pamvtml200.txt";
  ifstream inputStream(filename.c_str());
  if(inputStream.is_open())
  {
    string line;
    bool firstRead = true;                        // flag for catching the first line
    vector<string> firstLine;
    while(getline(inputStream, line))
    {
        istringstream buf(line);
        istream_iterator<string> beg(buf), end;
        vector<string> tokens(beg, end);
        if (firstRead) {
          firstLine = tokens;
          firstRead = false;
        }
        else {
          bool aminoChar = true;                   // flag for catching first char/string bit
          for(int i = 0; i < tokens.size(); i++){
            if (aminoChar) {
              aminoChar = false;
            }
            else {
              subRates[tokens[0]][firstLine[i-1]] = atoi(tokens[i].c_str()); //save the match to the right character character combo
            }
          }
        }
    }
  }
  else
  {
    cout << "200pamvtml.txt" << " has not been found!!" << endl;
  }
}
