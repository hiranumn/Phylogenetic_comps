/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* main.cpp is the entry point for use of the random data generator,
* the program for outputting a randomly generated phylogenetic tree with known relationships
* based on input from the options file and/or command line (mostly the command line because of bash scripts)
*/

#include "main.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]){
  srand (time(0));
  //srand(atoi(argv[5]));
  //calling rand to fix the non-random issue with the rand().
  for(int i = 0; i < 10; ++i){
    rand();
  }

  if (argc < 5) {
    cout << "Too few arguments passed to RandomDataGenerator." << endl;
  }
  else {
    Options myOptions; //load in options from phylo/util/options
    //input cleaner for bash scripts
    //Overwrite a few from terminal
    myOptions.numS = atoi(argv[2]); // number of species to construct tree for
    myOptions.randSpeciesDNALength = atoi(argv[3]); // length of the dna 
    myOptions.rndRate = atof(argv[4]); // global mutation rate
    //myOptions.insRate = argv[3]; // global insertion rate
    //myOptions.delRate = argv[4]; // global deletion rate

    //first argument is the round identifier for creating unique names
    string fileName = "../data/generated_data/random/" + (string) argv[1] + "-" + (string) argv[2] + "-" + (string) argv[3] + "-" + (string) argv[4] + ".txt"; // number of species to construct tree for
    RandomDataGenerator rnd = RandomDataGenerator(myOptions);

    //create the topology of a tree, generate a species list based on the tree.
    Tree *myTree = rnd.makeRandomTree();
    rnd.assignSpecies(myTree);

    //save the tree's information to be read later
    myTree->writeTree(fileName);    
  }
  return 0;
  
}
