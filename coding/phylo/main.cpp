//main.cpp
#include "main.h"

int main(int argc, char** argv)
{  
    Options myOptions;

    if(argc < 4){
      cout << endl << "USAGE: ./phylo <algorithmName> <sourcefile> <outfile> <optionfile (optional)>" << endl;
      cout << "  algorithmName can be.." << endl;
      cout << "     -CW for clustalW" << endl;
      cout << "     -CS for center star" << endl;
      cout << "     -MSC for muscle" << endl;
      cout << "     -NJ for neighbor joining" << endl;
      cout << "     -MPP for maximum parsimony (progressive)" << endl;
      cout << "     -MPH for maximum parsimony (hill-climbing)" << endl;
      cout << "     -MLP for maximum likelihood (progressive)" << endl;
      cout << "     -MLH for maximum likelihood (hill-climbing)" << endl;
      cout << "     -MC for markov chain monte carlo" << endl;
      cout << "  other options are" << endl;
      cout << "     -makeApeTree for making an ape tree" << endl;
      cout << "     -visualizeTree for converting tree.txt to tree.dot" << endl << endl;
      cout << "EXAMPLE: ./phylo -MLP clustalW.aln result.dot" << endl << endl;
      return 1;
    }

    if(argc == 5)
    {
      //reading in specified option file
      myOptions.readFromFile(argv[4]);
    }
    
    //inittializing random number generator
    srand (static_cast <unsigned> (time(0)));
    //calling rand to fix the non-random issue with the rand().
    for(int i = 0; i < 10; ++i){
      rand();
    }
    
    string alg = argv[1];
    string sourcefile = argv[2];
    string outfile = argv[3];
    
    SequenceData myData;


    clock_t t1, t2;
    t1 = clock();
    
    if(alg.compare("-makeApeTree")==0){

      cout << "Making an ape tree" << endl;
      ApeTreeBuilder* at = new ApeTreeBuilder();
      Tree apetree = at->makeTree();
      apetree.serializeTree(outfile);

    }else if(alg.compare("-visualizeTree")==0){

      cout << "Visualizing the tree file" << endl;
      Tree tree;
      tree.readTree(sourcefile);
      tree.serializeTree(outfile);

    }else  if(alg.compare("-CW")==0){

      cout << "Running clustalW..." << endl;
      myData.readFromFile(sourcefile);
      AbstractAligner* cw = new ClustalW();
      cw->align(myData,myOptions);
      myData.saveAlignment(outfile);

    }else if(alg.compare("-CS")==0){

      cout << "Running center star..." << endl;
      myData.readFromFile(sourcefile);
      AbstractAligner* cs = new CenterStar();
      cs->align(myData,myOptions);
      myData.saveAlignment(outfile);

    }else if(alg.compare("-MSC")==0){

      cout << "Running MUSCLE..." << endl;
      myData.readFromFile(sourcefile);
      AbstractAligner* msc = new Muscle();
      msc->align(myData,myOptions);
      myData.saveAlignment(outfile);

    }else if(alg.compare("-NJ")==0){

      cout << "Running Neighbor Joining..." << endl;
      myData.readAlignedFromFile(sourcefile); 
      AbstractTreeBuilder* NJ = new NeighborJoining();
      Tree mcTree = NJ->solve(myData, myOptions);
      mcTree.writeTree(outfile);

    }else if(alg.compare("-MPP")==0){

      cout << "Running maximum parsimony (progressive)..." << endl;
      myData.readAlignedFromFile(sourcefile); 
      AbstractTreeBuilder* MPP = new MaximumParsimonyProgressive();
      Tree mcTree = MPP->solve(myData, myOptions);
      mcTree.writeTree(outfile);

    }else if(alg.compare("-MPH")==0){

      cout << "Running maximum parsimony (hill climb)..." << endl;
      myData.readAlignedFromFile(sourcefile); 
      AbstractTreeBuilder* MPH = new MaximumParsimonyClimb();
      Tree mcTree = MPH->solve(myData, myOptions);
      mcTree.writeTree(outfile);

    }else if(alg.compare("-MLP")==0){

      cout << "Running maximum likelihood (progressive)..." << endl;
      myData.readAlignedFromFile(sourcefile); 
      AbstractTreeBuilder* MLP = new MaximumLikelihoodProgressive();
      Tree mcTree = MLP->solve(myData, myOptions);
      mcTree.writeTree(outfile);

    }else if(alg.compare("-MLH")==0){

      cout << "Running maximum likelihood (hill climb)..." << endl;
      myData.readAlignedFromFile(sourcefile); 
      AbstractTreeBuilder* MLH = new MaximumLikelihoodClimb();
      Tree mcTree = MLH->solve(myData, myOptions);
      mcTree.writeTree(outfile);

    }else if(alg.compare("-MC")==0){

      cout << "Running markov chain monte carlo..." << endl;
      myData.readAlignedFromFile(sourcefile); 
      AbstractTreeBuilder* monteCarlo = new MonteCarlo();
      Tree mcTree = monteCarlo->solve(myData, myOptions);
      mcTree.writeTree(outfile);

    }
    t2 = clock();

    double diff = (((double)t2 - (double)t1) / 1000000.0 ) * 1000;   
    ofstream myFile;
    myFile.open((outfile + to_string(".time")).c_str());
    myFile << diff << endl;
    myFile.close();
}
