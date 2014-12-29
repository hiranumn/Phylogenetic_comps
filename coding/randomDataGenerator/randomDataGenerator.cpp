/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* randomDataGenerator.cpp
* Program for outputting a randomly generated phylogenetic tree with known relationships
* based on input from the options.
*/
#include "randomDataGenerator.h"

using namespace std;

/* 
  constructer class to store all the options specified for the random data generator
  also allows the values to be altered as desired if the defaults set in the options file are not acceptable
*/
RandomDataGenerator::RandomDataGenerator(Options opt)
{
  this->numS = opt.numS;                                  // number of species to construct tree for
  this->maxBranchLength = opt.maxBranchLength;            // upper limit for the random branch length generating
  this->randSpeciesDNALength = opt.randSpeciesDNALength;  // length of the dna 
  this->nameBase = opt.nameBase;                          // naming convention for species
  this->scount = 1;                                       // keeps track of the current species name number
  this->uRate = opt.rndRate;                              // global mutation rate
  this->insRate = opt.insRate;                            // global insertion rate
  this->delRate = opt.delRate;                            // global deletion rate
  this->mMap = MutationMap();                             // houses the mutation rates across nucleotides
}

/* 
  generate a random tree given the desired attributes
  the tree will contain species nodes with no information, it is highly recommended this function call
  is immediately followed by the assignSpecies() function
*/
Tree* RandomDataGenerator::makeRandomTree(){
  // make an empty tree
  Species* newS;
  Tree *randTree = new Tree();
  
  // make the random number generator for an int
  rng = RNGType(time(0));
  uniform_int<> intgen( 1, 3628800 ); //3628800 = 10!
  variate_generator< RNGType, uniform_int<> > randint(rng, intgen);

  // make the random number generator for a float
  uniform_real<float> floatgen(0, this->maxBranchLength);
  variate_generator< RNGType, uniform_real<float> > randfloat(rng, floatgen);

  // if numS is 0 or less, return an empty tree
  if(this->numS < 1) 
    return randTree;
  else
  {
    vertex_descriptor newID;
    // adding the first species
    newS = new Species();
    newID = add_vertex(randTree->g);
    put(Vertex_t(), randTree->g, newID, newS);
    randTree->numV++;
    if(this->numS == 1)
      return randTree;   // stop if the tree is specified as 1 species
    
    // add the second species, connect the first and second with a random weight
    newS = new Species();
    newID = add_vertex(randTree->g);
    put(Vertex_t(), randTree->g, newID, newS);
    pair<edge_descriptor, bool> e = add_edge(0,1,randTree->g);
    put(edge_weight_t(), randTree->g, e.first, randfloat());

    randTree->numV++;
    if(this->numS == 2)
      return randTree;   // stop if the tree is specified as 2 species
    
    //Add species until the counter = numS
    int counter = 2;
    int edgeCount = 1;
    int targetEdge = 0;
    float newWeight1 = 0;
    float newWeight2 = 0;
    int midV = 0;
    int newV = 0;
    Species* midS;
    pair<edge_descriptor, bool> newEdge1;
    pair<edge_descriptor, bool> newEdge2;
    pair<edge_descriptor, bool> newEdge3;
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    
    while(counter < this->numS)
    {
      targetEdge = randint()%edgeCount; // pick a random edge
      tie(ei, ei_end) = edges(randTree->g); // get the iterators.
      for (int i = 0; i != targetEdge; ++i) 
      {
	      ei++; // get the nth edge.
      }
      
      // split the edge into two sub edges with one vertex between them.
      midS = new Species();
      midV = add_vertex(randTree->g);
      put(Vertex_t(), randTree->g, midV, midS);

      newEdge1 = add_edge(midV,source(*ei, randTree->g), randTree->g);
      newEdge2 = add_edge(midV,target(*ei, randTree->g), randTree->g);

      // fresh lengths for each branch
      newWeight1 = randfloat(); 
      newWeight2 = randfloat();

      // assigning weight to divided branches
      put(edge_weight_t(), randTree->g, newEdge1.first, newWeight1);
      put(edge_weight_t(), randTree->g, newEdge2.first, newWeight2);

      // remove the old edge
      remove_edge(source(*ei, randTree->g),target(*ei, randTree->g),randTree->g);
      // now make a new vertex
      newV = add_vertex(randTree->g);
      newS = new Species();
      put(Vertex_t(), randTree->g, newV, newS);      
      newEdge3 = add_edge(midV, newV, randTree->g);
      put(edge_weight_t(), randTree->g, newEdge3.first, newWeight2);

      //increment the edge count
      edgeCount += 2;
      randTree->numV += 2;

      counter++;
    };
    return randTree;
  }
};

/* 
  Names the species in the tree given the randomized tree topology and the mutation procedure creates their content
*/
void RandomDataGenerator::assignSpecies(Tree* t){
  if(t->numV){
    vd v = vertex(0,t->g);    // get the first node
    
    // make a starting species with a random DNA sequence of a specified length
    // random integer generator to be used for the entire recursive process
    rng = RNGType(time(0)); // TODO switch this to nonboost random? or does it really not matter?
    uniform_int<> intgen( 1, 2147483647 ); //2147483647 = 2^31 - 1 (max number for the boost rnd gen)
    variate_generator< RNGType, uniform_int<> > randint(rng, intgen);
    string dna = "";
    int counter = 0;
    while(counter < this->randSpeciesDNALength){
      switch(randint()%4){
      case 0:
      	dna += "A";
      	break;
      case 1:
      	dna += "T";
      	break;
      case 2:
      	dna += "C";
      	break;
      case 3:
      	dna += "G";
      	break;
      }
      counter++;
    }
    Species* source = new Species("source", dna);
    //now given v as a root, recursively fill in nodes with species using the mutate function
    map<vertex_descriptor, bool> visited;
    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
    for(tie(vi, vi_end) = vertices(t->g); vi != vi_end; ++vi)
    {
      visited[*vi] = false;
    }
    vector<Species*> leaves;
    recursivelyMutate(t, v, source, visited, &leaves);
  }
};

/* 
  given a vertex in a tree, recursively mutate across the topology to fill the tree 
*/
void RandomDataGenerator::recursivelyMutate(Tree* t, 
                                            vertex_descriptor cur,
                                            Species* source,
                                            map<vertex_descriptor, bool> visited,
                                            vector<Species*> *leaves)
{
  Species* curS = get(Vertex_t(), t->g, cur);

  // name the current vertex
  string nameS = this->nameBase + lexical_cast<string>(this->scount);
  this->scount++;
  curS->name = nameS;


  // get the vertex descriptor for source to match with cur
  int src = t->searchSpecies(source);
  // mark the current vertex visited
  visited[src] = true;

  // mutate the source to obtain the current vertex/species sequence
  float dist = 2.0;
  if(cur == 0){
    curS->rawSequence = source->rawSequence;
  }else{
    //dist is bewteen the node asociated with source, and curV. 
    curS->rawSequence = source->mutate(get(edge_weight_t(), t->g, edge(src,cur,t->g).first), this->uRate, this->insRate, this->delRate, this->mMap); 
  }
  // get the list of adjacent vertices
  // iterate through that and go to the vertex with unvisited s
  graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
  source = curS;
  int isThisLeaf = 0;
  for (tie(ai, ai_end) = adjacent_vertices(cur, t->g); ai != ai_end; ++ai) {
    if(! visited[*ai]){ //get(Vertex_t(), t->g, *ai)->visited){
      dist =  get(edge_weight_t(), t->g, edge(cur,*ai,t->g).first);
      recursivelyMutate(t, *ai, curS, visited, leaves);
    }
    isThisLeaf ++;
  }

  if(isThisLeaf == 0 || isThisLeaf == 1){ //identified as a leaf
    leaves->push_back(curS);
  }
  else { //identified as a "dummy" taxa
    curS->dummy = true;
  }
};
