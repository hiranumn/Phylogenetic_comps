//tree.cpp
//This file contains actual implementation of boostTree class

#include "tree.h"
#include <iostream>
#include "../util/species.h"
#include <iomanip>


Tree::Tree()
{
    this->directed = false;
    this->root = 0;  
    this->speciesmap = get(Vertex_t(), this->g);
    this->weights = get(edge_weight_t(), this->g);
    this->numV = 0;
    this->numE = 0; //not being used right now, fix remove vertex if you wanna make it work.
}

Tree::Tree(const Tree &treeIn)
{
    this->directed = treeIn.directed;
    this->root = treeIn.root;
    this->speciesmap = treeIn.speciesmap;
    this->weights = treeIn.weights;
    this->numV = treeIn.numV;
    this->numE = treeIn.numE;
    this->g = Tree_graph(treeIn.g);
}


//Needs to be updated to handle edge wethis->numV = 0;ights and labels.
void Tree::serializeTree(string outputFile)
{
    ofstream output;
    output.open(outputFile.c_str());
    output << "graph\n{";
    //Everything goes here.
    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    
    int numEdge = 0;
    for(tie(vi, vi_end) = vertices(this->g); vi != vi_end; ++vi)
    {
	string vertexName = get(Vertex_t(), this->g, *vi)->name;
	output << numEdge <<  "[label=" << vertexName << "];\n";
	numEdge++;
    }
    
    for(tie(ei, ei_end) = edges(this->g); ei != ei_end; ++ei)
    {
	int id1 = source(*ei, this->g);
	int id2 = target(*ei, this->g);
	double weight = get(edge_weight_t(), this->g, *ei);
	output << id1 << " -- " << id2 << "[label=\"" << weight << "\", weight=\""<< weight <<"\"];\n";
    }
    
    output << "}";
    output.close();
}


void Tree::printTree()
{
    graph_traits<Tree_graph>::vertex_iterator vi1, vi_end1;
    graph_traits<Tree_graph>::vertex_iterator vi2, vi_end2;
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    int n = 0;
    
    //Print key
    cout << "-----------------------" << endl;
    cout << "Key:" << endl;
    cout << "-----------------------" << endl;
    for (tie(vi1, vi_end1) = vertices(this->g); vi1 != vi_end1; ++vi1) {
      cout << n << "\t:" << get(Vertex_t(), this->g, *vi1)->name;
      cout << "\t:" << get(Vertex_t(), this->g, *vi1)->visited;
      cout << "\t:" << get(Vertex_t(), this->g, *vi1)->rawSequence << endl;
	++n; //saves space compared to printing get(Vertex_t(), this->g, *vi1)->name
    }
    cout << "-----------------------" << endl;
    n = 0; //reset counter
    
    //Print header row
    cout << "Tree Graph:" << endl;
    cout << "-----------------------" << endl;
    cout << "\t";
    for (tie(vi1, vi_end1) = vertices(this->g); vi1 != vi_end1; ++vi1) {
	cout << "    " << n << "\t";
	++n;
    }
    cout << endl;
    //Print table border-top
    cout << "\t";
    for (int m=0; m<n; ++m){
	cout << "--------";
    }
    cout << endl;
    n = 0; //reset counter
    
    //Print 2d matrix
    for (tie(vi1, vi_end1) = vertices(this->g); vi1 != vi_end1; ++vi1) {
	cout << n << "\t";
	++n;
	for (tie(vi2, vi_end2) = vertices(this->g); vi2 != vi_end2; ++vi2) {
	    if (*vi1 == *vi2){
		cout << "\\" << "\t";
	    }
	    else if(edge(*vi1,*vi2,this->g).second){
		cout << "|" << setprecision(2) << fixed << get(edge_weight_t(), this->g, edge(*vi1,*vi2,this->g).first) << "\t";
	    }
	    else{
		cout << "|0" << "\t";
	    }
	}
	cout << "|" << endl;
    }
    //Print table border-bottom
    cout << "\t";
    for (int m=0; m<n; ++m){
	cout << "--------";
    }
    cout << endl;
    
};

int Tree::searchSpecies(Species *s)
{
  graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
  for(tie(vi, vi_end) = vertices(this->g); vi != vi_end; ++vi){
    Species *temp = get(Vertex_t(), this->g, *vi);
    if(! s->name.compare(temp->name)){
      return *vi;
    }
  }
  return -1;
};

int Tree::addNode(Species *s){
  int nodeID = add_vertex(this->g);
  put(Vertex_t(), this->g, nodeID ,s);
  numV++;
  return nodeID;
};
  
bool Tree::addEdge(Species *s1, Species *s2, float weight){
  pair<edge_descriptor, bool> e = add_edge(this->searchSpecies(s1),this->searchSpecies(s2),this->g);
  if(e.second == true){
    put(edge_weight_t(), this->g, e.first, weight);
  }
  return e.second;
};

//not implemented yet. 
bool Tree::isFullyConnected(){
  cout << "i dont know if this graph fully connected lol" << endl;
  return false;
}

vertex_descriptor Tree::divideEdge(Species *dmy, edge_descriptor e){

  edge_descriptor newEdge1;
  edge_descriptor newEdge2;
  float oldWeight = 0;
  float newWeight1 = 0;
  float newWeight2 = 0;

  int midV = add_vertex(this->g);
  //////////THIS IS HOW U ACTIALLY WANNA DO IT////////////
  put(Vertex_t(), this->g, midV, dmy);

  ///GOD DAMN IT WHY IS THIS HAPPENING
  //cout << dmy << endl;
  //cout << get(Vertex_t(), this->g, midV) << endl;

  newEdge1 = add_edge(midV ,source(e, this->g), this->g).first;
  newEdge2 = add_edge(midV ,target(e, this->g), this->g).first;
 
  oldWeight = get(edge_weight_t(), this->g, e);
  newWeight1 = oldWeight/2; 
  newWeight2 = oldWeight - newWeight1;

  put(edge_weight_t(), this->g, newEdge1, newWeight1);
  put(edge_weight_t(), this->g, newEdge2, newWeight2);

  remove_edge(source(e, this->g),target(e, this->g), this->g);

  return midV;
}


void Tree::branchNodeFromEdge(Species* newS, edge_descriptor e){
  //make an empty species and insert that into a branch
  Species* dmy = new Species("dmy","");
  int midV = divideEdge(dmy, e);


  //make another vertex and associate that with a species. 
  int newV = add_vertex(this->g);
  put(Vertex_t(), this->g, newV, newS);

  //add an edge between the first and the second vertex
  pair<edge_descriptor, bool> newEdge;
  newEdge = add_edge(midV , newV, this->g);

  //set the weight of the edge.
  float arbitraryNewBranchWeight = 5.00; 
  ////// I NEED TO SEE IF I CAN ACTUALLY DO THIS ///////
  put(edge_weight_t(), this->g, newEdge.first, arbitraryNewBranchWeight);

  //not returning anything for now.
}

void Tree::removeNode(Species *s1)
{
  clear_vertex(this->searchSpecies(s1), this->g);
  remove_vertex(this->searchSpecies(s1), this->g);
  numV--;
};

void Tree::removeEdge(Species *s1, Species *s2)
{
  remove_edge(this->searchSpecies(s1), this->searchSpecies(s2), this->g);
};


//Call this after using the visited boolean.
void Tree::markAllUnvisited(){
  graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
  for(tie(vi, vi_end) = vertices(this->g); vi != vi_end; ++vi) {
    get(Vertex_t(), this->g, *vi)->visited = false;
  }
}
//int Tree::countUnvisitedNeighbors(vertex_descriptor v){
//
//}

//test function 
void Tree::test()
{
  //adding a vertex and mapping it to a species.
  Species *mys = new Species("cat", "CAGCC");
  put(this->speciesmap, add_vertex(this->g) , mys);
  
  //accessing to the associated species of the vertex 0.
  cout << speciesmap[0]->name << endl;
  //or you can also do 
  cout << get(Vertex_t(), this->g, 0)->name << endl;

  //adding another one
  Species *mys2 = new Species("dog", "TATGGG");
  put(Vertex_t(), this->g, add_vertex(this->g) , mys2);

  //adding another one
  Species *mys3 = new Species("pig", "AAAACC");
  put(Vertex_t(), this->g, add_vertex(this->g) , mys3);

  //adding an edge between two nodes
  pair<edge_descriptor, bool> e = add_edge(0,1,this->g);
  if(e.second == true){
    put(edge_weight_t(), this->g, e.first, 100);
  }
  
  //printing out the weight
  cout << get(edge_weight_t(), this->g, edge(0,1,this->g).first) << endl;
  
  //adding another edge
  pair<edge_descriptor, bool> e2 = add_edge(1,2,this->g);
  if(e2.second == true){
    put(edge_weight_t(), this->g, e2.first, 50);
  } 
}
