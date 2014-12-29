/*
 * tree.cpp
 * This is a tree class that holds any tree structure related information.
 * It is implemented using boost graphs. 
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#include "tree.h"
#include <iostream>
#include "../util/species.h"
#include <iomanip>

#include <string>
#include <sstream>
#include <iterator>


Tree::Tree()
{
    this->directed = false;
    this->root = 0;  
    this->numV = 0;
    this->numE = 0; //not being used right now, fix remove vertex if you wanna make it work.
}

Tree::Tree(const Tree &treeIn)
{
    this->directed = treeIn.directed;
    this->root = treeIn.root;
    this->numV = treeIn.numV;
    this->numE = treeIn.numE;
    this->g = Tree_graph(treeIn.g);
}

void Tree::serializeTree(string outputFile)
{
    ofstream output;
    output.open(outputFile.c_str());
    output << "graph\n{";
    //Everything goes here.
    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    
    int numEdge = 0;

    //make nodes
    for(tie(vi, vi_end) = vertices(this->g); vi != vi_end; ++vi)
    {
	string vertexName = get(Vertex_t(), this->g, *vi)->name;
	output << numEdge <<  "[label=" << vertexName << "];\n";
	numEdge++;
    }
    
    //make edges
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


//Modified version of serializeTree, writes the tree to a better save format including all important information
void Tree::writeTree(string outputFile)
{
  ofstream output;
  output.open(outputFile.c_str());
  output << "ADJACENCY MATRIX\n";

  graph_traits<Tree_graph>::edge_iterator ei, ei_end;
  graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
  //Construct the adjacency matrix
  int dimension = num_vertices(this->g);
  double adjMatrix[dimension][dimension];
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      adjMatrix[i][j] = 0.0; //Zero should be safe as it makes no sense to have a 0 weight
    }
  }

  //Fill the matrix
  for(tie(ei, ei_end) = edges(this->g); ei != ei_end; ++ei)
  {
    int id1 = source(*ei, this->g);
    int id2 = target(*ei, this->g);
    double weight = get(edge_weight_t(), this->g, *ei);
    adjMatrix[id1][id2] = weight;
    //adjMatrix[id2][id1] = weight; // fill in the edge weights at both sites USE THIS IF YOU WANT TO ENABLE LEAF/INTERNAL IDENTITIES
  }


  //Output the matrix
  //Header line
  for (int i = 0; i < dimension; i++) {
    output << "\t" << i;
  }
  output << "\n";
  //Matrix contents
  //Write the contents and identify which species are leaf or internal nodes by number of edges
  for (int i = 0; i < dimension; i++) {
    int numEdges = 0;
    output << i << "\t";
    for (int j = 0; j < dimension; j++) {
      if (i == j){
        output << "\\\t";   
      }
      else {
        if (adjMatrix[i][j] != 0) {
          numEdges++;
        }
        output << adjMatrix[i][j] << "\t"; //Some formatting could be done to make the matrix human readable as long doubles make things ugly
      }
    }
    output << "\n";
  }

  output << "\n\n\n";
  output << "SPECIES INFORMATION\n";
  //Write all the info from the species
  for(tie(vi, vi_end) = vertices(this->g); vi != vi_end; ++vi)
  {
    Species* tempSpecies = get(Vertex_t(), this->g, *vi);
    string vertexName = tempSpecies->name;
    string status = "UNKNOWN";
    if (!tempSpecies->dummy)
      status = "LEAF";
    else
      status = "INTERNAL";
    string rawSeq = tempSpecies->rawSequence;
    string aliSeq = tempSpecies->alignedSequence;

    output << "ID:\t" << *vi << "\n";
    output << "name:\t" << vertexName << "\n";
    output << "status:\t" << status << "\n"; //status << "\n"; ENABLE THE ABOVE LINE THAT MENTIONS STATUS, CURRENTLY THIS WON'T WORK. KEEPING THIS OFF HELPS EDGE SOURCE-TARGET IDENTIFYING
    output << "rawSeq:\t" << rawSeq << "\n";
    output << "aliSeq:\t" << aliSeq << "\n";
    output << "\n"; // makes the file slightly more legible
  }

    output.close();
}

void Tree::readTree(string inputFile)
{
  ifstream inputStream(inputFile.c_str());
  if(inputStream.is_open())
  {
    //File line
    string line;
    //Edge connection info
    vector<string> headerLine;
    //Species info
    int ID;
    string name;
    string status;
    string rawSeq = "";
    string aliSeq = "";
    //Species map id to species
    map<int, Species*> speciesMap;
    //Edge map id to another id
    map<int, map<int, double> > branchMap;
    int dimension = 0;
    //Iterate through the file, gather info on edges, add in species
    while(getline(inputStream, line))
    {
      //Breaks each line up (on whitespace) into distinct strings
      istringstream buf(line);
      istream_iterator<string> beg(buf), end;
      vector<string> tokens(beg, end);

      if (tokens.size() < 1){ //empty line
        continue;
      }

      else if (tokens[0] == "ADJACENCY") { //Header line to adjaceny matrix
        getline(inputStream, line); //move down a line
        istringstream buf2(line);
        istream_iterator<string> beg2(buf2), end2;
        vector<string> tokens2(beg2, end2);

        headerLine = tokens2;

        bool readingMatrix = true;
        while(readingMatrix) {
          getline(inputStream, line); //move down a line
          istringstream buf(line);
          istream_iterator<string> beg(buf), end;
          vector<string> tokens(beg, end);

          if (tokens.size() == 0) {
            readingMatrix = false;
          }
          else {
            dimension++;
            ID = atoi(tokens[0].c_str());
            for (int i = 1; i< tokens.size(); i++) { //skip first int (its the ID)
              branchMap[ID][atoi(headerLine[i-1].c_str())] = atof(tokens[i].c_str());
            }
          }
        }        
      }

      else if (tokens[0] == "ID:") { // start of a taxa info
        Species *newSpecies; //our species

        ID = atoi(tokens[1].c_str()); //corresponds to what is usually the *vi

        getline(inputStream, line); //move down a line
        istringstream buf1(line);
        istream_iterator<string> beg1(buf1), end1;
        vector<string> tokens1(beg1, end1);

        string name = tokens1[1];

        getline(inputStream, line); //move down a line
        istringstream buf2(line);
        istream_iterator<string> beg2(buf2), end2;
        vector<string> tokens2(beg2, end2);

        string status = tokens2[1];

        getline(inputStream, line); //move down a line
        istringstream buf3(line);
        istream_iterator<string> beg3(buf3), end3;
        vector<string> tokens3(beg3, end3);

        if (tokens3.size() > 1) {
          rawSeq = tokens3[1];
        }

        getline(inputStream, line); //move down a line
        istringstream buf4(line);
        istream_iterator<string> beg4(buf4), end4;
        vector<string> tokens4(beg4, end4);

        if (tokens4.size() > 1) {
          aliSeq = tokens4[1];
        }

        newSpecies = new Species(name, rawSeq, aliSeq);
        if (status == "INTERNAL") //identify leaves
          newSpecies->dummy = true;
        speciesMap[ID] = newSpecies; //save ID that matches it to the adjacency matrix
        this->addNode(newSpecies); //add it to the tree
      }

      else { //not worth picking whatever is written down
        continue;
      } 
    }

    //File is read in, add in edges now to correct species
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
        if (branchMap[i][j] > 0) {
          this->addEdge(speciesMap[i],speciesMap[j],branchMap[i][j]);
        }
      }
    }
  }
  else
  {
    cout << "Cannot load tree, " << inputFile << " has not been found!!" << endl;
  }
}

SequenceData* Tree::getData()
{
  SequenceData* sqData = new SequenceData();

  graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
  for(tie(vi, vi_end) = vertices(this->g); vi != vi_end; ++vi)
  {
    Species* tempSpecies = get(Vertex_t(), this->g, *vi);
    //Excluding dummy species. 
    if (!tempSpecies->dummy) {
      sqData->speciesMap[tempSpecies->name] = tempSpecies;
    }
  }
  return sqData;
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
      cout << n << ":\t" << get(Vertex_t(), this->g, *vi1)->name << "  ";
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

vertex_descriptor Tree::divideEdge(Species *dmy, edge_descriptor e)
{
    float oldWeight = get(edge_weight_t(), this->g, e);
    float newWeight1 = oldWeight/2;
    float newWeight2 = oldWeight-newWeight1;

    remove_edge(source(e, this->g),target(e, this->g), this->g);
    vertex_descriptor midV = add_vertex(this->g);

    edge_descriptor newEdge1 = add_edge(midV, source(e, this->g), this->g).first;
    edge_descriptor newEdge2 = add_edge(midV, target(e, this->g), this->g).first;
    
    put(Vertex_t(), this->g, midV, dmy);
    put(edge_weight_t(), this->g, newEdge1, newWeight1);
    put(edge_weight_t(), this->g, newEdge2, newWeight2);
    
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
  float arbitraryNewBranchWeight = 0.5; 
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

