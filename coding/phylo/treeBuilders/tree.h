/*
 * tree.h
 * This is a tree class that holds any tree structure related information.
 * It is implemented using boost graphs. 
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */
#ifndef TREE
#define TREE

#include <string>
#include "../util/species.h"
#include "../util/sequenceData.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graphviz.hpp"
#include <string.h>
#include <stdio.h>

using namespace std;
using namespace boost;

// defining edge weight property in float.  
typedef property<edge_weight_t, float> edgeWeight;

// defining custom vertex property with a species pointer
struct Vertex_t {
  typedef vertex_property_tag kind;
}; 
typedef property<Vertex_t, Species*> Vertex_prop;

// specifying the properties of adjacency list
typedef adjacency_list<vecS, vecS, undirectedS, Vertex_prop, edgeWeight> Tree_graph;

typedef graph_traits<Tree_graph>::vertex_descriptor vertex_descriptor;
typedef graph_traits<Tree_graph>::edge_descriptor edge_descriptor;
//--------------------------------------------

class Tree
{
 public:
  /*
   * A constructor that makes an empty tree.
   */
  Tree();

    
  /*
   * A constructor that copies a given tree.
   * args: a const refernce to the original tree.
   */    
  Tree(const Tree &treeIn);

    
  //attributes
  bool directed;
  int root;
  Tree_graph g;
  int numV;
  int numE;
    
  /*
   * Iterates through the graph and returns the integer ID of the node that contains the input species.  
   */
  int searchSpecies(Species *s);

  /*
   * converts the tree into .dot format and write to a specified file location. 
   * args: output file name. 
   */  
  void serializeTree(string filename);

  /*
   * converts the tree into .txt format and write to a specified file location. 
   * args: output file name. 
   */  
  void writeTree(string filename);

  /*
   * reads a .txt format tree and creates a tree object accordingly.
   * make sure to use this function on an empty tree. 
   * args: input file name. 
   */  
  void readTree(string filename);

  /*
   * Retrieves the sequence data stored in the tree.
   * return: SequenceData object that contains the names and sequences of leaf nodes (non-dummy species)
   */  
  SequenceData* getData();

  /*
   * Prints the tree structure to stdout.
   */  
  void printTree();
    
  /*
   * Add a node that contains a given species.
   * Returns integer ID of the added node as an interger
   */
  int addNode(Species *s);

  /*
   * Add an edge between two species and set its weight.
   * Return a boolean that tells if the addtion of the edge was successful
   */
  bool addEdge(Species *s1, Species *s2, float weight);

  /*
   * Remove the edge between given two species 
   */
  void removeEdge(Species *s1, Species *s2);

  /*
   * Remove the given species from the graph. 
   */
  void removeNode(Species *s1);

  /*
   * divides an edge by inserting a new species.
   * args: a species to insert, an edge to divide
   * return: returns an interger ID of the added species. 
   */
  vertex_descriptor divideEdge(Species *dmy,  edge_descriptor e);

  /* 
   *  Given an edge and new species newS, this function adds newS to the specififed edge
   */
  void branchNodeFromEdge(Species* newS, edge_descriptor e);

};

#endif

