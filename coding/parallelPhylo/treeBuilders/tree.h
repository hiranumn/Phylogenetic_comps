//tree.h
//This file contains the non-abstract class tree, which holds the topology and other attributes of a tree. 

#ifndef TREE
#define TREE

#include <string>
#include "../util/species.h"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graphviz.hpp"
#include <string.h>
#include <stdio.h>

using namespace std;
using namespace boost;

//MAY HAVE TO MOVE THIS SECTION INTO TREE CLASS
//--------------------------------------------
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
    // constructor
    Tree();

    //Copy constructor
    //makes a copy of tree in
    Tree(const Tree &treeIn);

    
    //attributes
    bool directed;
    int root;
    Tree_graph g;
    int numV;
    int numE;
    
    //propertymaps of tree_graph
    property_map<Tree_graph, Vertex_t>::type speciesmap;
    property_map<Tree_graph, edge_weight_t>::type weights;
    
    //functions
    void test();
    void printTree();
    bool isFullyConnected();
    int searchSpecies(Species *s);
    void serializeTree(string filename);
    void markAllUnvisited();
    
    
    //these functions might be too slow idk man we need to talk about this.
    int addNode(Species *s);
    bool addEdge(Species *s1, Species *s2, float weight);
    void removeEdge(Species *s1, Species *s2);
    void removeNode(Species *s1);
    vertex_descriptor divideEdge(Species *dmy,  edge_descriptor e);
    void branchNodeFromEdge(Species* newS, edge_descriptor e);
    //int Tree::countUnvisitedNeighbors(vertex_descriptor v);

};

#endif

