/* 
 * monteCarlo.cpp
 * Given alinged sequences and options, this class uses Monte Carlo Markov Chain reconstruction of a phylogenetic tree.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#include "monteCarlo.h"
#include <math.h>

Tree MonteCarlo::solve(const SequenceData &sd, Options opt)
{
  //put this into options later
  this->phi = opt.phi;
  this->HastingsRatio = opt.hastingsRatio;
  this->uRate = opt.uRate;
  this->nuc = opt.nuc;
  this->alphabetSize = nuc.size();
  this->freq = opt.freq;
  this->iterations = opt.iterations;

  cout << "Running for " << this->iterations << " iterations" << endl;
  cout << "phi: " << this->phi << endl;
  
  Tree curTree = getStartingTopology(sd, &(this->root));
  Tree newTree;
  double accProb = 0.0;
  for(int i = 0; i < this->iterations+1; i++){
      if(i%1 == 0){
	  cout << i << " iterations" << endl;
      }
      if(i%2 == 0){
	  newTree = globalWithMol(curTree);
      }else{
	  newTree = localWithMol(curTree);
      }
      int siteNum = sd.speciesMap.at(sd.getFirstSpeciesName())->alignedSequence.size();
      accProb = metropolisHastings(curTree, newTree, siteNum);
      if(randPerturb(0,1) < accProb){
	  curTree = newTree;
      }
  }
  return curTree;
}


Tree MonteCarlo::getStartingTopology(const SequenceData &sd, vertex_descriptor* root)
{
  Tree returnTree;
  typedef map<string, Species*>::const_iterator MapIterator;
  vector<vertex_descriptor> leafNodes;
  for (MapIterator iter = sd.speciesMap.begin();
       iter != sd.speciesMap.end();
       ++iter)
    {
      vertex_descriptor curLeaf = add_vertex(returnTree.g);
      put(Vertex_t(), returnTree.g, curLeaf, iter->second);
      leafNodes.push_back(curLeaf);
    }
    
  vertex_descriptor curJoin = leafNodes[0];
  int numDummy = 0;
  //this loop doesn't apply to the first or last leaf nodes
  for(int i = 1; i < leafNodes.size(); ++i)
    {
      //adding root
      Species* dummy;
      if(i == leafNodes.size()-1)
	{
	  dummy = new Species("root");
	  dummy->dummy = true;
	}
      else
	{
	  ostringstream convert;
	  convert << numDummy;
	  string dummyName = string("DUMMY_") + convert.str();
	  dummy = new Species(dummyName);
	  dummy->dummy = true;
	}
      this->dummies.push_back(dummy);

      vertex_descriptor curInternal = add_vertex(returnTree.g);
      if(i == leafNodes.size()-1)
	{
	  *root = curInternal;
	}
      put(Vertex_t(), returnTree.g, curInternal, dummy); 

      edge_descriptor firstEdge = add_edge(curInternal, leafNodes[i], returnTree.g).first;
      edge_descriptor secondEdge = add_edge(curInternal, curJoin, returnTree.g).first;

      put(edge_weight_t(), returnTree.g, firstEdge, i*1.0);
      put(edge_weight_t(), returnTree.g, secondEdge, 1.0);

      curJoin = curInternal;
      ++ numDummy;
    }
  return returnTree;
}

vector<edge_descriptor> MonteCarlo::getInternalEdges(const Tree &curTree)
{
  vector<edge_descriptor> internalEdges;
  graph_traits<Tree_graph>::edge_iterator ei, ei_end;
  for (tie(ei, ei_end) = edges(curTree.g); ei != ei_end; ++ei)
    {
      vertex_descriptor vertexOne = source(*ei, curTree.g);
      vertex_descriptor vertexTwo = target(*ei, curTree.g);
      graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
      int neighborCountOne = 0;
      int neighborCountTwo = 0;
      //"Internal edges" connect two nodes with exactly three neighbors each.
      for (tie(ai, ai_end) = adjacent_vertices(vertexOne, curTree.g); ai != ai_end; ++ai)
	neighborCountOne ++;
      for (tie(ai, ai_end) = adjacent_vertices(vertexTwo, curTree.g); ai != ai_end; ++ai)
	neighborCountTwo ++;
      if (neighborCountOne > 1 && neighborCountTwo > 1)
	internalEdges.push_back(*ei);
    }
  return internalEdges;
}

vector<vertex_descriptor> MonteCarlo::getInternalNodes(const Tree &curTree)
{
  vector<vertex_descriptor> internalNodes;
  graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
  graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
  for (tie(vi, vi_end) = vertices(curTree.g); vi != vi_end; ++vi)
    {
      int neighborCount = 0;
      //"Internal nodes" have more than one neighbor node.
      for (tie(ai, ai_end) = adjacent_vertices(*vi, curTree.g); ai != ai_end; ++ai)
	neighborCount ++;
      if (neighborCount > 1)
	internalNodes.push_back(*vi);
    }
  return internalNodes;
}

double MonteCarlo::randPerturb(double fMin, double fMax){
    return ((fMax-fMin)*((double)rand()/RAND_MAX))+fMin;
}

vector<vertex_descriptor> MonteCarlo::getPath(vertex_descriptor start,
					      vertex_descriptor goal,
					      const Tree &tree)
{
  vector<vertex_descriptor> predMap(num_vertices(tree.g));
  vector<vertex_descriptor> path;
    
  dijkstra_shortest_paths(tree.g, start, predecessor_map(&predMap[0]));
    
  vertex_descriptor current = goal;
  while(current!=start)
    {
      path.push_back(current);
      current=predMap[current];
    }
  path.push_back(start);
  reverse(path.begin(), path.end());

  return path;
} 
    
Tree MonteCarlo::globalWithMol(const Tree &curTree){
  //Setting the Hastings-ratio
  this->HastingsRatio = 1;


  Tree newTree(curTree);
  vector<vertex_descriptor> internalNodes = getInternalNodes(newTree);
  double randAdj;

  for (vector<vertex_descriptor>::iterator curNode = internalNodes.begin(); curNode != internalNodes.end(); ++curNode){

    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    vector<vertex_descriptor> children;
    int minpath = 999999;
    vertex_descriptor minv;
    for (tie(ai, ai_end) = adjacent_vertices(*curNode, newTree.g); ai != ai_end; ++ai)
      {
	children.push_back(*ai);
	int comppath = getPath(*ai, this->root, newTree).size();
	if (comppath < minpath){
	  minpath = comppath;
	  minv = *ai;
	}	
      }
	
    //the one closest to the root is the parent of curNode.
    vertex_descriptor parent = NULL;
    if (children.size() == 3)
      {
	for(vector<vertex_descriptor>::iterator curchild = children.begin(); curchild != children.end(); ++curchild)
	  {
	    if(*curchild == minv)
	      {
		parent = *curchild;
		children.erase(curchild);
		break;
	      }
	  }
      }
    double lowerBound = (-1)*this->phi;
    double upperBound = this->phi;
    if(parent){
      double pl = get(edge_weight_t(), newTree.g, edge(*curNode, parent, newTree.g).first); //parent b length
      if(pl < this->phi){
	lowerBound = (-1)*pl;
      }
    }

    double scl = get(edge_weight_t(), newTree.g, edge(*curNode, children[0], newTree.g).first);
    if(get(edge_weight_t(), newTree.g, edge(*curNode, children[1], newTree.g).first)< scl){
      scl = get(edge_weight_t(), newTree.g, edge(*curNode, children[1], newTree.g).first);
    }
    if(scl< this->phi){
      upperBound = scl;
    }



    randAdj = randPerturb(lowerBound, upperBound);

    for(vector<vertex_descriptor>::iterator curchild = children.begin(); curchild != children.end(); ++curchild)
      {
	edge_descriptor e = edge(*curNode, *curchild, newTree.g).first;
	double weight = get(edge_weight_t(), newTree.g, e);
	weight = weight-randAdj;
	put(edge_weight_t(), newTree.g, e, weight);
      }
    //Subtract adjustment from parent
    if (parent){
      edge_descriptor e = edge(*curNode, parent, newTree.g).first;
      double weight = get(edge_weight_t(), newTree.g, e);
      weight = weight+randAdj;
      put(edge_weight_t(), newTree.g, e, weight);
    }

  }
    
  return newTree;
}


Tree MonteCarlo::localWithMol(const Tree &curTree){

  edge_descriptor e1;
  edge_descriptor e2;

  //duplicate the input tree;
  Tree newTree(curTree);

  //Pick a random internal edge
  vector<edge_descriptor> internalEdges = getInternalEdges(newTree);
  int randNum = rand() % internalEdges.size();
  edge_descriptor randEdge = internalEdges[randNum];
  
  //find u and v
  //get two nodes joined by this edge
  vertex_descriptor v1 = source(randEdge, newTree.g);
  vertex_descriptor v2 = target(randEdge, newTree.g);
  int dist1 = getPath(v1, this->root, newTree).size();
  int dist2 = getPath(v2, this->root, newTree).size();
  
  vertex_descriptor u;
  vertex_descriptor v;
  if(dist1<dist2){
    v = v1;
    u = v2;
  }else{
    v = v2;
    u = v1;
  }
  
  vertex_descriptor a; //child of u; 
  vertex_descriptor b; //child of u;
  vertex_descriptor c; //child of v that is not u; 
  
  graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
  vector<vertex_descriptor> children;
  for (tie(ai, ai_end) = adjacent_vertices(u, newTree.g); ai != ai_end; ++ai){
    if(*ai != v){
      children.push_back(*ai);
    }
  }

  //we might want to randomize this assingment
  a = children[0];
  b = children[1];
  //case 1: v is not root

  if(dist1 > 1 && dist2 > 1){
    vertex_descriptor w; //parent of v;
    
    vector<vertex_descriptor> neighbors;
    for (tie(ai, ai_end) = adjacent_vertices(v, newTree.g); ai != ai_end; ++ai){
      if(*ai != u){
	neighbors.push_back(*ai);
      }
    }

    dist1 = getPath(neighbors[0], this->root, newTree).size();
    dist2 = getPath(neighbors[1], this->root, newTree).size();

    if(dist1 < dist2){
      w = neighbors[0];
      c = neighbors[1];
    }else{
      w = neighbors[1];
      c = neighbors[0];
    }
    //cout << "a is " << get(Vertex_t(), newTree.g, a)->name << endl;
    //cout << "b is " << get(Vertex_t(), newTree.g, b)->name << endl;
    //cout << "c is " << get(Vertex_t(), newTree.g, c)->name << endl;
    //cout << "w is " << get(Vertex_t(), newTree.g, w)->name << endl;
    
    //Calculate the h array and sort it
    vector<vertexAndPathLength> h;
    double distA = getPathLength(a, w, newTree);
    vertexAndPathLength vdA;
    vdA.vert = a;
    vdA.path = distA; 
    double distB = getPathLength(b, w, newTree);
    vertexAndPathLength vdB;
    vdB.vert = b;
    vdB.path = distB;
    double distC = getPathLength(c, w, newTree);
    vertexAndPathLength vdC;
    vdC.vert = c;
    vdC.path = distC;
    h.push_back(vdA);
    h.push_back(vdB);
    h.push_back(vdC);
    

    //sorting the vector. make sure you include <algorithm>
    sort(h.begin(), h.end(), cmp);
    
    //get x and y values;
    double x = randPerturb(0, h[1].path);
    double y = randPerturb(0, h[0].path);

    //get distances for new u and v 
    double newUdist = max(x,y);
    double newVdist = min(x,y);

    //get uv and cv for the HR calculation
    double uv = getPathLength(u, v, newTree);
    double cv = getPathLength(c, v, newTree);
        
    //change the weight of w-u and u-v accordingly
    put(edge_weight_t(), newTree.g, edge(v, w, newTree.g).first, newVdist);
    put(edge_weight_t(), newTree.g, edge(u, v, newTree.g).first, newUdist-newVdist);


    //two cases. Case1 : if max(x,y), or newUdist < h[0]. 
    if(newUdist < h[0].path){

      randNum = rand()%3;
      switch(randNum){
      case 0: // a is joined to v
	//remove unneccessary edges
	remove_edge(a, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	
	e1 = add_edge(a, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distA-newVdist < 0) cout << "something wrong0" << endl;
	if(distC-newUdist < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distA-newVdist);
	put(edge_weight_t(), newTree.g, e2, distC-newUdist);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-newUdist);
	break;
      case 1: // b is joined to v
	//remove unneccessary edges
	remove_edge(b, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	e1 = add_edge(b, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distB-newVdist < 0) cout << "something wrong0" << endl;
	if(distC-newUdist < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distB-newVdist);
	put(edge_weight_t(), newTree.g, e2, distC-newUdist);
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-newUdist);
	break;
      case 2: //topology does not change. c is joined with v
	//just update weight
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-newUdist);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-newUdist);
	put(edge_weight_t(), newTree.g, edge(c, v, newTree.g).first, distC-newVdist);
	break;
      default:
	cout << "should never reach here" << endl;
	break;
      }

      //Setting the Hastings ratio
      //IF original uv > cv then HR = 3.
      //Otherwise, it is 1.
      
      if(uv > cv){
	this->HastingsRatio = 3;
      }else{
	this->HastingsRatio = 1;
      }
      
    }else{
      //case 2: max(x,y) >h[0]
      
      //join a child with minimum height to v
      if(c == h[0].vert){
	//just update the weights
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-newUdist);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-newUdist);
	put(edge_weight_t(), newTree.g, edge(c, v, newTree.g).first, distC-newVdist);
      }else if(a == h[0].vert){
	//cout << "hi" << endl;
	//remove unneccessary edges
	remove_edge(a, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	e1 = add_edge(a, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distA-newVdist < 0) cout << "something wrong0" << endl;
	if(distC-newUdist < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distA-newVdist);
	put(edge_weight_t(), newTree.g, e2, distC-newUdist);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-newUdist);
      }else{
	//cout << "hi" << endl;
	//remove unneccessary edges
	remove_edge(b, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	e1 = add_edge(b, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distB-newVdist < 0) cout << "something wrong0" << endl;
	if(distC-newUdist < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distB-newVdist);
	put(edge_weight_t(), newTree.g, e2, distC-newUdist);
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-newUdist);
      }  

      //Setting the Hastings ratio
      //IF original uv > cv then HR = 3.
      //Otherwise, it is 1.
      
      if(uv < cv){
	this->HastingsRatio = 1.0/3.0;
      }else{
	this->HastingsRatio = 1;
      }
      
    }
  }else{

    for (tie(ai, ai_end) = adjacent_vertices(v, newTree.g); ai != ai_end; ++ai){
      if(*ai != u){
	c = *ai;
      }
    }   
      
    //Calculate the h array and sort it
    vector<vertexAndPathLength> h;
    double distA = getPathLength(a, v, newTree);
    vertexAndPathLength vdA;
    vdA.vert = a;
    vdA.path = distA; 
    double distB = getPathLength(b, v, newTree);
    vertexAndPathLength vdB;
    vdB.vert = b;
    vdB.path = distB; 
    double distC = getPathLength(c, v, newTree);
    vertexAndPathLength vdC;
    vdC.vert = c;
    vdC.path = distC; 
    h.push_back(vdA);
    h.push_back(vdB);
    h.push_back(vdC);

    //sorting the vector. make sure you include <algorithm>
    sort(h.begin(), h.end(), cmp);
      
    vector<vertexAndPathLength> hStar;
    //hi*=h1 x e dot gamma^(U-0.5)

    ///we wanna talk about this later 
    double U = randPerturb(0.8,1.2);
    //double lamda = 2.0*log(2);
    //cout << M_E*pow(lamda, U-0.5) << endl;
      
    vertexAndPathLength temp1;
    temp1.vert = h[0].vert;
    temp1.path = h[0].path*U;//*M_E*pow(lamda, U-0.5);
    hStar.push_back(temp1);

    vertexAndPathLength temp2;
    temp2.vert = h[1].vert;
    temp2.path = h[1].path+hStar[0].path-h[0].path;
    hStar.push_back(temp2);
      
    vertexAndPathLength temp3;
    temp3.vert = h[2].vert;
    temp3.path = h[2].path+hStar[0].path-h[0].path;
    hStar.push_back(temp3);

    //get uv and cv for the HR calculation
    double uv = getPathLength(u, v, newTree);
    double cv = getPathLength(c, v, newTree);
      
    //place u* at a height x above v*
    double x = randPerturb(0, hStar[1].path);
    put(edge_weight_t(), newTree.g, edge(u, v, newTree.g).first, x);

    for(vector<vertexAndPathLength>::iterator cur = hStar.begin(); cur != hStar.end(); cur++){
      if(cur->vert == a){
	distA = cur->path;
      }else if(cur->vert == b){
	distB = cur->path;
      }else{
	distC = cur->path;
      }
    }

    if(x < hStar[0].path){
      //any topopogical swapping can happen;
      randNum = rand()%3;
      switch(randNum){
      case 0: // a is joined to v*
	//remove unneccessary edges
	remove_edge(a, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	e1 = add_edge(a, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distC-x < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distA);
	put(edge_weight_t(), newTree.g, e2, distC-x);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-x);
	break;
      case 1: // b is joined to v
	//remove unneccessary edges
	remove_edge(b, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	e1 = add_edge(b, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distC-x < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distB);
	put(edge_weight_t(), newTree.g, e2, distC-x);
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-x);
	break;
      case 2: //topology does not change. c is joined with v
	//just update weight
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-x);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-x);
	put(edge_weight_t(), newTree.g, edge(c, v, newTree.g).first, distC);
	break;
      default:
	cout << "should never reach here" << endl;
	break;
      }

      if(uv > cv){
	this->HastingsRatio = 3;
      }else{
	this->HastingsRatio = 1;
      }
      
    }else{
      if(c == hStar[0].vert){
	//just update the weights
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-x);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-x);
	put(edge_weight_t(), newTree.g, edge(c, v, newTree.g).first, distC);
      }else if(a == hStar[0].vert){
	//cout << "hi" << endl;
	//remove unneccessary edges
	remove_edge(a, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	e1 = add_edge(a, v, newTree.g).first;
	e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distC-x < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distA);
	put(edge_weight_t(), newTree.g, e2, distC-x);
	put(edge_weight_t(), newTree.g, edge(b, u, newTree.g).first, distB-x);
      }else{
	//cout << "hi" << endl;
	//remove unneccessary edges
	remove_edge(b, u, newTree.g);
	remove_edge(c, v, newTree.g);

	//add new edges
	edge_descriptor e1 = add_edge(b, v, newTree.g).first;
	edge_descriptor e2 = add_edge(c, u, newTree.g).first;

	//update edge weights
	if(distC-x < 0) cout << "something wrong1" << endl;
	put(edge_weight_t(), newTree.g, e1, distB);
	put(edge_weight_t(), newTree.g, e2, distC-x);
	put(edge_weight_t(), newTree.g, edge(a, u, newTree.g).first, distA-x);
      }
      if(uv < cv){
	this->HastingsRatio = 1.0/3.0;
      }else{
	this->HastingsRatio = 1;
      }
    }
    
    this->HastingsRatio = this->HastingsRatio*(hStar[0].path/h[0].path);

  }
  
  return newTree;
}

double MonteCarlo::getPathLength(vertex_descriptor start, vertex_descriptor goal, const Tree &tree){
  vector<vertex_descriptor> path = getPath(start, goal, tree);
    
  double sum = 0;
  for(int i = 0; i < path.size()-1; i++){
    edge_descriptor e = edge(path[i], path[i+1], tree.g).first;
    double weight = get(edge_weight_t(), tree.g, e);
    sum += weight;
  }
  return sum;
}

double MonteCarlo::metropolisHastings(const Tree &curTree, const Tree &newTree, int sitenum){
    double logLCur = calculateLikelihood(curTree, sitenum);
    double logLNew = calculateLikelihood(newTree, sitenum);
    double LRatio = pow(M_E, logLNew-logLCur)*this->HastingsRatio;
    if (logLNew-logLCur > 10)
	return 1.0;
    else if (logLNew-logLCur < -10)
	return 0.0;
    else
	return min(1.0, LRatio);
}

double MonteCarlo::calculateLikelihood(const Tree& inputTree, int siteNum)
{
    //First, we need to root the topology by inserting a node ON an edge.
    //Get first edge and treat it as root

    Tree topologyIn(inputTree);
    
    double logLikelihood = 0;
    //This can be parallelized.
    #pragma omp parallel for reduction(+:logLikelihood)
    for(int i = 0; i < siteNum; ++i)
    {
	double perSiteL = 0;
	for(int j = 0; j < this->alphabetSize; j++)
	{
	    map<vertex_descriptor, bool> visited;
	    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
	    for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
	    {
		visited[*vi] = false;
	    }
	    perSiteL += this->freq[this->nuc[j]] * calculateLikelihoodRecursively(topologyIn,
										  this->root,
										  this->nuc[j],
										  visited,
										  i);
	}
	logLikelihood += log(perSiteL);
    }

    return logLikelihood;
}

double MonteCarlo::calculateLikelihoodRecursively(const Tree& t,
						  vertex_descriptor rootNode,
						  char state,
						  map<vertex_descriptor, bool> &visited,
						  int seqIndex)
{
    visited[rootNode] = true;
    
    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    vector <vertex_descriptor> unvisitedNeighbors;

    // get the vector of unvisited neighbors, which will either
    // have length of 2 or 0.
    for (tie(ai, ai_end) = adjacent_vertices(rootNode, t.g); ai != ai_end; ++ai)
    {
	if (! visited[*ai])
	    unvisitedNeighbors.push_back(*ai);
    }
        
    //Leaf case
    if(unvisitedNeighbors.size() == 0)
    {
	// Get sequence
	string seq = get(Vertex_t(), t.g, rootNode)->alignedSequence;
	char nucleotide = seq[seqIndex];
	if(nucleotide == state)
	{
	    //map clean up
	    visited[rootNode] = false;
	    return 1;
	}
	else
	{
	    //map clean up
	    visited[rootNode] = false;
	    return 0;
	}
    }
    else
    {
	vertex_descriptor childOne = unvisitedNeighbors[0];
	vertex_descriptor childTwo = unvisitedNeighbors[1];
	
	double childP1 = 0;
	double childP2 = 0;
	
	// Check all possible nucleotides that can appear in child nodes and how that will affect the likelihood; 
	for(int j = 0; j<this->alphabetSize; j++)
	{
	    childP1 += 
		mutationProbability(state, this->nuc[j], 
				    get(edge_weight_t(), t.g, edge(rootNode,childOne,t.g).first))
		* calculateLikelihoodRecursively(t, childOne, this->nuc[j], visited, seqIndex);
	    childP2 += 
		mutationProbability(state, this->nuc[j],
				    get(edge_weight_t(), t.g, edge(rootNode,childTwo,t.g).first))
		* calculateLikelihoodRecursively(t, childTwo, this->nuc[j], visited, seqIndex);
	}
	//reset the map
	visited[rootNode] = false;
	return childP1*childP2;
    }
}

double MonteCarlo::mutationProbability(char fromState, char toState, float branchLength)
{
    /* Calculates the probability of transitioning from a state to another in a single site, currently only accounts for freq */  
    double prob;
    if (fromState == toState)
    {
	prob = exp(-(uRate*branchLength)) + (1 - exp(-(uRate*branchLength)))*this->freq[toState];
    }
    else
    {
	prob = (1 - exp(-(uRate*branchLength)))*this->freq[toState];
    }
    return prob;
}
    


bool cmp(const vertexAndPathLength &a, const vertexAndPathLength &b)
{
  return a.path< b.path;
}
