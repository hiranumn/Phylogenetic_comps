/*
 * maximumLikelihoodProg.cpp
 * Given an aligned sequences and options, it creates a tree topology that explains the evolutionary history of sequences
 * using maximumLiokelihood method and progressive tree search approach.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#include "maximumLikelihoodProg.h"

Tree MaximumLikelihoodProgressive::solve(const SequenceData& sd, Options opt)
{
  //measure the length of the sequences
  int siteNum = sd.speciesMap.at(sd.getFirstSpeciesName())->alignedSequence.size();

  //get parameters from the option file
  this->uRate = opt.uRate;
  this->convergence = opt.convergence;
  this->nuc = opt.nuc;
  this->alphabetSize = nuc.size();
  this->freq = opt.freq;
    
    
  // incorrect usage check
  if(sd.speciesMap.size() < 2)
    {
      cout << "Please specify more sequences" << endl;
      Tree empty;
      return empty;
    }
    

  // Place all species into a vector for ease of iteration later
  map<string, Species*>::const_iterator iter;
  vector<Species*> species;
  for(iter = sd.speciesMap.begin(); iter != sd.speciesMap.end(); iter++)
    {
      species.push_back(iter->second);
    }
    
  Tree curTree;

  // Init base tree
  vertex_descriptor v1, v2;
  v1 = add_vertex(curTree.g);
  v2 = add_vertex(curTree.g);
  put(Vertex_t(), curTree.g, v1, species[0]);
  put(Vertex_t(), curTree.g, v2, species[1]);
  edge_descriptor firstEdge = add_edge(v1, v2, curTree.g).first;
  put(edge_weight_t(), curTree.g, firstEdge, 0.5);

  //build up the tree
  for(int i = 2; i!=species.size(); i++)
    {
      // Generate possible topologies of the size + 1, find the best choice
      vector<Tree> toplist = getPossibleTopology(curTree, species[i], i-2); 
      adjustTopology(toplist[0], siteNum);

      // start with the first possible tree topology
      int bestNeighborIndex = 0;
      double bestNeighborLikelihood = calculateLikelihood(toplist[0], siteNum);

      // iterate through the rest to see if they are any better. 
      for(int j=1; j!=toplist.size(); j++)
	{ 
	  adjustTopology(toplist[j], siteNum);
	  double compLikelihood = calculateLikelihood(toplist[j], siteNum);
	  if(compLikelihood > bestNeighborLikelihood)
	    {
	      bestNeighborIndex = j;
	      bestNeighborLikelihood = compLikelihood;
	    } 
	}
      curTree = toplist[bestNeighborIndex];
    }
  return curTree; // Returns the constructed tree with all species inserted and branches adjusted
}


vector<Tree> MaximumLikelihoodProgressive::getPossibleTopology(const Tree& baseTree,
							       Species* s,
							       int dummyNum)
{
  /* Generate a vector containing all possible new trees of size n+1 from the original structure */
  vector<Tree> ret;
    
  graph_traits<Tree_graph>::edge_iterator ei, ei_end;
  // For each branch in the given tree.
  for(tie(ei, ei_end) = edges(baseTree.g); ei != ei_end; ++ei)
    {
      // Duplicate the tree
      Tree t(baseTree);
	
      // Get the edge information from the original tree by vertex descriptor
      vertex_descriptor v1 = source(*ei, baseTree.g);
      vertex_descriptor v2 = target(*ei, baseTree.g);
	
      // Put a dummy node on a new edge
      edge_descriptor e = edge(v1,v2,t.g).first;
	
      float curEdgeLength = get(edge_weight_t(), t.g, e);
      float ew1 = curEdgeLength/2.0;
      float ew2 = 1 - ew1;
	
      remove_edge(v1, v2, t.g);
      vertex_descriptor dummyNodeId = add_vertex(t.g);

      string dummyName = "DUMMY_"; dummyName += to_string(dummyNum);
      Species* dummy = new Species(dummyName);
      dummy -> dummy = true;
      put(Vertex_t(), t.g, dummyNodeId, dummy);
      dummies.push_back(dummy);

      edge_descriptor newEdge1 = add_edge(dummyNodeId, v1, t.g).first;
      edge_descriptor newEdge2 = add_edge(dummyNodeId, v2, t.g).first;
	
      put(edge_weight_t(), t.g, newEdge1, ew1);
      put(edge_weight_t(), t.g, newEdge2, ew2);

      //Now attach current species to dummy node with edge weight 5
      vertex_descriptor nodeId = add_vertex(t.g);
      put(Vertex_t(), t.g, nodeId, s);
	
      edge_descriptor newEdgeReal = add_edge(nodeId, dummyNodeId, t.g).first;
      put(edge_weight_t(), t.g, newEdgeReal, .5);

      // Push t into the vector
      ret.push_back(t);
    }
  return ret;
}


void MaximumLikelihoodProgressive::adjustTopology(Tree &t, int siteNum)
{
  // Tweaks the topology by adjusting branch lengths one at a time to find improvements
  graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    
  for(tie(ei, ei_end) = edges(t.g); ei != ei_end; ++ei)
    {
      float weight = getOptimalBranchLength(*ei, t, siteNum);
      put(edge_weight_t(), t.g, *ei, weight);
    }
}


float MaximumLikelihoodProgressive::getOptimalBranchLength(edge_descriptor e, 
							   const Tree &t,
							   int siteNum)
{
  /* Finds the optimal length for a single branch */
  double weight = get(edge_weight_t(), t.g, e);
    
  double curP = 1 - exp(-weight);
  double oldP = 9999999; // Dramatic increase in likelihood of entering the improvement loop
    
  // Get the source and the target of an edge
  vertex_descriptor v1 = source(e, t.g);
  vertex_descriptor v2 = target(e, t.g);
    
  // Single time calculation of the A & B values which relate to subtree likelihood
  vector<double> vecA(siteNum);
  vector<double> vecB(siteNum);
  #pragma omp parallel for shared(vecA, vecB)
  for(int i = 0; i < siteNum; ++i)
  {
      vecA[i] = A(v1, v2, i, t);
      vecB[i] = B(v1, v2, i, t);
  }
    
  // While curP does not converge according to our set criteria
  while(abs(curP-oldP) > this->convergence)
  {
      oldP = curP;
      curP = 0;
      #pragma omp parallel for reduction(+:curP)
      for(int i = 0; i < siteNum; ++i)
	  curP += (vecB[i]*oldP)/(vecA[i]*(1-oldP)+vecB[i]*oldP);
      curP = (1.0/siteNum)*curP;
  }
  
  float newBranchLength = -log(1-curP);
  return newBranchLength;
}

double MaximumLikelihoodProgressive::A(vertex_descriptor v1, vertex_descriptor v2, int index,
				       const Tree& t)
{
  Tree copy(t);
  //temporarily remove the edge between v1 and v2, 
  //so that likelihood calculation does not recurse in a wrong direction 
  remove_edge(v1, v2, copy.g);

  double sum = 0;
  double subprod;
  Species* s1 = get(Vertex_t(), t.g, v1); // Vertex 1 info
  Species* s2 = get(Vertex_t(), t.g, v2); // Vertex 2 info
    
  for(int i=0; i<this->alphabetSize; i++)
    {
      double curTerm = 1;
      curTerm *= this->freq[this->nuc[i]];
      map<vertex_descriptor, bool> visited;
      graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
      //Make sure to set visited to false for all nodes. 
      for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	  visited[*vi] = false;
	}
      //calculate the likelihood of a subtree rooted at v1. 
      curTerm *= calculateLikelihoodRecursively(copy, v1, this->nuc[i],
						visited, index);

      //Again, make sure to set visited to false for all nodes.	
      for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	  visited[*vi] = false;
	}
      //calculate the likelihood of a subtree rooted at v2. 
      curTerm *= calculateLikelihoodRecursively(copy, v2, this->nuc[i],
						visited, index);
      sum += curTerm;
    }
  return sum;
}


double MaximumLikelihoodProgressive::B(vertex_descriptor v1, vertex_descriptor v2, int index,
				       const Tree& t)
{
  //copy the tree...
  Tree copy(t);
  //temporarily remove the edge between v1 and v2, 
  //so that likelihood calculation does not recurse in a wrong direction
  remove_edge(v1, v2, copy.g);

  /* Calculates the likelihood of a subtree given mismatched characters in a state */
  double sum1 = 0;
  double sum2 = 0;

  for(int i = 0; i < this->nuc.size(); ++i)
    {
      map<vertex_descriptor, bool> visited;
      graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
      //reset visited
      for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	  visited[*vi] = false;
	}
      //calculate the likelihood of a subtree rooted at v1.
      sum1 += freq[this->nuc[i]]*calculateLikelihoodRecursively(copy, v1, this->nuc[i],
								visited, index);
      //reset visited
      for(tie(vi, vi_end) = vertices(copy.g); vi != vi_end; ++vi)
	{
	  visited[*vi] = false;
	}
      //calculate the likelihood of a subtree rooted at v2.
      sum2 += freq[this->nuc[i]]*calculateLikelihoodRecursively(copy, v2, this->nuc[i],
								visited, index);
    }
  double ret = sum1*sum2;
  return ret;
}



double MaximumLikelihoodProgressive::mutationProbability(char fromState, 
							 char toState, 
							 float branchLength)
{
  //Calculates the probability of transitioning from a state to another
  //in a single site, currently only accounts for freq
  double prob;
  if (fromState == toState)
    {
      prob = exp(-(uRate*branchLength)) + (1 - exp(-(uRate*branchLength)))*this->freq[toState]; //e^-ut*1 +(1-e^-ut)*p_j
    }
  else 
    {
      prob = (1 - exp(-(uRate*branchLength)))*this->freq[toState]; //e^-ut*0 +(1-e^-ut)*p_j
    }
  return prob;
}

double MaximumLikelihoodProgressive::calculateLikelihood(const Tree& inputTree, int siteNum)
{
  //First, we need to root the topology by inserting a node ON an edge.
  //Get first edge and treat it as root

  Tree topologyIn(inputTree);

  graph_traits<Tree_graph>::edge_iterator ei, ei_end;
  tie(ei, ei_end) = edges(topologyIn.g);
  vertex_descriptor neighborOne = source(*ei, topologyIn.g);
  vertex_descriptor neighborTwo = target(*ei, topologyIn.g);
    
  float oldEdgeWeight = get(edge_weight_t(),
			    topologyIn.g,
			    edge(neighborOne, neighborTwo, topologyIn.g).first);

  float ew1 = oldEdgeWeight/2;
  float ew2 = oldEdgeWeight - ew1;

  remove_edge(neighborOne, neighborTwo, topologyIn.g);
  //Put the "dummy root" into the graph.
  Species* dummyRoot = new Species("DUMMY_ROOT");
  vertex_descriptor nodeId = add_vertex(topologyIn.g);
  put(Vertex_t(), topologyIn.g, nodeId, dummyRoot);
    
  //Add edges from dummy root to the neighbors.
  edge_descriptor firstEdge = add_edge(nodeId, neighborOne, topologyIn.g).first;
  edge_descriptor secondEdge = add_edge(nodeId, neighborTwo, topologyIn.g).first;

  put(edge_weight_t(), topologyIn.g, firstEdge, ew1);
  put(edge_weight_t(), topologyIn.g, secondEdge, ew2);
    
  double logLikelihood = 0;

  //Calculate likelihood of the tree at every DNA site and multiply them together.
  //We are using log likelihood so we actually log them and add them.
  #pragma omp parallel for reduction(+:logLikelihood)
  for(int i = 0; i < siteNum; ++i)
  {
      double perSiteL = 0;
      //The state (nucleotide) at the root is UNKOWN so, do it over every possible DNA assingment.
      //Add the resulting likelihood values together.
      for(int j = 0; j < this->alphabetSize; j++)
      {
	  map<vertex_descriptor, bool> visited;
	  graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
	  for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
	  {
	      visited[*vi] = false;
	  }
	  perSiteL += this->freq[this->nuc[j]] * calculateLikelihoodRecursively(topologyIn,
										nodeId,
										this->nuc[j],
										visited,
										i);
      }
      logLikelihood += log(perSiteL);
  }
  
  delete dummyRoot;
    
  return logLikelihood;
}

double MaximumLikelihoodProgressive::calculateLikelihoodRecursively(const Tree& t,
								    vertex_descriptor rootNode,
								    char state,
								    map<vertex_descriptor, bool> &visited,
								    int seqIndex)
{
  //Mark itself visited   
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
      // Get the nucleotide and see if matches or not with the assigned state
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

