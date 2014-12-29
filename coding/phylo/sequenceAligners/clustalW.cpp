/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* clustalW.cpp
* Contains an implementation of the clustal W multiple sequence alignment algorithm
*/

#include "clustalW.h"

/*
	The main entry point for the clustal W algorithm
*/
void ClustalW::align(SequenceData &data, const Options &options)
{
	// Construct a distance matrix by quickly aligning all sequences to one another and getting the distance
    map<string, map<string, double> > distances;
    vector<string> names = data.getAllSpeciesNames();
    for(int i = 0; i < names.size(); ++i)
    {
		for(int j = i + 1; j < names.size(); ++j)
		{
		    pair<string, string> alignment = pairwiseAlign(data.speciesMap[names[i]]->rawSequence,
								   data.speciesMap[names[j]]->rawSequence,
								   options);
		    int numNonGapPositions = 0;
		    int numMatches = 0;
		    for(int k = 0; k < alignment.first.length(); ++k)
		    {
				if(alignment.first.at(k) != '-' && alignment.second.at(k) != '-') numNonGapPositions++;
				if(alignment.first.at(k) == alignment.second.at(k)) numMatches++;
		    }
		    double distance = 1 - 1.0*numMatches/numNonGapPositions;

		    distances[names[i]][names[j]] = distance;
		    distances[names[j]][names[i]] = distance;
		}
    }

    // Construct the tree from the distance matrix
    DistanceMatrix distMatrix(distances);
    NeighborJoining neighborJoiner;
    Tree guideTree = neighborJoiner.solve(data, distMatrix, options);

    // Traverse this tree to progressively align all of the sequences
    alignWithTree(guideTree, options);
    
    data.aligned = true;
    data.alignedLength = data.speciesMap[data.getFirstSpeciesName()]->alignedSequence.size();
}

/*
	Given two strings, find their optimal alignment
*/
pair<string, string> ClustalW::pairwiseAlign(string s1, 
					     string s2, 
					     const Options &options)
{
    vector<vector<double> > alignMatrix(s1.length()+1,vector<double>(s2.length()+1,0));
    vector<vector<int> > pathMatrix(s1.length()+1,vector<int>(s2.length()+1,0));
    for(int i = 1; i < s1.length() + 1; ++i)
    {
		for(int j = 1; j < s2.length() + 1; ++j)
		{
		    //choiceVal represents which alignment choice we made
		    //this corresponds to where in the distance matrix we
		    //computed this value from, 0 represents inserting
		    //a space into string 1, 1 represents "accepting" the
		    //character match, and 2 represents inserting a space
		    //into string 2.
		    int choiceVal = -1;
		    double choiceOne = alignMatrix[i][j-1] + options.gapReward;
		    double choiceTwo = alignMatrix[i-1][j-1];
		    if(s1.at(i-1) == s2.at(j-1)) //if we have a character match
			choiceTwo += options.matchReward;
		    double choiceThree = alignMatrix[i-1][j] + options.gapReward;
		    
		    vector<double> choiceVec;
		    choiceVec.push_back(choiceOne);
		    choiceVec.push_back(choiceTwo);
		    choiceVec.push_back(choiceThree);
		    vector<double>::iterator max = max_element(choiceVec.begin(), choiceVec.end());
		    int maxIndex = distance(choiceVec.begin(), max);
		    switch(maxIndex)
		    {
			    case 0:
					alignMatrix[i][j] = choiceOne;
					pathMatrix[i-1][j-1] = 0;
					break;
			    case 1:
					alignMatrix[i][j] = choiceTwo;
					pathMatrix[i-1][j-1] = 1;
					break;
			    case 2:
					alignMatrix[i][j] = choiceThree;
					pathMatrix[i-1][j-1] = 2;
					break;
		    }
		}
    }

    string alignedOne = "";
    string alignedTwo = "";
    int curI = s1.length()-1;
    int curJ = s2.length()-1;
    while(true)
    {
		if(curI < 0)
		{
		    alignedOne += "-"; alignedTwo += s2.at(curJ);
		    curJ --;
		}
		else if(curJ < 0)
		{
		    alignedTwo += "-"; alignedOne += s1.at(curI);
		    curI --;
		}
		else
		{
		    switch(pathMatrix[curI][curJ])
		    {
			    case 0:
					alignedOne += "-";
					alignedTwo += s2.at(curJ);
					curJ--;
					break;
			    case 1:
					alignedOne += s1.at(curI);
					alignedTwo += s2.at(curJ);
					curI--;
					curJ--;
					break;
			    case 2:
					alignedOne += s1.at(curI);
					alignedTwo += "-";
					curI--;
				break;
		    }
		}
		if(curI < 0 && curJ < 0) break;
    }
    alignedOne = string(alignedOne.rbegin(), alignedOne.rend());
    alignedTwo = string(alignedTwo.rbegin(), alignedTwo.rend());
    
    pair<string, string> returnPair(alignedOne, alignedTwo);

    return returnPair;
}

/*
	A scoring value to choose whether to insert or pass
*/
double ClustalW::calcSP(vector <Species*> profile1, 
			vector <Species*> profile2, 
			int profile1Col, int profile2Col,
			const Options &options)
{    
    vector<double> profile1Freqs;
    vector<double> profile2Freqs;
    
    //map for tracking the number of occurrences of each character
    map<char, int> profile1CharCount;
    map<char, int> profile2CharCount;

    for (int i = 0; i < profile1.size(); ++i)
		profile1CharCount[profile1[i]->alignedSequence.at(profile1Col)]++;
    for (int i = 0; i < profile2.size(); ++i)
		profile2CharCount[profile2[i]->alignedSequence.at(profile2Col)]++;

    for (int i = 0; i < options.alphabet.size(); ++i)
    {
		profile1Freqs.push_back(1.0 * profile1CharCount[options.alphabet[i]] / profile1.size());
		profile2Freqs.push_back(1.0 * profile2CharCount[options.alphabet[i]] / profile2.size());
    }
    
    double runningSum = 0;
    for (int i = 0; i < options.alphabet.size(); ++i)
    {
		for (int j = 0; j < options.alphabet.size(); ++j)
		{
		    if (options.similarities.find(options.alphabet[i]) == options.similarities.end())
		    {
				continue;
		    }
		    else
		    {
				if(options.similarities.at(options.alphabet[i]).find(options.alphabet[j]) == options.similarities.at(options.alphabet[i]).end())
			    	continue;
		    }
		    double sim = options.similarities.at(options.alphabet[i]).at(options.alphabet[j]);
		    runningSum += profile1Freqs[i] * profile2Freqs[j] * sim;
		}
    }
    return runningSum;
}
    
/*
    Create the alignment profile from the given two profiles
*/
void ClustalW::alignProfiles(vector<Species*> profile1,
			     vector<Species*> profile2,
			     const Options &options)
{
    //length of strings in profile1.
    int n1 = profile1[0]->alignedSequence.length();
    //length of strings in profile2.
    int n2 = profile2[0]->alignedSequence.length();
    
    vector<vector<double> > alignMatrix(n1+1,vector<double>(n2+1,0));
    vector<vector<int> > pathMatrix(n1,vector<int>(n2,0));

    for(int i = 1; i < n1 + 1; ++i)
    {
		for(int j = 1; j < n2 + 1; ++j)
		{
		    //choiceVal represents which alignment choice we made
		    //this corresponds to where in the distance matrix we
		    //computed this value from
		    //0 is coming from the right
		    //1 is coming diagonally down/right
		    //2 is coming down
		    int choiceVal = -1;
		    double choiceOne = alignMatrix[i][j-1] + options.profileGapReward;
		    double choiceTwo = alignMatrix[i-1][j-1] + calcSP(profile1, profile2, i-1, j-1, options);
		    double choiceThree = alignMatrix[i-1][j] + options.profileGapReward;
		    
		    vector<double> choiceVec;
		    choiceVec.push_back(choiceOne);
		    choiceVec.push_back(choiceTwo);
		    choiceVec.push_back(choiceThree);
		    vector<double>::iterator max = max_element(choiceVec.begin(), choiceVec.end());
		    int maxIndex = distance(choiceVec.begin(), max);
		    switch(maxIndex)
		    {
			    case 0:
					alignMatrix[i][j] = choiceOne;
					pathMatrix[i-1][j-1] = 0;
					break;
			    case 1:
					alignMatrix[i][j] = choiceTwo;
					pathMatrix[i-1][j-1] = 1;
					break;
			    case 2:
					alignMatrix[i][j] = choiceThree;
					pathMatrix[i-1][j-1] = 2;
					break;
		    }
		}
    }

    int curI = n1-1;
    int curJ = n2-1;
    while(true)
    {
		if(curI < 0)
		{
		    for(int k = 0; k < profile1.size(); ++k)
				profile1[k]->alignedSequence = '-' + profile1[k]->alignedSequence;
		    curJ --;
		}
		else if(curJ < 0)
		{
		    for(int k = 0; k < profile2.size(); ++k)
				profile2[k]->alignedSequence = '-' + profile2[k]->alignedSequence;
		    curI --;
		}
		else
		{
		    switch(pathMatrix[curI][curJ])
		    {
			    case 0:
					for(int k = 0; k < profile1.size(); ++k)
					    profile1[k]->alignedSequence.insert(profile1[k]->alignedSequence.begin()+curI+1,'-');
					curJ--;
					break;
			    case 1:
					curI--;
					curJ--;
					break;
			    case 2:
					for(int k = 0; k < profile2.size(); ++k)
					    profile2[k]->alignedSequence.insert(profile2[k]->alignedSequence.begin()+curJ+1,'-');
					curI--;
					break;
			}
		}
		if(curI < 0 && curJ < 0) break;
    }
}

/*
    Insert a root, set up for alignment, call alignwithtreerecursive
*/
void ClustalW::alignWithTree(const Tree &constTopologyIn, const Options &options)
{
    //First, we need to root the topology by inserting a node ON an edge.
    //Get first edge and treat it as root
    Tree topologyIn(constTopologyIn);

    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    tie(ei, ei_end) = edges(topologyIn.g);
    vertex_descriptor neighborOne = source(*ei, topologyIn.g);
    vertex_descriptor neighborTwo = target(*ei, topologyIn.g);
    
    remove_edge(neighborOne, neighborTwo, topologyIn.g);
    //Put the "dummy root" into the graph.
    Species* dummyRoot = new Species("DUMMY_ROOT");
    vertex_descriptor nodeId = add_vertex(topologyIn.g);
    put(Vertex_t(), topologyIn.g, nodeId, dummyRoot);
    
    //Add edges from dummy root to the neighbors.
    edge_descriptor firstEdge = add_edge(nodeId, neighborOne, topologyIn.g).first;
    edge_descriptor secondEdge = add_edge(nodeId, neighborTwo, topologyIn.g).first;

    put(edge_weight_t(), topologyIn.g, firstEdge, 1.0);
    put(edge_weight_t(), topologyIn.g, secondEdge, 1.0);
    
    int totalLength = 0;
    //This can be parallelized.
    map<vertex_descriptor, bool> visited;
    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
    for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
    {
		visited[*vi] = false;
    }
    map<vertex_descriptor, vector<Species*> > nodeProfiles;
    alignWithTreeRecursive(topologyIn, nodeId, visited, nodeProfiles, options);

    delete dummyRoot;
        
}

/*
    Gather information on the tree's structure, call alignprofiles to make profiles and push that up
*/
void ClustalW::alignWithTreeRecursive(const Tree &topology,
				      vertex_descriptor rootNode,
				      map<vertex_descriptor, bool> &visited,
				      map<vertex_descriptor, vector<Species*> > &nodeProfiles,
				      const Options &options)
{
    visited[rootNode] = true;

    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    vector <vertex_descriptor> unvisited_neighbors;

    // get the vector of unvisited neighbors, which will either
    // have length of 2 or 0.
    for (tie(ai, ai_end) = adjacent_vertices(rootNode, topology.g); ai != ai_end; ++ai)
    {
		if (! visited[*ai])
		    unvisited_neighbors.push_back(*ai);
    }
    
    // base case when looking at a leaf node:
    // add the sequence character to charSet, 
    // and return length 0.
    if (unvisited_neighbors.size() == 0)
    {
		vector<Species*> curProfile;
		curProfile.push_back(get(Vertex_t(), topology.g, rootNode));
		nodeProfiles[rootNode] = curProfile;
		return;
    }
    
    vertex_descriptor childOne = unvisited_neighbors[0];
    vertex_descriptor childTwo = unvisited_neighbors[1];
    alignWithTreeRecursive(topology, childOne, visited, nodeProfiles, options);
    alignWithTreeRecursive(topology, childTwo, visited, nodeProfiles, options);
        
    vector<Species*> childOneProfile = nodeProfiles[childOne];
    vector<Species*> childTwoProfile = nodeProfiles[childTwo];

    alignProfiles(childOneProfile, childTwoProfile, options);

    vector<Species*> curProfile;
    for(int i = 0; i < childOneProfile.size(); ++i)
		curProfile.push_back(childOneProfile[i]);

    for(int i = 0; i < childTwoProfile.size(); ++i)
		curProfile.push_back(childTwoProfile[i]);

    nodeProfiles[rootNode] = curProfile;
    
}
