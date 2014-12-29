//clustalW.cpp
//contains implementations of the ClustalWParallel methods
#include "clustalWParallel.h"
    
    
void ClustalWParallel::align(SequenceData &data, const Options &options)
{
    map<string, map<string, double> > distances;
    vector<string> names = data.getAllSpeciesNames();
    int numThreads = names.size() * (names.size()-1) / 2;
    thread_group threads;
    vector<double> resultsArray(numThreads);
    
    int threadIndex1 = 0;
    int threadIndex2 = 1;
    for(int i = 0; i < numThreads; ++i)
    {
	thread* curThread = new thread(&ClustalWParallel::parallelDistanceCalc, this,
				       threadIndex1, threadIndex2, boost::ref(resultsArray[i]),
				       boost::ref(data), boost::ref(names), boost::ref(options));
    
	threads.add_thread(curThread);
	threadIndex2++;
	if(threadIndex2 == names.size())
	{
	    threadIndex1 ++;
	    threadIndex2 = threadIndex1 + 1;
	}
    }
    threads.join_all();
    
    //put paralell results in distance map correctly.
    int curIndex1 = 0;
    int curIndex2 = 1;
    for(int i = 0; i < numThreads; ++i)
    {
	distances[names[curIndex1]][names[curIndex2]] = resultsArray[i];
	distances[names[curIndex2]][names[curIndex1]] = resultsArray[i];
	curIndex2++;
	if(curIndex2 == names.size())
	{
	    curIndex1 ++;
	    curIndex2 = curIndex1 + 1;
	}
    }


    DistanceMatrix distMatrix(distances);
    NeighborJoining neighborJoiner;
    Tree guideTree = neighborJoiner.solve(data, distMatrix, options);

    //Now, we're going to traverse this tree to progressively align all of the sequences
    alignWithTree(guideTree, options);
    
    data.aligned = true;
    data.alignedLength = data.speciesMap[data.getFirstSpeciesName()]->alignedSequence.size();
}

void ClustalWParallel::parallelDistanceCalc(int index1,
					    int index2,
					    double &result,
					    const SequenceData &data,
					    const vector<string> &names,
					    const Options &options)
{
    string seq1 = data.speciesMap.at(names[index1])->rawSequence;
    string seq2 = data.speciesMap.at(names[index2])->rawSequence;
    
    pair<string, string> alignment = pairwiseAlign(seq1,
						   seq2,
						   options);
    
    int numNonGapPositions = 0;
    int numMatches = 0;
    for(int k = 0; k < alignment.first.length(); ++k)
    {
	if(alignment.first.at(k) != '-' && alignment.second.at(k) != '-') numNonGapPositions++;
	if(alignment.first.at(k) == alignment.second.at(k)) numMatches++;
    }
    double distance = 1 - 1.0*numMatches/numNonGapPositions;
    
    result = distance;
}



pair<string, string> ClustalWParallel::pairwiseAlign(string s1, 
						     string s2, 
						     const Options &options)
{
    vector<vector<double> > alignMatrix(s1.length()+1,vector<double>(s2.length()+1,0));
    vector<vector<int> > pathMatrix(s1.length(),vector<int>(s2.length(),0));

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

double ClustalWParallel::calcSP(vector <Species*> profile1, 
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
    
void ClustalWParallel::alignProfiles(vector<Species*> profile1,
				     vector<Species*> profile2,
				     const Options &options)
{
    //variable naming in accordance with book
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



void ClustalWParallel::alignWithTree(const Tree &constTopologyIn, const Options &options)
{
    //First, we need to root the topology by inserting a node ON an edge.
    //Get first edge and treat it as root

    Tree topologyIn(constTopologyIn);

    edge_descriptor rootEdge = ParallelUtils::pickRootEdge(topologyIn);
    vertex_descriptor neighborOne = source(rootEdge, topologyIn.g);
    vertex_descriptor neighborTwo = target(rootEdge, topologyIn.g);
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
    map<vertex_descriptor, bool> visited;
    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
    for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
    {
	visited[*vi] = false;
    }
    //Now, we need to make the neighbor dictionary, which maps from
    //vertex_descriptor -> vector<vertex_descriptor> where the key is the parent
    //and the vector is a vector of children.
    map<vertex_descriptor, vector<vertex_descriptor> > childMap;
    vector<vertex_descriptor> verticesToMap;
    verticesToMap.push_back(nodeId);
    while(verticesToMap.size() != 0)
    {
	vertex_descriptor curVertex = verticesToMap.back();
	verticesToMap.pop_back();
	visited[curVertex] = true;

	graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
	vector <vertex_descriptor> unvisitedNeighbors;
	
	// get the vector of unvisited neighbors, which will either
	// have length of 2 or 0.
	for(tie(ai, ai_end) = adjacent_vertices(curVertex, topologyIn.g); ai != ai_end; ++ai)
	{
	    if(! visited[*ai])
	    {
		unvisitedNeighbors.push_back(*ai);
		verticesToMap.push_back(*ai);
	    }
	}
	childMap[curVertex] = unvisitedNeighbors;
    }


    vector<Species*> rootProfile;

    alignWithTreeRecursive(topologyIn, nodeId, rootProfile, childMap, options);

    delete dummyRoot;
        
}

void ClustalWParallel::alignWithTreeRecursive(const Tree &topology,
					      const vertex_descriptor& rootNode,
					      vector<Species*> &myProfile,
					      const map<vertex_descriptor, vector<vertex_descriptor> >& childMap,
					      const Options &options)
{

    vector<vertex_descriptor> children = childMap.at(rootNode);

    if(children.size() == 0)
    {
	myProfile.push_back(topology.speciesmap[rootNode]);
	return;
    }

    vertex_descriptor childOne = children[0];
    vertex_descriptor childTwo = children[1];
    
    vector<Species*> childOneProfile, childTwoProfile;

    
    thread threadOne = thread(&ClustalWParallel::alignWithTreeRecursive, this, topology, childOne, boost::ref(childOneProfile), childMap, options);
    thread threadTwo = thread(&ClustalWParallel::alignWithTreeRecursive, this, topology, childTwo, boost::ref(childTwoProfile), childMap, options);
    
    threadOne.join();
    threadTwo.join();

    alignProfiles(childOneProfile, childTwoProfile, options);

    vector<Species*> curProfile;
    for(int i = 0; i < childOneProfile.size(); ++i)
	myProfile.push_back(childOneProfile[i]);

    for(int i = 0; i < childTwoProfile.size(); ++i)
	myProfile.push_back(childTwoProfile[i]);
}
