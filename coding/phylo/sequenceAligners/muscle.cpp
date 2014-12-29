/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* muscle.cpp 
* contains the implementation of the MUSCLE multiple sequence alignment algorithm
*/

#include "muscle.h"

/*
    Utility structure for splitting the tree
*/
struct edgeStruct 
{
    vertex_descriptor source;
    vertex_descriptor target;
};

/*
    Main entry function to the muscle algorithm
*/
void Muscle::align(SequenceData &data, const Options &options)
{
    NeighborJoining neighborJoiner;
    map<string, map<string, double> > distances;
    
    // DRAFT
    // KMER DISTANCE MEASUREMENT
    distances = kmer(data, options);
    DistanceMatrix distMatrixKMER(distances);
    // GUIDE TREE CREATION
    Tree guideTree1 = neighborJoiner.solve(data, distMatrixKMER, options);
    // PROGRESSIVE ALIGNMENT
    alignWithTree(guideTree1, options);
    data.aligned = true;
    data.alignedLength = data.speciesMap[data.getFirstSpeciesName()]->alignedSequence.size();
    // EVALUATION
    double SP1 = SP(data);
    cout << "SP distance - KMER : " << SP1 << endl;


    // IMPROVEMENT
    // KIMURA DISTANCE MEASUREMENT
    distances = kimura(data, options);
    DistanceMatrix distMatrixKIMURA(distances);
    // GUIDE TREE CREATION
    Tree guideTree2 = neighborJoiner.solve(data, distMatrixKIMURA, options);
    // PROGRESSIVE ALIGNMENT
    alignWithTree(guideTree2, options); 
    data.aligned = true; 
    data.alignedLength = data.speciesMap[data.getFirstSpeciesName()]->alignedSequence.size();
    // EVALUATION
    double SP2 = SP(data);
    cout << "SP distance - KIMURA : " << SP2 << endl;


    // REFINEMENT
    double oldSP = 0; // convergence trackers
    double newSP = 0;
    bool converged = false;
    vector<edgeStruct> edgeVector;
	map<vertex_descriptor, bool> visited;   
	map<vertex_descriptor, vector<Species*> > nodeProfiles; 

    // Iterate through edges of the tree and store them in a more convenient form
    graph_traits<Tree_graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(guideTree2.g); ei != ei_end; ++ei)
    {
        struct edgeStruct targetEdge;
        targetEdge.source = source(*ei, guideTree2.g);
        targetEdge.target = target(*ei, guideTree2.g);
        edgeVector.push_back(targetEdge);
    }

	// BEGIN REFINING
    while (!converged){
        cout << endl;
        for (int i = 0; i < edgeVector.size(); ++i)
        {
            cout << "\x1b[A" << "Refining, " << 100.0*(1.0*i/edgeVector.size()) << "\% done.               " << endl;
            // Helper map for aligning recursively
            graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
            for(tie(vi, vi_end) = vertices(guideTree2.g); vi != vi_end; ++vi)
            {
                visited[*vi] = false;
            }

            // Cut the tree into two
            remove_edge(edgeVector[i].source, edgeVector[i].target, guideTree2.g);

            // Try to improve the subtrees
            alignWithTreeRecursive(guideTree2,edgeVector[i].source,visited,nodeProfiles,options);
            alignWithTreeRecursive(guideTree2,edgeVector[i].target,visited,nodeProfiles,options);
            
            // Reconnect the tree
            pair<edge_descriptor, bool> reconnectingEdge = add_edge(edgeVector[i].source, edgeVector[i].target, guideTree2.g);
            put(edge_weight_t(), guideTree2.g, reconnectingEdge.first, 1.0); //Edge weights shouldn't matter

            // Realign to the whole tree
            alignWithTree(guideTree2, options);
        }
    	// Check for convergence in improvments
    	newSP = SP(data);
    	if (newSP-oldSP < options.convergenceValue) {
    		converged = true;
    	}
    	oldSP = newSP;
    	cout << "SP distance - REFINEMENT : " << oldSP << endl;
    }    
    cout << "MUSCLE done." << endl;
}

/*
    Evalutates total tree distance with aligned sequences
*/
double Muscle::SP(SequenceData &data) 
{
	double totalDistance = 0.0;
    vector<string> names = data.getAllSpeciesNames();
	for(int i = 0; i < names.size(); ++i)
    {
		for(int j = i + 1; j < names.size(); ++j)
		{
		    int numNonGapPositions = 0;
		    int numMatches = 0;
		    string seq1 = data.speciesMap[names[i]]->alignedSequence;
		    string seq2 = data.speciesMap[names[j]]->alignedSequence;
		    for(int k = 0; k < seq1.length(); ++k)
		    {
				if(seq1[k] != '-' && seq2[k] != '-') 
                {
                    numNonGapPositions++;
				    if(seq1[k] == seq2[k]) 
                        numMatches++;
                }
		    }
		    totalDistance += 1 - 1.0*numMatches/numNonGapPositions;
		}
    }
    return	totalDistance;
}

/*
    Metric for unaligned sequence relations
    Given some sequences, convert every sequence into it's amino acid form, and identify the frequencies of
    chains of amino acids of a given length, compared to other sequences using euclidean squared distance
    to compute their differences.
*/
map<string, map<string, double> > Muscle::kmer(SequenceData &data, 
						const Options &options)
{
    map<string, map<string, double> > distances;
    vector<string> names = data.getAllSpeciesNames();
    translate = AminoAlphabet();

    //TRANSLATE
    for(int i = 0; i < names.size(); ++i)
    {
      string raw = data.speciesMap[names[i]]->rawSequence;
      for (int k = 0; k < raw.size(); k+=3)
      {
		if (raw.size() - (k+3) > 0) // leftover characters aren't shoed into amino acids, probably sloppy
	  	{
		  data.speciesMap[names[i]]->aminoSequence += translate.amnAlpha[raw.substr(k,3)];
		}
      }
    }

    //COUNT "KMERS/WORDS/KTUPLES"
    typedef map<string, double>::iterator iter;
    for(int i = 0; i < names.size(); ++i)
    {
        // get the counts
    	int totalKmers = 0;
		for(int j = 0; j < data.speciesMap[names[i]]->aminoSequence.size()-options.kmerSize+1; ++j){
			data.speciesMap[names[i]]->kmerFreqs[data.speciesMap[names[i]]->aminoSequence.substr(j,options.kmerSize)] += 1;
			totalKmers += 1;
		}
		// alter value from raw count to frequency
		for(iter iterator = data.speciesMap[names[i]]->kmerFreqs.begin(); iterator != data.speciesMap[names[i]]->kmerFreqs.end(); iterator++){
			iterator->second = (float) data.speciesMap[names[i]]->kmerFreqs[iterator->first] / totalKmers;
		}
    }

    //COMPARE KMER FREQUENCIES FOR EUC DISTANCE
    double distance = 0;
    for(int i = 0; i < names.size(); ++i)
    {
		for(int j = i + 1; j < names.size(); ++j)
		{
			// for each entry in the name[i] species kmer freq map calculate the squared difference and add it to the distance total
			for(iter iterator = data.speciesMap[names[i]]->kmerFreqs.begin(); iterator != data.speciesMap[names[i]]->kmerFreqs.end(); iterator++){
				// if the entry matches to a frequency in the name[j] species kmer freq map
				if (data.speciesMap[names[j]]->kmerFreqs.find(iterator->first) != data.speciesMap[names[j]]->kmerFreqs.end()) {
					distance += pow(iterator->second - data.speciesMap[names[j]]->kmerFreqs[iterator->first], 2);
				}
				// else it doesnt match
				else {
					distance += pow(iterator->second, 2);
				}
			}  

		    distances[names[i]][names[j]] = distance;
		    distances[names[j]][names[i]] = distance;
		    distance = 0; //reset for next round
	   	}
    }
    return distances;
}


/*
    Metric for aligned sequence relations
    Calculates sequence difference using the proportion of mismatched bases (nucleotides here, could use amino acids)
    insertion characters are ignored if matched to an insertion, else counted as a mismatch.
*/
map<string, map<string, double> > Muscle::kimura(SequenceData &data, 
						const Options &options)
{
    map<string, map<string, double> > distances;
    vector<string> names = data.getAllSpeciesNames();

    //COMPARE
    for(int i = 0; i < names.size(); ++i)
    {
        string target = data.speciesMap[names[i]]->alignedSequence;
		for(int j = i + 1; j < names.size(); ++j)
		{
		    string compare = data.speciesMap[names[j]]->alignedSequence;
		    float D = 0;
		    for(int k = 0; k < target.size(); ++k) 
		    {
		    	if (target[k] == '-' && compare[k] == '-') { //insertion character ignoring
		    		continue;
		    	}
		    	else if (target[k] != compare[k]) {
		    		D++;
		    	}
		    	else {
		    		continue;
		    	}
			}
			D = D / target.size();

            double distance = -log(1-D);
		    distances[names[i]][names[j]] = distance;
		    distances[names[j]][names[i]] = distance;
	   	}
    }
    return distances;
}

/*
    Insert a root, set up for alignment, call alignwithtreerecursive
*/
void Muscle::alignWithTree(const Tree &constTopologyIn,
						const Options &options)
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

    put(edge_weight_t(), topologyIn.g, firstEdge, 1.0); //Edge weights shouldn't matter
    put(edge_weight_t(), topologyIn.g, secondEdge, 1.0);
    
    //Track what vertices have been visited
    map<vertex_descriptor, bool> visited;
    graph_traits<Tree_graph>::vertex_iterator vi, vi_end;
    for(tie(vi, vi_end) = vertices(topologyIn.g); vi != vi_end; ++vi)
    {
		visited[*vi] = false;
    }
    map<vertex_descriptor, vector<Species*> > nodeProfiles;

    //Proceed with alignment
    alignWithTreeRecursive(topologyIn, nodeId, visited, nodeProfiles, options);

    delete dummyRoot;    
}

/*
    Gather information on the tree's structure, call alignprofiles to make profiles and push that up
*/
void Muscle::alignWithTreeRecursive(const Tree &topology,
				      vertex_descriptor rootNode,
				      map<vertex_descriptor, bool> &visited,
				      map<vertex_descriptor, vector<Species*> > &nodeProfiles,
				      const Options &options)
{
    //Mark current 'root' visited
    visited[rootNode] = true;

    //Get the vector of unvisited neighbors, which will either
    //have length of 2 or 0.
    graph_traits<Tree_graph>::adjacency_iterator ai, ai_end;
    vector <vertex_descriptor> unvisited_neighbors;
    for (tie(ai, ai_end) = adjacent_vertices(rootNode, topology.g); ai != ai_end; ++ai)
    {
		if (! visited[*ai]) {
	    	unvisited_neighbors.push_back(*ai);
	    }
    }

    //Base case when looking at a leaf node:
    //add the sequence character to charSet
    //and return length 0.
    if (unvisited_neighbors.size() == 0)
    {
		vector<Species*> curProfile;
		curProfile.push_back(get(Vertex_t(), topology.g, rootNode));
		nodeProfiles[rootNode] = curProfile;
		return;
    }
    
    //Internal node, 2 neighbors found
    vertex_descriptor childOne = unvisited_neighbors[0];
    vertex_descriptor childTwo = unvisited_neighbors[1];
    alignWithTreeRecursive(topology, childOne, visited, nodeProfiles, options);
    alignWithTreeRecursive(topology, childTwo, visited, nodeProfiles, options);
    
    vector<Species*> childOneProfile = nodeProfiles[childOne];
    vector<Species*> childTwoProfile = nodeProfiles[childTwo];
    
    //Build the profile from leaves to dummy root
    alignProfiles(childOneProfile, childTwoProfile, options);

    vector<Species*> curProfile;
    for(int i = 0; i < childOneProfile.size(); ++i) curProfile.push_back(childOneProfile[i]);
    for(int i = 0; i < childTwoProfile.size(); ++i) curProfile.push_back(childTwoProfile[i]);

    nodeProfiles[rootNode] = curProfile;
}

/*
    Create the alignment profile from the given two profiles
*/
void Muscle::alignProfiles(vector<Species*> profile1,
			     vector<Species*> profile2,
			     const Options &options)
{
    //length of strings in profile1.
    int n1 = profile1[0]->alignedSequence.length();
    //length of strings in profile2.
    int n2 = profile2[0]->alignedSequence.length();
    
    //create the alignment matrix
    vector<vector<double> > alignMatrix(n1+1,vector<double>(n2+1,0));
    vector<vector<int> > pathMatrix(n1,vector<int>(n2,0));

    for(int i = 1; i < n1+1; ++i) //for every character in sequence 1 
    {
		for(int j = 1; j < n2+1; ++j) //for every character in sequence 2 
		{
		    //choiceVal represents which alignment choice we made
		    //this corresponds to where in the distance matrix we
		    //computed this value from
		    //0 is coming from the right
		    //1 is coming diagonally down/right
		    //2 is coming down
		    
			int choiceVal = -1;

            double choiceOne = alignMatrix[i][j-1] + options.mProfileGapReward;
            double choiceTwo = alignMatrix[i-1][j-1] + logExpecScore(profile1, profile2, i-1, j-1, options);
            double choiceThree = alignMatrix[i-1][j] + options.mProfileGapReward; 

		    vector<double> choiceVec;
		    choiceVec.push_back(choiceOne);
		    choiceVec.push_back(choiceTwo);
		    choiceVec.push_back(choiceThree);
		    vector<double>::iterator max = max_element(choiceVec.begin(), choiceVec.end()); //select greatest value of choice
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

    //align the sequences (from sequence tail to head)
    int curI = n1-1;
    int curJ = n2-1;
    while(true)
    {
		if(curI < 0) // no remaining characters in I
		{
		    for(int k = 0; k < profile1.size(); ++k)
		    {
				profile1[k]->alignedSequence = '-' + profile1[k]->alignedSequence;
		    	curJ --;
		    }
		}
		else if(curJ < 0) // no remaining characters in J
		{
		    for(int k = 0; k < profile2.size(); ++k)
		    {
				profile2[k]->alignedSequence = '-' + profile2[k]->alignedSequence;
		    	curI --;
			}		
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
    Scoring for profile construction, evaluated the probability of change to characters or insertions or no change
*/
double Muscle::logExpecScore(vector <Species*> profile1, 
			vector <Species*> profile2, 
			int profile1Col, 
			int profile2Col,
			const Options &options)
{    
    vector<double> profile1Freqs;
    vector<double> profile2Freqs;
    
    //map for tracking the number of occurrences of each character (INCLUDING GAPS)
    map<char, int> profile1CharCount;
    map<char, int> profile2CharCount;


    for (int i = 0; i < profile1.size(); ++i)
		profile1CharCount[profile1[i]->alignedSequence.at(profile1Col)]++;
    for (int i = 0; i < profile2.size(); ++i)
		profile2CharCount[profile2[i]->alignedSequence.at(profile2Col)]++;

	//Generate frequencies of each character in the presented alphabet
    for (int i = 0; i < options.mAlphabet.size(); ++i)
    {
		profile1Freqs.push_back(1.0 * profile1CharCount[options.mAlphabet[i]] / profile1.size());
		profile2Freqs.push_back(1.0 * profile2CharCount[options.mAlphabet[i]] / profile2.size());
    }

    
    double runningSum = 0;
    //for every combination of the mAlphabet
    for (int i = 0; i < options.mAlphabet.size()-1; ++i) //exclude the last one which should be the gaps
    {
		for (int j = 0; j < options.mAlphabet.size()-1; ++j) //exclude the last one which should be the gaps
		{
		    if (options.mSimilarities.find(options.mAlphabet[i]) == options.mSimilarities.end()) 
		    {
				continue;
		    }
		    else
		    {
				if(options.mSimilarities.at(options.mAlphabet[i]).find(options.mAlphabet[j]) == options.mSimilarities.at(options.mAlphabet[i]).end())
			    	continue;
		    }
		    double sim = options.mSimilarities.at(options.mAlphabet[i]).at(options.mAlphabet[j]);
		    runningSum += profile1Freqs[i] * profile2Freqs[j] * sim; // Honestly should be using the similarity between amino acids but that brings into scope the issue of aligning aminos not nucleotides
		}
    }
    return runningSum;     
}




