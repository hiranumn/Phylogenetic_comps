/*
* Authors: Evan Albright, Jack Hessel, Naozumi Hiranuma, Cody Wang
* Last modified: 2/16/2014
* centerStar.cpp
* Contains an implementation of the center star multiple sequence alignment algorithm
*/

#include "centerStar.h"


/*
    Main entry point to the center star algorithm
    Locates the sequence which minimizes the total distance between all other sequences, names it the center
    aligns everything to the center and progressively combines the 2 sequence alignments into a single multiple alignment
*/
void CenterStar::align(SequenceData &data, const Options &options)
{
    map<string, map<string, int> > distances;
    vector<string> names = data.getAllSpeciesNames();

    cout << "Identifying center star sequence..." << endl;

    //Compute all optimal distances between sequences
    for(int i = 0; i < names.size(); ++i)
    {
        for(int j = i + 1; j < names.size(); ++j)
        {
            pair<string, string> alignment = pairwiseAlign(data.speciesMap[names[i]]->rawSequence,
                                   data.speciesMap[names[j]]->rawSequence,
                                   false,
                                   options);

            int misMatchOrIndel = 0;
            for(int k = 0; k < alignment.first.length(); ++k)
            {
                if (alignment.first.at(k) != alignment.second.at(k))
                {
                    misMatchOrIndel++;
                }
            }
            int distance = misMatchOrIndel;

            distances[names[i]][names[j]] = distance;
            distances[names[j]][names[i]] = distance; 

        }
    }


    //Identify center string
    int total;
    string centerStar;
    int min = numeric_limits<int>::max(); //infinity value
    for(int i = 0; i < names.size(); ++i)
    {
        total = 0;
        for(int j = 0; j < names.size(); ++j)
        {
            if (i != j)
                total += distances[names[i]][names[j]];
        }
        if (min > total)
        {
            centerStar = names[i];
            min = total;
        }
    }

    cout << "Computing alignments to center star..." << endl;

    //Align everything to the center string
    vector<vector<string> > alignments;
    for(int i = 0; i < names.size(); ++i)
    {
        if (names[i] != centerStar)
        {
            pair<string, string> alignment = pairwiseAlign(data.speciesMap[centerStar]->rawSequence,
                                   data.speciesMap[names[i]]->rawSequence,
                                   false,
                                   options);
            vector<string> tempPair;
            tempPair.push_back(alignment.first);

            //store the alignment to the center star
            data.speciesMap[names[i]]->alignedSequence = alignment.second;
            tempPair.push_back(names[i]);
            alignments.push_back(tempPair);
        }
    }

    cout << "Combining alignments..." << endl;

    //combine them all
    while (alignments.size() > 1)
    {
        vector<string> firstAlignment = alignments[0]; //first in the list
        vector<string> secondAlignment = alignments[1]; //second in the list
        vector<string> combinedAlignments; //store the combinations

        //"Introduce" space in the center star
        pair<string, string> alignment = pairwiseAlign(firstAlignment[0],
                               secondAlignment[0],
                               false,
                               options);

        string combinedStar = alignment.first; //grab the combined star, should be the same as alignment.second
        combinedAlignments.push_back(combinedStar); //store the new star

        //pairwise align the sequences of the two groups to the combined center
        for (int i = 1; i < firstAlignment.size(); ++i)
        {
            pair<string, string> aligned = pairwiseAlign(combinedStar,
                                   data.speciesMap[firstAlignment[i]]->alignedSequence,
                                   true,
                                   options);

            //store the altered alignment to the center star
            data.speciesMap[firstAlignment[i]]->alignedSequence = aligned.second;
            combinedAlignments.push_back(firstAlignment[i]);
        }
        for (int i = 1; i < secondAlignment.size(); ++i)
        {
            pair<string, string> aligned = pairwiseAlign(combinedStar,
                                   data.speciesMap[secondAlignment[i]]->alignedSequence,
                                   true,
                                   options);
            //store the altered alignment to the center star
            data.speciesMap[secondAlignment[i]]->alignedSequence = aligned.second;
            combinedAlignments.push_back(secondAlignment[i]);
        }

        alignments.erase(alignments.begin(),alignments.begin()+2); //scrap the two used up alignments

        alignments.push_back(combinedAlignments); //add in the conglomerate
    }
    //store the center star's alignment
    data.speciesMap[centerStar]->alignedSequence = alignments[0][0];

    cout << "Center star aligment completed." << endl;

}

/*
    Given two strings finds the optimal alignment for them
    if the flag fitToStar is true it will restrict the alignment from altering the first string, used in
    matching strings to the conglomerate star sequences.
*/
pair<string, string> CenterStar::pairwiseAlign(string s1, 
                         string s2, 
                         bool fitToStar,
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
            double choiceOne;
            choiceOne = alignMatrix[i][j-1] + options.csGapReward;
            if (fitToStar)
                choiceOne = alignMatrix[i][j-1] - 100000; //penalty to ensure it never alters the first string
            double choiceTwo = alignMatrix[i-1][j-1];
            if(s1.at(i-1) == s2.at(j-1)) //if we have a character match
                choiceTwo += options.csMatchReward;
            double choiceThree = alignMatrix[i-1][j] + options.csGapReward;
            
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
        if(curI < 0 && curJ < 0) 
            break;
    }
    alignedOne = string(alignedOne.rbegin(), alignedOne.rend());
    alignedTwo = string(alignedTwo.rbegin(), alignedTwo.rend());
    
    pair<string, string> returnPair(alignedOne, alignedTwo);

    return returnPair;
}