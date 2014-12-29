//main.cpp
#include "main.h"

int main(int argc, char** argv)
{  
    Options myOptions;

    //string species[3] = {"5","8","12"};
    //string alignments[2] = {"CW","CS","DRIVE5"};
	string species[5] = {"4","5","6","7","8"};
    string alignments[3] = {"CW","CS","DRIVE5"};

	ofstream output;
    //output.open("realExperiment/alignmentEvaluation.txt");
    output.open("syntheticExperiment/alignmentEvaluation.txt");


    //for each alignment style
    for(int i = 0; i < 3; ++i)
    {
	    output << alignments[i] << "\n";

	    // for each value of species
		for (int k = 0; k < 5; ++k)
		{
			output << "\t" << species[k] << "\n\t";

	    	//for each entry of the style
			for(int j = 1; j < 15; ++j)
			{
				//string fileName = to_string("realExperiment/alignedSequences/random" + to_string(species[k]) + "COX1_" + to_string(j) + to_string(alignments[i]) + ".aln");
				string fileName = to_string("syntheticExperiment/alignedFiles/200high/" + to_string(species[k]) + "-" + to_string(j) + to_string(alignments[i]) + ".aln");
				SequenceData *sd = new SequenceData(fileName);
				double td = sd->totalDistance();
				output << "\t" << td;
			}

			output << "\n\n";
		}
    }
    output << "hi" << endl;
    output.close();
}
