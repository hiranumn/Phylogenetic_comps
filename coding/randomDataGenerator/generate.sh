#!/bin/bash

#To run this guy just enter "source generate.sh" while in the current directory

for i in {1..14} #number of runs for each mix of variables
	do
	for species in {4..8} #number of species
		do
		for mutation in {0.00005} #mutation rate, (should have like a low, medium, and high probability, don't exceed mutation rate = 0.1!!)
			do
			./randomDataGenerator $i $species 200 $mutation
     		sleep 2
		done
	done
done

