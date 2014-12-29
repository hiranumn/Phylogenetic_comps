from scipy import std, mean
from scipy.stats.mstats import mquantiles
import matplotlib.pyplot as plt
from pylab import polyfit
from math import sqrt
from scipy.stats import t
import numpy as np 
import pylab 
import scipy.stats as stats
import random

# For QQ plots
# stats.probplot(myThirdQ, dist="norm", plot=pylab)
# pylab.savefig("normalDistTest/"+str(upd)+".eps")
# pylab.figure()
numSpecies = range(4, 9)
algs = ["NJ", "MPP", "MPH", "MLP", "MLH", "MC"]


fontBig = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
        
fontMed = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
        
fontSm = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}


def bootstrapCI(data, N=10000):
	dist = []
	count = 0
	for i in range(N):
		sample = []
		#Sample data with replacement
		for j in range(len(data)):
			randIndex = random.randint(0, len(data)-1)
			sample.append(data[randIndex])
		tStat = (mean(sample) - mean(data))/((std(sample))/sqrt(len(sample)))
		dist.append(tStat)
	q2 = stats.scoreatpercentile(dist,97.5)
	q1 = stats.scoreatpercentile(dist,2.5)
	subVal = q2*std(data)/sqrt(len(data))
	addVal = -q1*std(data)/sqrt(len(data))
	print "done!"
	return subVal, addVal

def saveFig(statsFileName, title, outFileName):
	random.seed(8701)
	curFile = open(statsFileName)
        #maps from {"alg" -> {speciesNum -> (ave, -CI, +CI)}}
	dataDict = {}
	myLines = curFile.readlines()
	i = 0
	while 1:
		if i == len(myLines):
			break
		curAlg = myLines[i].strip()
                #maps from {speciesNum -> (ave, -CI, +CI)} 
		curDict = {}
		for j in range(4,9):
			i += 1
			curData = [float(x) for x in myLines[i].strip().split()]
			average = mean(curData)
			subVal, addVal = bootstrapCI(curData)
			curDict[j] = (average, subVal, addVal)
			dataDict[curAlg] = curDict
		i += 1
		
	plt.figure()
	plt.title(title, **fontBig)
	plt.xlabel("Number of Species", **fontMed)
	plt.ylabel("Runtime (clock cycles, in log scale)", **fontMed)
		
	for alg, speciesDict in dataDict.iteritems():
		myVals = [0,0,0,0,0]
		myPlusStd = [0,0,0,0,0]
		myMinusStd = [0,0,0,0,0]
		for speciesNum, stats in speciesDict.iteritems():
			myVals[speciesNum-4] = stats[0]
			myMinusStd[speciesNum-4] = stats[1]
			myPlusStd[speciesNum-4] = stats[2]
		plt.errorbar(numSpecies, myVals,
		     yerr=[myMinusStd, myPlusStd], label=alg)


	plt.legend(loc=2, prop={'size':11})
	plt.yscale('log')
	plt.xticks(numSpecies)
	plt.savefig(outFileName)

def main():
	algs = ["CS", "CW", "MSC"]
	fileNames = ["timeStats_" + alg for alg in algs]
	titles = ["Runtime vs. Number of Species for " + alg for alg in algs]
	outfiles = [alg + ".eps" for alg in algs]
	for i in range(len(algs)):
		saveFig(fileNames[i], titles[i], outfiles[i])

if __name__ == "__main__":
	main()
