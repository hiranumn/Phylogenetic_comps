from subprocess import Popen
# import scipy.stats as stats
# from scipy import std, mean
# import random

#generate all file identities
batches = [[str(i) + "-" + str(j) for j in range(1,15)] for i in range(4,9)]
fileFrom = "./syntheticExperiment/unalignedFiles/200high/"
fileTo = "./syntheticExperiment/drive5/"
filePath = "./syntheticExperiment/alignedFiles/200high/"


# def bootstrapCI(data, N=10000):
#     dist = []
#     count = 0
#     for i in range(N):
#         sample = []
#         #Sample data with replacement
#         for j in range(len(data)):
#             randIndex = random.randint(0, len(data)-1)
#             sample.append(data[randIndex])
#         tStat = (mean(sample) - mean(data))/((std(sample))/sqrt(len(sample)))
#         dist.append(tStat)
#     q2 = stats.scoreatpercentile(dist,97.5)
#     q1 = stats.scoreatpercentile(dist,2.5)
#     subVal = q2*std(data)/sqrt(len(data))
#     addVal = -q1*std(data)/sqrt(len(data))
#     return subVal, addVal

def cleanFormat(item,fileFrom,fileTo):
	inFile = fileFrom+item+".unaln"
	input_file = open(inFile, 'r')
	outFile = fileTo+item
	output_file = open(outFile, 'w')

	for line in input_file:
		if line[0:4] == "Taxa":
			output_file.write(">" + line)
		else:
			output_file.write(line)

	input_file.close()
	output_file.close()


def cleanOutput(item, fileTo, filePath):
	inFile = fileTo+item+"DRIVE5.afa"
	input_file = open(inFile, 'r')
	outFile = filePath+item+"MSC.aln"
	output_file = open(outFile, 'w')

	sequence = ""
	firstline = True
	for line in input_file:

		if line[0:1] == ">":
			if firstline:
				output_file.write(line[1:])
				firstline = False
			else:
				output_file.write(sequence+"\n")
				output_file.write(line[1:])
			sequence = ""
		else:
			sequence += line.rstrip()

	output_file.write(sequence)
	input_file.close()
	output_file.close()


def getCommand(item, fileFrom, fileTo, filePath):
	source = fileTo+item
	return "drive5MUSCLE -in " + source + " -out " + fileTo + item + "DRIVE5.afa"



for batch in batches:
	for item in batch:
		cleanFormat(item,fileFrom,fileTo)

for batch in batches:
    commands = [getCommand(item,fileFrom,fileTo,filePath) for item in batch]
    processes = [Popen(cmd, shell=True) for cmd in commands]
    for p in processes: 
    	p.wait()

for batch in batches:
	for item in batch:
		cleanOutput(item,fileTo,filePath)
