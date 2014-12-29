from subprocess import Popen
# import scipy.stats as stats
# from scipy import std, mean
# import random

#generate all file identities
batches = [["random" + str(i) + "COX1_" + str(j) for j in range(1,15)] for i in [5,8,12]]
fileFrom = "../../data/real_data/experiment1/"
fileTo = "./drive5/"
filePath = "./alignedSequences/"


def cleanFormat(item,fileFrom,fileTo):
	inFile = fileFrom+item+".txt"
	input_file = open(inFile, 'r')
	outFile = fileTo+item
	output_file = open(outFile, 'w')

	oddEven = 0
	for line in input_file:
		if oddEven == 0:
			output_file.write(">" + line)
			oddEven += 1
		else:
			output_file.write(line)
			oddEven -= 1

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



# for batch in batches:
# 	for item in batch:
# 		cleanFormat(item,fileFrom,fileTo)

# for batch in batches:
#     commands = [getCommand(item,fileFrom,fileTo,filePath) for item in batch]
#     processes = [Popen(cmd, shell=True) for cmd in commands]
#     for p in processes: 
#     	p.wait()

for batch in batches:
	for item in batch:
		cleanOutput(item,fileTo,filePath)
