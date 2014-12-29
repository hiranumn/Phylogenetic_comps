from subprocess import Popen

batches = [["random" + str(i) + "COX1_" + str(j) for j in range(1,15)] for i in [5,8]]

filepath = "realExperiment/alignedSequences/"
filepath2 = "realExperiment/outputTrees/"

def getCommand(filename, filepath, algprefix):
    return "./phylo -" + algprefix + " " + filepath + filename + "CW.aln " + filepath2 +  filename + algprefix + "_CW.txt options.txt" 

for algprefix in ["MC"]:
    for batch in batches:
        commands = [getCommand(x, filepath, algprefix) for x in batch]
        processes = [Popen(cmd, shell=True) for cmd in commands]
        for p in processes: p.wait()
