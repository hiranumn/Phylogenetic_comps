from subprocess import Popen

batches = [[str(i) + "-" + str(j) for j in range(1,15)] for i in range(4,9)]
inPath = "./syntheticExperiment/trueTrees/200high/"
outPath = "./syntheticExperiment/unalignedFiles/200high/"


def getCommand(filename, inPath, outPath, algprefix):
    return "./phylo -" + algprefix + " " + inPath + filename + ".tree" + " " + outPath + filename + ".unaln"

for algprefix in ["CS"]:
    for batch in batches:
        commands = [getCommand(x, inPath, outPath, algprefix) for x in batch]
        processes = [Popen(cmd, shell=True) for cmd in commands]
        for p in processes: p.wait()
