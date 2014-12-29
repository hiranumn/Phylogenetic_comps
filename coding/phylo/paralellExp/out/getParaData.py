#5MPP.tree.time

outfile = open("data.out", 'w')

for alg in ["MPP","MPH","MLP","MLH","MC"]:
    outfile.write(alg + "\n")
    for species in range(4,8):
        sfilename = "seq/" + str(species) + alg + ".tree.time"
        pfilename = "para/" + str(species) + alg + ".tree.time"
        sfile = open(sfilename)
        pfile = open(pfilename)
        sval = float(sfile.readlines()[0].strip())
        pval = float(pfile.readlines()[0].strip())
        rat = sval/pval
        outfile.write(str(rat) + "\n")
