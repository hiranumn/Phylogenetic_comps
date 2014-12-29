outFile=open("statsMSC","w")
for alg in ["NJ","MPP","MPH","MLP","MLH","MC"]:
    outFile.write(alg+'\n')
    for species in range(4,9):
        for trial in range(1,15):
            try:
                curFile=open(str(species)+'-'+str(trial)+alg+"_MSC.tree.time")
                curVal=float(curFile.readlines()[0].strip())
                outFile.write(str(curVal)+' ')
            except IOError:
                continue
        outFile.write('\n')
            
    
