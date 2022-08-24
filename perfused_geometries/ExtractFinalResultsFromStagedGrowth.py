import os
import sys


for rootdir in sys.argv[1:]:
    print ("Exctracting simulations from %s" %(rootdir))
    for subdir, dirs, files in os.walk(rootdir):
        
        eligibleFiles = []
        stages = []
        for file in files:
            if not file.endswith(('.vtp', '.root', '.conf', '.dat')):
                eligibleFiles.append(file[:-4])

                for i, token in enumerate(eligibleFiles[-1]):
                    if token == '_':
                        stages.append(int(eligibleFiles[-1][i+1:]))

        if stages:
            fileIndex = stages.index(max(stages))
            fileName = os.path.join(subdir, eligibleFiles[fileIndex])
            print("For sim", subdir, "copy files", fileName+'.cco/.vtp', "as", subdir+'.cco/.vtp')
            
            os.system("cp %s %s" %(fileName+'.cco', subdir+'.cco'))
            os.system("cp %s %s" %(fileName+'.vtp', subdir+'.vtp'))
