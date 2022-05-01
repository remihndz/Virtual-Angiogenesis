import sys
from utils import ReadTree, VBC
import numpy as np

vbc_means, vbc_stds = [], []
plot = len(sys.argv[1:])==1     # Only generate plot for a single file
for ccoFile in sys.argv[1:]:
    data, connectivity = ReadTree(ccoFile)
    vbc = VBC(data, connectivity, plot=plot)
    mean, std = vbc.mean(), vbc.std()
    vbc_means.append(mean)
    vbc_stds.append(std)
    
vbc_means = np.array(vbc_means)
vbc_stds  = np.array(vbc_stds)

print("Individual means:\n", vbc_means)
print("Individual standard deviations:\n", vbc_stds)
print("VBC of population = %f +- %f" %(vbc_means.mean(), vbc_means.std()))
    


