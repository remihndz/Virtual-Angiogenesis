import sys
from utils import IntercapillaryDistanceMap
from numpy import array, savetxt
import numpy as np


assert len(sys.argv)>1, "No files found. Use 'python IntercapillaryArea.py List-Of-PNG-Files'."

plotResults = len(sys.argv[1:])==1
count = 0
means, variances = [], []

f = open("ICD.dat", 'w')
f.write("File Mean Std\n")
for filename in sys.argv[1:]:
    
    D = IntercapillaryDistanceMap(filename, plotResult=plotResults)
    mean, std = D.mean(), D.std()   # Mean/std of Intercapillary distance, includes distance from FAZ to vessels
    print(f"Mean (std) for {filename}: {mean}, ({std})")
    count+=1
    means.append(mean)
    variances.append(std**2)

    f.write(f"{filename} {mean} {std}\n")
    
    

mean = array(means).mean()
# std  = sum(variances)**.5
std = array(means).std()
print(u'Mean\u00B1std of intercapillary distance out of', count, 'images (in microns):', mean, u"\u00B1", std)

