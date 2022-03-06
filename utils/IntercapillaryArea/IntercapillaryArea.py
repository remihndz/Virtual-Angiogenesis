import sys
from utils import IntercapillaryDistanceMap
from numpy import array

assert(len(sys.argv)>1)

plotResults = False
count = 0
means, variances = [], []
for filename in sys.argv[1:]:
    
    D = IntercapillaryDistanceMap(filename, plotResult=plotResults)
    mean, std = D.mean(), D.std()   # Mean/std of Intercapillary distance, includes distance from FAZ to vessels
    print(mean, std)
    count+=1
    means.append(mean)
    variances.append(std**2)

mean = array(means).mean()
std  = sum(variances)**.5
print(u'Mean\u00B1std of intercapillary distance out of', count, 'images:', mean, u"\u00B1", std)



