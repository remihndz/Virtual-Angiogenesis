import sys
from utils import *

assert(len(sys.argv) > 1)
listOfFiles = sys.argv[1:]
print("The mean and standard deviation of the vessel density is:", StatisticsMultipleTrees(listOfFiles))

    
