import sys
from utils import *
from math import pi

assert(len(sys.argv) > 1)
listOfFiles = sys.argv[1:]

domainAreaInCm = 0.09 - 0.04**2*pi

VAD, VPI, VDI, VCI, VSD = StatisticsMultipleTrees(listOfFiles, domainAreaInCm=domainAreaInCm)
print(f"Found {len(listOfFiles):d} files.")
print("Metric  |  Mean+-std")
print(f"VAD     |  {VAD.mean()} +- {VAD.std()}")
print(f"VPI     |  {VPI.mean()} +- {VPI.std()}")
print(f"VDI     |  {VDI.mean()} +- {VDI.std()}")
print(f"VCI     |  {VCI.mean()} +- {VCI.std()}")
print(f"VSD     |  {VSD.mean()} +- {VSD.std()}")

